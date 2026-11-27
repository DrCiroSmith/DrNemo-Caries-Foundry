#!/usr/bin/env python3
"""
Run the retargeted docking campaign against 8FJC (GtfB).

Design notes worth reading before you trust the output
------------------------------------------------------
1. Ligands are built by biophysics_engine.build_macrocycle_3d, imported unchanged
   from the repository root. The comparison against the old 3AIC numbers only
   means something if the molecules are built by the same code.

2. That function fixes the ETKDG conformer seed at 42 internally, so tier 4
   (which varies the conformer seed) cannot use it. Tier 4 uses a local copy with
   the seed parameterised and nothing else changed. Before any tier 4 run, the
   local copy is made to build at seed 42 and its output is compared byte for
   byte with the imported function's output. If they differ, the run stops. An
   unverified re-implementation is not worth having.

3. Every Vina call passes an explicit --seed. The existing pipeline passes none,
   so Vina seeds itself and the published scores are not reproducible runs of a
   deterministic protocol. Seed is treated here as a factor to measure.

4. Output paths are unique per ligand, receptor, box, docking seed and conformer
   seed. The existing pipeline reuses any pose file that already exists, which
   with a seed sweep would silently collapse every seed onto the first result.

5. The existing pipeline falls back to --autobox when fixed-box Vinardo rescoring
   fails to produce a number. That changes the scoring context. The fallback is
   kept, because dropping it would lose rows, but every use is recorded in the
   result so that affected rows can be excluded at analysis time.

Writes only into retarget_gtfb/out. Nothing outside this directory is touched.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
import time
import uuid
from pathlib import Path

HERE = Path(__file__).parent.resolve()
REPO_ROOT = HERE.parent
OUT = HERE / "out"
LIGAND_DIR = OUT / "ligands"
POSE_DIR = OUT / "poses"
RESULTS = OUT / "results.jsonl"
CONFIG = json.loads((HERE / "config.json").read_text())
SEQUENCES = json.loads((HERE / "sequences.json").read_text())

sys.path.insert(0, str(REPO_ROOT))


def md5_of(path) -> str:
    return hashlib.md5(Path(path).read_bytes()).hexdigest()


def sha256_of(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def file_tag(path, n: int = 8) -> str:
    """Short content tag for a file, for use inside a deterministic output path.

    Two different receptor or ligand files can therefore never write to the same
    pose path, whatever they are called.
    """
    return sha256_of(path)[:n]


def die(msg: str) -> None:
    print(f"\nSTOP: {msg}", file=sys.stderr)
    sys.exit(1)


# --------------------------------------------------------------------------- #
# Preflight
# --------------------------------------------------------------------------- #
def preflight() -> dict:
    env = {"python": sys.version.split()[0], "platform": platform.platform()}

    for exe in ("vina", "vina_split", "obabel"):
        path = shutil.which(exe)
        if not path:
            die(f"{exe} is not on PATH. The campaign cannot run without it.")
        env[f"{exe}_path"] = path

    res = subprocess.run(["vina", "--version"], capture_output=True, text=True)
    version_text = (res.stdout + res.stderr).strip()
    env["vina_version_text"] = version_text
    required = CONFIG["docking"]["required_version_substring"]
    if required not in version_text:
        die(f"this package expects AutoDock Vina {required} and found:\n  {version_text}\n"
            f"Scores are not comparable across Vina versions, so the run is refusing to "
            f"start rather than produce numbers that cannot be compared with the old ones. "
            f"Install {required}, or change required_version_substring in config.json "
            f"deliberately and note it in the write-up.")

    try:
        import rdkit, meeko, numpy  # noqa: F401
        env["rdkit"] = rdkit.__version__
        env["meeko"] = getattr(meeko, "__version__", "unknown")
        env["numpy"] = numpy.__version__
    except ImportError as exc:
        die(f"a required python package is missing: {exc}. "
            f"Activate the project virtualenv first.")

    stamp = OUT / "receptor_validation.json"
    boxes_path = OUT / "boxes.json"
    if not stamp.exists():
        die("out/receptor_validation.json is missing. Run validate_receptor.py first.")
    if not boxes_path.exists():
        die("out/boxes.json is missing. Run prepare_receptor.py first.")
    validation = json.loads(stamp.read_text())
    boxes = json.loads(boxes_path.read_text())
    if boxes.get("source_sha256") != validation.get("sha256"):
        die("out/boxes.json was built from a different receptor file than the one currently "
            "validated. Re-run prepare_receptor.py.")
    for key in ("receptor_pdbqt", "native_ligand_pdbqt"):
        if not Path(boxes[key]).exists():
            die(f"{boxes[key]} is missing. Re-run prepare_receptor.py.")

    env["receptor_sha256"] = validation["sha256"]
    env["receptor_title"] = validation["title"]
    return {"env": env, "boxes": boxes, "validation": validation}


# --------------------------------------------------------------------------- #
# Ligand construction
# --------------------------------------------------------------------------- #
def build_local(seq: str, out_pdbqt: str, out_pdb: str | None,
                target_center, conf_seed: int) -> int:
    """
    Local copy of biophysics_engine.build_macrocycle_3d with the ETKDG seed
    parameterised. Nothing else is changed. Verified against the imported
    function at seed 42 by verify_builder_equivalence() before use.
    """
    import numpy as np
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from meeko import MoleculePreparation, PDBQTWriterLegacy

    mol = Chem.MolFromSequence(seq)
    if mol is None:
        raise ValueError(f"Invalid peptide sequence: {seq}")
    mol = Chem.AddHs(mol)

    N = len(seq)

    def find_sg(m):
        a1 = aN = None
        for a in m.GetAtoms():
            info = a.GetPDBResidueInfo()
            if info and info.GetResidueName().strip() == "CYS" and info.GetName().strip() == "SG":
                if info.GetResidueNumber() == 1:
                    a1 = a.GetIdx()
                elif info.GetResidueNumber() == N:
                    aN = a.GetIdx()
        return a1, aN

    cys1_sg, cysN_sg = find_sg(mol)
    if cys1_sg is not None and cysN_sg is not None:
        atoms_to_remove = []
        for s_idx in (cys1_sg, cysN_sg):
            for nbr in mol.GetAtomWithIdx(s_idx).GetNeighbors():
                if nbr.GetSymbol() == "H":
                    atoms_to_remove.append(nbr.GetIdx())
                    break
        emol = Chem.EditableMol(mol)
        for h_idx in sorted(atoms_to_remove, reverse=True):
            emol.RemoveAtom(h_idx)
        mol = emol.GetMol()
        cys1_sg, cysN_sg = find_sg(mol)
        emol = Chem.EditableMol(mol)
        emol.AddBond(cys1_sg, cysN_sg, Chem.BondType.SINGLE)
        mol = emol.GetMol()

    Chem.SanitizeMol(mol)

    params = AllChem.ETKDGv3()
    params.randomSeed = conf_seed                      # the only change
    embed_res = AllChem.EmbedMolecule(mol, params)
    if embed_res != 0:
        embed_res = AllChem.EmbedMolecule(mol, useRandomCoords=True,
                                          maxAttempts=5000, randomSeed=conf_seed)
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
    except Exception:
        try:
            AllChem.UFFOptimizeMolecule(mol, maxIters=500)
        except Exception:
            pass

    if target_center is not None:
        conf = mol.GetConformer()
        positions = conf.GetPositions()
        centroid = np.mean(positions, axis=0)
        translation = np.array(target_center) - centroid
        for i in range(mol.GetNumAtoms()):
            conf.SetAtomPosition(i, positions[i] + translation)

    if out_pdb:
        os.makedirs(os.path.dirname(os.path.abspath(out_pdb)), exist_ok=True)
        Chem.MolToPDBFile(Chem.RemoveHs(mol), out_pdb)

    os.makedirs(os.path.dirname(os.path.abspath(out_pdbqt)), exist_ok=True)
    preparator = MoleculePreparation()
    preparator.prepare(mol)
    pdbqt_str = PDBQTWriterLegacy().write_string(preparator.setup)[0]
    Path(out_pdbqt).write_text(pdbqt_str, encoding="utf-8")
    return mol.GetNumHeavyAtoms()


_equivalence_checked = False


def verify_builder_equivalence(centre) -> None:
    """Build the same molecule both ways at seed 42 and require identical bytes."""
    global _equivalence_checked
    if _equivalence_checked:
        return
    from biophysics_engine import build_macrocycle_3d

    probe = "CAHWRPWHC"
    tmp = OUT / "_builder_check"
    tmp.mkdir(parents=True, exist_ok=True)
    a, b = tmp / "imported.pdbqt", tmp / "local.pdbqt"
    build_macrocycle_3d(probe, str(a), None, tuple(centre))
    build_local(probe, str(b), None, tuple(centre), 42)
    ha, hb = (hashlib.sha256(p.read_bytes()).hexdigest() for p in (a, b))
    if ha != hb:
        die("the local conformer builder does not reproduce "
            "biophysics_engine.build_macrocycle_3d at seed 42.\n"
            f"  imported: {ha}\n  local:    {hb}\n"
            "Tier 4 would therefore be measuring a difference between two builders rather "
            "than a difference between conformer seeds. Refusing to continue.")
    (tmp / "equivalence.json").write_text(json.dumps(
        {"probe_sequence": probe, "sha256_both": ha, "verified": True}, indent=2) + "\n")
    print(f"  builder equivalence verified at seed 42 (sha256 {ha[:16]})")
    _equivalence_checked = True


def torsion_count(pdbqt: Path) -> int:
    return sum(1 for l in pdbqt.read_text().splitlines() if l.startswith("BRANCH"))


GLUE_TYPES = {"CG0", "CG1", "CG2", "CG3", "G0", "G1", "G2", "G3"}


def has_glue_atoms(pdbqt: Path) -> bool:
    """True when Meeko opened a macrocycle ring for sampling.

    Meeko marks an opened ring by inserting a pair of pseudo-atoms whose AutoDock
    types are CGn and Gn. Read that type as the last whitespace-separated token on
    the ATOM line, not by a fixed column slice: the type sits at the end of the
    record after the partial charge, its column varies between Meeko versions, and
    a two-character slice silently truncates CG0 to CG. An earlier version of this
    function used a [77:79] slice and I reported the opposite of the truth from it.
    """
    for l in pdbqt.read_text(errors="replace").splitlines():
        if l.startswith(("ATOM", "HETATM")):
            tokens = l.split()
            if tokens and tokens[-1].upper() in GLUE_TYPES:
                return True
    return False


def terminal_state(pdbqt: Path) -> dict:
    """
    Record whether the built ligand has free termini, read off the file Vina is given.

    This exists because the repository's two models of the same molecule disagree.
    Chem.MolFromSequence, which the ligand builder calls, produces a peptide with a free
    N-terminal amine and a free C-terminal carboxylate: the macrocycle is closed by the
    Cys1-Cys9 disulfide, not head to tail. solubility_shield.calc_net_charge_at_ph takes
    the opposite view, skipping both terminal ionisations on the stated assumption that
    they are consumed by cyclization. One of those is wrong about this molecule.

    Nothing here changes the molecule. Changing it would make the new numbers
    incomparable with the old ones, which is the one thing this package is for. What it
    does is write the capping state into the per-ligand record, so that the ambiguity
    cannot be inherited silently by whatever reads these results next.
    """
    try:
        names, charges = [], []
        for line in pdbqt.read_text(errors="replace").splitlines():
            if line[:6] in ("ATOM  ", "HETATM"):
                names.append(line[12:16].strip().upper())
                try:
                    charges.append(float(line[70:76]))
                except ValueError:
                    pass
        return {
            "free_c_terminal_carboxylate_OXT_present": "OXT" in names,
            "atom_name_count": len(names),
            "summed_partial_charge": round(sum(charges), 3) if charges else None,
            "capping": "uncapped, free termini" if "OXT" in names else "no OXT atom found",
            "note": "Built by Chem.MolFromSequence, which does not cap termini. "
                    "solubility_shield.calc_net_charge_at_ph models this molecule as if the "
                    "termini were consumed by cyclization; on a disulfide macrocycle they are "
                    "not. The two disagree, and this field records which molecule was docked.",
        }
    except Exception as exc:  # a diagnostic must never stop a run
        return {"error": f"{type(exc).__name__}: {exc}"}


def ligand_for(seq: str, centre, conf_seed: int) -> tuple[Path, dict]:
    LIGAND_DIR.mkdir(parents=True, exist_ok=True)
    pdbqt = LIGAND_DIR / f"{seq}_conf{conf_seed}.pdbqt"
    meta_path = LIGAND_DIR / f"{seq}_conf{conf_seed}.meta.json"
    if pdbqt.exists() and meta_path.exists():
        return pdbqt, json.loads(meta_path.read_text())

    if conf_seed == 42:
        from biophysics_engine import build_macrocycle_3d
        heavy = build_macrocycle_3d(seq, str(pdbqt),
                                    str(LIGAND_DIR / f"{seq}_conf42.pdb"), tuple(centre))
        builder = "biophysics_engine.build_macrocycle_3d (imported unchanged)"
    else:
        verify_builder_equivalence(centre)
        heavy = build_local(seq, str(pdbqt), None, tuple(centre), conf_seed)
        builder = "retarget_gtfb local copy, ETKDG seed parameterised only"

    meta = {"sequence": seq, "conformer_seed": conf_seed, "heavy_atoms": heavy,
            "builder": builder,
            "pdbqt_sha256": hashlib.sha256(pdbqt.read_bytes()).hexdigest(),
            "rotatable_branches": torsion_count(pdbqt),
            "macrocycle_opened_by_meeko": has_glue_atoms(pdbqt),
            "terminal_state": terminal_state(pdbqt)}
    meta_path.write_text(json.dumps(meta, indent=2) + "\n")
    return pdbqt, meta


# --------------------------------------------------------------------------- #
# Docking
# --------------------------------------------------------------------------- #
def parse_vina_modes(pose: Path) -> list[dict]:
    modes = []
    for line in pose.read_text().splitlines():
        if line.startswith("REMARK VINA RESULT:"):
            parts = line.split()
            modes.append({"affinity_kcal": float(parts[3]),
                          "rmsd_lb": float(parts[4]), "rmsd_ub": float(parts[5])})
    return modes


def dock_one(ligand_pdbqt: Path, receptor: str, centre, size,
             dock_seed: int, tag: str, cpu: int, exhaustiveness: int,
             vina_version_text: str = "", allow_overwrite: bool = False) -> dict:
    POSE_DIR.mkdir(parents=True, exist_ok=True)
    pose = POSE_DIR / f"{tag}.pdbqt"
    log = POSE_DIR / f"{tag}.vina.log"

    # A pose file already sitting at a deterministic path, with no recorded result to
    # match it, is an orphan from an interrupted or edited run. Overwriting it silently
    # is how one result ends up standing in for another.
    if pose.exists() and not allow_overwrite:
        return {"ok": False, "seconds": 0.0,
                "error": f"refusing to overwrite an existing result at {pose.name}. "
                         f"Nothing in out/results.jsonl claims it, so it cannot be "
                         f"attributed to a run. Delete it, or re-run with --overwrite.",
                "overwrite_refused": True}

    run_id = uuid.uuid4().hex
    cx, cy, cz = centre
    sx, sy, sz = size

    cmd = ["vina", "--receptor", receptor, "--ligand", str(ligand_pdbqt),
           "--center_x", str(cx), "--center_y", str(cy), "--center_z", str(cz),
           "--size_x", str(sx), "--size_y", str(sy), "--size_z", str(sz),
           "--cpu", str(cpu), "--exhaustiveness", str(exhaustiveness),
           "--seed", str(dock_seed), "--out", str(pose)]
    # Provenance travels inside the result record, because the two files he sends back
    # are all I will have. A record that cannot be traced to a command line, a seed and
    # a pair of input hashes is not evidence of anything.
    prov = {
        "run_id": run_id,
        "vina_command": cmd,
        "vina_command_string": " ".join(cmd),
        "vina_version_text": vina_version_text,
        "exhaustiveness": exhaustiveness,
        "dock_seed": dock_seed,
        "cpu": cpu,
        "box_centre": [cx, cy, cz],
        "box_size": [sx, sy, sz],
        "receptor_path": str(receptor),
        "receptor_md5": md5_of(receptor),
        "receptor_sha256": sha256_of(receptor),
        "ligand_path": str(ligand_pdbqt),
        "ligand_md5": md5_of(ligand_pdbqt),
        "ligand_sha256": sha256_of(ligand_pdbqt),
        "pose_path_planned": str(pose),
        "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }

    t0 = time.time()
    res = subprocess.run(cmd, capture_output=True, text=True)
    elapsed = time.time() - t0
    prov["finished_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    log.write_text(" ".join(cmd) + f"\n# run_id {run_id}\n\n"
                   + res.stdout + "\n--- stderr ---\n" + res.stderr)
    if res.returncode != 0:
        return {"ok": False, "error": f"vina exit {res.returncode}: {res.stderr.strip()[:400]}",
                "seconds": round(elapsed, 1), **prov}

    modes = parse_vina_modes(pose)
    if not modes:
        return {"ok": False, "error": "no REMARK VINA RESULT lines in output",
                "seconds": round(elapsed, 1), **prov}

    # split to rescore mode 1 with Vinardo, same as the existing pipeline
    split_prefix = str(pose).replace(".pdbqt", "_split_")
    subprocess.run(["vina_split", "--input", str(pose), "--ligand", split_prefix],
                   capture_output=True, text=True)
    mode1 = Path(f"{split_prefix}1.pdbqt")

    vinardo, vinardo_mode = None, None
    if mode1.exists():
        fixed = ["vina", "--receptor", receptor, "--ligand", str(mode1),
                 "--center_x", str(cx), "--center_y", str(cy), "--center_z", str(cz),
                 "--size_x", str(sx), "--size_y", str(sy), "--size_z", str(sz),
                 "--scoring", "vinardo", "--score_only"]
        out = subprocess.run(fixed, capture_output=True, text=True).stdout
        for line in out.splitlines():
            if "Estimated Free Energy of Binding" in line:
                try:
                    vinardo = float(line.split(":")[1].split()[0])
                    vinardo_mode = "fixed_box"
                except (ValueError, IndexError):
                    pass
                break
        if vinardo is None:
            auto = ["vina", "--receptor", receptor, "--ligand", str(mode1),
                    "--scoring", "vinardo", "--score_only", "--autobox"]
            out = subprocess.run(auto, capture_output=True, text=True).stdout
            for line in out.splitlines():
                if "Estimated Free Energy of Binding" in line:
                    try:
                        vinardo = float(line.split(":")[1].split()[0])
                        vinardo_mode = "autobox_fallback"
                    except (ValueError, IndexError):
                        pass
                    break

    from biophysics_engine import measure_disulfide_bond
    return {
        "ok": True,
        **prov,
        "vinardo_command_string": " ".join(fixed) if mode1.exists() else None,
        "seconds": round(elapsed, 1),
        "vina_kcal": modes[0]["affinity_kcal"],
        "vina_modes": modes,
        "mode_gap_kcal": (round(modes[1]["affinity_kcal"] - modes[0]["affinity_kcal"], 3)
                          if len(modes) > 1 else None),
        "n_modes": len(modes),
        "vinardo_kcal": vinardo,
        "vinardo_scoring_context": vinardo_mode,
        "ss_bond_angstrom": measure_disulfide_bond(str(pose)),
        "pose_pdbqt": str(pose),
        "pose_sha256": hashlib.sha256(pose.read_bytes()).hexdigest(),
    }


# --------------------------------------------------------------------------- #
# Job list
# --------------------------------------------------------------------------- #
def build_jobs(boxes: dict, tiers: list[int]) -> list[dict]:
    d = CONFIG["docking"]
    wide = boxes["boxes"]["wide_28"]["size"]
    tight = boxes["boxes"]["tight_auto"]["size"]
    centre = boxes["centre"]
    gtfb = boxes["receptor_pdbqt"]
    native = boxes["native_ligand_pdbqt"]
    contenders = SEQUENCES["contenders"]
    null_a = SEQUENCES["null_a_wt_permutations"]
    null_b = SEQUENCES["null_b_lead_permutations"]
    off = {"1SMD": str(REPO_ROOT / "demo_data" / "receptors" / "1SMD_receptor.pdbqt"),
           "1LZ1": str(REPO_ROOT / "demo_data" / "receptors" / "1LZ1_receptor.pdbqt")}
    off_centres = {"1SMD": (8.778, 49.724, 17.278), "1LZ1": (10.319, 8.503, 25.960)}

    jobs: list[dict] = []

    exhaustiveness = d["exhaustiveness"]
    receptor_tags: dict[str, str] = {}

    def add(kind, ident, seq, receptor_id, receptor, rcentre, box_id, size, dseed, cseed,
            prebuilt=None, tier=None):
        # Protonation is named explicitly rather than left implicit. There is exactly one
        # state in this package and this is what it actually is: build_macrocycle_3d goes
        # through Chem.MolFromSequence and AddHs, which applies no pH model at all and
        # leaves the termini free. Naming it means a second state added later cannot
        # collide with this one in a path.
        protonation = ("native_cocrystal_as_deposited" if prebuilt
                       else "as_built_no_pH_model_free_termini")
        if receptor not in receptor_tags:
            try:
                receptor_tags[receptor] = file_tag(receptor)
            except OSError:
                receptor_tags[receptor] = "nofile00"
        jobs.append({"tier": tier, "kind": kind, "id": ident, "sequence": seq,
                     "receptor_id": receptor_id, "receptor": receptor,
                     "receptor_tag": receptor_tags[receptor],
                     "centre": list(rcentre), "box_id": box_id, "size": list(size),
                     "dock_seed": dseed, "conformer_seed": cseed,
                     "exhaustiveness": exhaustiveness, "protonation": protonation,
                     "prebuilt_ligand": prebuilt})

    if 1 in tiers:
        for s in d["seeds"]:
            add("calibration", "acarbose_native", None, "8FJC", gtfb, centre,
                "wide_28", wide, s, None, prebuilt=native, tier=1)
            for c in contenders:
                add("contender", c["id"], c["sequence"], "8FJC", gtfb, centre,
                    "wide_28", wide, s, 42, tier=1)
    if 2 in tiers:
        for s in d["seeds_for_null_sets"]:
            for n in null_a:
                add("null_a", n["id"], n["sequence"], "8FJC", gtfb, centre,
                    "wide_28", wide, s, 42, tier=2)
            for n in null_b:
                add("null_b", n["id"], n["sequence"], "8FJC", gtfb, centre,
                    "wide_28", wide, s, 42, tier=2)
    if 3 in tiers:
        for s in d["seeds_for_null_sets"]:
            add("calibration", "acarbose_native", None, "8FJC", gtfb, centre,
                "tight_auto", tight, s, None, prebuilt=native, tier=3)
            for c in contenders:
                add("contender", c["id"], c["sequence"], "8FJC", gtfb, centre,
                    "tight_auto", tight, s, 42, tier=3)
            for n in null_a[:12]:
                add("null_a", n["id"], n["sequence"], "8FJC", gtfb, centre,
                    "tight_auto", tight, s, 42, tier=3)
            for n in null_b[:12]:
                add("null_b", n["id"], n["sequence"], "8FJC", gtfb, centre,
                    "tight_auto", tight, s, 42, tier=3)
    if 4 in tiers:
        subset = [c for c in contenders if c["role"] in
                  ("reference_scaffold_P2A", "wild_type_lead")]
        for cseed in CONFIG["ligand_preparation"]["conformer_seeds"]:
            if cseed == 42:
                continue
            for c in subset:
                add("contender", c["id"], c["sequence"], "8FJC", gtfb, centre,
                    "wide_28", wide, 1, cseed, tier=4)
            for n in null_a[:10]:
                add("null_a", n["id"], n["sequence"], "8FJC", gtfb, centre,
                    "wide_28", wide, 1, cseed, tier=4)
            for n in null_b[:10]:
                add("null_b", n["id"], n["sequence"], "8FJC", gtfb, centre,
                    "wide_28", wide, 1, cseed, tier=4)
    if 6 in tiers:
        # The calcium control. Same ligands, same box, same seeds as tier 3, but against
        # the receptor that keeps the conserved structural Ca2+ instead of stripping it.
        # Skipped without comment here if prepare_receptor.py could not build or validate
        # that receptor; it writes the reason into out/boxes.json.
        ca = boxes.get("calcium_condition") or {}
        if ca.get("usable") and ca.get("receptor_pdbqt"):
            gtfb_ca = ca["receptor_pdbqt"]
            for s in d["seeds_for_null_sets"]:
                add("calibration", "acarbose_native", None, "8FJC_CA", gtfb_ca, centre,
                    "wide_28", wide, s, None, prebuilt=native, tier=6)
                for c in contenders:
                    add("contender", c["id"], c["sequence"], "8FJC_CA", gtfb_ca, centre,
                        "wide_28", wide, s, 42, tier=6)
                for n in null_a[:12]:
                    add("null_a", n["id"], n["sequence"], "8FJC_CA", gtfb_ca, centre,
                        "wide_28", wide, s, 42, tier=6)
                for n in null_b[:12]:
                    add("null_b", n["id"], n["sequence"], "8FJC_CA", gtfb_ca, centre,
                        "wide_28", wide, s, 42, tier=6)
    if 5 in tiers:
        for s in d["seeds_for_null_sets"]:
            for rid, rpath in off.items():
                if not Path(rpath).exists():
                    continue
                for c in contenders:
                    add("off_target", c["id"], c["sequence"], rid, rpath,
                        off_centres[rid], "wide_28", wide, s, 42, tier=5)
    return jobs


# --------------------------------------------------------------------------- #
# Naming
#
# In the repository this package corrects, two docked results appear under five
# filenames: the 3AIC result at -8.191 kcal/mol under three names, and the 1SMD
# result at -7.343 under two. Nothing inside any of those files records a seed, a
# command line, a run id or a time, so from the files alone there is no way to say
# whether that was one run saved several times or several runs that collapsed onto
# one. It took the repository's version history to settle that, and it turned out to
# be one result under three names; the cause is unexplained and is not asserted. The
# -8.191 is the lead's published headline number. (A third duplicated
# content in that directory is the input ligand conformer under three names, which
# is expected from one builder at one seed and is not a defect.) The two functions
# below are why none of this can happen here.
#
# identity()  the scientific identity of a run: every parameter that can change
#             the number. Resume and deduplication both key on this, so the same
#             run requested by two tiers is docked once, not twice.
# job_key()   identity plus the tier that first asked for it. This is the stem of
#             every file the run writes. Nothing is hand-labelled.
#
# The pose path additionally carries a short content tag of the ligand file, which
# is only known once the ligand is built, so two different ligand files cannot
# land on one pose path even under identical parameters.
# --------------------------------------------------------------------------- #
IDENTITY_FIELDS = ("kind", "id", "receptor_id", "protonation", "box_id",
                   "dock_seed", "conformer_seed", "exhaustiveness", "receptor_tag")


def identity(j: dict) -> str:
    return "|".join([str(j["kind"]), str(j["id"]), str(j["receptor_id"]),
                     str(j["protonation"]), str(j["box_id"]),
                     f"d{j['dock_seed']}", f"c{j['conformer_seed']}",
                     f"ex{j['exhaustiveness']}", f"r{j['receptor_tag']}"])


def job_key(j: dict) -> str:
    return f"t{j['tier']}|" + identity(j)


def path_stem(j: dict, ligand_tag: str) -> str:
    """Deterministic filesystem stem for one run. Same parameters, same stem, always."""
    stem = (job_key(j) + f"|L{ligand_tag}").replace("|", "__")
    stem = "".join(ch if (ch.isalnum() or ch in "_-.") else "-" for ch in stem)
    if len(stem) > 150:  # keep it deterministic and inside every filesystem limit
        stem = stem[:120] + "__" + hashlib.sha256(stem.encode()).hexdigest()[:16]
    return stem


def dedupe_jobs(jobs: list[dict]) -> list[dict]:
    """One run per scientific identity, tagged with every tier that asked for it."""
    seen: dict[str, dict] = {}
    for j in jobs:
        ident = identity(j)
        if ident in seen:
            seen[ident]["also_requested_by_tiers"].append(j["tier"])
        else:
            j = dict(j)
            j["also_requested_by_tiers"] = []
            seen[ident] = j
    return list(seen.values())


# --------------------------------------------------------------------------- #
# Naming self-check, with no test-runner dependency
#
# The full suite is test_provenance.py and needs pytest. This is the subset that
# gates the docking, written so that a machine without pytest still cannot start a
# nine-hour run on a naming scheme that has not been checked.
# --------------------------------------------------------------------------- #
def naming_self_check() -> int:
    checks: list[tuple[str, bool, str]] = []

    def base(**over):
        j = {"tier": 2, "kind": "contender", "id": "mutant_P2A_CAHWRPWHC",
             "receptor_id": "8FJC", "receptor": "/r/8FJC.pdbqt",
             "receptor_tag": "1111aaaa", "centre": [1.0, 2.0, 3.0], "box_id": "wide_28",
             "size": [28.0] * 3, "dock_seed": 1, "conformer_seed": 42,
             "exhaustiveness": 16, "protonation": "as_built_no_pH_model_free_termini",
             "prebuilt_ligand": None}
        j.update(over)
        return j

    checks.append(("the same parameters always give the same path",
                   path_stem(base(), "lig1") == path_stem(base(), "lig1"), ""))

    for field, value in (("tier", 3), ("id", "wild_type"), ("receptor_id", "8FJC_CA"),
                         ("protonation", "explicit_pH_5_2"), ("box_id", "tight_auto"),
                         ("dock_seed", 2), ("conformer_seed", 101),
                         ("exhaustiveness", 32), ("receptor_tag", "2222bbbb")):
        changed = path_stem(base(**{field: value}), "lig1") != path_stem(base(), "lig1")
        checks.append((f"changing {field} changes the output path", changed, ""))

    checks.append(("changing the ligand file changes the output path",
                   path_stem(base(), "lig1") != path_stem(base(), "lig2"), ""))
    checks.append(("the path is bounded and filesystem safe",
                   len(path_stem(base(id="X" * 400), "lig1")) <= 150
                   and all(c.isalnum() or c in "_-." for c in path_stem(base(), "lig1")),
                   ""))
    checks.append(("long names still give distinct paths",
                   path_stem(base(id="Y" * 400), "lig1")
                   != path_stem(base(id="Y" * 399 + "Z"), "lig1"), ""))
    checks.append(("identity ignores the tier, so resume works across tiers",
                   identity(base(tier=2)) == identity(base(tier=6))
                   and job_key(base(tier=2)) != job_key(base(tier=6)), ""))
    deduped = dedupe_jobs([base(tier=2), base(tier=6),
                           base(tier=3, box_id="tight_auto")])
    checks.append(("one run requested by two tiers is docked once",
                   len(deduped) == 2
                   and any(j["also_requested_by_tiers"] == [6] for j in deduped), ""))

    # the glue-atom read, which an earlier version of this file got wrong
    tmp = OUT / ".naming_self_check"
    tmp.mkdir(parents=True, exist_ok=True)
    opened = tmp / "opened.pdbqt"
    opened.write_text("ATOM      1  N   CYS A   1       0.0  0.0  0.0  0.00  0.00  -0.3 N\n"
                      "ATOM      2  CG  UNL     1       1.0  1.0  1.0  0.00  0.00  +0.0 CG0\n")
    rigid = tmp / "rigid.pdbqt"
    rigid.write_text("ATOM      1  N   CYS A   1       0.0  0.0  0.0  0.00  0.00  -0.3 N\n"
                     "ATOM      5  CG  LEU A   2       1.0  1.0  1.0  0.00  0.00  +0.0 C\n")
    checks.append(("an opened macrocycle ring is detected", has_glue_atoms(opened), ""))
    checks.append(("a side-chain CG atom is not mistaken for ring glue",
                   not has_glue_atoms(rigid), ""))
    for f in (opened, rigid):
        f.unlink()
    tmp.rmdir()

    # the failure this is all for: one content under several names
    poses = REPO_ROOT / "data" / "poses"
    if poses.is_dir():
        by: dict[str, list[str]] = {}
        for f in sorted(poses.glob("*.pdbqt")):
            by.setdefault(hashlib.sha256(f.read_bytes()).hexdigest(), []).append(f.name)
        dups = {h: v for h, v in by.items() if len(v) > 1}
        note = (f"{len(dups)} contents under more than one name in data/poses, "
                f"worst case {max((len(v) for v in dups.values()), default=0)} names "
                f"for one result" if dups else "data/poses holds no duplicates")
        checks.append(("the duplicate-content check runs on the existing repository",
                       True, note))

    print("Naming and provenance self-check")
    failed = 0
    for name, ok, note in checks:
        print(f"  {'PASS' if ok else 'FAIL'}  {name}" + (f"  [{note}]" if note else ""))
        failed += 0 if ok else 1
    if failed:
        print(f"\n{failed} naming checks failed. Not docking.")
        return 1
    print(f"\nAll {len(checks)} naming checks passed. "
          f"Run test_provenance.py for the full provenance suite if pytest is available.")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tiers", default="1,2,3,4,5,6",
                    help="comma-separated tiers to run, in order. Default all. Tier 6 is "
                         "the calcium control and is skipped automatically if the "
                         "Ca-retained receptor could not be built or this Vina build "
                         "refused it; out/boxes.json records why.")
    ap.add_argument("--cpu", type=int, default=CONFIG["docking"]["cpu"])
    ap.add_argument("--check-naming", action="store_true",
                    help="run the naming and provenance checks and stop. No Vina, no "
                         "pytest, about a second. RUN_ME.sh runs this before docking.")
    ap.add_argument("--dry-run", action="store_true",
                    help="list the work and stop without docking")
    ap.add_argument("--overwrite", action="store_true",
                    help="allow a run to overwrite an existing pose file that no recorded "
                         "result claims. Off by default: an unclaimed file at a "
                         "deterministic path means an interrupted or edited run, and "
                         "overwriting it silently is how one result comes to stand in for "
                         "another. Nothing in the analysis needs this flag.")
    args = ap.parse_args()
    if args.check_naming:
        return naming_self_check()
    tiers = [int(t) for t in args.tiers.split(",") if t.strip()]

    print("Preflight")
    pre = preflight()
    print(f"  vina: {pre['env']['vina_version_text']}")
    print(f"  receptor: {pre['env']['receptor_title']}")
    print(f"  receptor sha256: {pre['env']['receptor_sha256'][:16]}")

    raw_jobs = build_jobs(pre["boxes"], tiers)
    jobs = dedupe_jobs(raw_jobs)
    shared = sum(1 for j in jobs if j["also_requested_by_tiers"])
    if len(raw_jobs) != len(jobs):
        print(f"\n  {len(raw_jobs) - len(jobs)} requested runs are the same run under "
              f"two tiers; each is docked once and reported in both.")

    done = set()
    if RESULTS.exists():
        for line in RESULTS.read_text().splitlines():
            if line.strip():
                try:
                    rec = json.loads(line)
                except json.JSONDecodeError:
                    continue
                # Resume keys on scientific identity, not on the tier label, so a run
                # already done under one tier is not repeated under another.
                if "identity" in rec:
                    done.add(rec["identity"])
                elif "key" in rec:  # a file written before identity existed
                    done.add(rec["key"].split("|", 1)[-1])
    todo = [j for j in jobs if identity(j) not in done]

    print(f"\n{len(jobs)} distinct docking runs in tiers {tiers}; "
          f"{len(done)} already recorded; {len(todo)} to run.")
    if shared:
        print(f"  {shared} of them serve more than one tier.")
    if args.dry_run:
        for j in todo[:20]:
            print("  " + job_key(j))
        if len(todo) > 20:
            print(f"  ... and {len(todo) - 20} more")
        return 0
    if not todo:
        print("Nothing to do. out/results.jsonl is already complete for these tiers.")
        return 0

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "manifest.json").write_text(json.dumps({
        "generated_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "config": CONFIG, "environment": pre["env"], "boxes": pre["boxes"],
        "sequence_provenance": SEQUENCES["provenance"],
        "tiers_requested": tiers, "total_runs_planned": len(jobs),
        "requested_runs_before_deduplication": len(raw_jobs),
        "naming": {
            "identity_fields": list(IDENTITY_FIELDS),
            "path_stem_rule": "t<tier>__<identity fields>__L<ligand content tag>, "
                              "non-alphanumerics replaced, hashed past 150 characters",
            "why": "In the repository this package corrects, two docked results appear "
                   "under five filenames and nothing inside those files can attribute "
                   "them. Every path here is derived from the run's own parameters and "
                   "the content of its two input files, so no two runs can share one.",
        },
    }, indent=2) + "\n")

    t_start = time.time()
    completed = failed = 0
    paths_written: dict[str, str] = {}
    with RESULTS.open("a") as fh:
        for i, j in enumerate(todo, 1):
            key = job_key(j)
            if j["prebuilt_ligand"]:
                lig, lmeta = Path(j["prebuilt_ligand"]), {"source": "native co-crystal ligand"}
            else:
                try:
                    lig, lmeta = ligand_for(j["sequence"], pre["boxes"]["centre"],
                                            j["conformer_seed"])
                except Exception as exc:
                    rec = {"key": key, "identity": identity(j), **j,
                           "ok": False, "error": f"ligand build failed: {exc}"}
                    fh.write(json.dumps(rec) + "\n")
                    fh.flush()
                    failed += 1
                    print(f"[{i}/{len(todo)}] {key} LIGAND FAILED: {exc}")
                    continue

            ligand_tag = file_tag(lig)
            stem = path_stem(j, ligand_tag)
            if stem in paths_written:
                die(f"two distinct runs derived the same output path {stem}.\n"
                    f"  first:  {paths_written[stem]}\n  second: {key}\n"
                    f"This is a bug in the naming scheme, not in your machine. Stopping "
                    f"rather than letting one run overwrite another.")
            paths_written[stem] = key

            out = dock_one(lig, j["receptor"], j["centre"], j["size"], j["dock_seed"],
                           tag=stem, cpu=args.cpu,
                           exhaustiveness=j["exhaustiveness"],
                           vina_version_text=pre["env"]["vina_version_text"],
                           allow_overwrite=args.overwrite)
            rec = {"key": key, "identity": identity(j), "path_stem": stem,
                   "ligand_tag": ligand_tag, **j, "ligand_meta": lmeta,
                   "recorded_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), **out}
            fh.write(json.dumps(rec) + "\n")
            fh.flush()

            if out["ok"]:
                completed += 1
                print(f"[{i}/{len(todo)}] {key}  vina={out['vina_kcal']:+.3f} "
                      f"vinardo={out['vinardo_kcal']} ({out['seconds']}s)")
            else:
                failed += 1
                print(f"[{i}/{len(todo)}] {key}  FAILED: {out['error']}")

            if i == 3:
                per = (time.time() - t_start) / 3
                left = per * (len(todo) - 3)
                print(f"\n  Timing: about {per:.0f}s per run, so roughly "
                      f"{left / 3600:.1f} h for the remaining {len(todo) - 3}.")
                print("  Safe to stop with Ctrl-C at any point. Re-running resumes "
                      "from out/results.jsonl.\n")

    print(f"\nDone. {completed} succeeded, {failed} failed, "
          f"{(time.time() - t_start) / 3600:.2f} h elapsed.")
    print(f"Results: {RESULTS}")
    if failed:
        print("Some runs failed. They are recorded with ok=false and will be retried on the "
              "next run only if you delete their lines; otherwise they stay as failures.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
