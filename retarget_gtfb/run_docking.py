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


def has_glue_atoms(pdbqt: Path) -> bool:
    """Meeko opens a macrocycle by inserting CG/G pseudo-atoms. Absence means a rigid ring."""
    for l in pdbqt.read_text().splitlines():
        if l.startswith(("ATOM", "HETATM")) and l[77:79].strip() in ("CG", "G"):
            return True
        if l.startswith(("ATOM", "HETATM")) and l[12:16].strip().startswith(("CG", "G")) \
                and l[17:20].strip() == "UNL":
            return True
    return False


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
            "macrocycle_opened_by_meeko": has_glue_atoms(pdbqt)}
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
             dock_seed: int, tag: str, cpu: int, exhaustiveness: int) -> dict:
    POSE_DIR.mkdir(parents=True, exist_ok=True)
    pose = POSE_DIR / f"{tag}.pdbqt"
    log = POSE_DIR / f"{tag}.vina.log"
    cx, cy, cz = centre
    sx, sy, sz = size

    cmd = ["vina", "--receptor", receptor, "--ligand", str(ligand_pdbqt),
           "--center_x", str(cx), "--center_y", str(cy), "--center_z", str(cz),
           "--size_x", str(sx), "--size_y", str(sy), "--size_z", str(sz),
           "--cpu", str(cpu), "--exhaustiveness", str(exhaustiveness),
           "--seed", str(dock_seed), "--out", str(pose)]
    t0 = time.time()
    res = subprocess.run(cmd, capture_output=True, text=True)
    elapsed = time.time() - t0
    log.write_text(res.stdout + "\n--- stderr ---\n" + res.stderr)
    if res.returncode != 0:
        return {"ok": False, "error": f"vina exit {res.returncode}: {res.stderr.strip()[:400]}",
                "seconds": round(elapsed, 1)}

    modes = parse_vina_modes(pose)
    if not modes:
        return {"ok": False, "error": "no REMARK VINA RESULT lines in output",
                "seconds": round(elapsed, 1)}

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

    def add(kind, ident, seq, receptor_id, receptor, rcentre, box_id, size, dseed, cseed,
            prebuilt=None):
        jobs.append({"kind": kind, "id": ident, "sequence": seq,
                     "receptor_id": receptor_id, "receptor": receptor,
                     "centre": list(rcentre), "box_id": box_id, "size": list(size),
                     "dock_seed": dseed, "conformer_seed": cseed, "prebuilt_ligand": prebuilt})

    if 1 in tiers:
        for s in d["seeds"]:
            add("calibration", "acarbose_native", None, "8FJC", gtfb, centre,
                "wide_28", wide, s, None, prebuilt=native)
            for c in contenders:
                add("contender", c["id"], c["sequence"], "8FJC", gtfb, centre,
                    "wide_28", wide, s, 42)
    if 2 in tiers:
        for s in d["seeds_for_null_sets"]:
            for n in null_a:
                add("null_a", n["id"], n["sequence"], "8FJC", gtfb, centre,
                    "wide_28", wide, s, 42)
            for n in null_b:
                add("null_b", n["id"], n["sequence"], "8FJC", gtfb, centre,
                    "wide_28", wide, s, 42)
    if 3 in tiers:
        for s in d["seeds_for_null_sets"]:
            add("calibration", "acarbose_native", None, "8FJC", gtfb, centre,
                "tight_auto", tight, s, None, prebuilt=native)
            for c in contenders:
                add("contender", c["id"], c["sequence"], "8FJC", gtfb, centre,
                    "tight_auto", tight, s, 42)
            for n in null_a[:12]:
                add("null_a", n["id"], n["sequence"], "8FJC", gtfb, centre,
                    "tight_auto", tight, s, 42)
            for n in null_b[:12]:
                add("null_b", n["id"], n["sequence"], "8FJC", gtfb, centre,
                    "tight_auto", tight, s, 42)
    if 4 in tiers:
        subset = [c for c in contenders if c["role"] in
                  ("reference_scaffold_P2A", "wild_type_lead")]
        for cseed in CONFIG["ligand_preparation"]["conformer_seeds"]:
            if cseed == 42:
                continue
            for c in subset:
                add("contender", c["id"], c["sequence"], "8FJC", gtfb, centre,
                    "wide_28", wide, 1, cseed)
            for n in null_a[:10]:
                add("null_a", n["id"], n["sequence"], "8FJC", gtfb, centre,
                    "wide_28", wide, 1, cseed)
            for n in null_b[:10]:
                add("null_b", n["id"], n["sequence"], "8FJC", gtfb, centre,
                    "wide_28", wide, 1, cseed)
    if 5 in tiers:
        for s in d["seeds_for_null_sets"]:
            for rid, rpath in off.items():
                if not Path(rpath).exists():
                    continue
                for c in contenders:
                    add("off_target", c["id"], c["sequence"], rid, rpath,
                        off_centres[rid], "wide_28", wide, s, 42)
    return jobs


def job_key(j: dict) -> str:
    return "|".join([j["kind"], j["id"], j["receptor_id"], j["box_id"],
                     f"d{j['dock_seed']}", f"c{j['conformer_seed']}"])


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tiers", default="1,2,3,4,5",
                    help="comma-separated tiers to run, in order. Default all.")
    ap.add_argument("--cpu", type=int, default=CONFIG["docking"]["cpu"])
    ap.add_argument("--dry-run", action="store_true",
                    help="list the work and stop without docking")
    args = ap.parse_args()
    tiers = [int(t) for t in args.tiers.split(",") if t.strip()]

    print("Preflight")
    pre = preflight()
    print(f"  vina: {pre['env']['vina_version_text']}")
    print(f"  receptor: {pre['env']['receptor_title']}")
    print(f"  receptor sha256: {pre['env']['receptor_sha256'][:16]}")

    jobs = build_jobs(pre["boxes"], tiers)
    done = set()
    if RESULTS.exists():
        for line in RESULTS.read_text().splitlines():
            if line.strip():
                try:
                    done.add(json.loads(line)["key"])
                except (json.JSONDecodeError, KeyError):
                    continue
    todo = [j for j in jobs if job_key(j) not in done]

    print(f"\n{len(jobs)} docking runs in tiers {tiers}; "
          f"{len(done)} already recorded; {len(todo)} to run.")
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
    }, indent=2) + "\n")

    t_start = time.time()
    completed = failed = 0
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
                    rec = {"key": key, **j, "ok": False, "error": f"ligand build failed: {exc}"}
                    fh.write(json.dumps(rec) + "\n")
                    fh.flush()
                    failed += 1
                    print(f"[{i}/{len(todo)}] {key} LIGAND FAILED: {exc}")
                    continue

            out = dock_one(lig, j["receptor"], j["centre"], j["size"], j["dock_seed"],
                           tag=key.replace("|", "__"), cpu=args.cpu,
                           exhaustiveness=CONFIG["docking"]["exhaustiveness"])
            rec = {"key": key, **j, "ligand_meta": lmeta,
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
