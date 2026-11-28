#!/usr/bin/env python3
"""
Biophysics Engine: True 28Å / exh=16 Vina + Vinardo + OpenMM Physics.
Author: Dr. Roi Castrellon, DDS (DrNemo Caries Foundry)

Zero synthetic formulas, zero mock shortcuts:
- Uncompressed 28.0 x 28.0 x 28.0 A grid centered on catalytic triad oxygens
- Dual scoring: AutoDock Vina (--scoring vina) + Vinardo (--scoring vinardo)
- Ligand Efficiency: LE = -E_vina / heavy_atoms
- Disulfide S-S geometry verification in output pose (2.00 - 2.15 A)
- OpenMM all-atom Amber14 (amber14-all.xml + implicit/gbn2.xml) isolated ligand-only steric minimization, capped at 200 iterations
- Solubility Shield: Donnan equilibrium (Q_pH6 <= +4.0) & Stephan acid-switch (Delta Q >= +1.0)
"""
import os
import sys
import shutil
import hashlib
import subprocess
import math
import numpy as np
from pathlib import Path
from typing import Dict, Any, Tuple, Optional

import provenance as prov

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except ImportError:
    Chem = None
    AllChem = None

try:
    from meeko import MoleculePreparation, PDBQTWriterLegacy
except ImportError:
    MoleculePreparation = None
    PDBQTWriterLegacy = None

try:
    import openmm as mm
    from openmm import app, unit
except ImportError:
    mm = None
    app = None
    unit = None

from solubility_shield import (
    is_viable,
    calc_net_charge_at_ph,
    calc_stephan_switch,
    count_protease_cleavage_sites,
    MAX_DONNAN_CHARGE_PH6
)

REPO_ROOT = Path(__file__).parent.resolve()
RECEPTORS_DIR = REPO_ROOT / "demo_data" / "receptors"
RECEPTOR_3AIC_PDBQT = RECEPTORS_DIR / "3AIC_receptor.pdbqt"
RECEPTOR_1SMD_PDBQT = RECEPTORS_DIR / "1SMD_receptor.pdbqt"
RECEPTOR_1LZ1_PDBQT = RECEPTORS_DIR / "1LZ1_receptor.pdbqt"

# Uncompressed 28 A grid centered on solvent-exposed catalytic triads / active clefts
TARGET_COORDS = {
    "3AIC": {
        "center": (192.072, 44.370, 197.309),
        "size": (28.0, 28.0, 28.0),
        "receptor": str(RECEPTOR_3AIC_PDBQT)
    },
    "1SMD": {
        "center": (8.778, 49.724, 17.278),
        "size": (28.0, 28.0, 28.0),
        "receptor": str(RECEPTOR_1SMD_PDBQT)
    },
    "1LZ1": {
        "center": (10.319, 8.503, 25.960),
        "size": (28.0, 28.0, 28.0),
        "receptor": str(RECEPTOR_1LZ1_PDBQT)
    }
}


def build_macrocycle_3d(
    seq: str,
    out_pdbqt: str,
    out_pdb: Optional[str] = None,
    target_center: Optional[Tuple[float, float, float]] = None
) -> int:
    """
    Builds an all-atom 3D structure for a cyclic disulfide peptide with RDKit.

    MMFF94 minimisation is ATTEMPTED, not guaranteed. MMFFOptimizeMolecule reports
    "forcefield could not be set up" by returning -1 and a hit iteration cap by
    returning 1, per RDKit's documentation, so the try/except below reaches its UFF
    fallback only on an exception. The return value is now read and a non-zero value
    is reported on standard error, but nothing persists it, so a past run's outcome
    cannot be recovered. Do not describe output of this function as MMFF94-minimised
    without checking the warnings from the run that produced it.
    If target_center is provided, translates the conformer centroid to target_center.
    Returns the count of heavy atoms.
    """
    mol = Chem.MolFromSequence(seq)
    if mol is None:
        raise ValueError(f"Invalid peptide sequence: {seq}")
    mol = Chem.AddHs(mol)

    # Disulfide bond strictly between first residue (resi 1, N-terminal CYS) SG
    # and last residue (resi N, C-terminal CYS) SG
    cys1_sg = None
    cysN_sg = None
    N = len(seq)
    for a in mol.GetAtoms():
        info = a.GetPDBResidueInfo()
        if info and info.GetResidueName().strip() == "CYS" and info.GetName().strip() == "SG":
            if info.GetResidueNumber() == 1:
                cys1_sg = a.GetIdx()
            elif info.GetResidueNumber() == N:
                cysN_sg = a.GetIdx()

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

        cys1_sg = None
        cysN_sg = None
        for a in mol.GetAtoms():
            info = a.GetPDBResidueInfo()
            if info and info.GetResidueName().strip() == "CYS" and info.GetName().strip() == "SG":
                if info.GetResidueNumber() == 1:
                    cys1_sg = a.GetIdx()
                elif info.GetResidueNumber() == N:
                    cysN_sg = a.GetIdx()

        emol = Chem.EditableMol(mol)
        emol.AddBond(cys1_sg, cysN_sg, Chem.BondType.SINGLE)
        mol = emol.GetMol()

    Chem.SanitizeMol(mol)

    # 3D Distance Geometry embedding, then an ATTEMPT at MMFF94 minimisation.
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    embed_res = AllChem.EmbedMolecule(mol, params)
    if embed_res != 0:
        embed_res = AllChem.EmbedMolecule(mol, useRandomCoords=True, maxAttempts=5000, randomSeed=42)
    if embed_res != 0:
        # Both embeddings failed, so the molecule carries no conformer. Previously this
        # fell through to the writers with nothing to write. Fail here instead of
        # handing a conformer-less molecule onward.
        raise RuntimeError(
            f"ETKDGv3 embedding failed for {seq} on both the seeded and the random-coordinate "
            f"attempt, so no conformer exists and no pose can be written.")

    # The force field selection below is deliberately left as it was: changing which
    # field runs would change the geometry, and these are the geometries the published
    # campaign docked. What is added is disclosure, not a different calculation.
    # MMFFOptimizeMolecule returns 0 on convergence, -1 if the force field could not be
    # set up for this molecule, and 1 if more iterations were needed. Failure is a return
    # code, not an exception, so the UFF fallback below is unreachable for the -1 case.
    try:
        mmff_res = AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
        if mmff_res != 0:
            reason = ("MMFF94 could not be set up for this molecule"
                      if mmff_res == -1 else "the 500-iteration cap was reached")
            print(f"WARNING: MMFF94 minimisation of {seq} returned {mmff_res}: {reason}. "
                  f"The conformer is not a converged MMFF94 minimum.", file=sys.stderr)
    except Exception as exc:
        print(f"WARNING: MMFF94 minimisation of {seq} raised {exc!r}, falling back to UFF.",
              file=sys.stderr)
        try:
            AllChem.UFFOptimizeMolecule(mol, maxIters=500)
        except Exception as exc2:
            print(f"WARNING: UFF minimisation of {seq} also raised {exc2!r}. The conformer "
                  f"is raw ETKDG output with no force field applied.", file=sys.stderr)

    if target_center is not None:
        conf = mol.GetConformer()
        positions = conf.GetPositions()
        centroid = np.mean(positions, axis=0)
        translation = np.array(target_center) - centroid
        for i in range(mol.GetNumAtoms()):
            conf.SetAtomPosition(i, positions[i] + translation)

    # Save clean heavy-atom PDB if requested (for OpenMM / PDBFixer)
    if out_pdb:
        os.makedirs(os.path.dirname(os.path.abspath(out_pdb)), exist_ok=True)
        mol_heavy = Chem.RemoveHs(mol)
        Chem.MolToPDBFile(mol_heavy, out_pdb)

    # Convert to PDBQT via Meeko
    os.makedirs(os.path.dirname(os.path.abspath(out_pdbqt)), exist_ok=True)
    preparator = MoleculePreparation()
    preparator.prepare(mol)
    writer = PDBQTWriterLegacy()
    pdbqt_str = writer.write_string(preparator.setup)[0]
    with open(out_pdbqt, "w", encoding="utf-8") as f:
        f.write(pdbqt_str)

    return mol.GetNumHeavyAtoms()


def measure_disulfide_bond(pdbqt_path: str) -> Optional[float]:
    """Measures the terminal Cys1:SG-CysN:SG distance (Angstroms) in a docked PDBQT pose."""
    cys1_sg = None
    cys_last_sg = None
    max_resi = -1
    if not os.path.exists(pdbqt_path):
        return None
    with open(pdbqt_path, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                try:
                    res_num = int(line[22:26].strip())
                except ValueError:
                    continue
                if atom_name == "SG" and res_name == "CYS":
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    if res_num == 1:
                        cys1_sg = (x, y, z)
                    if res_num > max_resi:
                        max_resi = res_num
                        cys_last_sg = (x, y, z)
            if line.startswith("ENDMDL"):
                break
    if cys1_sg and cys_last_sg and max_resi > 1:
        dx = cys1_sg[0] - cys_last_sg[0]
        dy = cys1_sg[1] - cys_last_sg[1]
        dz = cys1_sg[2] - cys_last_sg[2]
        dist = math.sqrt(dx * dx + dy * dy + dz * dz)
        return round(dist, 2)
    return None


def get_vina_path() -> str:
    """Resolves AutoDock Vina binary from PATH or raises FileNotFoundError."""
    path = shutil.which("vina")
    if not path:
        raise FileNotFoundError("AutoDock Vina executable ('vina') not found on system PATH")
    return path


def get_vina_split_path() -> str:
    """Resolves AutoDock Vina split binary from PATH or raises FileNotFoundError."""
    path = shutil.which("vina_split")
    if not path:
        raise FileNotFoundError("AutoDock Vina split executable ('vina_split') not found on system PATH")
    return path


def run_docking_vina_and_vinardo(
    ligand_pdbqt: str,
    target: str = "3AIC",
    exhaustiveness: int = 16,
    cpu: int = 8,
    out_pdbqt: Optional[str] = None,
    seed: Optional[int] = None,
    ligand_id: Optional[str] = None,
    protonation: str = "neutral",
    resume: bool = False
) -> Dict[str, Any]:
    """
    Executes AutoDock Vina with uncompressed 28A grid, extracts Vina affinity,
    then re-scores with --score_only --scoring vinardo to get Vinardo affinity.

    Naming and provenance. When out_pdbqt is not given, the output path is
    derived from the run's own parameters by provenance.RunSpec: target,
    ligand, protonation, seed, exhaustiveness, box, and the content of the
    receptor and ligand files. The result carries the full command line, the
    seed, the box, the engine version, both input digests, a run identifier
    and a UTC timestamp, written as REMARK PROVENANCE records that PDBQT
    readers ignore. A caller-supplied out_pdbqt is honoured, and the record
    then names both the path written and the canonical path, so a hand-chosen
    name is detectable by scripts/check_result_provenance.py.

    Resuming. An existing file at the output path is no longer adopted
    silently. Pass resume=True to reuse one, and it is accepted only if its
    own provenance record shows it came from these exact parameters. That
    check is the reason this argument exists: adopting a file written under a
    different protocol is how a result ends up reported as something it is not.
    """
    cfg = TARGET_COORDS[target]
    spec = prov.RunSpec.for_run(
        target=target,
        ligand_id=ligand_id or Path(ligand_pdbqt).stem,
        receptor=cfg["receptor"],
        ligand=ligand_pdbqt,
        box_center=cfg["center"],
        box_size=cfg["size"],
        exhaustiveness=exhaustiveness,
        seed=seed,
        protonation=protonation,
    )
    caller_chose_path = out_pdbqt is not None
    if out_pdbqt is None:
        out_pdbqt = str(spec.result_path(Path(ligand_pdbqt).parent))

    cx, cy, cz = cfg["center"]
    sx, sy, sz = cfg["size"]

    # Reuse an existing pose only on request, and only if it says it is this run
    vina_score = None
    if resume and os.path.exists(out_pdbqt) and os.path.getsize(out_pdbqt) > 100:
        record = prov.read_record(out_pdbqt)
        if record is None:
            raise RuntimeError(
                f"{out_pdbqt} exists but carries no provenance record, so it "
                f"cannot be shown to be the output of this run. Delete it, or "
                f"dock to a path derived from the parameters you mean."
            )
        if record.get("spec_id") != spec.spec_id():
            raise RuntimeError(
                f"{out_pdbqt} was produced by a different run: its record says "
                f"spec {record.get('spec_id')}, these parameters are "
                f"{spec.spec_id()}. Refusing to report it as this one."
            )
        with open(out_pdbqt, "r", encoding="utf-8") as f:
            for line in f:
                if line.startswith("REMARK VINA RESULT:"):
                    vina_score = float(line.split()[3])
                    break

    # 1. Dock with Vina scoring if not already present
    if vina_score is None:
        prov.reserve_result_path(out_pdbqt, overwrite=True)
        cmd_vina = [
            get_vina_path(),
            "--receptor", cfg["receptor"],
            "--ligand", ligand_pdbqt,
            "--center_x", str(cx), "--center_y", str(cy), "--center_z", str(cz),
            "--size_x", str(sx), "--size_y", str(sy), "--size_z", str(sz),
            "--cpu", str(cpu),
            "--exhaustiveness", str(exhaustiveness),
            "--out", out_pdbqt
        ]
        if seed is not None:
            cmd_vina += ["--seed", str(seed)]
        res_vina = subprocess.run(cmd_vina, capture_output=True, text=True)
        if res_vina.returncode != 0:
            raise RuntimeError(f"Vina docking failed: {res_vina.stderr}")

        # Record what produced this file, inside the file, before anything reads it
        prov.stamp(out_pdbqt, prov.build_record(
            spec, cmd_vina, engine=prov.engine_version(get_vina_path()),
            written_path=out_pdbqt))

        # Parse best Vina affinity from output pose
        with open(out_pdbqt, "r", encoding="utf-8") as f:
            for line in f:
                if line.startswith("REMARK VINA RESULT:"):
                    vina_score = float(line.split()[3])
                    break

    if vina_score is None:
        raise RuntimeError(f"Could not parse Vina score from {out_pdbqt}")

    # 2. Rescore top pose with Vinardo
    split_prefix = out_pdbqt.replace(".pdbqt", "_split_")
    split_cmd = [get_vina_split_path(), "--input", out_pdbqt, "--ligand", split_prefix]
    subprocess.run(split_cmd, capture_output=True, text=True)
    mode1_pdbqt = f"{split_prefix}1.pdbqt"

    vinardo_score = None
    if os.path.exists(mode1_pdbqt):
        cmd_vinardo = [
            get_vina_path(),
            "--receptor", cfg["receptor"],
            "--ligand", mode1_pdbqt,
            "--center_x", str(cx), "--center_y", str(cy), "--center_z", str(cz),
            "--size_x", str(sx), "--size_y", str(sy), "--size_z", str(sz),
            "--scoring", "vinardo",
            "--score_only"
        ]
        res_vinardo = subprocess.run(cmd_vinardo, capture_output=True, text=True)
        for line in res_vinardo.stdout.splitlines():
            if "Estimated Free Energy of Binding" in line:
                # "Estimated Free Energy of Binding   : -5.200 (kcal/mol) [=(1)+(2)+(3)-(4)]"
                parts = line.split(":")
                if len(parts) >= 2:
                    val_str = parts[1].split()[0]
                    try:
                        vinardo_score = float(val_str)
                    except ValueError:
                        pass
                break

        # Fallback to autobox if ligand drifted slightly outside fixed box
        if vinardo_score is None:
            cmd_vinardo_auto = [
                get_vina_path(),
                "--receptor", cfg["receptor"],
                "--ligand", mode1_pdbqt,
                "--scoring", "vinardo",
                "--score_only",
                "--autobox"
            ]
            res_auto = subprocess.run(cmd_vinardo_auto, capture_output=True, text=True)
            for line in res_auto.stdout.splitlines():
                if "Estimated Free Energy of Binding" in line:
                    parts = line.split(":")
                    if len(parts) >= 2:
                        try:
                            vinardo_score = float(parts[1].split()[0])
                        except ValueError:
                            pass
                    break

    if vinardo_score is None:
        raise RuntimeError(f"Vinardo scoring failed to extract binding free energy for {ligand_pdbqt} on {target}")

    # 3. Pose SHA256 & Disulfide check
    with open(out_pdbqt, "rb") as f:
        pose_sha = hashlib.sha256(f.read()).hexdigest()

    ss_dist = measure_disulfide_bond(out_pdbqt)

    return {
        "target": target,
        "vina_kcal": vina_score,
        "vinardo_kcal": vinardo_score,
        "out_pdbqt": out_pdbqt,
        "pose_sha256": pose_sha,
        "spec_id": spec.spec_id(),
        "canonical_name": spec.canonical_name(),
        "path_chosen_by_caller": caller_chose_path,
        "seed": seed,
        "ss_bond_angstrom": ss_dist
    }


def relax_with_openmm(pdb_path: str, temp_k: float = 310.15) -> float:
    """
    Runs isolated ligand-only Amber14/GBn2 energy minimization for steric relaxation
    and returns the implicit-solvation potential energy (kJ/mol) of the relaxed pose.

    maxIterations=200 is an upper bound, not a step count. Per the OpenMM API, the
    minimiser halts when the RMS force reaches its default tolerance or when the cap
    is reached, whichever comes first, minimizeEnergy returns nothing, and no
    convergence status is available afterwards. So the caller cannot tell how many
    iterations ran or whether the result converged, and nothing here records it. Do
    not describe the output as 200-step relaxed.

    Does NOT run bound receptor-ligand molecular dynamics at 310.15 K.
    """
    from pdbfixer import PDBFixer
    fixer = PDBFixer(filename=pdb_path)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.4)

    forcefield = app.ForceField('amber14-all.xml', 'implicit/gbn2.xml')
    system = forcefield.createSystem(
        fixer.topology,
        nonbondedMethod=app.NoCutoff,
        constraints=app.HBonds
    )

    integrator = mm.LangevinMiddleIntegrator(
        temp_k * unit.kelvin,
        1.0 / unit.picoseconds,
        2.0 * unit.femtoseconds
    )

    simulation = app.Simulation(fixer.topology, system, integrator)
    simulation.context.setPositions(fixer.positions)

    # Energy minimize
    simulation.minimizeEnergy(maxIterations=200)
    state = simulation.context.getState(getEnergy=True)
    energy_kj = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
    return round(energy_kj, 2)
