#!/usr/bin/env python3
"""Recompute the Angle 1 co-crystal redocking control from committed coordinates.

The repository reported a 4.51 Angstrom "Hungarian RMSD" for redocking the native
acarbose ligand of PDB 3AIC, but no code in the repository produced that number:
it was a literal in scripts/run_exhaustive_disproof.py, and the helper named
compute_heavy_atom_rmsd is a naive ordered RMSD that assumes the two files list
their atoms in the same order, which PDB and PDBQT do not guarantee.

This script computes the control from the files actually committed:

  reference : demo_data/receptors/3AIC_native_ligand.pdb   (crystal, chain I)
  redocked  : data/poses/3AIC_native_acarbose_docked_split_1.pdbqt  (Vina mode 1)

It reports three numbers, all heavy atom only:

  ordered RMSD    atoms paired in file order, which is what the old helper did
  assignment RMSD atoms paired optimally within each element, the standard
                  symmetry-corrected measure for a molecule with equivalent atoms
  centroid shift  distance between the two centres of geometry

Usage:  python3 scripts/compute_redocking_rmsd.py [--json]

No third-party dependency beyond numpy; the assignment solver is included below
so the result can be reproduced on a bare checkout.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import List, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
REFERENCE = ROOT / "demo_data" / "receptors" / "3AIC_native_ligand.pdb"
REDOCKED = ROOT / "data" / "poses" / "3AIC_native_acarbose_docked_split_1.pdbqt"


def _element(line: str) -> str:
    """Element symbol for a PDB or PDBQT atom record.

    PDBQT carries an AutoDock type in columns 78 to 79 (A for aromatic carbon,
    OA for an acceptor oxygen, and so on), so the leading letter is taken and
    aromatic carbon is folded back onto carbon.
    """
    tail = line[76:79].strip()
    if not tail:
        tail = line[12:16].strip()
    sym = tail[0].upper()
    return "C" if sym == "A" else sym


def read_atoms(path: Path) -> Tuple[np.ndarray, List[str]]:
    coords: List[List[float]] = []
    elements: List[str] = []
    for line in path.read_text().splitlines():
        if line.startswith("ENDMDL"):
            break
        if not line.startswith(("ATOM", "HETATM")):
            continue
        el = _element(line)
        if el == "H":
            continue
        # Meeko glue pseudo-atoms (type G, G0, G1) close an opened macrocycle
        # and are not atoms. Acarbose is not a macrocycle, so neither file in
        # this comparison carries one and the reported figure is unchanged by
        # this line; it is here so the routine stays correct if it is ever
        # pointed at a macrocyclic ligand.
        if line[76:79].strip().upper() in ("G", "G0", "G1", "G2", "G3"):
            continue
        coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        elements.append(el)
    if not coords:
        raise SystemExit(f"no heavy atoms parsed from {path}")
    return np.asarray(coords, dtype=float), elements


def linear_sum_assignment(cost: np.ndarray) -> List[Tuple[int, int]]:
    """Minimum-cost assignment, Jonker-Volgenant with potentials. Rows <= columns."""
    n, m = cost.shape
    if n > m:
        raise ValueError("rows must not exceed columns")
    INF = float("inf")
    u = [0.0] * (n + 1)
    v = [0.0] * (m + 1)
    p = [0] * (m + 1)
    way = [0] * (m + 1)
    for i in range(1, n + 1):
        p[0] = i
        j0 = 0
        minv = [INF] * (m + 1)
        used = [False] * (m + 1)
        while True:
            used[j0] = True
            i0 = p[j0]
            delta = INF
            j1 = -1
            for j in range(1, m + 1):
                if used[j]:
                    continue
                cur = cost[i0 - 1][j - 1] - u[i0] - v[j]
                if cur < minv[j]:
                    minv[j] = cur
                    way[j] = j0
                if minv[j] < delta:
                    delta = minv[j]
                    j1 = j
            for j in range(m + 1):
                if used[j]:
                    u[p[j]] += delta
                    v[j] -= delta
                else:
                    minv[j] -= delta
            j0 = j1
            if p[j0] == 0:
                break
        while True:
            j1 = way[j0]
            p[j0] = p[j1]
            j0 = j1
            if j0 == 0:
                break
    out = [0] * (n + 1)
    for j in range(1, m + 1):
        if p[j]:
            out[p[j]] = j
    return [(i - 1, out[i] - 1) for i in range(1, n + 1)]


def ordered_rmsd(a: np.ndarray, b: np.ndarray) -> float:
    n = min(len(a), len(b))
    d = a[:n] - b[:n]
    return float(np.sqrt((d * d).sum() / n))


def assignment_rmsd(a: np.ndarray, ea: List[str], b: np.ndarray, eb: List[str]) -> Tuple[float, int]:
    """Optimal pairing within each element, which is the symmetry-corrected RMSD."""
    total = 0.0
    paired = 0
    for el in sorted(set(ea) | set(eb)):
        ia = [i for i, e in enumerate(ea) if e == el]
        ib = [i for i, e in enumerate(eb) if e == el]
        if not ia or not ib:
            continue
        flip = len(ia) > len(ib)
        rows, cols = (ib, ia) if flip else (ia, ib)
        cost = np.zeros((len(rows), len(cols)))
        for r, i in enumerate(rows):
            pa = b[i] if flip else a[i]
            for c, j in enumerate(cols):
                pb = a[j] if flip else b[j]
                cost[r, c] = float(((pa - pb) ** 2).sum())
        for r, c in linear_sum_assignment(cost):
            total += cost[r, c]
            paired += 1
    return float(math.sqrt(total / paired)), paired


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", action="store_true", help="emit a JSON object")
    args = ap.parse_args()

    ref, eref = read_atoms(REFERENCE)
    dock, edock = read_atoms(REDOCKED)
    a_rmsd, paired = assignment_rmsd(ref, eref, dock, edock)
    o_rmsd = ordered_rmsd(ref, dock)
    centroid = float(np.linalg.norm(ref.mean(axis=0) - dock.mean(axis=0)))

    result = {
        "reference_file": str(REFERENCE.relative_to(ROOT)),
        "redocked_file": str(REDOCKED.relative_to(ROOT)),
        "reference_heavy_atoms": len(ref),
        "redocked_heavy_atoms": len(dock),
        "atoms_paired": paired,
        "assignment_rmsd_angstrom": round(a_rmsd, 3),
        "ordered_rmsd_angstrom": round(o_rmsd, 3),
        "centroid_displacement_angstrom": round(centroid, 3),
        "passes_2A_redocking_standard": bool(a_rmsd <= 2.0),
    }
    if args.json:
        print(json.dumps(result, indent=2))
        return 0
    print(f"reference        {result['reference_file']}  ({len(ref)} heavy atoms)")
    print(f"redocked         {result['redocked_file']}  ({len(dock)} heavy atoms)")
    print(f"atoms paired     {paired}")
    print(f"assignment RMSD  {a_rmsd:.3f} A   <- symmetry corrected, report this one")
    print(f"ordered RMSD     {o_rmsd:.3f} A   <- file-order pairing, not meaningful here")
    print(f"centroid shift   {centroid:.3f} A")
    print(f"2.0 A redocking standard: {'PASS' if a_rmsd <= 2.0 else 'FAIL'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
