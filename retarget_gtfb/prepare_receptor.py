#!/usr/bin/env python3
"""
Turn the validated 8FJC coordinate file into a docking receptor, extract the
native acarbose for the redocking calibration control, and derive the grid boxes.

Preparation deliberately mirrors scripts/prepare_receptors.py, which produced the
3AIC receptor the old campaign used: chain A standard residues only, then
obabel -xr -p 7.4 --partialcharge gasteiger. Same flags, so the only thing that
differs between the old numbers and the new ones is the protein itself.

Writes only into retarget_gtfb/receptor and retarget_gtfb/out.
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent.resolve()
RECEPTOR_DIR = HERE / "receptor"
OUT = HERE / "out"
CONFIG = json.loads((HERE / "config.json").read_text())

TRIAD_OXYGENS = {
    (451, "ASP"): ("OD1", "OD2"),
    (489, "GLU"): ("OE1", "OE2"),
    (562, "ASP"): ("OD1", "OD2"),
}
LIGAND_RESNAMES = {"GLC", "AC1", "ACR"}
TIGHT_MARGIN_A = 8.0
CA_KEEP_RADIUS_A = 15.0

STANDARD_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "MSE",
}


def die(msg: str) -> None:
    print(f"STOP: {msg}", file=sys.stderr)
    sys.exit(1)


def xyz(line: str) -> np.ndarray:
    return np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])


def main() -> int:
    stamp = OUT / "receptor_validation.json"
    if not stamp.exists():
        die("out/receptor_validation.json is missing. Run validate_receptor.py first; "
            "nothing here will touch an unvalidated coordinate file.")
    validation = json.loads(stamp.read_text())
    if not validation.get("validated"):
        die("the validation stamp does not record a pass.")

    src = Path(validation["path"])
    lines = src.read_text(errors="replace").splitlines()
    RECEPTOR_DIR.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)

    # --- chain A protein -----------------------------------------------------
    receptor_lines, het_lines = [], []
    for raw in lines:
        rec, chain = raw[:6], raw[21:22]
        resname = raw[17:20].strip().upper()
        altloc = raw[16:17]
        if rec == "ATOM  " and chain == "A" and resname in STANDARD_AA:
            if altloc in (" ", "A"):
                receptor_lines.append(raw)
        elif rec == "HETATM":
            het_lines.append(raw)

    if len(receptor_lines) < 5000:
        die(f"only {len(receptor_lines)} chain A protein atoms survived filtering; expected "
            f"roughly 6500. Refusing to build a receptor from a partial model.")

    receptor_pdb = RECEPTOR_DIR / "8FJC_receptor.pdb"
    receptor_pdb.write_text("\n".join(receptor_lines) + "\nTER\nEND\n")

    # --- catalytic triad centroid -------------------------------------------
    triad_pts, triad_found = [], {}
    for raw in receptor_lines:
        resnum = int(raw[22:26])
        resname = raw[17:20].strip().upper()
        atom = raw[12:16].strip()
        wanted = TRIAD_OXYGENS.get((resnum, resname))
        if wanted and atom in wanted:
            triad_pts.append(xyz(raw))
            triad_found.setdefault(f"{resname}{resnum}", []).append(atom)
    if len(triad_pts) < 6:
        die(f"found only {len(triad_pts)} of the 6 catalytic carboxylate oxygens "
            f"({triad_found}). The site cannot be located reliably.")
    triad_centre = np.mean(triad_pts, axis=0)

    # --- native acarbose, discovered rather than assumed --------------------
    # The ligand chain id is not assumed: every acarbose-like residue in the file is
    # collected, grouped by chain, and the group nearest the chain A triad is the one
    # that belongs to chain A's active site.
    by_chain: dict[str, list[str]] = {}
    for raw in het_lines:
        if raw[17:20].strip().upper() in LIGAND_RESNAMES and raw[16:17] in (" ", "A"):
            by_chain.setdefault(raw[21:22], []).append(raw)
    if not by_chain:
        die("no acarbose atoms (GLC/AC1/ACR) found anywhere in the file.")

    best_chain, best_dist = None, None
    for chain, group in by_chain.items():
        centre = np.mean([xyz(l) for l in group], axis=0)
        d = float(np.linalg.norm(centre - triad_centre))
        if best_dist is None or d < best_dist:
            best_chain, best_dist = chain, d

    ligand_lines = by_chain[best_chain]
    ligand_pts = np.array([xyz(l) for l in ligand_lines])
    ligand_centre = np.mean(ligand_pts, axis=0)

    if best_dist > 12.0:
        die(f"the nearest acarbose copy (chain {best_chain}) sits {best_dist:.2f} A from the "
            f"chain A catalytic triad. That is too far to be the bound inhibitor of this "
            f"copy of the protein, so the site cannot be assigned.")

    native_ligand_pdb = RECEPTOR_DIR / "8FJC_native_ligand.pdb"
    native_ligand_pdb.write_text("\n".join(ligand_lines) + "\nEND\n")

    # --- box definitions ----------------------------------------------------
    centre_separation = float(np.linalg.norm(ligand_centre - triad_centre))
    if centre_separation > 4.0:
        die(f"the acarbose centroid and the catalytic triad centroid are {centre_separation:.2f} A "
            f"apart, over the 4.0 A limit set in config.json. Two independent definitions of the "
            f"site disagreeing by that much means the site is not what this package assumes.")

    extent = ligand_pts.max(axis=0) - ligand_pts.min(axis=0)
    tight = np.ceil((extent + TIGHT_MARGIN_A) * 2) / 2
    wide = np.array(CONFIG["docking"]["boxes"]["wide_28"]["size"], dtype=float)

    # --- heteroatoms near the site, reported not hidden ---------------------
    near_site = []
    for raw in het_lines:
        resname = raw[17:20].strip().upper()
        if resname in LIGAND_RESNAMES or resname == "HOH":
            continue
        d = float(np.linalg.norm(xyz(raw) - ligand_centre))
        if d <= 6.0:
            near_site.append({"resname": resname, "chain": raw[21:22],
                              "resnum": raw[22:26].strip(), "atom": raw[12:16].strip(),
                              "distance_to_box_centre_A": round(d, 2)})

    # --- the conserved structural Ca2+, kept as an explicit condition -------
    # Every heteroatom is stripped from the main receptor, which is what the original
    # campaign did and what keeps the two comparable. But the GtfB catalytic domain
    # carries a conserved structural calcium next to the triad, and one of its
    # coordinating residues (Asn455) also hydrogen bonds acarbose, so it is not the same
    # kind of heteroatom as a crystallisation sulfate. Removing it leaves open volume
    # that Vina's steric terms can fill; keeping it gives Vina a metal that it types as a
    # hydrogen-bond donor. Neither is obviously right, so both receptors are built and
    # tier 6 docks against both. Nothing here assumes which answer is correct.
    ca_lines, ca_report = [], []
    for raw in het_lines:
        if raw[17:20].strip().upper() != "CA":
            continue
        d = float(np.linalg.norm(xyz(raw) - ligand_centre))
        dt = float(np.linalg.norm(xyz(raw) - triad_centre))
        ca_report.append({"chain": raw[21:22], "resnum": raw[22:26].strip(),
                          "distance_to_box_centre_A": round(d, 2),
                          "distance_to_triad_centroid_A": round(dt, 2),
                          "kept": d <= CA_KEEP_RADIUS_A})
        if d <= CA_KEEP_RADIUS_A:
            ca_lines.append(raw)

    receptor_ca_pdb = RECEPTOR_DIR / "8FJC_receptor_with_ca.pdb"
    if ca_lines:
        receptor_ca_pdb.write_text("\n".join(receptor_lines + ca_lines) + "\nTER\nEND\n")

    # --- obabel conversion --------------------------------------------------
    receptor_pdbqt = RECEPTOR_DIR / "8FJC_receptor.pdbqt"
    ligand_pdbqt = RECEPTOR_DIR / "8FJC_native_ligand.pdbqt"
    for src_pdb, dst, extra in [(receptor_pdb, receptor_pdbqt, ["-xr"]),
                                (native_ligand_pdb, ligand_pdbqt, [])]:
        cmd = ["obabel", "-ipdb", str(src_pdb), "-opdbqt", "-O", str(dst),
               *extra, "-p", "7.4", "--partialcharge", "gasteiger"]
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode != 0 or not dst.exists() or dst.stat().st_size < 1000:
            die(f"obabel failed on {src_pdb.name}:\n{res.stderr.strip()}")
        print(f"  wrote {dst.name} ({dst.stat().st_size} bytes)")

    # --- convert and smoke-test the Ca-retained receptor --------------------
    # Done here, in seconds, rather than discovered as 87 failed docking runs overnight.
    # Two things can go wrong and both are checked: Open Babel can silently drop the
    # HETATM, and a given Vina build can refuse the resulting atom type.
    receptor_ca_pdbqt = RECEPTOR_DIR / "8FJC_receptor_with_ca.pdbqt"
    ca_status: dict = {"ca_ions_in_file": len(ca_report), "ca_ions_kept": len(ca_lines),
                       "ions": ca_report, "keep_radius_A": CA_KEEP_RADIUS_A,
                       "usable": False, "receptor_pdbqt": None}
    if ca_lines:
        cmd = ["obabel", "-ipdb", str(receptor_ca_pdb), "-opdbqt", "-O", str(receptor_ca_pdbqt),
               "-xr", "-p", "7.4", "--partialcharge", "gasteiger"]
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode != 0 or not receptor_ca_pdbqt.exists():
            ca_status["problem"] = f"obabel failed: {res.stderr.strip()[:400]}"
        else:
            kept = [l for l in receptor_ca_pdbqt.read_text(errors="replace").splitlines()
                    if l[:6] in ("ATOM  ", "HETATM") and l.rstrip().split()[-1].upper() == "CA"]
            if not kept:
                ca_status["problem"] = (
                    "obabel produced a PDBQT with no Ca atom type in it, so the ion did not "
                    "survive conversion. The with-calcium condition cannot be run.")
            else:
                probe = subprocess.run(
                    ["vina", "--score_only", "--receptor", str(receptor_ca_pdbqt),
                     "--ligand", str(ligand_pdbqt)], capture_output=True, text=True)
                if probe.returncode != 0:
                    ca_status["problem"] = (
                        f"this Vina build refused the Ca-containing receptor: "
                        f"{probe.stderr.strip()[:400]}")
                else:
                    ca_status["usable"] = True
                    ca_status["receptor_pdbqt"] = str(receptor_ca_pdbqt)
                    ca_status["ca_atoms_in_pdbqt"] = len(kept)
                    ca_status["score_only_probe_ok"] = True
    else:
        ca_status["problem"] = (
            f"no Ca2+ ion within {CA_KEEP_RADIUS_A} A of the box centre. Every result will "
            f"come from the metal-free receptor, which is a stated limitation.")

    boxes = {
        "receptor_pdbqt": str(receptor_pdbqt),
        "calcium_condition": ca_status,
        "receptor_pdb": str(receptor_pdb),
        "native_ligand_pdbqt": str(ligand_pdbqt),
        "native_ligand_pdb": str(native_ligand_pdb),
        "native_ligand_chain_in_8fjc": best_chain,
        "native_ligand_atoms": len(ligand_lines),
        "centre": [round(float(v), 3) for v in ligand_centre],
        "centre_definition": "centroid of the heavy atoms of the bound acarbose nearest chain A",
        "triad_centre": [round(float(v), 3) for v in triad_centre],
        "triad_oxygens_used": triad_found,
        "centre_separation_A": round(centre_separation, 3),
        "acarbose_extent_A": [round(float(v), 2) for v in extent],
        "boxes": {
            "wide_28": {"size": [float(v) for v in wide]},
            "tight_auto": {"size": [float(v) for v in tight],
                           "margin_A": TIGHT_MARGIN_A},
        },
        "heteroatoms_within_6A_of_centre": near_site,
        "receptor_atoms": len(receptor_lines),
        "source_sha256": validation["sha256"],
    }
    (OUT / "boxes.json").write_text(json.dumps(boxes, indent=2) + "\n")

    print("\nReceptor prepared.")
    print(f"  chain A atoms kept:      {len(receptor_lines)}")
    print(f"  acarbose found in chain:  {best_chain} ({len(ligand_lines)} atoms)")
    print(f"  box centre (acarbose):    {boxes['centre']}")
    print(f"  catalytic triad centroid: {boxes['triad_centre']}")
    print(f"  the two centres differ by {centre_separation:.2f} A (limit 4.00)")
    print(f"  wide box:  {boxes['boxes']['wide_28']['size']}")
    print(f"  tight box: {boxes['boxes']['tight_auto']['size']} "
          f"(acarbose extent {boxes['acarbose_extent_A']} plus {TIGHT_MARGIN_A} A)")
    if near_site:
        print(f"  NOTE: {len(near_site)} crystallisation heteroatom(s) within 6 A of the box "
              f"centre, listed in out/boxes.json. 8FJC was grown from 2 M ammonium sulfate.")
    if ca_status["usable"]:
        d = min(i["distance_to_triad_centroid_A"] for i in ca_report if i["kept"])
        print(f"  Ca2+ condition: usable. {ca_status['ca_ions_kept']} ion(s) retained, "
              f"nearest {d} A from the catalytic triad centroid. Tier 6 will dock against "
              f"both receptors.")
    else:
        print(f"  Ca2+ condition: not available. {ca_status['problem']}")
    print(f"  wrote {OUT / 'boxes.json'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
