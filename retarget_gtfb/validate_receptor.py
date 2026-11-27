#!/usr/bin/env python3
"""
Validate the downloaded 8FJC coordinate file before anything is docked against it.

Why this file exists
--------------------
This package was prepared in an environment with no access to the RCSB or PDBe
download endpoints, so the 8FJC coordinates were never in hand when the package
was written. That has one consequence worth stating plainly: there is no
pre-computed SHA-256 of the correct file to compare against, because inventing
one would be worse than having none. A checksum you cannot verify is not
provenance, it is decoration.

So the check here is on content rather than on a hash. The file has to actually
be GtfB, actually be the tetragonal acarbose complex, actually contain the three
catalytic residues at the residue numbers the primary literature gives for GtfB,
and actually contain the bound inhibitor. Each of those is independently
checkable against Schormann et al. 2023 (DOI 10.1107/S2053230X23003199), and all
of them together are far harder to satisfy by accident than a hash is.

The SHA-256 of whatever file you supply is computed and recorded in
out/receptor_validation.json, so every later result is traceable to the exact
bytes that produced it. From the second run onwards the hash is compared with
the recorded one and any change stops the run.

Exit status is 0 only if every check passes. Any failure stops the pipeline.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

HERE = Path(__file__).parent.resolve()
OUT = HERE / "out"
CONFIG = json.loads((HERE / "config.json").read_text())

# From Schormann et al. 2023. Checked, not guessed.
EXPECTED_CATALYTIC = {451: "ASP", 489: "GLU", 562: "ASP", 491: "TRP"}
EXPECTED_CONTACTS = {449: "ARG", 455: "ASN", 511: "ASN", 514: "ARG",
                     561: "HIS", 883: "ASP", 934: "GLN"}
LIGAND_RESNAMES = {"GLC", "AC1", "ACR"}

# The conserved structural Ca2+ site, which sits adjacent to the catalytic triad.
# Schormann et al. 2023 give the GtfC coordinating residues as Glu431, Asp437, Asn481
# and Asp959 (from Ito et al. 2011), and state directly that GtfC Asp437 corresponds to
# GtfB Asp411. Every GtfC-to-GtfB pair in that paper's Table 4 differs by exactly 26
# (475/449, 477/451, 481/455, 515/489, 540/514, 587/561, 588/562, 909/883, 960/934), and
# 437 - 26 = 411 reproduces the one correspondence the paper states outright. The other
# three GtfB numbers below are therefore derived at that offset, not quoted, so they are
# reported rather than enforced: if one of them comes back as the wrong amino acid, the
# offset derivation is wrong and that is worth knowing before anything is docked.
CALCIUM_SITE_STATED = {411: "ASP"}
CALCIUM_SITE_DERIVED = {405: "GLU", 455: "ASN", 933: "ASP"}
METAL_RESNAMES = {"CA", "MG", "MN", "ZN", "NA", "K"}


class Failed(Exception):
    pass


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def parse(path: Path):
    """Minimal fixed-column PDB reader. No external dependencies on purpose."""
    protein: dict[tuple[str, int], str] = {}
    hetatms: dict[tuple[str, str, int], int] = {}
    metals: list[dict] = []
    header, title, remarks = "", [], []
    n_atom = 0
    for raw in path.read_text(errors="replace").splitlines():
        rec = raw[:6]
        if rec == "HEADER":
            header = raw[10:].strip()
        elif rec == "TITLE ":
            title.append(raw[10:].strip())
        elif rec.startswith("REMARK"):
            remarks.append(raw.strip())
        elif rec == "ATOM  ":
            n_atom += 1
            chain = raw[21]
            try:
                resnum = int(raw[22:26])
            except ValueError:
                continue
            protein[(chain, resnum)] = raw[17:20].strip().upper()
        elif rec == "HETATM":
            chain = raw[21]
            resname = raw[17:20].strip().upper()
            try:
                resnum = int(raw[22:26])
            except ValueError:
                continue
            key = (chain, resname, resnum)
            hetatms[key] = hetatms.get(key, 0) + 1
            if resname in METAL_RESNAMES:
                try:
                    pos = (float(raw[30:38]), float(raw[38:46]), float(raw[46:54]))
                except ValueError:
                    pos = None
                metals.append({"resname": resname, "chain": chain, "resnum": resnum,
                               "atom": raw[12:16].strip(), "xyz": pos})
    return {"header": header, "title": " ".join(title), "remarks": remarks,
            "protein": protein, "hetatms": hetatms, "n_atom": n_atom,
            "metals": metals}


def check(path: Path) -> dict:
    if not path.exists():
        raise Failed(f"{path} does not exist. Run fetch_receptor.sh first.")
    if path.stat().st_size < 200_000:
        raise Failed(
            f"{path} is only {path.stat().st_size} bytes. A two-chain 2.5 A structure "
            f"is around 1 MB. This is most likely a truncated download or an error page, "
            f"not a coordinate file. Delete it and fetch it again."
        )

    s = parse(path)
    results, problems = {}, []

    # 1. Is it 8FJC at all?
    if "8FJC" not in s["header"].upper():
        problems.append(f"HEADER does not mention 8FJC: {s['header']!r}")
    results["header"] = s["header"]

    # 2. Is it GtfB, and is it the tetragonal acarbose complex?
    title = s["title"].upper()
    results["title"] = s["title"]
    for token, why in [("GTFB", "the title must name GtfB"),
                       ("ACARBOSE", "the title must name the bound acarbose"),
                       ("P4322", "the title must name the tetragonal space group P4322")]:
        if token not in title.replace("(", "").replace(")", ""):
            problems.append(f"{why}; token {token} absent from TITLE")

    # A deposited GtfC structure must never be accepted here again.
    for wrong in ("GTFC", "GTF-SI", "3AIC"):
        if wrong in title:
            problems.append(f"TITLE mentions {wrong}. This is the error this whole "
                            f"exercise exists to correct. Refusing to continue.")

    # 3. Resolution, read from REMARK 2 rather than trusted from the config.
    res = None
    for r in s["remarks"]:
        if r.startswith("REMARK   2 RESOLUTION."):
            for tok in r.split():
                try:
                    res = float(tok)
                    break
                except ValueError:
                    continue
            break
    results["resolution_angstrom_from_file"] = res
    expected_res = CONFIG["receptor"]["resolution_angstrom"]
    if res is None:
        problems.append("could not read REMARK 2 RESOLUTION from the file")
    elif abs(res - expected_res) > 0.01:
        problems.append(f"file says {res} A, package expects {expected_res} A")

    # 4. Chain A present and the right size?
    chain_a = {rn: nm for (ch, rn), nm in s["protein"].items() if ch == "A"}
    results["chain_a_residues_observed"] = len(chain_a)
    if not chain_a:
        problems.append("chain A contains no ATOM records")
    else:
        lo, hi = min(chain_a), max(chain_a)
        results["chain_a_residue_range"] = [lo, hi]
        # Construct is 191-1051; disordered termini are normal, so this is a loose sanity bound.
        if lo < 150 or hi > 1100:
            problems.append(f"chain A spans {lo}-{hi}, outside the 191-1051 construct")
        if len(chain_a) < 700:
            problems.append(f"chain A has only {len(chain_a)} residues; expected roughly 860")

    # 5. The catalytic machinery, at GtfB numbering. This is the strongest single check:
    #    the same residue numbers in a GtfC structure would be different amino acids.
    seen = {}
    for resnum, expected_name in {**EXPECTED_CATALYTIC, **EXPECTED_CONTACTS}.items():
        actual = chain_a.get(resnum)
        seen[f"{expected_name}{resnum}"] = actual
        if actual is None:
            problems.append(f"residue {resnum} absent from chain A (expected {expected_name})")
        elif actual != expected_name:
            problems.append(
                f"residue {resnum} is {actual}, expected {expected_name}. "
                f"Either this is not GtfB or the numbering is not the deposited numbering."
            )
    results["catalytic_and_contact_residues"] = seen

    # 6. Is the inhibitor actually there?
    lig = {f"{ch}/{nm}/{rn}": n for (ch, nm, rn), n in s["hetatms"].items()
           if nm in LIGAND_RESNAMES}
    results["inhibitor_residues_found"] = lig
    if not lig:
        problems.append(
            "no GLC, AC1 or ACR residues found. Acarbose is deposited in 8FJC as two GLC "
            "moieties plus one AC1 moiety; without it there is no native ligand for the "
            "redocking calibration control."
        )

    # 7. Sulfate, worth knowing about rather than failing on.
    sulfates = [f"{ch}/{rn}" for (ch, nm, rn) in s["hetatms"] if nm == "SO4"]
    results["sulfate_sites"] = sulfates
    results["sulfate_note"] = (
        f"{len(sulfates)} sulfate ions present. 8FJC was grown from 2 M ammonium sulfate, "
        f"so a sulfate may sit in the negatively charged subsites. prepare_receptor.py "
        f"removes all heteroatoms and reports any that lay within 6 A of the box centre, "
        f"because an arginine-rich ligand docking into a site that a sulfate occupies in "
        f"the crystal is a result that needs stating, not hiding."
    )

    # 8. The conserved structural Ca2+ site. Reported, never a reason to stop.
    #    This exists because prepare_receptor.py strips every heteroatom, and a
    #    structural metal adjacent to the catalytic triad is not the same kind of
    #    heteroatom as a crystallisation sulfate. Stripping it leaves open volume where
    #    an ion sat, and AutoDock Vina types a retained metal as a hydrogen-bond donor,
    #    so the two receptors are not interchangeable in either direction. The package
    #    therefore docks against both and reports the difference instead of assuming one.
    ca_site_seen = {}
    for resnum, expected_name in CALCIUM_SITE_STATED.items():
        actual = chain_a.get(resnum)
        ca_site_seen[f"{expected_name}{resnum}"] = actual
        if actual != expected_name:
            problems.append(
                f"residue {resnum} is {actual}, expected {expected_name}. Schormann et al. "
                f"2023 state that GtfB Asp411 is one of the two Ca2+-coordinating "
                f"aspartates, so this is a direct check on the numbering, not a derived one."
            )
    ca_site_derived_seen, ca_site_derived_mismatches = {}, []
    for resnum, expected_name in CALCIUM_SITE_DERIVED.items():
        actual = chain_a.get(resnum)
        ca_site_derived_seen[f"{expected_name}{resnum}"] = actual
        if actual != expected_name:
            ca_site_derived_mismatches.append(f"{resnum}: expected {expected_name}, got {actual}")
    results["calcium_site_residues_stated_in_paper"] = ca_site_seen
    results["calcium_site_residues_derived_at_offset_26"] = ca_site_derived_seen
    results["calcium_site_derived_mismatches"] = ca_site_derived_mismatches

    results["metal_ions_in_file"] = s["metals"]
    ca_ions = [m for m in s["metals"] if m["resname"] == "CA"]
    results["calcium_ions_found"] = len(ca_ions)
    if ca_ions:
        results["calcium_note"] = (
            f"{len(ca_ions)} Ca2+ ion(s) present. Schormann et al. 2023 report a calcium ion "
            f"in the intact GtfB catalytic domain, at a site whose coordinating residues are "
            f"conserved across GtfB, GtfC and GtfD, and one of those residues (Asn455 in GtfB "
            f"numbering) also hydrogen bonds acarbose, so the metal sits immediately beside the "
            f"inhibitor subsites. prepare_receptor.py builds one receptor without it and one "
            f"with it, and tier 6 docks against both."
        )
    else:
        results["calcium_note"] = (
            "No Ca2+ ion found in this file. Either the deposited model does not include it or "
            "the record is named differently. Tier 6 will be skipped and every result will come "
            "from the heteroatom-stripped receptor, which is a limitation to state rather than "
            "a failure: an apo-metal receptor is not the physiological enzyme."
        )

    results["n_protein_atoms"] = s["n_atom"]
    results["sha256"] = sha256(path)
    results["path"] = str(path)

    if problems:
        raise Failed("receptor validation failed:\n  - " + "\n  - ".join(problems))
    return results


def main() -> int:
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / "receptor" / "8FJC.pdb"
    OUT.mkdir(parents=True, exist_ok=True)
    stamp_path = OUT / "receptor_validation.json"
    try:
        results = check(path)
    except Failed as exc:
        print(f"STOP: {exc}", file=sys.stderr)
        return 1

    if stamp_path.exists():
        previous = json.loads(stamp_path.read_text()).get("sha256")
        if previous and previous != results["sha256"]:
            print(
                "STOP: the receptor file has changed since the first validated run.\n"
                f"  recorded: {previous}\n  now:      {results['sha256']}\n"
                "Results already in out/results.jsonl were produced against the recorded "
                "file. Mixing two receptors in one result set would be unrecoverable, so "
                "either restore the original file or move out/ aside and start again.",
                file=sys.stderr,
            )
            return 1

    results["validated"] = True
    stamp_path.write_text(json.dumps(results, indent=2) + "\n")

    print("Receptor validation passed.")
    print(f"  {results['title']}")
    print(f"  resolution read from file: {results['resolution_angstrom_from_file']} A")
    print(f"  chain A residues: {results['chain_a_residues_observed']} "
          f"spanning {results.get('chain_a_residue_range')}")
    print(f"  catalytic residues confirmed at GtfB numbering: "
          f"{', '.join(k for k, v in results['catalytic_and_contact_residues'].items() if v)}")
    print(f"  inhibitor residues: {results['inhibitor_residues_found']}")
    print(f"  Ca2+ ions found: {results['calcium_ions_found']}")
    if results["calcium_site_derived_mismatches"]:
        print("  NOTE: a derived Ca2+-site residue number did not match: "
              + "; ".join(results["calcium_site_derived_mismatches"])
              + ". The stated Asp411 check passed, so the structure is GtfB; the offset "
                "derivation for the other three is what is in doubt. Reported, not fatal.")
    print(f"  sha256: {results['sha256']}")
    print(f"  recorded in {stamp_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
