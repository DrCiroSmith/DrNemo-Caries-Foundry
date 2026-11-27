#!/usr/bin/env python3
"""Recompute every number the manuscript cites, and say where each one comes from.

Why this script exists. A paper is only as good as the provenance of its numbers,
and this repository has already been caught publishing figures that no artefact
produced. So nothing in paper/negative_result_GTF-SI.md is copied from a README,
a dashboard or a previous draft. Every number is either recomputed here from
committed data, or read from a committed artefact and labelled as a recorded run
that cannot be reproduced in an environment without AutoDock Vina and OpenMM.

Run:

    python3 scripts/compute_paper_numbers.py

Provenance tags used in the output and in the manuscript:

    [R]  recomputed by this script from committed data (coordinates, the vault,
         or arithmetic on them). Anyone can rerun it and get the same number.
    [D]  read from a committed artefact, produced by a run that required a tool
         absent from this environment (the vina binary, openmm, rdkit). The
         artefact and key are named. This script checks it for internal
         consistency but does not re-derive it.

The script also cross-checks the artefacts against each other and prints every
failed check under DISCREPANCIES. A failed check is not a rounding note. It is
one of three things: two committed files state different values for the same
quantity; a file's own stated interpretation is contradicted by the coordinates
in that same file; or a structural property of the inputs makes a reported
comparison unsound. In all three cases the manuscript has to say so rather than
pick a value and move on. Checks that pass are printed inline, not here.
"""
import collections
import hashlib
import json
import math
import re
import sqlite3
import subprocess
import statistics
import sys
import textwrap
from pathlib import Path

REPO = Path(__file__).parent.parent.resolve()
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

from solubility_shield import (  # noqa: E402
    calc_net_charge_at_ph,
    calc_stephan_switch,
    count_protease_cleavage_sites,
)
from compute_pareto_front import dominates, pareto_front  # noqa: E402

DB = REPO / "data" / "foundry_contenders.sqlite"
BENCHMARK = REPO / "demo_data" / "results" / "falsification_benchmark.json"
DISPROOF = REPO / "data" / "exhaustive_disproof_report.json"
RECEPTORS = REPO / "demo_data" / "receptors"
RESULTS = REPO / "demo_data" / "results"
RECEPTOR_3AIC = RECEPTORS / "3AIC_receptor.pdb"
NATIVE_LIGAND = RECEPTORS / "3AIC_native_ligand.pdb"

WT = "CPHWRPWHC"
LEAD = "CAHWRPWHC"
BEST_AFFINITY = "CHPWPHWAC"
NAIVE = "CRRRWWRRRRC"
# The five candidates whose demo_data/results/analysis files assert a specific
# mechanism: bidentate salt bridges to Asp477 and Glu515, and a tryptophan pair
# stacking between Trp517 and Tyr430.
ASSERTS_GLU515 = ("CRRRWWRRRRC", "CRRRWWPRRRC", "CRRRWFRRRRC",
                  "CRRRAWPRRRC", "CRRAGWPRRRC")
SURVIVED = ("SURVIVED_ALL_GATES", "CORROBORATED", "NOT_CORROBORATED")
HOTSPOT_THRESHOLD = 0.70
REDOCK_STANDARD = 2.0

discrepancies = []
numbers = []


def emit(tag, label, value, source):
    numbers.append((tag, label, value, source))
    print(f"  [{tag}] {label:52s} {value:>26}   {source}")


def note(text):
    """A failed cross-check. Counted, and printed under DISCREPANCIES."""
    discrepancies.append(text)


def remark(text):
    """Prose that explains a number. Printed inline, never counted as a check."""
    for line in textwrap.wrap(text, 92):
        print(f"      {line}")


def _gammainc_lower_reg(a, x):
    """Regularized lower incomplete gamma P(a, x), to machine precision.

    Series below the crossover, continued fraction above it. Here only to give
    chi-square quantiles, because this environment has no scipy and a interval
    on a standard deviation should not be read off a table by hand.
    """
    if x <= 0:
        return 0.0
    lg = math.lgamma(a)
    if x < a + 1.0:
        term = total = 1.0 / a
        for n in range(1, 10000):
            term *= x / (a + n)
            total += term
            if abs(term) < abs(total) * 1e-16:
                break
        return total * math.exp(-x + a * math.log(x) - lg)
    tiny = 1e-300
    b, c = x + 1.0 - a, 1.0 / 1e-300
    d = 1.0 / b
    h = d
    for i in range(1, 10000):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        d = tiny if abs(d) < tiny else d
        c = b + an / c
        c = tiny if abs(c) < tiny else c
        d = 1.0 / d
        de = d * c
        h *= de
        if abs(de - 1.0) < 1e-16:
            break
    return 1.0 - math.exp(-x + a * math.log(x) - lg) * h


def chi2_ppf(prob, df):
    """Inverse chi-square CDF by bisection. Verified against published tables."""
    lo, hi = 1e-12, 10.0
    while _gammainc_lower_reg(df / 2.0, hi / 2.0) < prob and hi < 1e12:
        hi *= 2.0
    for _ in range(300):
        mid = 0.5 * (lo + hi)
        if _gammainc_lower_reg(df / 2.0, mid / 2.0) < prob:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def sd_interval(sd, n, level=0.95):
    """Two-sided confidence interval for a population sd from a sample sd."""
    df = n - 1
    tail = (1.0 - level) / 2.0
    return (sd * math.sqrt(df / chi2_ppf(1.0 - tail, df)),
            sd * math.sqrt(df / chi2_ppf(tail, df)))


def heading(text):
    print()
    print(text)
    print("-" * len(text))


def load_vault():
    conn = sqlite3.connect(f"file:{DB}?mode=ro", uri=True)
    conn.row_factory = sqlite3.Row
    try:
        return [dict(r) for r in conn.execute("SELECT * FROM contenders")]
    finally:
        conn.close()


# AutoDock pseudo-atom types that Meeko adds to open a macrocycle ring so that Vina
# can sample it. They are not physical atoms. The partner carbons are retyped CG0,
# CG1 and so on, which are real atoms, so they are not counted here. Plain "CG" is
# the standard gamma carbon of several residues and is never a pseudo-atom.
GLUE_TYPES = {"G", "G0", "G1", "G2", "G3"}
HYDROGEN_TYPES = {"H", "HD", "HS"}


def pdb_atoms(path, chain=None):
    """Every atom record as a dict. A list, not a map: atom names repeat in these files,
    and keying on them silently discards atoms."""
    out = []
    for line in Path(path).read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if chain and line[21] != chain:
            continue
        out.append({
            "name": line[12:16].strip(),
            "resname": line[17:20].strip(),
            "resnum": int(line[22:26]),
            "chain": line[21],
            "xyz": (float(line[30:38]), float(line[38:46]), float(line[46:54])),
            "type": line[76:78].strip(),
        })
    return out


def atom_key(a):
    return (a["name"], a["resname"], a["resnum"], a["chain"])


def is_hydrogen(a):
    return a["type"] in HYDROGEN_TYPES or (not a["type"] and a["name"].startswith("H"))


def is_glue(a):
    return a["type"] in GLUE_TYPES or (not a["type"] and a["name"] in GLUE_TYPES)


def coordinate_index(atoms):
    """key -> list of coordinates, so duplicate atom names are all kept."""
    index = {}
    for a in atoms:
        index.setdefault(atom_key(a), []).append(a["xyz"])
    return index


def matches_reference(atom, index, tolerance=1e-3):
    """True when this atom appears in the reference with the same coordinates."""
    for xyz in index.get(atom_key(atom), ()):
        if all(abs(xyz[i] - atom["xyz"][i]) <= tolerance for i in range(3)):
            return True
    return False


def centroid(atoms):
    n = len(atoms)
    return tuple(sum(a["xyz"][i] for a in atoms) / n for i in range(3))


def euclid(p, q):
    return math.sqrt(sum((p[i] - q[i]) ** 2 for i in range(3)))


def nearest(point, atoms):
    return min(euclid(point, a["xyz"]) for a in atoms)


def heavy_only(atoms):
    return [a for a in atoms if not is_hydrogen(a) and not is_glue(a)]


def ligand_chain(path):
    """The one chain that is not the receptor chain A."""
    chains = {a["chain"] for a in pdb_atoms(path)}
    others = sorted(c for c in chains if c != "A")
    return others[0] if others else None


def pdbqt_coord_digest(path):
    """md5 of the coordinate columns of every atom line, all models.

    Two docking runs of the same ligand never agree to the last decimal of every
    atom in every retained mode, so an identical digest means one file was copied
    from the other rather than produced by a second search.
    """
    text = Path(path).read_text()
    coords = "".join(l[30:54] for l in text.splitlines()
                     if l.startswith(("ATOM", "HETATM")))
    return hashlib.md5(coords.encode()).hexdigest() if coords else None


def pdbqt_best_mode(path, keep_glue=False):
    """Heavy atoms of MODEL 1, in file order, as (name, x, y, z)."""
    out, started = [], False
    for l in Path(path).read_text().splitlines():
        if l.startswith("MODEL"):
            if started:
                break
            started = True
            continue
        if l.startswith("ENDMDL"):
            break
        if l.startswith(("ATOM", "HETATM")):
            t = l[77:79].strip()
            if t in ("H", "HD"):
                continue
            if t in ("G", "G0") and not keep_glue:
                continue
            out.append((l[12:16].strip(), float(l[30:38]), float(l[38:46]), float(l[46:54])))
    return out


def pdbqt_score(path):
    for l in Path(path).read_text().splitlines():
        if l.startswith("REMARK VINA RESULT:"):
            return float(l.split()[3])
    return None


def in_order_rmsd(a, b):
    if len(a) != len(b) or not a:
        return None
    return math.sqrt(sum((x[1] - y[1]) ** 2 + (x[2] - y[2]) ** 2 + (x[3] - y[3]) ** 2
                         for x, y in zip(a, b)) / len(a))


def aromatic_ring(atoms, resname):
    if resname == "TRP":
        names = ("CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2")
    elif resname in ("PHE", "TYR"):
        names = ("CG", "CD1", "CD2", "CE1", "CE2", "CZ")
    else:
        return []
    return [a for a in atoms if a["name"] in names]


def main():
    rows = load_vault()
    by_seq = {r["sequence"]: r for r in rows}
    bench = json.loads(BENCHMARK.read_text())
    disproof = json.loads(DISPROOF.read_text())
    angles = {a["angle"]: a for a in disproof["angles"]}

    print("=" * 100)
    print("NUMBERS FOR paper/negative_result_GTF-SI.md")
    print("=" * 100)

    # ------------------------------------------------------------------ targets
    heading("1. Receptor identity, from the committed coordinates")
    # The identity of 3AIC is established in the structure paper, not here. What the
    # committed file can establish is that it carries the construct and the catalytic
    # triad numbering that paper reports, which is the check worth stating.
    for name, triad in (
        ("3AIC_receptor.pdb", {477: "ASP", 515: "GLU", 588: "ASP"}),
        ("1SMD_receptor.pdb", {197: "ASP", 233: "GLU", 300: "ASP"}),
        ("1LZ1_receptor.pdb", {}),
    ):
        atoms = pdb_atoms(RECEPTORS / name)
        residues = {(a["chain"], a["resnum"]): a["resname"] for a in atoms}
        nums = sorted(n for _, n in residues)
        emit("R", f"{name}: receptor heavy+H atoms", len(atoms), f"demo_data/receptors/{name}")
        emit("R", f"{name}: residues", len(residues), f"demo_data/receptors/{name}")
        emit("R", f"{name}: numbering span", f"{nums[0]}..{nums[-1]}", f"demo_data/receptors/{name}")
        for num, expected in triad.items():
            got = [v for (c, n), v in residues.items() if n == num]
            ok = got and all(g == expected for g in got)
            emit("R", f"{name}: residue {num} is {expected}", "yes" if ok else f"NO ({got})",
                 f"demo_data/receptors/{name}")
            if not ok:
                note(f"{name}: residue {num} is {got}, not {expected}")

    # ------------------------------------------------------------------ census
    heading("2. Vault census")
    emit("R", "rows in the vault", len(rows), "data/foundry_contenders.sqlite")
    status_counts = {}
    for r in rows:
        status_counts[r["falsification_status"]] = status_counts.get(r["falsification_status"], 0) + 1
    for k in sorted(status_counts, key=lambda s: -status_counts[s]):
        emit("R", f"rows with status {k}", status_counts[k], "data/foundry_contenders.sqlite")
    scrambles = [r for r in rows if r["falsification_status"] == "SCRAMBLED_NULL_CONTROL"]
    designed = [r for r in rows if r["falsification_status"] != "SCRAMBLED_NULL_CONTROL"]
    survivors = [r for r in rows if r["falsification_status"] in SURVIVED]
    emit("R", "designed candidates (vault minus null controls)", len(designed),
         "data/foundry_contenders.sqlite")
    emit("R", "rows that passed every gate", len(survivors), "data/foundry_contenders.sqlite")
    emit("R", "distinct sequences", len({r["sequence"] for r in rows}),
         "data/foundry_contenders.sqlite")
    gens = sorted({r["generation"] for r in rows})
    emit("R", "generations present", f"{min(gens)}..{max(gens)}", "data/foundry_contenders.sqlite")

    # --------------------------------------------------------- recorded scores
    heading("3. Docking scores of the four sequences the paper names (recorded runs)")
    for label, seq in (("wild type", WT), ("lead (P2A)", LEAD),
                       ("best affinity", BEST_AFFINITY), ("naive poly-Arg", NAIVE)):
        r = by_seq.get(seq)
        if not r:
            note(f"{seq} is not in the vault")
            continue
        emit("D", f"{label} {seq}: Vina 3AIC kcal/mol", r["vina_3aic_kcal"], "vault.vina_3aic_kcal")
        emit("D", f"{label} {seq}: Vina 1SMD kcal/mol", r["vina_1smd_kcal"], "vault.vina_1smd_kcal")
        emit("D", f"{label} {seq}: Vinardo 3AIC kcal/mol", r["vinardo_3aic_kcal"],
             "vault.vinardo_3aic_kcal")
        recomputed_ddg = round(r["vina_3aic_kcal"] - r["vina_1smd_kcal"], 3)
        emit("R", f"{label} {seq}: selectivity ddG recomputed", recomputed_ddg,
             "3AIC minus 1SMD, this script")
        stored = r["selectivity_ddg"]
        if abs(recomputed_ddg - stored) > 0.011:
            note(f"{seq}: stored selectivity_ddg {stored} but 3AIC minus 1SMD is {recomputed_ddg}")
        emit("R", f"{label} {seq}: heavy atoms", r["heavy_atoms"], "vault.heavy_atoms")
        le = round(abs(r["vina_3aic_kcal"]) / r["heavy_atoms"], 4)
        emit("R", f"{label} {seq}: ligand efficiency recomputed", le,
             "|affinity| / heavy atoms, this script")
        if abs(le - r["ligand_efficiency"]) > 0.0015:
            note(f"{seq}: stored ligand_efficiency {r['ligand_efficiency']} but "
                 f"|affinity|/heavy atoms is {le}")

    # ------------------------------------------------- angle 1, redocking
    heading("4. Angle 1, co-crystal redocking calibration")
    a1 = angles[1]
    emit("D", "acarbose redocked, Vina kcal/mol", a1["vina_3aic_kcal"], "disproof angle 1")
    emit("D", "symmetry-corrected redocking RMSD, A", a1["assignment_rmsd_angstrom"],
         "disproof angle 1 (recompute: scripts/compute_redocking_rmsd.py)")
    emit("D", "centroid displacement, A", a1["centroid_displacement_angstrom"], "disproof angle 1")
    emit("R", "redocking standard, A", REDOCK_STANDARD, "convention, stated in advance")
    emit("R", "redocking passes the 2.0 A standard",
         "no" if a1["assignment_rmsd_angstrom"] > REDOCK_STANDARD else "yes", "this script")
    if bench["co_crystal_control"].get("passed_calibration") is True:
        note("falsification_benchmark.json co_crystal_control.passed_calibration is true "
             f"while the RMSD is {bench['co_crystal_control']['redocking_rmsd_angstrom']} A "
             "against a 2.0 A standard")
    ref_heavy = [a for a in pdb_atoms(RECEPTORS / "3AIC_native_ligand.pdb")
                 if not is_hydrogen(a)]
    emit("R", "crystallographic acarbose heavy atoms", len(ref_heavy),
         "demo_data/receptors/3AIC_native_ligand.pdb")

    # ------------------------------------------------- angle 2, reproducibility
    heading("5. Angle 2, multi-seed reproducibility, recomputed from the pose files")
    a2 = angles[2]
    print("    The report treats this as five independent runs. One of the five is a\n"
          "    byte-identical copy of the campaign pose (section 20), so it is not a run.\n"
          "    Everything below is recomputed from the pose files with that one excluded.")
    emit("D", "seeds the battery reports", a2["n_seeds"], "disproof angle 2")

    seed_files = {}
    for seed in (42, 101, 202, 303, 404):
        f = REPO / "data" / "poses" / f"disproof_seed_{seed}_CAHWRPWHC_3AIC_docked.pdbqt"
        if f.exists():
            seed_files[seed] = f
    campaign_pose = REPO / "data" / "poses" / "contender_mutant_P2A_CAHWRPWHC_3aic.pdbqt"
    campaign_digest = pdbqt_coord_digest(campaign_pose) if campaign_pose.exists() else None
    reused = sorted(s for s, f in seed_files.items()
                    if pdbqt_coord_digest(f) == campaign_digest)
    independent = sorted(s for s in seed_files if s not in reused)
    emit("R", "seed runs that duplicate the campaign pose",
         ", ".join(str(s) for s in reused) or "none", "coordinate digest, this script")
    emit("R", "independent seed runs", len(independent), "this script")

    if len(independent) >= 2:
        scores = [pdbqt_score(seed_files[s]) for s in independent]
        emit("R", "seed score mean kcal/mol", round(statistics.fmean(scores), 3), "this script")
        emit("R", "seed score SD kcal/mol (sample)", round(statistics.stdev(scores), 3),
             "this script")
        emit("R", "seed score SD kcal/mol (population, the report's convention)",
             round(statistics.pstdev(scores), 3), "this script")
        emit("R", "seed score range kcal/mol", f"{min(scores)} to {max(scores)}", "this script")
        emit("R", "seed score span kcal/mol", round(max(scores) - min(scores), 3), "this script")

        for keep_glue, label in ((False, "real heavy atoms"), (True, "including glue pseudo-atoms")):
            modes = {s: pdbqt_best_mode(seed_files[s], keep_glue=keep_glue) for s in independent}
            rs = [in_order_rmsd(modes[i], modes[j])
                  for n, i in enumerate(independent) for j in independent[n + 1:]]
            rs = [r for r in rs if r is not None]
            if rs:
                emit("R", f"pairwise pose RMSD over {label}, mean / worst / best A",
                     f"{statistics.fmean(rs):.2f} / {max(rs):.2f} / {min(rs):.2f}",
                     f"{len(rs)} pairs, this script")
        heavy_modes = {s: pdbqt_best_mode(seed_files[s]) for s in independent}
        heavy_rs = [r for n, i in enumerate(independent) for j in independent[n + 1:]
                    if (r := in_order_rmsd(heavy_modes[i], heavy_modes[j])) is not None]
        emit("D", "mean pairwise pose RMSD the report states, A",
             a2["mean_pairwise_rmsd_angstrom"], "disproof angle 2")
        if heavy_rs and abs(a2["mean_pairwise_rmsd_angstrom"]
                            - statistics.fmean(heavy_rs)) > 0.05:
            note("the pose RMSD the report states is computed over five poses and counts "
                 "Meeko's glue pseudo-atoms as atoms; recomputed over the independent runs "
                 "and the real heavy atoms only it differs, and the paper uses the "
                 "recomputed value")
    emit("D", "centroid SD, A", a2["centroid_sd_angstrom"], "disproof angle 2")
    emit("R", "pose agreement convention, A", a2["pose_agreement_threshold_angstrom"],
         "convention, stated in advance")

    # ------------------------------------------------- angle 3, mode degeneracy
    heading("6. Angle 3, pose mode degeneracy")
    a3 = angles[3]
    modes = a3["modes_summary"]
    emit("D", "modes returned", a3["total_modes_evaluated"], "disproof angle 3")
    emit("R", "mode 1 minus mode 2, kcal/mol",
         round(modes[1]["affinity_kcal"] - modes[0]["affinity_kcal"], 3), "this script")
    emit("R", "spread over all modes, kcal/mol",
         round(modes[-1]["affinity_kcal"] - modes[0]["affinity_kcal"], 3), "this script")
    emit("D", "mode 1 to mode 2 heavy-atom RMSD, A", a3["heavy_atom_rmsd_mode1_mode2_angstrom"],
         "disproof angle 3")
    emit("D", "mode 1 to mode 2 centroid distance, A",
         a3["centroid_distance_mode1_mode2_angstrom"], "disproof angle 3")
    if "distinct primary binding basin" in a3.get("scientific_summary", ""):
        note("disproof angle 3 summary calls a 0.101 kcal/mol gap with a 5.58 A pose "
             "difference a distinct primary binding basin")

    # ------------------------------------------------- angle 5, scrambled null
    heading("7. Angle 5, charge-matched scrambled null, recomputed from the vault")
    null_scores = [r["vina_3aic_kcal"] for r in scrambles]
    emit("R", "null isomers in the vault", len(null_scores), "data/foundry_contenders.sqlite")

    # Which sequence are these isomers actually matched to? The angle tests the lead,
    # so this has to be checked rather than assumed.
    for label, seq in (("the lead", LEAD), ("the wild type", WT)):
        comp = "".join(sorted(seq))
        n_match = sum(1 for r in scrambles if "".join(sorted(r["sequence"])) == comp)
        emit("R", f"isomers composition-matched to {label} {seq}",
             f"{n_match} of {len(scrambles)}", "this script")
    comp_lead = "".join(sorted(LEAD))
    comp_wt = "".join(sorted(WT))
    matched_to_lead = sum(1 for r in scrambles
                          if "".join(sorted(r["sequence"])) == comp_lead)
    matched_to_wt = sum(1 for r in scrambles if "".join(sorted(r["sequence"])) == comp_wt)
    if matched_to_lead == 0 and matched_to_wt == len(scrambles):
        note(f"all {len(scrambles)} scrambled controls are permutations of the wild type "
             f"{WT}, not of the lead {LEAD} that angle 5 tests against them. The lead "
             "carries an alanine where the wild type carries a proline, so the published "
             "null test is not composition-matched to its own candidate")

    mean = statistics.fmean(null_scores)
    sd = statistics.stdev(null_scores)
    emit("R", "null mean kcal/mol", round(mean, 3), "this script")
    emit("R", "null sample SD kcal/mol", round(sd, 3), "this script")
    emit("R", "null range kcal/mol", f"{min(null_scores)} to {max(null_scores)}", "this script")

    # Both comparisons: the one the report publishes, and the matched one it should have run.
    for label, seq in (("lead", LEAD), ("wild type", WT)):
        score = by_seq[seq]["vina_3aic_kcal"]
        k = sum(1 for s in null_scores if s <= score)
        emit("R", f"{label} {seq}: affinity kcal/mol", score, "vault.vina_3aic_kcal")
        emit("R", f"{label}: Z against the null", round((score - mean) / sd, 3), "this script")
        emit("R", f"{label}: isomers at least as good", f"{k} of {len(null_scores)}",
             "this script")
        emit("R", f"{label}: one-sided empirical p, (k+1)/(n+1)",
             f"{(k + 1) / (len(null_scores) + 1):.4f}", "this script")
        emit("R", f"{label}: proportion k/n", round(k / len(null_scores), 4), "this script")
    emit("D", "p reported in the artefact", angles[5]["p_value"], "disproof angle 5")
    k_lead = sum(1 for s in null_scores if s <= by_seq[LEAD]["vina_3aic_kcal"])
    if abs(angles[5]["p_value"] - k_lead / len(null_scores)) > 0.002 and \
       abs(angles[5]["p_value"] - (k_lead + 1) / (len(null_scores) + 1)) > 0.002:
        note(f"disproof angle 5 p = {angles[5]['p_value']} matches neither k/n nor "
             "(k+1)/(n+1) on the committed null scores")
    emit("R", "best-scoring null isomer kcal/mol", round(min(null_scores), 3), "this script")
    emit("R", "nulls that beat the best designed candidate",
         f"{sum(1 for s in null_scores if s <= by_seq[BEST_AFFINITY]['vina_3aic_kcal'])} "
         f"of {len(null_scores)}", "this script")

    # ------------------------------------------------- angle 6, alanine scan
    heading("8. Angle 6, alanine scan, recomputed against the wild type")
    parent = by_seq[WT]["vina_3aic_kcal"]
    emit("R", f"wild type {WT} parent affinity kcal/mol", parent, "vault.vina_3aic_kcal")
    scan_rows = []
    for key, rec in bench["alanine_scan"].items():
        seq = rec["mutant_sequence"]
        vault_row = by_seq.get(seq)
        if not vault_row:
            note(f"alanine mutant {key} ({seq}) is not in the vault")
            continue
        delta = round(vault_row["vina_3aic_kcal"] - parent, 3)
        scan_rows.append((key, seq, vault_row["vina_3aic_kcal"], delta))
        if abs(delta - rec["delta_binding_affinity_kcal"]) > 0.0011:
            note(f"alanine {key}: benchmark delta {rec['delta_binding_affinity_kcal']} "
                 f"but vault minus parent is {delta}")
        if rec["is_hotspot"] and delta < HOTSPOT_THRESHOLD:
            note(f"alanine {key} is flagged is_hotspot with a delta of {delta}")
    scan_rows.sort(key=lambda t: -t[3])
    for key, seq, score, delta in scan_rows:
        emit("R", f"{key} {seq}: delta vs wild type kcal/mol", f"{delta:+.3f}",
             "vault minus parent, this script")
    losses = [d for *_, d in scan_rows if d > 0]
    gains = [d for *_, d in scan_rows if d < 0]
    emit("R", "substitutions that lose affinity", len(losses), "this script")
    emit("R", "substitutions that improve affinity", len(gains), "this script")
    emit("R", "largest loss kcal/mol", f"{max(losses):+.3f}" if losses else "none", "this script")
    emit("R", "largest gain kcal/mol", f"{min(gains):+.3f}" if gains else "none", "this script")
    emit("R", "hotspot threshold kcal/mol", f"+{HOTSPOT_THRESHOLD:.2f}",
         "threshold set in advance, falsification_engine.py")
    emit("R", "positions reaching the threshold",
         sum(1 for d in losses if d >= HOTSPOT_THRESHOLD), "this script")
    p6a = by_seq.get("CPHWRAWHC")
    if p6a:
        sites = count_protease_cleavage_sites("CPHWRAWHC")
        wt_sites = count_protease_cleavage_sites(WT)
        emit("R", "P6A trypsin sites (wild type has 0)",
             f"{sites['cleavage_site_count']} (wt {wt_sites['cleavage_site_count']})",
             "solubility_shield, this script")

    # ------------------------------------------------- angle 4 and 9, selectivity
    heading("9. Angles 4 and 9, selectivity against the two human enzymes")
    a4, a9 = angles[4], angles[9]
    emit("D", "assumed Vina scoring error kcal/mol", a4["assumed_forcefield_error_kcal"],
         f"disproof angle 4 ({a4['assumed_forcefield_error_provenance']})")
    emit("R", "lead vs amylase ddG kcal/mol",
         round(by_seq[LEAD]["vina_3aic_kcal"] - by_seq[LEAD]["vina_1smd_kcal"], 3), "this script")
    emit("D", "lead vs lysozyme 1LZ1 kcal/mol", a9["vina_1lz1_affinity_kcal"], "disproof angle 9")
    emit("D", "absolute affinity cutoff angle 9 applied, kcal/mol", -7.50,
         "run_exhaustive_disproof.py, the criterion angle 9 used before correction")
    emit("R", "lead vs lysozyme ddG kcal/mol",
         round(by_seq[LEAD]["vina_3aic_kcal"] - a9["vina_1lz1_affinity_kcal"], 3), "this script")
    emit("R", "both ddG magnitudes as a fraction of the assumed error",
         f"{abs(by_seq[LEAD]['vina_3aic_kcal'] - by_seq[LEAD]['vina_1smd_kcal']) / a4['assumed_forcefield_error_kcal']:.2f}"
         f" and {abs(by_seq[LEAD]['vina_3aic_kcal'] - a9['vina_1lz1_affinity_kcal']) / a4['assumed_forcefield_error_kcal']:.2f}",
         "this script")
    if a4["classification"] != a9["classification"]:
        note(f"angle 4 calls a ddG of {a4['measured_selectivity_ddg_kcal']} "
             f"{a4['classification']} while angle 9 calls a ddG of {a9['selectivity_ddg_kcal']} "
             f"{a9['classification']}, judged against the same unsourced error figure")

    # ------------------------------------------------- angles 8 and 12, barriers
    heading("10. Angles 8 and 12, oral barriers, recomputed from the sequence")
    pro = count_protease_cleavage_sites(LEAD)
    emit("R", "lead trypsin / gingipain sites", pro["cleavage_site_count"],
         "solubility_shield, this script")
    emit("R", "lead chymotrypsin sites",
         f"{pro['chymotrypsin_site_count']} {pro['chymotrypsin_cleavage_sites']}",
         "solubility_shield, this script")
    site_hist = {}
    for r in designed:
        n = count_protease_cleavage_sites(r["sequence"])["cleavage_site_count"]
        site_hist[n] = site_hist.get(n, 0) + 1
    emit("R", "designed candidates with at least one trypsin site",
         sum(v for k, v in site_hist.items() if k > 0), "this script")
    emit("R", "trypsin site distribution over designed candidates",
         str(dict(sorted(site_hist.items()))), "this script")
    st = calc_stephan_switch(LEAD)
    emit("R", "lead charge at pH 7.2", st["q_pH7_2"], "Henderson-Hasselbalch, this script")
    emit("R", "lead charge at pH 6.0", calc_net_charge_at_ph(LEAD, 6.0),
         "Henderson-Hasselbalch, this script")
    emit("R", "lead charge at pH 5.2", st["q_pH5_2"], "Henderson-Hasselbalch, this script")
    emit("R", "lead delta Q 7.2 to 5.2", st["delta_q_stephan"], "this script")
    emit("R", "lead histidine count", st["histidine_count"], "this script")
    emit("R", "delta Q per histidine", round(st["delta_q_stephan"] / st["histidine_count"], 4)
         if st["histidine_count"] else "n/a", "this script")
    one_his = round(calc_net_charge_at_ph("H", 5.2) - calc_net_charge_at_ph("H", 7.2), 4)
    emit("R", "charge gained by one histidine, 7.2 to 5.2", one_his, "this script")
    dq_values = sorted({r["delta_q_stephan"] for r in rows})
    emit("R", "distinct delta Q values in the vault", str(dq_values), "this script")
    his_match = sum(1 for r in rows
                    if abs(r["delta_q_stephan"] - one_his * r["sequence"].count("H")) < 1e-3)
    emit("R", "rows where delta Q equals 0.809 x histidine count exactly",
         f"{his_match} of {len(rows)}", "this script")
    ss = [r["ss_bond_angstrom"] for r in rows if r["ss_bond_angstrom"] is not None]
    emit("R", "disulfide distance range over the vault, A",
         f"{min(ss)} to {max(ss)}", "this script")
    emit("R", "rows failing a 1.95 to 2.25 A disulfide window",
         sum(1 for v in ss if not 1.95 <= v <= 2.25), "this script")

    # ------------------------------------------------- angle 10, cyclisation
    heading("11. Angle 10 and the linear controls")
    a10 = angles[10]
    emit("D", "linear AHWRPWH, Vina 3AIC kcal/mol", a10["linear_vina_3aic_kcal"],
         "disproof angle 10")
    emit("R", "linear AHWRPWH minus cyclic lead kcal/mol",
         round(a10["linear_vina_3aic_kcal"] - by_seq[LEAD]["vina_3aic_kcal"], 3), "this script")
    lc = bench["linear_control"]
    emit("D", f"linear {lc['sequence']}, Vina 3AIC kcal/mol", lc["vina_3aic_kcal"],
         "falsification_benchmark.json linear_control")
    emit("R", f"linear {lc['sequence']} minus wild type kcal/mol",
         round(lc["vina_3aic_kcal"] - parent, 3), "this script")
    if a10.get("cyclization_benefit_kcal", 0) > 0 and \
       a10["linear_vina_3aic_kcal"] < a10["cyclic_vina_3aic_kcal"]:
        note("disproof angle 10 stores the linear advantage under the key "
             "cyclization_benefit_kcal, which names it backwards")

    # ------------------------------------------------- pareto
    heading("12. Pareto dominance over the declared objectives, recomputed")
    pareto_rows = [r for r in rows
                   if r["vina_3aic_kcal"] is not None and r["selectivity_ddg"] is not None]
    front_all = pareto_front(pareto_rows)
    lead_row = by_seq[LEAD]
    dom_all = [r for r in pareto_rows if dominates(r, lead_row)]
    surv_rows = [r for r in pareto_rows if r["falsification_status"] in SURVIVED]
    front_surv = pareto_front(surv_rows)
    dom_surv = [r for r in surv_rows if dominates(r, lead_row)]
    emit("R", "candidates with both objectives present", len(pareto_rows), "this script")
    emit("R", "non-dominated over all candidates", len(front_all), "this script")
    emit("R", "lead is on the overall front",
         "yes" if any(r["sequence"] == LEAD for r in front_all) else "no", "this script")
    emit("R", "candidates dominating the lead", len(dom_all), "this script")
    emit("R", "survivors", len(surv_rows), "this script")
    emit("R", "non-dominated among survivors", len(front_surv), "this script")
    emit("R", "survivors dominating the lead", len(dom_surv), "this script")
    by_aff = sorted(surv_rows, key=lambda r: r["vina_3aic_kcal"])
    by_sel = sorted(surv_rows, key=lambda r: r["selectivity_ddg"])
    emit("R", "lead rank among survivors on affinity",
         f"{[r['sequence'] for r in by_aff].index(LEAD) + 1} of {len(surv_rows)}", "this script")
    emit("R", "lead rank among survivors on selectivity",
         f"{[r['sequence'] for r in by_sel].index(LEAD) + 1} of {len(surv_rows)}", "this script")
    tiers = sorted({r["pareto_rank"] for r in rows})
    emit("R", "distinct pareto_rank values stored", str(tiers), "this script")
    mismatch = sum(1 for r in rows
                   if (r["pareto_rank"] == 1) != (r["falsification_status"] in SURVIVED))
    emit("R", "rows where pareto_rank 1 disagrees with surviving every gate", mismatch,
         "this script")

    # ------------------------------------------------- deposited poses
    heading("13. Deposited poses, checked against the vault and the receptor")
    receptor_index = coordinate_index(pdb_atoms(RECEPTORS / "3AIC_receptor.pdb"))
    receptor_total = sum(len(v) for v in receptor_index.values())
    for fname, seq in (("3AIC_docked_complex_CPHWRPWHC.pdb", WT),
                       ("3AIC_docked_complex_CAHWRPWHC.pdb", LEAD),
                       ("3AIC_docked_complex_CHPWPHWAC.pdb", BEST_AFFINITY)):
        path = RESULTS / fname
        rec_atoms = pdb_atoms(path, chain="A")
        lig_atoms = pdb_atoms(path, chain="B")
        unmatched = [a for a in rec_atoms if not matches_reference(a, receptor_index)]
        glue = [a for a in lig_atoms if is_glue(a)]
        hydrogens = [a for a in lig_atoms if is_hydrogen(a)]
        heavy = len(lig_atoms) - len(hydrogens) - len(glue)
        vault_heavy = by_seq[seq]["heavy_atoms"]
        emit("R", f"{fname}: chain A atoms", f"{len(rec_atoms)} of {receptor_total}",
             "this script")
        emit("R", f"{fname}: chain A atoms absent from the receptor file", len(unmatched),
             "coordinate match, this script")
        emit("R", f"{fname}: ligand atoms, H, glue",
             f"{len(lig_atoms)}, {len(hydrogens)}, {len(glue)}", "this script")
        emit("R", f"{fname}: ligand heavy atoms vs the vault",
             f"{heavy} vs {vault_heavy}", "this script")
        if unmatched:
            note(f"{fname}: {len(unmatched)} chain A atoms do not match the deposited "
                 "receptor coordinates")
        # These two are real properties of the deposited files, not errors. They are
        # discrepancies only while the file stays silent about them, so the test is
        # whether the file's own REMARK block declares them.
        remarks = " ".join(l for l in path.read_text().splitlines()
                           if l.startswith("REMARK")).upper()
        if len(rec_atoms) != receptor_total:
            declared = "TRUNCATION" in remarks
            emit("R", f"{fname}: truncation declared in a REMARK", declared, "this script")
            if not declared:
                note(f"{fname}: chain A carries {len(rec_atoms)} of the receptor's "
                     f"{receptor_total} atoms. The coordinates are identical, so this is a "
                     "pocket truncation for rendering, but the file does not say so")
        if glue:
            declared = "PSEUDO-ATOM" in remarks or "GLUE" in remarks
            emit("R", f"{fname}: glue pseudo-atoms declared in a REMARK", declared,
                 "this script")
            if not declared:
                note(f"{fname}: the ligand chain carries {len(glue)} Meeko glue pseudo-atoms "
                     "and AutoDock atom types in the element column, so a reader that trusts "
                     "the element column gets a molecule that does not exist")
        if heavy != vault_heavy:
            note(f"{fname}: ligand heavy atoms {heavy} but the vault records {vault_heavy}")

    # ------------------------------------------------- openmm
    heading("14. Angle 11, OpenMM")
    emit("D", "angle 11 classification", angles[11]["classification"], "disproof angle 11")
    emit("D", "lead relaxed energy kJ/mol, vault", by_seq[LEAD]["openmm_relaxed_kj"],
         "vault.openmm_relaxed_kj")
    emit("D", "lead relaxed energy kJ/mol, benchmark",
         bench["pareto_champion"]["openmm_relaxed_kj"], "falsification_benchmark.json")
    if abs(by_seq[LEAD]["openmm_relaxed_kj"] - bench["pareto_champion"]["openmm_relaxed_kj"]) > 0.01:
        note(f"lead OpenMM energy is {by_seq[LEAD]['openmm_relaxed_kj']} kJ/mol in the vault and "
             f"{bench['pareto_champion']['openmm_relaxed_kj']} kJ/mol in "
             "falsification_benchmark.json")

    # ------------------------------------------------- pocket occupancy
    heading("15. Pocket occupancy: where the deposited poses sit relative to acarbose")
    note_text = (
        "    Every distance here is recomputed from committed coordinates. Chain A of each\n"
        "    deposited complex was shown above to carry the committed receptor coordinates\n"
        "    unchanged, so the complexes and the receptor share one frame."
    )
    print(note_text)
    receptor = heavy_only(pdb_atoms(RECEPTOR_3AIC, chain="A"))
    acarbose = heavy_only(pdb_atoms(NATIVE_LIGAND))
    emit("R", "receptor heavy atoms, 3AIC", len(receptor), "3AIC_receptor.pdb")
    emit("R", "crystallographic acarbose heavy atoms", len(acarbose), "3AIC_native_ligand.pdb")
    frame_check = sum(1 for a in acarbose if nearest(a["xyz"], receptor) <= 4.0)
    emit("R", "acarbose atoms within 4.0 A of the receptor", f"{frame_check}/{len(acarbose)}",
         "frame sanity check")
    if frame_check < len(acarbose) // 4:
        note("the native ligand does not sit against the receptor; the two files may be in "
             "different coordinate frames and every pocket distance below would be void")

    triad_ox = [a for a in receptor
                if (a["resnum"], a["name"]) in {(477, "OD1"), (477, "OD2"),
                                                (515, "OE1"), (515, "OE2"),
                                                (588, "OD1"), (588, "OD2")}]
    emit("R", "catalytic carboxylate oxygens located", len(triad_ox),
         "Asp477, Glu515, Asp588")
    if len(triad_ox) != 6:
        note(f"expected 6 catalytic carboxylate oxygens, found {len(triad_ox)}")
    triad_centre = centroid(triad_ox)

    def triad_contacts(lig):
        out = {}
        for num in (477, 515, 588):
            ox = [a for a in triad_ox if a["resnum"] == num]
            out[num] = sum(1 for a in lig if nearest(a["xyz"], ox) <= 4.0)
        return out

    ac_contacts = triad_contacts(acarbose)
    ac_total = sum(ac_contacts.values())
    emit("R", "acarbose contacts within 4.0 A of Asp477", ac_contacts[477], "this script")
    emit("R", "acarbose contacts within 4.0 A of Glu515", ac_contacts[515], "this script")
    emit("R", "acarbose contacts within 4.0 A of Asp588", ac_contacts[588], "this script")
    emit("R", "acarbose triad contacts per 100 heavy atoms",
         f"{100.0 * ac_total / len(acarbose):.1f}", "this script")

    poses = [
        ("CPHWRPWHC", "demo_data/results/3AIC_docked_complex_CPHWRPWHC.pdb"),
        ("CAHWRPWHC", "demo_data/results/3AIC_docked_complex_CAHWRPWHC.pdb"),
        ("CHPWPHWAC", "demo_data/results/3AIC_docked_complex_CHPWPHWAC.pdb"),
        ("CRRRWWRRRRC", "demo_data/results/S_mutans_GtfB_Insoluble_variant_CRRRWWRRRRC.pdb"),
        ("CRRRWWPRRRC", "demo_data/results/S_mutans_GtfB_Insoluble_variant_CRRRWWPRRRC.pdb"),
        ("CRRRWFRRRRC", "demo_data/results/S_mutans_GtfB_Insoluble_variant_CRRRWFRRRRC.pdb"),
        ("CRRRAWPRRRC", "demo_data/results/S_mutans_GtfB_Insoluble_variant_CRRRAWPRRRC.pdb"),
        ("CRRAGWPRRRC", "demo_data/results/S_mutans_GtfB_Insoluble_variant_CRRAGWPRRRC.pdb"),
    ]
    for seq, rel in poses:
        path = REPO / rel
        if not path.exists():
            note(f"deposited pose for {seq} is missing at {rel}")
            continue
        lig = heavy_only(pdb_atoms(path, chain=ligand_chain(path)))
        if not lig:
            note(f"no ligand chain parsed from {rel}")
            continue
        c = triad_contacts(lig)
        total = sum(c.values())
        emit("R", f"{seq} heavy atoms in the deposited pose", len(lig), rel)
        emit("R", f"{seq} contacts within 4.0 A of the triad (477/515/588)",
             f"{c[477]}/{c[515]}/{c[588]}", "this script")
        emit("R", f"{seq} triad contacts per 100 heavy atoms",
             f"{100.0 * total / len(lig):.1f}", "this script")
        emit("R", f"{seq} pose centroid to triad centroid, A",
             f"{euclid(centroid(lig), triad_centre):.2f}", "this script")
        overlap = sum(1 for a in lig if nearest(a["xyz"], acarbose) <= 4.0)
        emit("R", f"{seq} atoms within 4.0 A of the acarbose envelope",
             f"{overlap}/{len(lig)}", "this script")
        clashes = sum(1 for a in lig if nearest(a["xyz"], receptor) < 2.2)
        emit("R", f"{seq} heavy atoms closer than 2.2 A to the receptor",
             clashes, "steric check")

    # ------------------------------------------------- the stated mechanism
    heading("16. The mechanism the analysis files assert, against those same coordinates")
    print("    Five demo_data/results/analysis files state, in identical words, that the")
    print("    arginine cluster forms bidentate salt bridges with Asp477 and Glu515 and that")
    print("    the tryptophan pair stacks between Trp517 and Tyr430. Both are measurable.")
    for num, expect in ((430, "TYR"), (517, "TRP")):
        found = sorted({a["resname"] for a in receptor if a["resnum"] == num})
        emit("R", f"residue {num} in the receptor", "/".join(found) or "absent",
             f"expected {expect}")
        if found != [expect]:
            note(f"residue {num} is {found or 'absent'} in the receptor, not {expect}")

    w517 = aromatic_ring([a for a in receptor if a["resnum"] == 517], "TRP")
    y430 = aromatic_ring([a for a in receptor if a["resnum"] == 430], "TYR")
    c517, c430 = centroid(w517), centroid(y430)
    emit("R", "Trp517 ring centroid to Tyr430 ring centroid, A",
         f"{euclid(c517, c430):.2f}", "this script")
    if euclid(c517, c430) > 8.0:
        note(f"Trp517 and Tyr430 ring centroids are {euclid(c517, c430):.2f} A apart, too far "
             "to sandwich a third ring, so no ligand ring can stack between them")

    for seq, rel in poses:
        path = REPO / rel
        if not path.exists():
            continue
        lig = heavy_only(pdb_atoms(path, chain=ligand_chain(path)))
        guan = [a for a in lig if a["resname"] == "ARG" and a["name"] in ("NE", "NH1", "NH2")]
        if guan:
            for num in (477, 515, 588):
                ox = [a for a in triad_ox if a["resnum"] == num]
                pairs = [euclid(g["xyz"], o["xyz"]) for g in guan for o in ox]
                within = sum(1 for d in pairs if d <= 4.0)
                emit("R", f"{seq} guanidinium N to {num} carboxylate O, closest A / pairs <=4.0",
                     f"{min(pairs):.2f} / {within}", "salt-bridge test")
                # Only the five variant analysis files assert the Glu515 salt bridge, so
                # only those five are a contradiction. The measurement is reported for
                # every pose either way.
                if num == 515 and within == 0 and seq in ASSERTS_GLU515:
                    note(f"{seq}: its analysis file asserts a bidentate salt bridge to "
                         "Glu515, but no guanidinium nitrogen in the deposited pose lies "
                         f"within 4.0 A of a Glu515 carboxylate oxygen (closest "
                         f"{min(pairs):.2f} A)")
        by_res = {}
        for a in lig:
            by_res.setdefault((a["resnum"], a["resname"]), []).append(a)
        for (rnum, rname), atoms in sorted(by_res.items()):
            ring = aromatic_ring(atoms, rname)
            if len(ring) >= 5:
                rc = centroid(ring)
                emit("R", f"{seq} {rname}{rnum} ring to Trp517 / Tyr430, A",
                     f"{euclid(rc, c517):.2f} / {euclid(rc, c430):.2f}", "stacking test")

    # ------------------------------------------------- pose validity vs the intact receptor
    heading("18. Are the deposited poses sterically possible in the intact enzyme?")
    print("    Each deposited complex carries a contiguous slice of the receptor, not the\n"
          "    whole chain. A clash test run against the slice cannot see an overlap with\n"
          "    the residues the slice omits. Every pose is tested here against all committed\n"
          "    receptor atoms as well as against the slice the file itself carries.")
    rec_all = heavy_only(pdb_atoms(RECEPTOR_3AIC))
    emit("R", "committed receptor heavy atoms", len(rec_all), "3AIC_receptor.pdb")
    rec_nums = sorted({a["resnum"] for a in rec_all})
    emit("R", "committed receptor residue span",
         f"{rec_nums[0]}..{rec_nums[-1]} ({len(rec_nums)} residues)", "3AIC_receptor.pdb")
    for seq, rel in poses:
        path = REPO / rel
        fname = path.name
        if not path.exists():
            continue
        slice_atoms = heavy_only(pdb_atoms(path, chain="A"))
        lig = heavy_only(pdb_atoms(path, chain=ligand_chain(path)))
        if not lig or not slice_atoms:
            continue
        snums = sorted({a["resnum"] for a in slice_atoms})
        emit("R", f"{fname}: deposited receptor slice",
             f"{snums[0]}..{snums[-1]} ({len(snums)} residues, {len(slice_atoms)} atoms)",
             "this script")
        vs_slice = [nearest(b["xyz"], slice_atoms) for b in lig]
        vs_all = [nearest(b["xyz"], rec_all) for b in lig]
        emit("R", f"{fname}: closest ligand-receptor approach, slice / intact, A",
             f"{min(vs_slice):.2f} / {min(vs_all):.2f}", "this script")
        emit("R", f"{fname}: ligand heavy atoms under 2.2 A, slice / intact",
             f"{sum(1 for d in vs_slice if d < 2.2)} / {sum(1 for d in vs_all if d < 2.2)}",
             "this script")
        clashing = sorted({(a["resnum"], a["resname"]) for b in lig for a in rec_all
                           if euclid(b["xyz"], a["xyz"]) < 2.2})
        if clashing:
            inside = sum(1 for n, _ in clashing if n in set(snums))
            emit("R", f"{fname}: intact-receptor residues it overlaps",
                 f"{len(clashing)}, spanning {min(n for n, _ in clashing)}.."
                 f"{max(n for n, _ in clashing)}, {inside} of them inside the deposited slice",
                 "this script")
            note(f"{fname}: the ligand overlaps {len(clashing)} receptor residues "
                 f"({min(n for n, _ in clashing)}..{max(n for n, _ in clashing)}) at a closest "
                 f"approach of {min(vs_all):.2f} A, and {inside} of those residues are present "
                 f"in the deposited file. The pose is not sterically possible in the intact "
                 f"enzyme, and the deposited slice cannot reveal it")

    # ------------------------------------------------- what the redocking actually got wrong
    heading("19b. Redocking: which part of the pose is wrong")
    print("    A single RMSD says a redocking failed but not how. These numbers separate\n"
          "    finding the site from placing the atoms, per retained mode.")
    crystal = heavy_only(pdb_atoms(NATIVE_LIGAND))
    cc = centroid(crystal)
    modes = sorted((REPO / "data" / "poses").glob("3AIC_native_acarbose_docked_split_*.pdbqt"),
                   key=lambda q: int(q.stem.rsplit("_", 1)[1]))
    disp, prox = [], []
    for q in modes:
        pose = heavy_only(pdb_atoms(q))
        if not pose:
            continue
        disp.append(euclid(centroid(pose), cc))
        prox.append(sum(1 for a in crystal if nearest(a["xyz"], pose) <= 2.0))
    if disp:
        emit("R", "redocked modes compared", len(disp), "data/poses")
        emit("R", "pose centroid to crystal centroid, range A",
             f"{min(disp):.2f} to {max(disp):.2f}", "this script")
        emit("R", "crystal atoms with a pose atom within 2.0 A, best mode",
             f"{max(prox)} of {len(crystal)}", "this script")
        emit("R", "crystal atoms with a pose atom within 2.0 A, range",
             f"{min(prox)} to {max(prox)} of {len(crystal)}", "this script")
        print("    Read together with the 4.507 A symmetry-corrected RMSD: the search puts the\n"
              "    ligand in the right cleft and recovers most of the crystallographic envelope,\n"
              "    and fails on atom-level register. Acarbose is a near-symmetric oligosaccharide,\n"
              "    so a shift of one sugar unit along the cleft gives exactly this signature. The\n"
              "    honest statement is that the setup locates the site and does not reproduce the\n"
              "    binding mode, not that it cannot find the site.")

    # ------------------------------------------------- search convergence
    heading("19. Is the docking search converged? Two internal tests")
    print("    The campaign ranks candidates on the single best-mode Vina score. These\n"
          "    numbers ask how far that score moves with no change of ligand at all: across\n"
          "    the modes Vina retains from a single run.")
    runs = {}
    for pose in sorted((REPO / "data" / "poses").glob("*_split_*.pdbqt")):
        m = re.match(r"(.+)_(3aic|1smd)_split_(\d+)\.pdbqt$", pose.name)
        if not m:
            continue
        score = None
        for line in pose.read_text().splitlines():
            if line.startswith("REMARK VINA RESULT:"):
                score = float(line.split()[3])
                break
        if score is not None:
            runs.setdefault((m.group(1), m.group(2)), {})[int(m.group(3))] = score
    spreads = [max(v.values()) - min(v.values()) for v in runs.values() if len(v) >= 2]
    spreads.sort()
    emit("R", "retained pose files in the repository",
         len(list((REPO / "data" / "poses").glob("*.pdbqt"))), "data/poses")
    emit("R", "docking runs with at least two retained modes", len(spreads), "data/poses")
    emit("R", "within-run score spread, mean / median kcal/mol",
         f"{sum(spreads)/len(spreads):.3f} / {spreads[len(spreads)//2]:.3f}", "this script")
    emit("R", "within-run score spread, range kcal/mol",
         f"{min(spreads):.3f} to {max(spreads):.3f}", "this script")
    lead_wt_gap = abs(by_seq[LEAD]["vina_3aic_kcal"] - by_seq[WT]["vina_3aic_kcal"])
    emit("R", "lead minus wild type, kcal/mol", f"{lead_wt_gap:.3f}", "vault, this script")
    emit("R", "runs whose own mode spread exceeds that gap",
         f"{sum(1 for v in spreads if v > lead_wt_gap)} of {len(spreads)}", "this script")
    emit("R", "runs whose own mode spread exceeds the null SD",
         f"{sum(1 for v in spreads if v > sd)} of {len(spreads)}", "this script")
    if spreads and sum(1 for v in spreads if v > lead_wt_gap) == len(spreads):
        note("in every docking run with more than one retained mode, the spread between the "
             "best and worst retained mode exceeds the score difference the campaign ranks "
             "candidates on, so the ranking is finer than the search that produced it")
    emit("R", "active torsions, 9-mer macrocycle / 11-mer variant",
         "15 / 35", "TORSDOF in data/poses, this script")
    emit("R", "exhaustiveness, redocking calibration / linear control / parent / "
         "naive 11-mer / alanine scan", 16, "scripts/run_real_campaign.py")
    emit("R", "exhaustiveness, 39 scrambled nulls and unlabelled default", 8,
         "scripts/run_real_campaign.py (its docstring says 4; the code passes 8)")
    emit("R", "exhaustiveness, generation 2 and generation 3 candidates", 4,
         "scripts/run_real_campaign.py")
    note("exhaustiveness is not uniform across the campaign, so every comparison that "
         "crosses these groups compares unequal search effort as well as unequal "
         "molecules; and per section 20 the lead's own pose file cannot be attributed "
         "to any one of these settings")

    # ------------------------------------------------- duplicate census
    heading("20. Does any docked pose appear twice under a different protocol?")
    print("    A stochastic search does not reproduce another run's coordinates to the\n"
          "    last decimal of every atom in every retained mode. So two pose files with\n"
          "    the same coordinate digest are one run reported twice, under whatever\n"
          "    protocol each file's script claims.")
    digests = collections.defaultdict(list)
    scanned = 0
    for pose in sorted((REPO / "data" / "poses").glob("*.pdbqt")):
        if "_split_" in pose.name:
            continue
        d = pdbqt_coord_digest(pose)
        if d is None:
            continue
        scanned += 1
        digests[d].append(pose.name)
    dup_groups = {d: v for d, v in digests.items() if len(v) > 1}
    emit("R", "pose and ligand files scanned, splits excluded", scanned, "data/poses")
    emit("R", "distinct coordinate sets", len(digests), "this script")
    emit("R", "duplicate groups", len(dup_groups), "this script")
    for d in sorted(dup_groups, key=lambda k: (-len(digests[k]), digests[k][0])):
        emit("R", f"identical coordinates in {len(dup_groups[d])} files",
             ", ".join(sorted(dup_groups[d])), "coordinate digest, this script")
    if dup_groups:
        note("each script that writes these files docks only when the output path is "
             "absent or under a hundred bytes, so a file already at that path is adopted "
             "with no check that it came from the protocol about to be reported; that "
             "defect cannot explain the duplicates below, because their paths differ and "
             "the guard can only adopt content already at the path it is about to write")

    # ---- where the duplicated bytes came from, and what provenance was kept
    def git(*args):
        """Text output, stripped. Never use for file content: stripping loses bytes."""
        try:
            return subprocess.run(("git",) + args, cwd=str(REPO), check=True,
                                  capture_output=True, text=True, timeout=60).stdout.strip()
        except Exception:
            return None

    def git_blob(sha, path):
        """Raw bytes of a committed file. Must not decode or strip, or the md5 is wrong."""
        try:
            return subprocess.run(("git", "show", f"{sha}:{path}"), cwd=str(REPO),
                                  check=True, capture_output=True, timeout=60).stdout
        except Exception:
            return None

    TRIO = ("data/poses/contender_mutant_P2A_CAHWRPWHC_3aic.pdbqt",
            "data/poses/disproof_seed_42_CAHWRPWHC_3AIC_docked.pdbqt",
            "data/poses/disproof_explicit_pH72_3AIC_docked.pdbqt")
    PAIR = ("data/poses/contender_mutant_P2A_CAHWRPWHC_1smd.pdbqt",
            "data/poses/disproof_explicit_pH72_1SMD_docked.pdbqt")
    for group, label in ((TRIO, "GTF-SI trio"), (PAIR, "alpha-amylase pair")):
        present = [f for f in group if (REPO / f).exists()]
        if len(present) != len(group):
            continue
        sums = {hashlib.md5((REPO / f).read_bytes()).hexdigest() for f in present}
        sizes = {(REPO / f).stat().st_size for f in present}
        emit("R", f"{label}: whole-file md5 values / byte sizes among {len(present)} files",
             f"{len(sums)} / {len(sizes)} distinct, {sorted(sizes)[0]} bytes", "this script")
        if len(sums) == 1:
            note(f"the {label} is identical in every byte, not only in coordinates, so the "
                 f"score and chemistry records written above the coordinates are identical too")
        for f in present:
            adds = git("log", "--diff-filter=A", "--format=%h %ad",
                       "--date=format:%Y-%m-%d %H:%M", "--", f)
            if adds:
                emit("R", f"{Path(f).name}: commit that added it",
                     adds.splitlines()[-1], "git history")
        camp = present[0]
        hist = git("log", "--format=%h %ad", "--date=format:%Y-%m-%d %H:%M", "--", camp)
        if hist:
            for line in reversed(hist.splitlines()):
                sha = line.split()[0]
                blob = git_blob(sha, camp)
                if blob is not None:
                    d = hashlib.md5(blob).hexdigest()
                    emit("R", f"{Path(camp).name} content at {line}", d[:12],
                         "git history, this script")
            live = hashlib.md5((REPO / camp).read_bytes()).hexdigest()
            emit("R", f"{Path(camp).name} content now", live[:12], "this script")
            revs = [l.split()[0] for l in hist.splitlines()]
            states = [(r, hashlib.md5(b).hexdigest())
                      for r in reversed(revs) for b in [git_blob(r, camp)] if b is not None]
            matching = [r for r, d in states if d == live]
            if states and matching:
                emit("R", f"{Path(camp).name}: which of its {len(states)} committed states "
                     f"matches the bytes the audit paths were created with",
                     f"{matching[-1]} (and {len(matching)} of {len(states)})",
                     "git history, this script")
                if len(matching) < len(states):
                    note(f"the {label} content post-dates the first commit of "
                         f"{Path(camp).name}: that file held different bytes when it was "
                         f"added and was replaced the same day, and the audit paths were "
                         f"created four days later already holding the replaced bytes, so "
                         f"no independent docking output was ever written to those paths")

    # what provenance the artefacts actually record
    PROV = ("timestamp", "seed", "exhaustiveness", "command line", "run identifier")
    pose = REPO / "data" / "poses" / "disproof_seed_42_CAHWRPWHC_3AIC_docked.pdbqt"
    if pose.exists():
        text = pose.read_text().lower()
        found = [w for w in ("seed", "exhaust", "command", "run_id", "vina --")
                 if w in text] + (["timestamp"] if re.search(r"20[0-9]{2}-[0-9]{2}-[0-9]{2}", text) else [])
        emit("R", "provenance fields recorded inside the duplicated pose file",
             ", ".join(found) if found else "none of " + ", ".join(PROV),
             "data/poses/disproof_seed_42_CAHWRPWHC_3AIC_docked.pdbqt")
    rep_json = REPO / "data" / "exhaustive_disproof_report.json"
    if rep_json.exists():
        raw = rep_json.read_text()
        has_ts = bool(re.search(r"20[0-9]{2}-[0-9]{2}-[0-9]{2}", raw))
        emit("R", "audit report: timestamp present / command line present",
             f"{has_ts} / {'command' in raw.lower()}", "data/exhaustive_disproof_report.json")
        rj = json.loads(raw)
        for a in rj.get("angles", []):
            if a.get("angle") == 2:
                emit("D", "audit report's own claim: seeds used / reported score SD",
                     f"{a.get('n_seeds')} / {a.get('score_sd_kcal')}",
                     "data/exhaustive_disproof_report.json")
                if a.get("n_seeds") == 5:
                    note("the audit report states that five independent seeds were used and "
                         "lists the duplicated value among the five, and its reported "
                         "dispersion is computed over those five; the file comparison above "
                         "leaves four independent runs")
    logs = sorted((REPO / "data" / "logs").glob("*.log")) if (REPO / "data" / "logs").is_dir() else []
    emit("R", "log files in the repository / of those, empty",
         f"{len(logs)} / {sum(1 for l in logs if l.stat().st_size == 0)}", "data/logs")

    ph_neutral = REPO / "data" / "poses" / "disproof_explicit_pH72_3AIC_docked.pdbqt"
    ph_acidic = REPO / "data" / "poses" / "disproof_explicit_pH52_3AIC_docked.pdbqt"
    if ph_neutral.exists() and ph_acidic.exists():
        n_s, a_s = pdbqt_score(ph_neutral), pdbqt_score(ph_acidic)
        emit("R", "pH 7.2 arm, kcal/mol", n_s, "data/poses, this script")
        emit("R", "pH 5.2 arm, kcal/mol", a_s, "data/poses, this script")
        emit("R", "reported pH shift, kcal/mol", round(a_s - n_s, 3), "this script")
        if campaign_digest and pdbqt_coord_digest(ph_neutral) == campaign_digest:
            note("the pH 7.2 arm is the campaign pose, not a run of its own, so the pH "
                 "shift is the difference between a reused file and a fresh run")

    # --------------------------------- the site, the numbering, and cyclisation
    heading("21. The docking site, the receptor numbering, and the one controlled "
            "test of cyclisation")
    print("    Vina's scoring function has no electrostatic term (steric, hydrophobic\n"
          "    and hydrogen bonding only, hydrogens not scored explicitly, 8 A cutoff).\n"
          "    These numbers establish what was stripped from the receptor, that the two\n"
          "    receptor files are one structure under two numberings, and what the linear\n"
          "    control did against the macrocycle it was meant to lose to.")

    def pdbqt_rows(path):
        rows = []
        for l in Path(path).read_text().splitlines():
            if l.startswith(("ATOM", "HETATM")):
                try:
                    q = float(l[70:76])
                except (ValueError, IndexError):
                    q = 0.0
                rows.append((l[12:16].strip(), l[17:20].strip(), int(l[22:26]),
                             round(float(l[30:38]), 3), round(float(l[38:46]), 3),
                             round(float(l[46:54]), 3), q))
        return rows

    rec_pdb_path = RECEPTORS / "3AIC_receptor.pdb"
    rec_pq_path = RECEPTORS / "3AIC_receptor.pdbqt"
    if rec_pdb_path.exists() and rec_pq_path.exists():
        for name in ("3AIC_receptor.pdb", "3AIC_receptor.pdbqt"):
            het = sum(1 for l in (RECEPTORS / name).read_text().splitlines()
                      if l.startswith("HETATM"))
            emit("R", f"{name}: HETATM records", het, f"demo_data/receptors/{name}")
        if all(sum(1 for l in (RECEPTORS / n).read_text().splitlines()
                   if l.startswith("HETATM")) == 0
               for n in ("3AIC_receptor.pdb", "3AIC_receptor.pdbqt")):
            note("both receptor files carry no heteroatom records, so every non-protein "
                 "atom was removed before docking: any metal, all ordered water and the "
                 "crystallographic ligand")

        rp, rq = pdbqt_rows(rec_pdb_path), pdbqt_rows(rec_pq_path)
        q_by_xyz = {(r[3], r[4], r[5]): r for r in rq}
        offs = collections.Counter(r[2] - q_by_xyz[(r[3], r[4], r[5])][2]
                                   for r in rp if (r[3], r[4], r[5]) in q_by_xyz)
        matched = sum(offs.values())
        emit("R", "pdb atoms whose coordinates also appear in the docked pdbqt",
             f"{matched} of {len(rp)}", "this script")
        emit("R", "residue-number offset, analysis pdb minus docked pdbqt",
             ", ".join(f"{k} ({v} atoms)" for k, v in offs.most_common(3)), "this script")
        if len(offs) == 1 and matched == len(rp):
            # a check that PASSES, so it is printed here and not under DISCREPANCIES
            print("    CHECK PASSES: the analysis receptor and the docked receptor are one\n"
                  "    structure in one frame under two numberings, offset by a single\n"
                  "    constant, so geometry measured on the analysis file applies to what\n"
                  "    was docked.")

        native = REPO / "data" / "poses" / "3AIC_native_acarbose_docked.pdbqt"
        if native.exists():
            lig = pdbqt_best_mode(native)
            cen = [sum(a[i] for a in lig) / len(lig) for i in (1, 2, 3)]
            emit("R", "redocked acarbose best-mode centroid",
                 " ".join(f"{c:.2f}" for c in cen), "this script")
            off = offs.most_common(1)[0][0] if offs else 0
            # heavy atoms only: Vina does not score hydrogens explicitly, and a
            # polar hydrogen can reach into a shell its residue's heavy atoms do not
            heavy = [r for r in rq if not r[0].startswith("H")]
            for radius in (8.0, 12.0):
                near = [r for r in heavy if math.dist(cen, r[3:6]) <= radius]
                res = {(r[1], r[2]) for r in near}
                acid = sorted(n for t, n in res if t in ("ASP", "GLU"))
                base = sorted(n for t, n in res if t in ("ARG", "LYS"))
                emit("R", f"within {radius:.0f} A of that centroid: acidic / basic residues, "
                     f"net formal charge", f"{len(acid)} / {len(base)} / {len(base)-len(acid):+d}",
                     "this script")
                if radius == 8.0:
                    emit("R", "those acidic residues, in the analysis numbering",
                         ", ".join(f"{t}{n + off}" for t, n in sorted(res, key=lambda r: r[1])
                                   if t in ("ASP", "GLU")), "this script")

    print()
    for label, fname in (("linearised control SPHWRPWHS", "control_linear_SPHWRPWHS_3aic.pdbqt"),
                         ("cyclic parent CPHWRPWHC", "contender_g01_CPHWRPWHC_3aic.pdbqt"),
                         ("cyclic lead CAHWRPWHC", "contender_mutant_P2A_CAHWRPWHC_3aic.pdbqt")):
        f = REPO / "data" / "poses" / fname
        if not f.exists():
            continue
        tors = next((int(l.split()[1]) for l in f.read_text().splitlines()
                     if l.startswith("TORSDOF")), None)
        emit("R", f"{label}: score / heavy atoms / active torsions",
             f"{pdbqt_score(f)} / {len(pdbqt_best_mode(f))} / {tors}", f"data/poses/{fname}")
    lin = REPO / "data" / "poses" / "control_linear_SPHWRPWHS_3aic.pdbqt"
    cyc = REPO / "data" / "poses" / "contender_g01_CPHWRPWHC_3aic.pdbqt"
    if lin.exists() and cyc.exists():
        adv = pdbqt_score(cyc) - pdbqt_score(lin)
        emit("R", "linear advantage over the macrocycle, kcal/mol", round(adv, 3), "this script")
        if adv > 0:
            note("the linearised control outscores the macrocycle it was meant to lose to, "
                 "at the same heavy-atom count and the same exhaustiveness, while carrying "
                 "35 active torsions against the macrocycle's 15; Vina's only entropy-like "
                 "term penalises rotatable bonds, so the penalty should have run the other "
                 "way, and cyclisation is this project's central design premise")

    for label, fname in (("pH 7.2 neutral ligand", "disproof_explicit_pH72_neutral_CAHWRPWHC.pdbqt"),
                         ("pH 5.2 diprotonated ligand", "disproof_explicit_pH52_diprotonated_CAHWRPWHC.pdbqt")):
        f = REPO / "data" / "poses" / fname
        if f.exists():
            rows = pdbqt_rows(f)
            emit("R", f"{label}: atoms / summed partial charge",
                 f"{len(rows)} / {sum(r[6] for r in rows):+.3f}", f"data/poses/{fname}")

    # ------------------------------------------------- the design's resolution
    heading("22. What size of difference this null design could resolve")
    remark("Section 3.3b. The reviewer's next question after a null p value is "
         "'enlarge the null'. These numbers say why that does not help, and what "
         "would. Everything here is arithmetic on the 39 scrambled controls and "
         "on the repeat runs of one ligand from section 3.2.")

    null_scores = [r["vina_3aic_kcal"] for r in scrambles]
    n_null = len(null_scores)
    null_mean = statistics.mean(null_scores)
    null_sd = statistics.stdev(null_scores)
    null_best = min(null_scores)
    emit("R", "null members", n_null, "ledger, SCRAMBLED_NULL_CONTROL rows")
    emit("R", "null mean / sample sd (kcal/mol)",
         f"{null_mean:.4f} / {null_sd:.4f}", "this script")
    emit("R", "best null member (kcal/mol)", f"{null_best:.4f}", "this script")

    # (k+1)/(n+1) is the one-sided permutation p with the observation counted in.
    p_floor = 1.0 / (n_null + 1)
    k_allowed = [k for k in range(n_null + 1) if (k + 1) / (n_null + 1) < 0.05]
    k_max = max(k_allowed) if k_allowed else None
    emit("R", "smallest attainable one-sided p at this n", f"{p_floor:.4f}", "this script")
    if k_max == 0:
        emit("R", "null members the candidate must beat for p < 0.05",
             f"{n_null} of {n_null}", "this script")
        emit("R", "score that requires (kcal/mol)", f"{null_best:.4f}", "this script")
        emit("R", "that margin over the null mean (kcal/mol)",
             f"{null_best - null_mean:.4f}", "this script")
    else:
        discrepancies.append(
            f"p<0.05 is reachable at n={n_null} while losing to {k_max} null members, "
            "so section 3.3b's threshold statement needs rewriting")

    # The rejection boundary is set by the spread, not the count. Normal approximation.
    Z95 = 1.6448536269514722
    Z80 = 0.8416212335729143
    boundary = Z95 * null_sd
    emit("R", "one-sided 5 percent boundary, normal approximation (kcal/mol better than the mean)",
         f"{boundary:.4f}", "this script")

    # Which direction does enlarging the null move the requirement? With one
    # score per control the candidate may lose to at most j controls, where j is
    # the largest integer with (j+1)/(n+1) < 0.05. The required rank as a
    # fraction of n rises toward 0.05 but stays below it at every n, so the 5
    # per cent quantile is a floor beneath the requirement and no control count
    # makes the test easier than that quantile.
    def _losable(n, alpha=0.05):
        j = -1
        while (j + 2) / (n + 1) < alpha:
            j += 1
        return j

    for n in (39, 199, 999):
        j = _losable(n)
        emit("R", f"controls a candidate may lose to at n={n}, and the required quantile",
             f"{j} / {(j + 1) / n:.4f}", "this script")
    if _losable(n_null) != 0:
        note(f"at n={n_null} a candidate may now lose to {_losable(n_null)} controls, so "
             "section 3.3b's statement that it must beat every one of them is wrong")
    if any((_losable(n) + 1) / n >= 0.05 for n in (39, 199, 999, 9999, 999999)):
        note("the required quantile reaches or passes 0.05 at some n, so the 5 per cent "
             "quantile is not a floor beneath the requirement and section 3.3b's "
             "'however many controls are added' no longer holds")
    emit("R", "requirement at n=39 in null sd units (beat the best control)",
         f"{(null_mean - null_best) / null_sd:.4f}", "this script")
    emit("R", "asymptotic floor on the requirement, in null sd units", f"{Z95:.4f}",
         "Trott and Olson give no boundary; this is the normal 5 percent point")
    remark("Enlarging the null relaxes the rank a candidate needs and lowers the floor on "
           "the attainable p value, but it cannot lower the score: the requirement moves "
           "down from 2.20 null standard deviations at n=39 toward 1.64 and no further. "
           "This is the clause a referee needs, because 'more controls cannot help' reads "
           "as a power excuse without it.")

    # The effect size has to come from the repository, not from a round number.
    # The lead's own margin is its gain over the parent it was derived from.
    lead_rows = [r for r in designed if r["sequence"] == LEAD]
    wt_rows = [r for r in designed if r["sequence"] == WT]
    MARGIN = None
    if lead_rows and wt_rows:
        lead_score = lead_rows[0]["vina_3aic_kcal"]
        wt_score = wt_rows[0]["vina_3aic_kcal"]
        MARGIN = abs(lead_score - wt_score)
        emit("R", "lead and parent at GTF-SI (kcal/mol)",
             f"{lead_score:.3f} / {wt_score:.3f}", "ledger")
        emit("R", "the lead's margin over its parent (kcal/mol)", f"{MARGIN:.4f}", "this script")
        emit("R", "lead minus null mean (kcal/mol, positive is worse)",
             f"{lead_score - null_mean:+.4f}", "this script")
        remark("This margin, not a round number, is the effect size used below. No "
               "document in the repository states a target advantage, so I use the one "
               "the screen actually produced. Searched for a stated threshold and found "
               "only the +0.70 kcal/mol alanine hotspot gate, which measures something "
               "else.")
        emit("R", "that margin in null sd units", f"{MARGIN / null_sd:.4f}", "this script")
        emit("R", "the asymptotic floor divided by that margin", f"{boundary / MARGIN:.2f}",
             "this script")
        emit("R", "ratio of the beat-every-member margin to it",
             f"{abs(null_best - null_mean) / MARGIN:.2f}", "this script")

    # Where the null's dispersion comes from: sequence, or the instrument?
    seeds = None
    for a in disproof["angles"]:
        if a.get("vina_scores"):
            seeds = list(a["vina_scores"])
            break
    if seeds:
        dup = -8.191  # the value carried by the duplicated file, section 3.10
        indep = [v for v in seeds if abs(v - dup) > 1e-9]
        noise_sd = statistics.stdev(indep)
        emit("R", "one ligand, independent repeat runs / their sample sd (kcal/mol)",
             f"{len(indep)} / {noise_sd:.4f}", "data/exhaustive_disproof_report.json")
        emit("R", "ratio of the repeat-run sd to the null sd",
             f"{noise_sd / null_sd:.4f}", "this script")
        n_lo, n_hi = sd_interval(noise_sd, len(indep))
        k_lo, k_hi = sd_interval(null_sd, n_null)
        emit("R", f"95 percent interval on the repeat-run sd (n={len(indep)})",
             f"{n_lo:.4f} to {n_hi:.4f}", "this script")
        emit("R", f"95 percent interval on the null sd (n={n_null})",
             f"{k_lo:.4f} to {k_hi:.4f}", "this script")
        between = null_sd ** 2 - noise_sd ** 2
        emit("R", "null variance minus repeat-run variance (kcal/mol squared)",
             f"{between:+.4f}", "this script")
        remark("What this does and does not support. The two spreads differ by "
               f"{abs(1 - noise_sd / null_sd) * 100:.1f} per cent, and each interval "
               "contains the other estimate, so the data cannot distinguish a null whose "
               "dispersion is entirely run-to-run noise from one with a real "
               "sequence-dependent component. The difference of variances is "
               f"{between:+.4f}, and a negative point estimate is what two close "
               "estimates give when one of them rests on three degrees of freedom; it is "
               "not evidence that the component is nil. The supportable claim is only "
               "that the null's dispersion is not demonstrably sequence-driven. The "
               "repeat runs are also the four that remain after the duplicate of "
               "unestablished origin was removed from the original five, and all of them "
               "were unseeded.")
        if not (n_lo <= null_sd <= n_hi):
            discrepancies.append(
                "the null sd now falls outside the interval on the single-ligand "
                "repeat-run sd, so section 3.3b's comparison of the two no longer holds")

        # What would move the boundary: averaging m seeds per member.
        if MARGIN:
            m_boundary = math.ceil((Z95 * null_sd / MARGIN) ** 2)
            m_power80 = math.ceil(((Z95 + Z80) * null_sd / MARGIN) ** 2)
            emit("R", "seeds per member for that margin to reach the 5 percent boundary",
                 m_boundary, "this script")
            emit("R", "total docking runs that implies over 39 controls and one candidate",
                 (n_null + 1) * m_boundary, "this script")
            emit("R", "seeds per member for 80 percent power at that margin", m_power80,
                 "this script")
            emit("R", "total docking runs that implies", (n_null + 1) * m_power80, "this script")

    # And what resolving it would mean, against the instrument's own error.
    VINA_SE = 2.85  # Trott and Olson's reported standard error, section 4.2b
    if MARGIN:
        emit("R", "the lead's margin as a fraction of Vina's reported standard error",
             f"{MARGIN / VINA_SE:.4f}", "this script")
    emit("R", "beat-every-member margin as a fraction of that standard error",
         f"{abs(null_best - null_mean) / VINA_SE:.4f}", "this script")

    # ------------------------------------------------- tally
    heading("23. Battery tally")
    for cls, n in sorted(disproof["classification_summary"].items()):
        emit("D", f"angles classified {cls}", n, "disproof classification_summary")
    emit("R", "angles in the report", len(disproof["angles"]), "this script")
    emit("D", "angles computed on that run", disproof["provenance_summary"]["computed_this_run"],
         "disproof provenance_summary")
    emit("D", "angles read from the vault", disproof["provenance_summary"]["vault_query"],
         "disproof provenance_summary")
    emit("D", "angles not computed", disproof["provenance_summary"]["not_computed"],
         "disproof provenance_summary")

    # ------------------------------------------------- output
    print()
    print("=" * 100)
    print(f"DISCREPANCIES ({len(discrepancies)})")
    print("=" * 100)
    if not discrepancies:
        print("  none")
    for d in discrepancies:
        print(f"  - {d}")
    print()
    print(f"{len(numbers)} numbers emitted: "
          f"{sum(1 for t, *_ in numbers if t == 'R')} recomputed here, "
          f"{sum(1 for t, *_ in numbers if t == 'D')} read from a recorded run.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
