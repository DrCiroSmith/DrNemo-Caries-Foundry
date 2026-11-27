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

The script also cross-checks the artefacts against each other and prints any
disagreement under DISCREPANCIES. A disagreement is not a rounding note: it
means two committed files state different values for the same quantity, and the
manuscript has to say so rather than pick one.
"""
import json
import math
import sqlite3
import statistics
import sys
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
    discrepancies.append(text)


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
    heading("5. Angle 2, multi-seed reproducibility")
    a2 = angles[2]
    emit("D", "independent seeds", a2["n_seeds"], "disproof angle 2")
    emit("R", "seed score mean kcal/mol", round(statistics.fmean(a2["vina_scores"]), 3),
         "mean of disproof angle 2 vina_scores, this script")
    emit("R", "seed score SD kcal/mol", round(statistics.stdev(a2["vina_scores"]), 3),
         "sample SD of the same five scores, this script")
    emit("R", "seed score range kcal/mol",
         f"{min(a2['vina_scores'])} to {max(a2['vina_scores'])}", "this script")
    emit("D", "mean pairwise pose RMSD, A", a2["mean_pairwise_rmsd_angstrom"], "disproof angle 2")
    emit("D", "worst pairwise pose RMSD, A", a2["max_pairwise_rmsd_angstrom"], "disproof angle 2")
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
             round((k + 1) / (len(null_scores) + 1), 4), "this script")
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

    # ------------------------------------------------- tally
    heading("17. Battery tally")
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
