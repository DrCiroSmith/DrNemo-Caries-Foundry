#!/usr/bin/env python3
"""
Analyse out/results.jsonl and run the disproof battery.

Everything here is recomputed from the raw docking records. No number is carried
over from the existing repository artefacts. Where a test cannot be run because
the relevant rows are absent, it reports that rather than substituting a value.

Dependencies: numpy only. scipy is deliberately avoided so that the same code
runs unchanged wherever the results are analysed. The assignment and permutation
routines are therefore implemented here and checked by --self-test.

Usage:
    python3 analyze.py --self-test      # verify the statistics code first
    python3 analyze.py                  # analyse out/results.jsonl
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).parent.resolve()
OUT = HERE / "out"
RESULTS = OUT / "results.jsonl"
CONFIG = json.loads((HERE / "config.json").read_text())
SEQUENCES = json.loads((HERE / "sequences.json").read_text())

# AutoDock Vina's own reported standard error on its training/test set. Used as the
# yardstick for whether a difference is interpretable at all, not as a p-value.
VINA_SE_KCAL = 2.85
VINA_SE_SOURCE = ("the error figure the existing campaign already uses for Vina's "
                  "scoring function; it is a property of the scoring function, not of "
                  "this run, and it is quoted here only to keep effect sizes in "
                  "perspective")


# --------------------------------------------------------------------------- #
# Assignment problem, for symmetry-corrected RMSD
# --------------------------------------------------------------------------- #
def hungarian(cost: np.ndarray) -> np.ndarray:
    """
    Minimum-cost assignment of rows to columns. Returns an array of length n_rows
    giving the column chosen for each row. Requires n_rows <= n_cols.
    Classic O(n^2 m) potentials-and-augmenting-path method.
    """
    n, m = cost.shape
    if n > m:
        raise ValueError("hungarian requires n_rows <= n_cols")
    INF = float("inf")
    u = np.zeros(n + 1)
    v = np.zeros(m + 1)
    p = np.zeros(m + 1, dtype=int)
    way = np.zeros(m + 1, dtype=int)

    for i in range(1, n + 1):
        p[0] = i
        j0 = 0
        minv = np.full(m + 1, INF)
        used = np.zeros(m + 1, dtype=bool)
        while True:
            used[j0] = True
            i0 = p[j0]
            delta, j1 = INF, -1
            for j in range(1, m + 1):
                if not used[j]:
                    cur = cost[i0 - 1, j - 1] - u[i0] - v[j]
                    if cur < minv[j]:
                        minv[j], way[j] = cur, j0
                    if minv[j] < delta:
                        delta, j1 = minv[j], j
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

    out = np.full(n, -1, dtype=int)
    for j in range(1, m + 1):
        if p[j] > 0:
            out[p[j] - 1] = j - 1
    return out


def read_pdbqt_atoms(path: Path, first_model_only: bool = True):
    """Return (coords, elements) for a PDBQT or PDB file, heavy atoms only."""
    coords, elems = [], []
    for line in path.read_text().splitlines():
        if line.startswith("ENDMDL") and first_model_only and coords:
            break
        if line.startswith(("ATOM", "HETATM")):
            el = line[76:79].strip() or line[12:16].strip()[:1]
            el = el.upper()
            if el in ("H", "HD", "HS"):
                continue
            # AutoDock types map onto elements
            el = {"A": "C", "NA": "N", "OA": "O", "SA": "S"}.get(el, el)
            coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            elems.append(el[0] if el not in ("CL", "BR") else el)
    return np.array(coords), elems


def symmetry_corrected_rmsd(a_path: Path, b_path: Path) -> dict:
    """
    RMSD between two poses of the same molecule, allowing like atoms to be matched
    to each other. Carbohydrate rings are near-symmetric, so a fixed atom order
    exaggerates the deviation; this is the standard correction and it is what the
    'Hungarian RMSD' in the existing report means.
    """
    ca, ea = read_pdbqt_atoms(a_path)
    cb, eb = read_pdbqt_atoms(b_path)
    if len(ca) == 0 or len(cb) == 0:
        return {"ok": False, "error": "one of the poses has no heavy atoms"}

    total_sq, matched = 0.0, 0
    per_element = {}
    for el in sorted(set(ea) | set(eb)):
        ia = [i for i, e in enumerate(ea) if e == el]
        ib = [i for i, e in enumerate(eb) if e == el]
        if not ia or not ib:
            continue
        A, B = ca[ia], cb[ib]
        if len(ia) <= len(ib):
            cost = ((A[:, None, :] - B[None, :, :]) ** 2).sum(-1)
            assign = hungarian(cost)
            sq = float(sum(cost[r, assign[r]] for r in range(len(ia))))
            k = len(ia)
        else:
            cost = ((B[:, None, :] - A[None, :, :]) ** 2).sum(-1)
            assign = hungarian(cost)
            sq = float(sum(cost[r, assign[r]] for r in range(len(ib))))
            k = len(ib)
        total_sq += sq
        matched += k
        per_element[el] = {"n_matched": k, "rmsd_A": round(math.sqrt(sq / k), 3)}

    if matched == 0:
        return {"ok": False, "error": "no atoms of matching element between the two poses"}
    return {"ok": True, "rmsd_A": round(math.sqrt(total_sq / matched), 3),
            "atoms_matched": matched, "atoms_a": len(ca), "atoms_b": len(cb),
            "per_element": per_element,
            "note": "symmetry-corrected: like atoms matched by minimum-cost assignment"}


def plain_rmsd(a: np.ndarray, b: np.ndarray) -> float:
    """Direct RMSD, valid only when atom order is identical."""
    return float(math.sqrt(((a - b) ** 2).sum(-1).mean()))


# --------------------------------------------------------------------------- #
# Statistics
# --------------------------------------------------------------------------- #
def permutation_p_one_sided(observed: float, null: list[float]) -> dict:
    """
    Exact empirical one-sided p for 'the lead binds better than chance'.
    Scores are negative-is-better, so a null member at least as good as the
    observed value counts against the lead. The +1 in numerator and denominator
    is the standard conservative correction; it keeps p from ever being 0, which
    a finite null cannot justify.
    """
    arr = np.asarray(null, dtype=float)
    at_least_as_good = int((arr <= observed).sum())
    p = (at_least_as_good + 1) / (len(arr) + 1)
    mean, sd = float(arr.mean()), float(arr.std(ddof=1)) if len(arr) > 1 else float("nan")
    z = (observed - mean) / sd if sd and not math.isnan(sd) and sd > 0 else float("nan")
    return {"observed": round(observed, 3), "n_null": len(arr),
            "null_mean": round(mean, 3), "null_sd": round(sd, 3) if sd == sd else None,
            "null_best": round(float(arr.min()), 3), "null_worst": round(float(arr.max()), 3),
            "null_members_at_least_as_good": at_least_as_good,
            "empirical_p_one_sided": round(p, 4),
            "z_score": round(z, 3) if z == z else None,
            "rejects_null_at_0.05": bool(p < 0.05),
            "p_definition": "(k+1)/(n+1) where k is the number of null members scoring at "
                            "least as well as the candidate; exact and distribution-free"}


def agg(values: list[float]) -> dict:
    a = np.asarray(values, dtype=float)
    return {"n": int(a.size), "mean": round(float(a.mean()), 3),
            "sd": round(float(a.std(ddof=1)), 3) if a.size > 1 else None,
            "best": round(float(a.min()), 3), "worst": round(float(a.max()), 3),
            "range": round(float(a.max() - a.min()), 3)}


# --------------------------------------------------------------------------- #
# Self-test
# --------------------------------------------------------------------------- #
def self_test() -> int:
    ok = True

    def chk(name, cond, detail=""):
        nonlocal ok
        print(f"  {'PASS' if cond else 'FAIL'}  {name}{'  ' + detail if detail else ''}")
        ok = ok and cond

    print("Hungarian assignment")
    c = np.array([[4.0, 1.0, 3.0], [2.0, 0.0, 5.0], [3.0, 2.0, 2.0]])
    a = hungarian(c)
    chk("known 3x3 optimum is 5", abs(sum(c[i, a[i]] for i in range(3)) - 5.0) < 1e-9,
        f"got {sum(c[i, a[i]] for i in range(3))}")
    chk("assignment is a permutation", sorted(a.tolist()) == [0, 1, 2], f"{a.tolist()}")
    rng = np.random.default_rng(7)
    worse = 0
    for _ in range(200):
        n = rng.integers(2, 7)
        m = np.asarray(rng.random((n, n)) * 10)
        got = sum(m[i, hungarian(m)[i]] for i in range(n))
        import itertools
        brute = min(sum(m[i, pp[i]] for i in range(n)) for pp in itertools.permutations(range(n)))
        if got > brute + 1e-9:
            worse += 1
    chk("matches brute force on 200 random matrices", worse == 0, f"{worse} mismatches")

    print("Identity and translation RMSD")
    pts = rng.random((20, 3)) * 10
    chk("RMSD of a pose with itself is 0", abs(plain_rmsd(pts, pts)) < 1e-12)
    chk("RMSD under a 3 A shift along x is 3", abs(plain_rmsd(pts, pts + [3, 0, 0]) - 3.0) < 1e-9)

    print("Permutation p-value")
    r = permutation_p_one_sided(-10.0, [-1.0] * 9)
    chk("candidate better than all 9 gives p = 0.1", abs(r["empirical_p_one_sided"] - 0.1) < 1e-9,
        f"p={r['empirical_p_one_sided']}")
    r = permutation_p_one_sided(0.0, [-1.0] * 9)
    chk("candidate worse than all 9 gives p = 1.0", abs(r["empirical_p_one_sided"] - 1.0) < 1e-9,
        f"p={r['empirical_p_one_sided']}")
    r = permutation_p_one_sided(-5.0, [-9.0, -7.0, -5.0, -3.0, -1.0])
    chk("ties count against the candidate", r["null_members_at_least_as_good"] == 3,
        f"k={r['null_members_at_least_as_good']}")
    chk("p is never 0", permutation_p_one_sided(-99.0, [0.0] * 1000)["empirical_p_one_sided"] > 0)

    print("\nSelf-test " + ("passed." if ok else "FAILED. Do not trust the analysis."))
    return 0 if ok else 1


# --------------------------------------------------------------------------- #
# Main analysis
# --------------------------------------------------------------------------- #
def load() -> list[dict]:
    if not RESULTS.exists():
        print(f"STOP: {RESULTS} does not exist. Run run_docking.py first.", file=sys.stderr)
        sys.exit(1)
    rows = []
    for line in RESULTS.read_text().splitlines():
        if line.strip():
            rows.append(json.loads(line))
    return rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    if args.self_test:
        return self_test()

    if self_test() != 0:
        return 1
    print()

    rows = load()
    good = [r for r in rows if r.get("ok")]
    report: dict = {
        "rows_total": len(rows), "rows_usable": len(good),
        "rows_failed": len(rows) - len(good),
        "receptor": json.loads((OUT / "receptor_validation.json").read_text())["title"]
        if (OUT / "receptor_validation.json").exists() else None,
        "caveats": [],
    }

    autobox = [r["key"] for r in good if r.get("vinardo_scoring_context") == "autobox_fallback"]
    if autobox:
        report["caveats"].append(
            f"{len(autobox)} runs had their Vinardo rescoring fall back to --autobox, which "
            f"scores in a different box than the docking used. Their Vinardo values are not "
            f"comparable with the rest and are excluded from Vinardo aggregates.")
    report["autobox_fallback_keys"] = autobox

    # ---- per-entity aggregation across docking seeds, wide box, conformer 42 ----
    def subset(kind=None, box="wide_28", cseed=42, receptor="8FJC"):
        return [r for r in good
                if (kind is None or r["kind"] == kind)
                and r["box_id"] == box and r["receptor_id"] == receptor
                and (r["conformer_seed"] == cseed or r["kind"] == "calibration")]

    by_id: dict[str, list[dict]] = defaultdict(list)
    for r in subset():
        by_id[r["id"]].append(r)

    primary = {}
    for ident, rs in sorted(by_id.items()):
        vina = [r["vina_kcal"] for r in rs]
        vinardo = [r["vinardo_kcal"] for r in rs
                   if r.get("vinardo_kcal") is not None
                   and r.get("vinardo_scoring_context") == "fixed_box"]
        gaps = [r["mode_gap_kcal"] for r in rs if r.get("mode_gap_kcal") is not None]
        ss = [r["ss_bond_angstrom"] for r in rs if r.get("ss_bond_angstrom") is not None]
        entry = {"kind": rs[0]["kind"], "sequence": rs[0].get("sequence"),
                 "seeds": sorted(r["dock_seed"] for r in rs),
                 "vina": agg(vina),
                 "vinardo": agg(vinardo) if vinardo else None,
                 "mode1_to_mode2_gap": agg(gaps) if gaps else None,
                 "disulfide_SS_A": agg(ss) if ss else None}
        # pose convergence across seeds, which is not the same thing as score convergence
        poses = [Path(r["pose_pdbqt"]) for r in rs if Path(r["pose_pdbqt"]).exists()]
        if len(poses) > 1:
            mats = []
            for p in poses:
                c, _ = read_pdbqt_atoms(p)
                mats.append(c)
            if len({m.shape for m in mats}) == 1:
                pw = [plain_rmsd(mats[i], mats[j])
                      for i in range(len(mats)) for j in range(i + 1, len(mats))]
                entry["pose_rmsd_across_seeds_A"] = agg(pw)
        primary[ident] = entry
    report["per_entity_wide_box_conformer42"] = primary

    # ---- redocking calibration ----
    boxes_path = OUT / "boxes.json"
    if boxes_path.exists():
        boxes = json.loads(boxes_path.read_text())
        crystal = Path(boxes["native_ligand_pdbqt"])
        cal = {}
        for box in ("wide_28", "tight_auto"):
            rs = [r for r in good if r["kind"] == "calibration" and r["box_id"] == box]
            per_seed = []
            for r in rs:
                pose = Path(str(r["pose_pdbqt"]).replace(".pdbqt", "_split_1.pdbqt"))
                if not pose.exists():
                    pose = Path(r["pose_pdbqt"])
                if pose.exists() and crystal.exists():
                    m = symmetry_corrected_rmsd(pose, crystal)
                    if m.get("ok"):
                        per_seed.append({"dock_seed": r["dock_seed"], "rmsd_A": m["rmsd_A"],
                                         "vina_kcal": r["vina_kcal"]})
            if per_seed:
                cal[box] = {
                    "per_seed": per_seed,
                    "rmsd_A": agg([d["rmsd_A"] for d in per_seed]),
                    "vina_kcal": agg([d["vina_kcal"] for d in per_seed]),
                    "passes_2A_standard": bool(max(d["rmsd_A"] for d in per_seed) < 2.0),
                    "interpretation":
                        "A redocking control reproduces the crystallographic pose of the ligand "
                        "that is actually bound in this structure. Under 2 A is the usual bar. "
                        "Above it, pose-level claims about any other ligand in this box are not "
                        "supportable, whatever the scores say.",
                }
        report["redocking_calibration"] = cal or {
            "status": "not run", "reason": "no calibration rows in results"}

    # ---- null tests ----
    # The candidates are docked with 5 seeds and the null members with 3. Taking the
    # best score across seeds would therefore give every candidate two extra draws at
    # the minimum, which flatters it for no reason connected to its sequence. So the
    # null tests are run on the seed set the two sides have in common, and on both the
    # best and the mean across those seeds, because "best of n" is itself a choice.
    NULL_SEEDS = set(CONFIG["docking"]["seeds_for_null_sets"])

    def scores_for(kind: str, box="wide_28", cseed=42, seeds=None,
                   stat="best") -> dict[str, float]:
        d: dict[str, list[float]] = defaultdict(list)
        for r in subset(kind, box, cseed):
            if seeds is None or r["dock_seed"] in seeds:
                d[r["id"]].append(r["vina_kcal"])
        if stat == "best":
            return {k: min(v) for k, v in d.items()}
        return {k: float(np.mean(v)) for k, v in d.items()}

    nulls = {}
    seq_of = {c["id"]: c["sequence"] for c in SEQUENCES["contenders"]}
    lead_id = next(c["id"] for c in SEQUENCES["contenders"]
                   if c["role"] == "reference_scaffold_P2A")
    wt_id = next(c["id"] for c in SEQUENCES["contenders"] if c["role"] == "wild_type_lead")

    # Unrestricted best-of-all-seeds, used only for the descriptive sections below.
    cand = scores_for("contender")

    # Seed-matched, used for every null test.
    cand_m = {s: scores_for("contender", seeds=NULL_SEEDS, stat=s) for s in ("best", "mean")}
    a_m = {s: scores_for("null_a", seeds=NULL_SEEDS, stat=s) for s in ("best", "mean")}
    b_m = {s: scores_for("null_b", seeds=NULL_SEEDS, stat=s) for s in ("best", "mean")}
    report["null_test_seed_matching"] = {
        "seeds_used_on_both_sides": sorted(NULL_SEEDS),
        "why": "Candidates are docked with more seeds than null members. Comparing the best "
               "of 5 against the best of 3 would favour the candidate through the extra "
               "draws alone, so both sides are restricted to the seeds they share.",
        "statistics_reported": ["best across the shared seeds", "mean across the shared seeds"],
    }

    pairings = [
        ("reference_scaffold_vs_null_B_composition_matched", lead_id, b_m,
         "PRIMARY TEST, pre-registered. The reference scaffold CAHWRPWHC against "
         "permutations of its own residues, so composition and net charge are held fixed "
         "and only the sequence order varies."),
        ("wild_type_vs_null_A_composition_matched", wt_id, a_m,
         "SECONDARY TEST, pre-registered. The wild-type lead CPHWRPWHC against null A, "
         "which is composed of permutations of the wild-type sequence."),
        ("reference_scaffold_vs_null_A_as_originally_reported", lead_id, a_m,
         "NOT A VALID TEST, reported for comparison only. This is the pairing the existing "
         "campaign reports. Null A members carry 86 heavy atoms against the reference "
         "scaffold's 84, so it confounds sequence order with a composition difference. "
         "Included to show what the original number becomes on the new receptor."),
    ]
    for name, ident, null_by_stat, note in pairings:
        entry = {"candidate_id": ident, "candidate_sequence": seq_of.get(ident), "note": note}
        for stat in ("best", "mean"):
            c, null = cand_m[stat], null_by_stat[stat]
            if ident in c and len(null) >= 5:
                entry[f"by_{stat}_across_shared_seeds"] = permutation_p_one_sided(
                    c[ident], list(null.values()))
            else:
                entry[f"by_{stat}_across_shared_seeds"] = {
                    "status": "not run",
                    "reason": f"candidate present: {ident in c}; null members: {len(null)}"}
        bests = [entry[f"by_{s}_across_shared_seeds"].get("rejects_null_at_0.05")
                 for s in ("best", "mean")]
        if None not in bests:
            entry["both_statistics_agree"] = bests[0] == bests[1]
            if not entry["both_statistics_agree"]:
                entry["disagreement_warning"] = (
                    "The best-of-seeds and mean-of-seeds statistics reach opposite verdicts. "
                    "A result that depends on which summary statistic is chosen is not a "
                    "result. Report it as inconclusive.")
        nulls[name] = entry
    report["null_tests"] = nulls

    # ---- box-size effect ----
    # Seed-matched for the same reason as the null tests: the tight box is run with 3
    # seeds and the wide box with 5.
    wide_c = scores_for("contender", "wide_28", seeds=NULL_SEEDS)
    tight_c = scores_for("contender", "tight_auto", seeds=NULL_SEEDS)
    shared = sorted(set(wide_c) & set(tight_c))
    if shared:
        deltas = [tight_c[k] - wide_c[k] for k in shared]
        report["box_size_effect"] = {
            "per_entity": {k: {"wide_28": round(wide_c[k], 3),
                               "tight_auto": round(tight_c[k], 3),
                               "delta_tight_minus_wide": round(tight_c[k] - wide_c[k], 3)}
                           for k in shared},
            "delta": agg(deltas),
            "interpretation":
                "If shrinking the box changes scores by more than it changes the ranking, the "
                "28 A cube was adding search noise rather than exploring real binding modes. "
                "A reordering of lead against controls between boxes would mean the original "
                "ranking was an artefact of box size.",
        }
    else:
        report["box_size_effect"] = {"status": "not run", "reason": "tier 3 rows absent"}

    # ---- conformer-seed effect ----
    conf: dict[str, dict[int, float]] = defaultdict(dict)
    for r in good:
        if r["box_id"] == "wide_28" and r["kind"] != "calibration" \
                and r["conformer_seed"] is not None and r["dock_seed"] == 1:
            conf[r["id"]][r["conformer_seed"]] = r["vina_kcal"]
    multi = {k: v for k, v in conf.items() if len(v) > 1}
    if multi:
        spreads = [max(v.values()) - min(v.values()) for v in multi.values()]
        report["conformer_seed_effect"] = {
            "per_entity": {k: {"by_conformer_seed": {str(s): round(x, 3) for s, x in sorted(v.items())},
                               "spread": round(max(v.values()) - min(v.values()), 3)}
                           for k, v in sorted(multi.items())},
            "spread": agg(spreads),
            "interpretation":
                "Meeko leaves the disulfide ring rigid, so each conformer seed gives one fixed "
                "backbone. If the spread across conformer seeds is comparable to the lead-minus-"
                "null difference, then that difference is substantially conformer-generation "
                "luck and no claim about sequence-specific binding survives it.",
        }
    else:
        report["conformer_seed_effect"] = {"status": "not run", "reason": "tier 4 rows absent"}

    # ---- off-target selectivity ----
    # Seed-matched: the counter-screens use 3 seeds, the on-target set 5.
    on_target = scores_for("contender", seeds=NULL_SEEDS)
    off = {}
    for rid in ("1SMD", "1LZ1"):
        s: dict[str, list[float]] = defaultdict(list)
        for r in good:
            if r["receptor_id"] == rid and r["box_id"] == "wide_28" \
                    and r["dock_seed"] in NULL_SEEDS:
                s[r["id"]].append(r["vina_kcal"])
        if s:
            off[rid] = {k: round(min(v), 3) for k, v in sorted(s.items())}
    if off:
        sel = {}
        for rid, scores in off.items():
            for ident, v in scores.items():
                if ident in on_target:
                    sel.setdefault(ident, {})[f"ddg_vs_{rid}"] = round(on_target[ident] - v, 3)
        report["off_target_selectivity"] = {
            "seeds_used_on_both_sides": sorted(NULL_SEEDS),
            "on_target_best": {k: round(v, 3) for k, v in sorted(on_target.items())},
            "off_target_best": off, "delta_delta_g": sel,
            "vina_standard_error_kcal": VINA_SE_KCAL,
            "vina_standard_error_source": VINA_SE_SOURCE,
            "interpretation":
                "Selectivity differences of under about 1 kcal/mol sit far inside the scoring "
                "function's own error and cannot be called selectivity on docking evidence alone.",
        }
    else:
        report["off_target_selectivity"] = {"status": "not run", "reason": "tier 5 rows absent"}

    # ---- determinism and geometry integrity ----
    nondet = {k: v["vina"]["range"] for k, v in primary.items()
              if v["vina"]["n"] > 1 and v["vina"]["range"] is not None}
    if nondet:
        report["seed_sensitivity"] = {
            "per_entity_score_range_kcal": {k: v for k, v in sorted(nondet.items())},
            "median_range": round(float(np.median(list(nondet.values()))), 3),
            "max_range": round(max(nondet.values()), 3),
            "interpretation":
                "This is the spread the existing pipeline could not see, because it never passed "
                "a seed to Vina. Any score difference smaller than this spread is not a result.",
        }
    broken_ss = {k: v["disulfide_SS_A"] for k, v in primary.items()
                 if v.get("disulfide_SS_A") and not (2.00 <= v["disulfide_SS_A"]["mean"] <= 2.15)}
    report["disulfide_integrity_failures"] = broken_ss or "none outside 2.00 to 2.15 A"

    path = OUT / "analysis.json"
    path.write_text(json.dumps(report, indent=2) + "\n")

    print(f"Analysed {len(good)} usable runs of {len(rows)} recorded.")
    for caveat in report["caveats"]:
        print(f"  CAVEAT: {caveat}")
    if isinstance(report.get("redocking_calibration"), dict):
        for box, c in report["redocking_calibration"].items():
            if isinstance(c, dict) and "rmsd_A" in c:
                print(f"  redocking {box}: RMSD {c['rmsd_A']['mean']} A "
                      f"(range {c['rmsd_A']['best']} to {c['rmsd_A']['worst']}), "
                      f"under 2 A: {c['passes_2A_standard']}")
    for name, t in report["null_tests"].items():
        print(f"  {name}:")
        for stat in ("best", "mean"):
            s = t.get(f"by_{stat}_across_shared_seeds", {})
            if "empirical_p_one_sided" in s:
                print(f"      by {stat} across seeds: observed {s['observed']}, null "
                      f"{s['null_mean']} +/- {s['null_sd']}, "
                      f"{s['null_members_at_least_as_good']}/{s['n_null']} at least as good, "
                      f"p = {s['empirical_p_one_sided']}")
            else:
                print(f"      by {stat} across seeds: {s.get('reason', 'not run')}")
        if t.get("disagreement_warning"):
            print(f"      WARNING: {t['disagreement_warning']}")
    print(f"\nFull report: {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
