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
import hashlib
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


def spearman(x: list[float], y: list[float]) -> float | None:
    """Rank correlation, implemented here because scipy is not a dependency."""
    if len(x) < 3 or len(x) != len(y):
        return None

    def ranks(v):
        order = np.argsort(np.asarray(v, dtype=float))
        r = np.empty(len(v), dtype=float)
        r[order] = np.arange(1, len(v) + 1, dtype=float)
        # average ties
        vals = np.asarray(v, dtype=float)
        for u in np.unique(vals):
            m = vals == u
            if m.sum() > 1:
                r[m] = r[m].mean()
        return r

    rx, ry = ranks(x), ranks(y)
    if rx.std() == 0 or ry.std() == 0:
        return None
    return round(float(np.corrcoef(rx, ry)[0, 1]), 4)


def minimum_detectable_effect(null: list[float], alpha: float = 0.05) -> dict:
    """What size of advantage this null could have detected at all.

    An exact one-sided permutation test on n nulls returns p = (k+1)/(n+1), where k
    is the number of nulls at least as good as the observation. The smallest p it can
    return is 1/(n+1), so n bounds the question before any score is computed. This is
    the measurement that says whether a null result means "no effect" or "no power",
    and the existing campaign reports a p value without it.
    """
    n = len(null)
    if n < 2:
        return {"n_null": n, "reason": "too few null members to characterise"}
    arr = np.sort(np.asarray(null, dtype=float))  # ascending: more negative is better
    smallest_p = 1.0 / (n + 1)
    # largest k with (k+1)/(n+1) < alpha
    k_max = -1
    for k in range(n + 1):
        if (k + 1) / (n + 1) < alpha:
            k_max = k
        else:
            break
    out = {
        "n_null": n,
        "alpha": alpha,
        "smallest_attainable_p": round(smallest_p, 4),
        "null_mean_kcal": round(float(arr.mean()), 3),
        "null_sd_kcal": round(float(arr.std(ddof=1)), 3),
        "null_best_kcal": round(float(arr[0]), 3),
    }
    if smallest_p >= alpha:
        out["verdict"] = (f"no observation could reach p < {alpha} with n = {n}. The test "
                          f"cannot return a significant result whatever the lead scores.")
        return out
    # the observation must beat all but k_max of the nulls
    threshold = float(arr[k_max]) if k_max >= 0 else float(arr[0])
    out["nulls_the_observation_may_not_beat"] = max(k_max, 0)
    out["score_threshold_kcal"] = round(threshold, 3)
    out["minimum_detectable_effect_kcal"] = round(float(arr.mean()) - threshold, 3)
    out["minimum_detectable_effect_in_null_sd"] = (
        round((float(arr.mean()) - threshold) / float(arr.std(ddof=1)), 2)
        if arr.std(ddof=1) > 0 else None)
    # What a larger null would and would not buy. This test asks whether the candidate
    # is extreme within its own permutation distribution, so the requirement converges
    # on the alpha-quantile of that distribution and goes no further. Raising n tightens
    # the quantile estimate; it does not move the quantile.
    sd = float(arr.std(ddof=1))
    z_alpha = 1.6449  # one-sided normal quantile at 0.05, for the large-n limit
    out["asymptotic_requirement_in_null_sd"] = z_alpha
    out["asymptotic_requirement_kcal"] = round(z_alpha * sd, 3) if sd > 0 else None
    out["what_a_larger_null_buys"] = (
        f"At any n this test asks whether the candidate falls in the most extreme "
        f"{alpha:.0%} of its own permutation distribution. With a null spread of "
        f"{sd:.3f} kcal/mol that boundary sits about {z_alpha * sd:.3f} kcal/mol better "
        f"than the null mean. A larger null sharpens the estimate of that boundary; it "
        f"does not move it. An advantage below roughly {z_alpha * sd:.2f} kcal/mol "
        f"therefore cannot be made significant by this test at any sample size, and "
        f"enlarging the null is not the remedy.")
    out["verdict"] = (
        f"with n = {n} at alpha = {alpha}, the observation must score at or below "
        f"{threshold:.3f} kcal/mol, which is {float(arr.mean()) - threshold:.3f} kcal/mol "
        f"better than the null mean. An advantage smaller than that cannot be detected "
        f"by this test at this n, however real it is.")
    return out


def effect_in_null_terms(observed: float, null: list[float]) -> dict:
    """Where an observed score sits inside the null distribution, in SD and percentile.

    Reported because an effect that sits in the middle of its own null is not a small
    significant effect waiting for a bigger sample; it is not an effect this design can
    ever call one.
    """
    if len(null) < 3:
        return {"reason": "too few null members"}
    arr = np.asarray(null, dtype=float)
    mean, sd = float(arr.mean()), float(arr.std(ddof=1))
    advantage = mean - observed  # positive means the candidate scores better
    pct = 100.0 * float((arr <= observed).sum()) / arr.size
    out = {
        "observed_kcal": round(observed, 3),
        "null_mean_kcal": round(mean, 3),
        "null_sd_kcal": round(sd, 3),
        "advantage_over_null_mean_kcal": round(advantage, 3),
        "advantage_in_null_sd": round(advantage / sd, 2) if sd > 0 else None,
        "percentile_of_the_null_the_observation_sits_at": round(pct, 1),
        "percent_of_null_members_scoring_at_least_as_well": round(pct, 1),
    }
    if sd > 0 and advantage / sd < 1.6449:
        out["reading"] = (
            f"The observation sits {advantage / sd:.2f} null standard deviations from the "
            f"null mean, inside the bulk of its own permutation distribution. The 5 percent "
            f"boundary is 1.64 SD. No sample size moves that boundary, so this is not a "
            f"small effect awaiting a larger null: it is an effect this design cannot "
            f"distinguish from sequence order noise at all.")
    else:
        out["reading"] = (
            f"The observation sits {advantage / sd:.2f} null standard deviations from the "
            f"null mean, outside the 1.64 SD boundary that a one-sided test at 0.05 "
            f"requires in the large-sample limit.")
    return out


def size_confound(entries: list[dict]) -> dict:
    """Regress score on heavy-atom count, because Vina's score grows with ligand size.

    Relevant here because the existing campaign's 39-member null is built from
    permutations of the wild type at 86 heavy atoms while the lead has 84, so the
    comparison carries a size difference as well as a sequence difference.
    """
    pts = [(e["heavy_atoms"], e["score"]) for e in entries
           if e.get("heavy_atoms") and e.get("score") is not None]
    if len(pts) < 5:
        return {"reason": f"only {len(pts)} entries carry both a score and an atom count"}
    ha = np.asarray([a for a, _ in pts], dtype=float)
    sc = np.asarray([b for _, b in pts], dtype=float)
    span = float(ha.max() - ha.min())
    if span < 5.0:
        # A slope fitted over a two-atom span is a number that looks authoritative and
        # is not. Refuse it rather than print it.
        return {"n": len(pts), "heavy_atom_range": [float(ha.min()), float(ha.max())],
                "heavy_atom_span": span,
                "reason": f"the heavy-atom counts in this set span only {span:g} atoms, "
                          f"which is too narrow to fit a slope worth reporting. Size "
                          f"cannot be ruled out as a confound from within a set whose "
                          f"members are all nearly one size. Use the slope measured over "
                          f"the full 260-contender ledger instead, recorded in "
                          f"PROTOCOL.md section 4.",
                "slope_withheld": True}
    A = np.vstack([ha, np.ones_like(ha)]).T
    slope, intercept = np.linalg.lstsq(A, sc, rcond=None)[0]
    r = float(np.corrcoef(ha, sc)[0, 1])
    return {
        "n": len(pts),
        "heavy_atom_range": [float(ha.min()), float(ha.max())],
        "slope_kcal_per_heavy_atom": round(float(slope), 4),
        "intercept_kcal": round(float(intercept), 3),
        "pearson_r": round(r, 3),
        "r_squared": round(r * r, 3),
        "interpretation":
            "Multiply the slope by the heavy-atom difference between two molecules to get "
            "the score difference attributable to size alone, before any sequence effect. "
            "Compare that with the difference being claimed.",
    }


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
# --------------------------------------------------------------------------- #
# Provenance audit
#
# In the repository this package corrects, two docked results appear under five
# filenames: the 3AIC result at -8.191 kcal/mol under three names, and the 1SMD
# result at -7.343 under two. The -8.191 is the lead's published headline number.
# Nothing inside those files records a seed, a command line, a run id or a time, so
# from the files alone there is no way to tell whether that was one run saved
# several times or several runs that collapsed onto one. A third duplicated content
# in that directory is the input ligand conformer under three names, which is
# expected from one builder at one seed and is not a defect; the distinction matters
# and the files cannot make it, which is the point. This audit exists so that the
# question can always be answered here, from the two files he sends back and
# nothing else.
# --------------------------------------------------------------------------- #
REQUIRED_PROVENANCE = ("run_id", "identity", "path_stem", "vina_command_string",
                       "dock_seed", "exhaustiveness", "box_centre", "box_size",
                       "vina_version_text", "receptor_md5", "ligand_md5",
                       "started_utc", "finished_utc", "pose_sha256")


def provenance_audit(rows: list[dict]) -> dict:
    """Return findings on whether every recorded result can be attributed to one run.

    Hard failures are the ones that mean two results cannot be told apart. They are
    reported and they make the analysis exit non-zero; they do not stop the rest of
    the report being written, because a nine-hour run should not be thrown away over
    a check that can be read afterwards.
    """
    ok = [r for r in rows if r.get("ok")]
    out: dict = {"records_audited": len(ok), "hard_failures": [], "notes": []}

    missing: dict[str, list[str]] = {}
    for r in ok:
        absent = [f for f in REQUIRED_PROVENANCE if r.get(f) in (None, "", [])]
        if absent:
            missing[r.get("identity") or r.get("key") or "<unkeyed record>"] = absent
    out["records_missing_provenance_fields"] = len(missing)
    if missing:
        out["records_missing_provenance_detail"] = dict(list(missing.items())[:10])
        out["notes"].append(
            f"{len(missing)} of {len(ok)} records lack at least one provenance field. "
            f"Those rows can be read but not independently attributed.")

    # 1. one content under more than one name: the failure in the existing repository
    by_content: dict[str, list[dict]] = {}
    for r in ok:
        h = r.get("pose_sha256")
        if h:
            by_content.setdefault(h, []).append(r)
    shared = []
    for h, group in by_content.items():
        idents = {r.get("identity") or r.get("key") for r in group}
        if len(idents) > 1:
            shared.append({
                "pose_sha256": h,
                "distinct_runs_sharing_this_content": len(idents),
                "identities": sorted(str(i) for i in idents),
                "path_stems": sorted({str(r.get("path_stem")) for r in group}),
                "run_ids": sorted({str(r.get("run_id")) for r in group}),
                "scores_kcal": sorted({r.get("vina_kcal") for r in group
                                       if r.get("vina_kcal") is not None}),
            })
    out["identical_pose_content_under_distinct_runs"] = shared or "none"
    if shared:
        out["hard_failures"].append(
            f"{len(shared)} pose file contents each appear under more than one distinct "
            f"run. Byte-identical coordinates from runs that differ in seed, box, "
            f"receptor or ligand mean the runs were not independent, whatever the record "
            f"says. This is the failure pattern already present in data/poses/ in this "
            f"repository, where two docked results sit under five filenames, and it must "
            f"be resolved before any of these rows is used.")

    # 2. one name under more than one content: the resume logic failed
    by_stem: dict[str, set] = {}
    for r in ok:
        st = r.get("path_stem")
        if st:
            by_stem.setdefault(st, set()).add(r.get("pose_sha256"))
    stem_conflicts = {k: sorted(str(x) for x in v) for k, v in by_stem.items() if len(v) > 1}
    out["one_path_holding_more_than_one_result"] = stem_conflicts or "none"
    if stem_conflicts:
        out["hard_failures"].append(
            f"{len(stem_conflicts)} output paths were written with more than one distinct "
            f"content, so an earlier result at that path was overwritten by a later one.")

    # 3. a run id reused
    run_ids = [r.get("run_id") for r in ok if r.get("run_id")]
    dup_ids = sorted({i for i in run_ids if run_ids.count(i) > 1})
    out["reused_run_ids"] = dup_ids or "none"
    if dup_ids:
        out["hard_failures"].append(
            f"{len(dup_ids)} run ids appear on more than one record, so those records "
            f"cannot be separated.")

    # 4. an identity recorded twice
    idents = [r.get("identity") for r in ok if r.get("identity")]
    dup_idents = sorted({i for i in idents if idents.count(i) > 1})
    out["identities_recorded_more_than_once"] = dup_idents[:20] or "none"
    if dup_idents:
        out["notes"].append(
            f"{len(dup_idents)} run identities appear on more than one record. Resume is "
            f"meant to prevent this. Later rows are duplicates of earlier ones and "
            f"inflate every count they enter.")

    out["passed"] = not out["hard_failures"]
    return out


def audit_pose_directory(directory: Path) -> dict:
    """Content audit of any directory of pose files, by content alone.

    This is what catches the failure in the existing repository, where the files
    carry no provenance at all and only their bytes can be compared.
    """
    by: dict[str, list[str]] = {}
    for f in sorted(Path(directory).rglob("*.pdbqt")):
        if f.is_file():
            by.setdefault(hashlib.sha256(f.read_bytes()).hexdigest(), []).append(f.name)
    dups = {h: names for h, names in by.items() if len(names) > 1}
    return {
        "directory": str(directory),
        "files_hashed": sum(len(v) for v in by.values()),
        "distinct_contents": len(by),
        "contents_under_more_than_one_name": len(dups),
        "detail": {h[:16]: names for h, names in sorted(dups.items())},
        "passed": not dups,
    }


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

    report["provenance_audit"] = provenance_audit(rows)
    if not report["provenance_audit"]["passed"]:
        report["caveats"].extend(report["provenance_audit"]["hard_failures"])

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
    # boxes is initialised here because later blocks read it. He is told that partial
    # output is useful, so the analysis has to survive a results file written before
    # boxes.json exists rather than crash on it.
    boxes = None
    boxes_path = OUT / "boxes.json"
    if boxes_path.exists():
        boxes = json.loads(boxes_path.read_text())
        crystal = Path(boxes["native_ligand_pdbqt"])
        cal = {}
        for box in ("wide_28", "tight_auto"):
            rs = [r for r in good if r["kind"] == "calibration" and r["box_id"] == box
                  and r["receptor_id"] == "8FJC"]
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
                   stat="best", receptor="8FJC") -> dict[str, float]:
        d: dict[str, list[float]] = defaultdict(list)
        for r in subset(kind, box, cseed, receptor):
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
        # What this null could have detected at all. Without it a p value above 0.05
        # cannot be read, because "no effect" and "no power" look identical.
        for stat in ("best", "mean"):
            null = null_by_stat[stat]
            if len(null) >= 2:
                entry[f"power_by_{stat}"] = minimum_detectable_effect(list(null.values()))
        obs_mde = entry.get("power_by_best", {})
        obs = entry.get("by_best_across_shared_seeds", {}).get("observed")
        if obs is not None and len(null_by_stat["best"]) >= 3:
            entry["effect_in_null_terms"] = effect_in_null_terms(
                obs, list(null_by_stat["best"].values()))
        if obs is not None and "minimum_detectable_effect_kcal" in obs_mde:
            observed_effect = obs_mde["null_mean_kcal"] - obs
            entry["observed_effect_kcal"] = round(observed_effect, 3)
            entry["observed_effect_against_minimum_detectable"] = {
                "observed_advantage_over_null_mean_kcal": round(observed_effect, 3),
                "minimum_detectable_at_this_n_kcal":
                    obs_mde["minimum_detectable_effect_kcal"],
                "ratio": (round(observed_effect / obs_mde["minimum_detectable_effect_kcal"], 2)
                          if obs_mde["minimum_detectable_effect_kcal"] else None),
                "reading":
                    "If the observed advantage is well below the minimum detectable effect, "
                    "a non-significant p value says the test was too small to see an effect "
                    "of that size, not that no effect exists. Neither direction is then "
                    "supportable from this test and the honest word is uninformative.",
            }

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

    # ---- calcium condition ----
    # prepare_receptor.py strips every heteroatom from the main receptor, which removes
    # the conserved structural Ca2+ that sits next to the catalytic triad. Tier 6 repeats
    # the same ligands in the same box against a receptor that keeps it. Neither receptor
    # is treated as correct here: the point is to measure how much the choice moves the
    # answer, because a result that survives only one of the two is not a result.
    #
    # Vina has no electrostatic term and does not read PDBQT partial charges, so a
    # retained metal cannot act on a charged ligand by Coulomb attraction. What it does
    # change is real but narrower: Vina types a metal as a hydrogen-bond donor, and an
    # absent ion leaves open volume that the steric terms can fill. Both effects are
    # local to the site, which is why this is run rather than argued about.
    ca_cfg = (boxes.get("calcium_condition") or {}) if boxes else {}
    ca_cont = scores_for("contender", seeds=NULL_SEEDS, receptor="8FJC_CA")
    if ca_cont:
        base_cont = scores_for("contender", seeds=NULL_SEEDS)
        ca_b = scores_for("null_b", seeds=NULL_SEEDS, receptor="8FJC_CA")
        base_b = {k: v for k, v in b_m["best"].items() if k in ca_b}
        ca_a = scores_for("null_a", seeds=NULL_SEEDS, receptor="8FJC_CA")
        base_a = {k: v for k, v in a_m["best"].items() if k in ca_a}
        shared_c = sorted(set(ca_cont) & set(base_cont))
        deltas = [ca_cont[k] - base_cont[k] for k in shared_c]

        # Does the lead keep the same standing against its own null when the ion is kept?
        def standing(cand: dict, null: dict):
            if lead_id not in cand or not null:
                return None
            obs = cand[lead_id]
            pt = permutation_p_one_sided(obs, list(null.values()))
            return {"lead_kcal": round(obs, 3),
                    "null_members": pt["n_null"],
                    "null_members_beating_lead": pt["null_members_at_least_as_good"],
                    "null_mean_kcal": pt["null_mean"],
                    "exact_p_one_sided": pt["empirical_p_one_sided"]}

        with_ca, without_ca = standing(ca_cont, ca_b), standing(base_cont, base_b)
        with_ca_a, without_ca_a = standing(ca_cont, ca_a), standing(base_cont, base_a)
        flipped = None
        if with_ca and without_ca:
            flipped = ((with_ca["exact_p_one_sided"] < 0.05)
                       != (without_ca["exact_p_one_sided"] < 0.05))
        flipped_a = None
        if with_ca_a and without_ca_a:
            flipped_a = ((with_ca_a["exact_p_one_sided"] < 0.05)
                         != (without_ca_a["exact_p_one_sided"] < 0.05))

        # Redocking, which is the check that decides whether any pose claim holds at all.
        ca_cal = []
        crystal = Path(boxes["native_ligand_pdbqt"]) if boxes else None
        for r in good:
            if r["kind"] != "calibration" or r["receptor_id"] != "8FJC_CA":
                continue
            pose = Path(str(r["pose_pdbqt"]).replace(".pdbqt", "_split_1.pdbqt"))
            if not pose.exists():
                pose = Path(r["pose_pdbqt"])
            if pose.exists() and crystal and crystal.exists():
                m = symmetry_corrected_rmsd(pose, crystal)
                if m.get("ok"):
                    ca_cal.append({"dock_seed": r["dock_seed"], "rmsd_A": m["rmsd_A"],
                                   "vina_kcal": r["vina_kcal"]})

        report["calcium_condition"] = {
            "receptor_built_by_prepare_step": {k: ca_cfg.get(k) for k in
                                               ("ca_ions_in_file", "ca_ions_kept",
                                                "ca_atoms_in_pdbqt", "keep_radius_A", "ions")},
            "per_contender": {k: {"no_ca": round(base_cont[k], 3),
                                  "with_ca": round(ca_cont[k], 3),
                                  "delta_with_minus_without": round(ca_cont[k] - base_cont[k], 3)}
                              for k in shared_c},
            "delta_kcal": agg(deltas),
            "lead_against_null_b_with_ca": with_ca,
            "lead_against_null_b_without_ca": without_ca,
            "verdict_flipped_by_the_ion": flipped,
            "lead_against_null_a_with_ca": with_ca_a,
            "lead_against_null_a_without_ca": without_ca_a,
            "verdict_flipped_by_the_ion_null_a": flipped_a,
            "redocking_with_ca": {"per_seed": ca_cal,
                                  "rmsd_A": agg([d["rmsd_A"] for d in ca_cal]) if ca_cal else None,
                                  "passes_2A_standard": bool(ca_cal) and
                                  max(d["rmsd_A"] for d in ca_cal) < 2.0} if ca_cal else None,
            "null_b_subset_note":
                f"Tier 6 runs {len(ca_b)} of the {len(b_m['best'])} null B members, so the "
                f"smallest exact p attainable in this condition is "
                f"{round(1.0 / (len(ca_b) + 1), 4) if ca_b else 'n/a'}. This condition is a "
                f"sensitivity check on the receptor, not an independent significance test.",
            "interpretation":
                "If the mean delta is small and the lead keeps its standing against null B "
                "either way, the metal is not driving the answer and the stripped receptor is "
                "an acceptable simplification to state. If the verdict flips, then no "
                "conclusion from the stripped receptor stands on its own, and the campaign has "
                "to be reported against the metal-containing receptor or not at all.",
        }
    else:
        report["calcium_condition"] = {
            "status": "not run",
            "reason": ca_cfg.get("problem", "tier 6 rows absent"),
            "consequence":
                "Every score in this report comes from a receptor with the conserved "
                "structural Ca2+ removed. That is a stated limitation, not a neutral choice: "
                "the metal sits beside the inhibitor subsites and one of its coordinating "
                "residues also hydrogen bonds acarbose.",
        }

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
                "Meeko opens the disulfide macrocycle at one bond, so Vina samples the "
                "torsions that opening exposes, but bond lengths, bond angles and the "
                "ring-closure geometry are fixed by the ETKDGv3 embedding and Vina varies "
                "none of them. Each embedding seed therefore gives one fixed reference "
                "geometry. If the spread across embedding seeds is comparable to the "
                "lead-minus-null difference, that difference is substantially "
                "conformer-generation luck and no claim about sequence-specific binding "
                "survives it. An earlier revision of this file said the ring was held "
                "rigid; that was wrong and is corrected in PROTOCOL.md section 2.3.",
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

    # ---- size confound, ligand efficiency, and scoring-function agreement ----
    # None of these costs a docking run. All three are measurements the existing
    # campaign could have made on data it already had and did not.
    ha_of: dict[str, int] = {}
    for r in good:
        lm = r.get("ligand_meta") or {}
        if lm.get("heavy_atoms"):
            ha_of.setdefault(r["id"], lm["heavy_atoms"])

    entries = []
    for ident, e in primary.items():
        if e["kind"] == "calibration":
            continue
        best = (e.get("vina") or {}).get("best")
        if best is None:
            continue
        entries.append({"id": ident, "kind": e["kind"], "heavy_atoms": ha_of.get(ident),
                        "score": best})

    report["size_confound"] = {
        "all_docked_entities": size_confound(entries),
        "within_null_a": size_confound([e for e in entries if e["kind"] == "null_a"]),
        "within_null_b": size_confound([e for e in entries if e["kind"] == "null_b"]),
        "why_this_is_here":
            "Vina's score grows with ligand size. The null distribution the existing "
            "campaign reports is built from permutations of the wild type at 86 heavy "
            "atoms, while the lead it is compared against has 84, so that comparison "
            "carries a size difference as well as a sequence difference. Null B in this "
            "package fixes the pairing; the slope is reported so the residual size term "
            "can be subtracted rather than assumed away.",
    }

    le = {}
    for e in entries:
        if e["heavy_atoms"]:
            le[e["id"]] = {"vina_kcal": round(e["score"], 3),
                           "heavy_atoms": e["heavy_atoms"],
                           "kcal_per_heavy_atom": round(e["score"] / e["heavy_atoms"], 4)}
    report["ligand_efficiency"] = {
        "per_entity": dict(sorted(le.items())),
        "note": "Score divided by heavy-atom count, reported beside the raw score rather "
                "than instead of it. A ranking that reverses between the two is a ranking "
                "by molecular weight.",
    }
    if le:
        by_raw = sorted(le, key=lambda k: le[k]["vina_kcal"])
        by_eff = sorted(le, key=lambda k: le[k]["kcal_per_heavy_atom"])
        report["ligand_efficiency"]["rank_agreement_raw_vs_efficiency"] = {
            "spearman": spearman([le[k]["vina_kcal"] for k in by_raw],
                                 [le[k]["kcal_per_heavy_atom"] for k in by_raw]),
            "best_by_raw_score": by_raw[0],
            "best_by_efficiency": by_eff[0],
            "same_winner": by_raw[0] == by_eff[0],
        }

    pairs = [(e["vina"]["best"], e["vinardo"]["best"]) for e in primary.values()
             if e.get("vina") and e.get("vinardo")
             and e["vina"].get("best") is not None and e["vinardo"].get("best") is not None]
    report["vina_vinardo_agreement"] = {
        "n_entities": len(pairs),
        "spearman_rank_correlation": spearman([a for a, _ in pairs], [b for _, b in pairs]),
        "mean_absolute_difference_kcal": (round(float(np.mean([abs(a - b) for a, b in pairs])), 3)
                                          if pairs else None),
        "what_it_is_not":
            "Vinardo is a reparameterisation of the same functional form as Vina, scored "
            "on the same pose by the same program. High agreement between them is not "
            "independent confirmation of anything; it is one method agreeing with itself. "
            "Reported because disagreement would be informative, not because agreement is.",
    }

    # ---- the noise floor, in one place ----
    # Tests 4, 5 and 6 in section 5 of PROTOCOL.md cannot be adjudicated unless the
    # claimed difference and every measured source of spread are set side by side.
    seed_ranges = [e["vina"]["range"] for e in primary.values()
                   if e.get("vina") and e["vina"].get("range") is not None
                   and e["kind"] != "calibration"]
    pose_spreads = [e["pose_rmsd_across_seeds_A"]["mean"] for e in primary.values()
                    if e.get("pose_rmsd_across_seeds_A")]
    conf_block = report.get("conformer_seed_effect") or {}
    conf_spread = conf_block.get("spread") if isinstance(conf_block, dict) else None
    primary_test = report["null_tests"].get(
        "reference_scaffold_vs_null_B_composition_matched", {})
    claimed = primary_test.get("observed_effect_kcal")
    floor = {
        "claimed_advantage_of_lead_over_null_mean_kcal": claimed,
        "median_score_range_across_docking_seeds_kcal":
            round(float(np.median(seed_ranges)), 3) if seed_ranges else None,
        "max_score_range_across_docking_seeds_kcal":
            round(max(seed_ranges), 3) if seed_ranges else None,
        "median_pose_rmsd_across_docking_seeds_A":
            round(float(np.median(pose_spreads)), 3) if pose_spreads else None,
        "conformer_seed_spread_kcal": (conf_spread.get("worst")
                                       if isinstance(conf_spread, dict) else None),
        "vina_reference_standard_error_kcal": VINA_SE_KCAL,
        "vina_reference_standard_error_source": VINA_SE_SOURCE,
        "minimum_detectable_effect_kcal":
            (primary_test.get("power_by_best") or {}).get("minimum_detectable_effect_kcal"),
        "boundary_a_larger_null_would_not_move_kcal":
            (primary_test.get("power_by_best") or {}).get("asymptotic_requirement_kcal"),
        "advantage_in_null_sd":
            (primary_test.get("effect_in_null_terms") or {}).get("advantage_in_null_sd"),
    }
    exceeded = []
    if claimed is not None and claimed <= 0:
        floor["claimed_advantage_is_smaller_than"] = (
            "not applicable: the candidate did not score better than the null mean at all, "
            "so there is no advantage to compare against the noise.")
        floor["reading"] = (
            "The candidate scored no better than the mean of its own matched permutations. "
            "Read together with the minimum detectable effect, which bounds what this test "
            "could have found, the test is uninformative rather than negative: it supports "
            "neither the candidate nor its absence.")
        report["noise_floor"] = floor
        exceeded = None
    if claimed is not None and exceeded is not None:
        for label, key in (("docking-seed spread", "max_score_range_across_docking_seeds_kcal"),
                           ("conformer-seed spread", "conformer_seed_spread_kcal"),
                           ("the test's minimum detectable effect",
                            "minimum_detectable_effect_kcal")):
            v = floor.get(key)
            if v is not None and claimed < v:
                exceeded.append(f"{label} ({v} kcal/mol)")
    if exceeded is None:
        exceeded = []
    else:
        floor["claimed_advantage_is_smaller_than"] = (
            exceeded or "none of the measured sources")
        floor["reading"] = (
            "A score difference smaller than any one of these is not a finding. Where "
            "the claimed advantage is smaller than the minimum detectable effect, the "
            "test is uninformative in both directions: it supports neither the lead nor "
            "its absence."
            if exceeded else
            "The claimed advantage exceeds every measured source of spread listed here, "
            "which is a necessary condition for reporting it and not a sufficient one.")
    floor["sign_convention"] = ("a more negative Vina score is better, so an advantage is "
                                "null mean minus candidate score and a negative advantage "
                                "means the candidate scored worse than its own null")
    report["noise_floor"] = floor

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

    pa = report["provenance_audit"]
    if not pa["passed"]:
        print("\n" + "=" * 72)
        print("PROVENANCE CHECK FAILED. The report above is written but must not be used")
        print("until this is resolved, because two or more results cannot be told apart.")
        for f in pa["hard_failures"]:
            print(f"\n  {f}")
        print("\nDetail is under provenance_audit in " + str(path))
        print("=" * 72)
        return 2
    if pa.get("records_missing_provenance_fields"):
        print(f"  note: {pa['records_missing_provenance_fields']} records lack a "
              f"provenance field; see provenance_audit in the report.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
