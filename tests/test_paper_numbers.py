"""Pins the figures the manuscript's resolution-limit argument rests on.

Section 3.3b of paper/negative_result_GTF-SI.md states that the margin this
screen produced cannot be brought inside a 5 per cent rejection region by
adding control sequences, because the null's dispersion is run-to-run noise
rather than sequence. That argument is arithmetic on two committed artefacts.
If either changes, the argument has to be rewritten rather than quietly
inherited, so these tests recompute every number in it and also check that the
prose still carries the recomputed values.

The load-bearing one is test_the_nulls_spread_is_not_sequence_driven. If it
ever fails, Section 3.3b is wrong and the paper must not be published as it
stands.
"""

import csv
import json
import math
import statistics
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
PAPER = REPO_ROOT / "paper" / "negative_result_GTF-SI.md"
LEDGER = REPO_ROOT / "data" / "contenders_ledger.csv"
REPORT = REPO_ROOT / "data" / "exhaustive_disproof_report.json"

Z95 = 1.6448536269514722
Z80 = 0.8416212335729143
DUPLICATED_SCORE = -8.191  # the value carried by the duplicated file, Section 3.10


def _ledger():
    return list(csv.DictReader(LEDGER.open()))


def _null_scores():
    return [float(r["vina_3aic_kcal"]) for r in _ledger()
            if r["falsification_status"] == "SCRAMBLED_NULL_CONTROL"]


def _repeat_runs():
    """The independent repeat runs of one ligand, the duplicate excluded."""
    report = json.loads(REPORT.read_text())
    seeds = next(a["vina_scores"] for a in report["angles"] if a.get("vina_scores"))
    return [v for v in seeds if abs(v - DUPLICATED_SCORE) > 1e-9]


def test_the_null_is_thirty_nine_members_with_the_stated_statistics():
    ns = _null_scores()
    assert len(ns) == 39
    assert round(statistics.mean(ns), 3) == -8.343
    assert round(statistics.stdev(ns), 3) == 0.470
    assert round(min(ns), 3) == -9.378


def test_significance_requires_beating_every_control():
    """At n = 39 the permutation test has no room: p < 0.05 means k = 0."""
    n = len(_null_scores())
    assert round(1 / (n + 1), 4) == 0.0250
    reachable = [k for k in range(n + 1) if (k + 1) / (n + 1) < 0.05]
    assert reachable == [0], "the paper's threshold statement assumes k must be zero"
    ns = _null_scores()
    assert round(min(ns) - statistics.mean(ns), 3) == -1.035


def test_the_leads_margin_over_its_parent():
    rows = {r["sequence"]: float(r["vina_3aic_kcal"]) for r in _ledger()
            if r["falsification_status"] != "SCRAMBLED_NULL_CONTROL"}
    assert round(rows["CAHWRPWHC"] - rows["CPHWRPWHC"], 3) == -0.129
    assert round(rows["CAHWRPWHC"] - statistics.mean(_null_scores()), 3) == 0.152


def test_the_nulls_spread_is_not_sequence_driven():
    """The claim Section 3.3b rests on.

    Treating a control's score as a sequence effect plus run-to-run noise, the
    between-sequence variance is the null variance minus the repeat-run
    variance. The point estimate is zero or negative. If new data ever makes it
    clearly positive, the section is wrong.
    """
    null_var = statistics.variance(_null_scores())
    repeat = _repeat_runs()
    assert len(repeat) == 4, "Section 3.2 reports four independent runs, not five"
    repeat_var = statistics.variance(repeat)
    assert repeat_var >= null_var, (
        f"between-sequence variance is now {null_var - repeat_var:+.4f}, so the null's "
        "spread is no longer consistent with pure run-to-run noise and Section 3.3b "
        "must be rewritten")
    assert round(null_var - repeat_var, 4) == -0.0053


@pytest.mark.parametrize("z, expected_m, expected_runs", [
    (Z95, 36, 1440),
    (Z95 + Z80, 82, 3280),
])
def test_seeds_per_control_needed_to_resolve_that_margin(z, expected_m, expected_runs):
    """Averaging seeds moves the boundary; adding controls does not."""
    ns = _null_scores()
    rows = {r["sequence"]: float(r["vina_3aic_kcal"]) for r in _ledger()
            if r["falsification_status"] != "SCRAMBLED_NULL_CONTROL"}
    margin = abs(rows["CAHWRPWHC"] - rows["CPHWRPWHC"])
    m = math.ceil((z * statistics.stdev(ns) / margin) ** 2)
    assert m == expected_m
    assert (len(ns) + 1) * m == expected_runs


def test_the_five_percent_boundary_is_set_by_the_spread():
    ns = _null_scores()
    assert round(Z95 * statistics.stdev(ns), 3) == 0.773
    rows = {r["sequence"]: float(r["vina_3aic_kcal"]) for r in _ledger()
            if r["falsification_status"] != "SCRAMBLED_NULL_CONTROL"}
    margin = abs(rows["CAHWRPWHC"] - rows["CPHWRPWHC"])
    assert round(margin / statistics.stdev(ns), 2) == 0.27


def test_the_manuscript_still_carries_the_recomputed_figures():
    """Paper and data must not drift apart silently."""
    text = PAPER.read_text()
    assert "## 3.3b" in text, "Section 3.3b is missing from the manuscript"
    for figure in ("-9.378", "1.035", "0.129", "0.773", "0.470", "0.475",
                   "0.2207", "0.2259", "m = 36", "1,440", "m = 82", "3,280"):
        assert figure in text, f"the manuscript no longer states {figure}"
