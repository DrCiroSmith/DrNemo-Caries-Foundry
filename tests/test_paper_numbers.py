"""Pins the figures the manuscript's resolution-limit argument rests on.

Section 3.3b of paper/negative_result_GTF-SI.md states that the margin this
screen produced cannot be brought inside a 5 per cent rejection region by
adding control sequences, because the rejection boundary is set by the null's
dispersion rather than by the number of controls. That argument is arithmetic
on two committed artefacts. If either changes, the argument has to be rewritten
rather than quietly inherited, so these tests recompute every number in it and
also check that the prose still carries the recomputed values.

The load-bearing one is test_the_nulls_spread_is_not_demonstrably_sequence_driven.
It pins the comparison the section is allowed to make, which is that the two
standard deviations are of the same order and that neither is resolved well
enough to separate them. It deliberately does not pin the sign of the variance
difference, which sits inside the noise and would flip on a rerun for no
scientific reason.
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

# Chi-square quantiles for a two-sided 95 per cent interval on a standard
# deviation from four observations (df = 3). Published table values, and
# reproduced to four decimals by chi2_ppf in scripts/compute_paper_numbers.py.
CHI2_DF3_LOWER_TAIL = 0.2158
CHI2_DF3_UPPER_TAIL = 9.3484
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


def _all_recorded_runs():
    """All five values the record states, including the duplicated one."""
    report = json.loads(REPORT.read_text())
    return next(a["vina_scores"] for a in report["angles"] if a.get("vina_scores"))


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


# Chi-square quantiles for a two-sided 95 per cent interval at df = 4, for the
# five-run reading. Published table values.
CHI2_DF4_LOWER_TAIL = 0.4844
CHI2_DF4_UPPER_TAIL = 11.1433


def _sd_interval(data):
    """Two-sided 95 per cent interval on the sd of four or five observations."""
    df = len(data) - 1
    lower, upper = {3: (CHI2_DF3_LOWER_TAIL, CHI2_DF3_UPPER_TAIL),
                    4: (CHI2_DF4_LOWER_TAIL, CHI2_DF4_UPPER_TAIL)}[df]
    sd = statistics.stdev(data)
    return sd * math.sqrt(df / upper), sd * math.sqrt(df / lower)


def _repeat_sd_interval():
    """The four-run interval, which is the pair Section 3.3b prints first."""
    return _sd_interval(_repeat_runs())


def _betacf(a, b, x):
    """Continued fraction for the incomplete beta, modified Lentz."""
    tiny, eps = 1e-300, 3e-16
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c, d = 1.0, 1.0 - qab * x / qap
    d = 1.0 / (tiny if abs(d) < tiny else d)
    h = d
    for m in range(1, 300):
        m2 = 2 * m
        for aa in (m * (b - m) * x / ((qam + m2) * (a + m2)),
                   -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))):
            d = 1.0 + aa * d
            d = tiny if abs(d) < tiny else d
            c = 1.0 + aa / c
            c = tiny if abs(c) < tiny else c
            d = 1.0 / d
            delta = d * c
            h *= delta
        if abs(delta - 1.0) < eps:
            return h
    raise RuntimeError("incomplete beta did not converge")


def _f_sf(x, df1, df2):
    """P(F > x). Verified against published F tables in test_the_f_test_machinery."""
    a, b = df1 / 2.0, df2 / 2.0
    z = df1 * x / (df1 * x + df2)
    lead = math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
    if z < (a + 1.0) / (a + b + 2.0):
        cdf = math.exp(lead + a * math.log(z) + b * math.log1p(-z)) * _betacf(a, b, z) / a
    else:
        cdf = 1.0 - math.exp(lead + b * math.log1p(-z) + a * math.log(z)) \
            * _betacf(b, a, 1.0 - z) / b
    return 1.0 - cdf


@pytest.mark.parametrize("df1, df2, critical", [
    (3, 3, 9.2766), (5, 5, 5.0503), (10, 10, 2.9782),
    (2, 10, 4.1028), (10, 2, 19.3959), (20, 20, 2.1242), (1, 30, 4.1709),
])
def test_the_f_test_machinery_reproduces_published_tables(df1, df2, critical):
    """No scipy here, so the F distribution is implemented and must be checked."""
    assert round(_f_sf(critical, df1, df2), 4) == 0.0500


def test_the_f_test_machinery_matches_the_closed_form_at_one_and_one():
    for x in (0.25, 1.0, 4.0, 10.0):
        closed = 1.0 - (2.0 / math.pi) * math.atan(math.sqrt(x))
        assert abs(_f_sf(x, 1, 1) - closed) < 1e-12


def test_the_interval_the_manuscript_quotes_is_the_one_the_data_gives():
    """Anchors the two bounds Section 3.3b prints.

    This is a prose check, not the scientific claim. If the committed runs ever
    change, this test should fail so the manuscript gets updated, while
    test_the_nulls_spread_is_not_demonstrably_sequence_driven keeps holding as
    long as the claim itself survives.
    """
    low, high = _repeat_sd_interval()
    assert (round(low, 3), round(high, 3)) == (0.269, 1.772), (
        f"Section 3.3b prints 0.269 to 1.772; the data now gives {low:.3f} to {high:.3f}")


@pytest.mark.parametrize("reading", ["four", "five"])
def test_the_nulls_spread_is_not_demonstrably_sequence_driven(reading):
    """The claim Section 3.3b is allowed to make, and no more than that.

    The section says the spread of one ligand docked repeatedly is of the same
    order as the spread across 39 different sequences, so the null's dispersion
    is not demonstrably sequence-driven. Two things have to hold for that: the
    two standard deviations must be close, and the repeat-run estimate must be
    too imprecise to separate them. Both are checked here.

    What is deliberately not checked is the sign of the null variance minus the
    repeat-run variance. That difference is -0.005 kcal squared per mol squared
    against standard deviations near 0.47, so it is inside the noise, and an
    assertion on its sign would fail on a rerun for no scientific reason.
    """
    null_scores = _null_scores()
    null_sd = statistics.stdev(null_scores)
    repeat = _repeat_runs() if reading == "four" else _all_recorded_runs()
    assert len(repeat) == (4 if reading == "four" else 5), (
        "Section 3.2 excludes the duplicate and reports four runs; the record states five")
    repeat_sd = statistics.stdev(repeat)

    ratio = repeat_sd / null_sd
    assert 0.8 <= ratio <= 1.25, (
        f"on the {reading}-run reading the single-ligand sd is now {ratio:.3f} times the "
        "null sd, so the two are no longer of the same order and Section 3.3b must be "
        "rewritten")

    low, high = _sd_interval(repeat)
    assert low <= null_sd <= high, (
        f"on the {reading}-run reading the null sd {null_sd:.4f} falls outside the "
        f"interval {low:.3f} to {high:.3f} on the single-ligand sd, so the two estimates "
        "are distinguishable and Section 3.3b's comparison no longer holds")

    # The conventional statement of the same thing, and the one a referee will
    # look for: an F test for equal variances must fail to reject.
    f_stat = statistics.variance(null_scores) / statistics.variance(repeat)
    p_value = _f_sf(f_stat, len(null_scores) - 1, len(repeat) - 1)
    assert p_value > 0.05, (
        f"on the {reading}-run reading an F test now separates the two variances "
        f"(F = {f_stat:.3f}, p = {p_value:.4f}), so Section 3.3b's claim that the data "
        "cannot resolve a sequence-dependent component is wrong")


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


def _losable(n, alpha=0.05):
    """Controls a candidate may lose to and still reach p < alpha, at one score each."""
    j = -1
    while (j + 2) / (n + 1) < alpha:
        j += 1
    return j


def test_enlarging_the_null_cannot_take_the_requirement_below_the_five_percent_point():
    """Why "however many controls are added" is not a power excuse.

    Section 3.3b claims the 5 per cent quantile is a floor beneath the score a
    candidate must reach, so no control count makes the test easier than that
    quantile. The required rank as a fraction of n must therefore stay below
    0.05 at every n. If it ever reached 0.05 the claim would be wrong, and the
    paragraph would have to be rewritten.
    """
    assert _losable(19) == -1, "p < 0.05 should be unattainable with 19 controls"
    assert _losable(39) == 0, "at n = 39 the candidate must beat every control"
    assert (_losable(199), _losable(999)) == (8, 48)
    previous = 0.0
    for n in (39, 59, 99, 199, 499, 999, 9999, 999999):
        quantile = (_losable(n) + 1) / n
        assert quantile < 0.05, (
            f"at n = {n} the required quantile is {quantile:.5f}, so the 5 per cent point "
            "is not a floor beneath the requirement and Section 3.3b is wrong")
        assert quantile > previous, "the requirement should relax monotonically with n"
        previous = quantile


def test_the_requirement_at_thirty_nine_is_stricter_than_the_asymptotic_floor():
    """The direction matters: more controls make the test easier, not harder."""
    ns = _null_scores()
    sd = statistics.stdev(ns)
    at_39 = (statistics.fmean(ns) - min(ns)) / sd
    assert round(at_39, 2) == 2.20
    assert at_39 > Z95, (
        "beating the best of 39 should ask more than beating the 5 per cent point; "
        "if it asks less, Section 3.3b's direction of travel is inverted")
    rows = {r["sequence"]: float(r["vina_3aic_kcal"]) for r in _ledger()
            if r["falsification_status"] != "SCRAMBLED_NULL_CONTROL"}
    margin = abs(rows["CAHWRPWHC"] - rows["CPHWRPWHC"])
    assert round(Z95 * sd / margin, 2) == 5.99


def test_the_five_percent_boundary_is_set_by_the_spread():
    ns = _null_scores()
    assert round(Z95 * statistics.stdev(ns), 3) == 0.773
    rows = {r["sequence"]: float(r["vina_3aic_kcal"]) for r in _ledger()
            if r["falsification_status"] != "SCRAMBLED_NULL_CONTROL"}
    margin = abs(rows["CAHWRPWHC"] - rows["CPHWRPWHC"])
    assert round(margin / statistics.stdev(ns), 2) == 0.27


def _section_3_3b():
    text = PAPER.read_text()
    start = text.index("## 3.3b")
    return text[start:text.index("## 3.4", start)]


def test_section_3_3b_warns_the_reader_about_the_five_run_reading():
    """A reader shown two readings will reach for the one with more data.

    So the paragraph that offers both must say, in the same place, that the
    five-run reading includes a file whose coordinates match the campaign pose
    across all nine modes and is therefore not an independent replicate. If
    that warning is ever dropped, the weaker reading goes to a referee with
    nothing attached to it.
    """
    section = _section_3_3b()
    assert "seed 42" in section, "Section 3.3b no longer names the duplicated run"
    assert "nine" in section, "Section 3.3b no longer says the match covers all nine modes"
    for phrase in ("0.415", "0.475", "both ways"):
        assert phrase in section, f"Section 3.3b no longer carries {phrase}"
    assert "upper bound" in section, (
        "Section 3.3b no longer calls the single-ligand figure an upper bound, which is "
        "the whole of the conformer-reseeding caveat")


def test_the_manuscript_still_carries_the_recomputed_figures():
    """Paper and data must not drift apart silently."""
    text = PAPER.read_text()
    assert "## 3.3b" in text, "Section 3.3b is missing from the manuscript"
    for figure in ("-9.378", "1.035", "0.129", "0.152", "0.773", "0.470", "0.475",
                   "2.20", "1.64", "0.026", "0.045", "0.049",
                   "0.415", "0.248", "1.192", "0.88",
                   "0.269", "1.772", "0.384", "0.605",
                   "m = 36", "1,440", "m = 82", "3,280"):
        assert figure in text, f"the manuscript no longer states {figure}"


# ---------------------------------------------- the reference pose's own spread

NATIVE_LIGAND = REPO_ROOT / "demo_data" / "receptors" / "3AIC_native_ligand.pdb"


def _native_ligand_b_factors():
    out = []
    for line in NATIVE_LIGAND.read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            out.append((float(line[54:60]), float(line[60:66])))
    return out


def test_the_reference_ligand_displacement_figures_match_the_coordinates():
    """Section 3.1's counter-argument is computed from the committed file.

    The paper now states the reference acarbose's displacement parameters,
    because they are the strongest argument against its own redocking
    conclusion and a reader is entitled to the number. Those figures are
    recomputed here rather than trusted.
    """
    if not NATIVE_LIGAND.exists():
        pytest.skip("the 3AIC reference ligand is no longer committed")
    rows = _native_ligand_b_factors()
    assert len(rows) == 44, (
        f"the reference ligand now has {len(rows)} atoms, not the 44 the paper "
        f"states; every figure in Section 3.1 that divides by it needs redoing")
    occupancies = {o for o, _ in rows}
    assert occupancies == {1.0}, (
        f"the reference ligand is no longer at full occupancy ({occupancies}), "
        f"so the paper's 'at full occupancy' is wrong")

    b = [x for _, x in rows]
    mean_b = statistics.fmean(b)
    u = math.sqrt(mean_b / (8 * math.pi ** 2))
    paper = PAPER.read_text()

    for value in (f"{min(b):.2f}", f"{max(b):.2f}", f"{mean_b:.2f}"):
        assert value in paper, (
            f"the paper does not carry the recomputed displacement parameter "
            f"{value} from the committed reference ligand")
    assert f"{u:.2f}" in paper, (
        f"the paper does not carry the RMS displacement {u:.2f} A implied by "
        f"those parameters")
    assert f"{4.507 / u:.1f} times" in paper, (
        f"the paper does not state the redocking RMSD as "
        f"{4.507 / u:.1f} times the reference ligand's mean displacement, "
        f"which is the bound that keeps the objection from closing the gap")


def test_the_crystallization_reference_is_cited_in_the_body():
    """A listed-but-uncited reference is a loose end in a paper about rigour."""
    paper = PAPER.read_text()
    body, _, refs = paper.partition("# References")
    assert "[6]" in body, (
        "reference 6 is listed but never cited in the body; it is the source "
        "of the 2.1 A figure that the superseded preprint attached to 3AIC")
    assert "Crystallization and preliminary X-ray analysis" in refs
