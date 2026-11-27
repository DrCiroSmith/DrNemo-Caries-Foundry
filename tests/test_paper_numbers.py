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

import re
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


def test_nothing_still_calls_the_scoring_error_uncited():
    """That claim was withdrawn, so it must not survive anywhere else.

    The paper traced the 2.85 kcal/mol constant to Trott and Olson's reported
    standard error over their 190-complex test set and withdrew the earlier
    "untraceable" note. The README, the project page and the battery's own
    summary line went on calling it assumed and uncited, which left the
    repository asserting two different things about the same number. A line may
    mention both only if it is the withdrawal itself.
    """
    paths = ["README.md", "docs/index.html", "paper/negative_result_GTF-SI.md",
             "scripts/run_exhaustive_disproof.py", "dashboard.py"]
    withdrawal = ("and it is not", "I withdraw", "was recorded here as uncited")
    offenders = []
    for rel in paths:
        f = REPO_ROOT / rel
        if not f.exists():
            continue
        for n, line in enumerate(f.read_text().splitlines(), 1):
            if "2.85" not in line:
                continue
            if ("uncited" in line or "untraceable" in line) and not any(
                    w in line for w in withdrawal):
                offenders.append(f"{rel}:{n}")
    assert not offenders, (
        f"these lines still call the 2.85 kcal/mol constant uncited or "
        f"untraceable after that claim was withdrawn: {offenders}")


def _reference_list():
    """The numbered references, as {number: line}."""
    text = PAPER.read_text()
    start = text.index("# References")
    refs = {}
    for line in text[start:].splitlines():
        m = re.match(r"^(\d+)\.\s+(.*)$", line.strip())
        if m:
            refs[int(m.group(1))] = m.group(2)
    return refs


def _body_citations():
    """Every citation marker in the body, excluding the reference list.

    Handles [n], [n,m] and [n-m]. A bracket holding digits in any other shape
    is a reported failure rather than a skipped token: an earlier version of
    this helper matched only [n], so the grouped citation [2,3] was invisible
    and it reported reference 3 as uncited when the paper cites it correctly.
    A citation parser that silently drops a form it does not know produces
    exactly the false accusation this paper exists to avoid making.
    """
    text = PAPER.read_text()
    body = text[: text.index("# References")]
    cited, parsed_spans = set(), []
    for m in re.finditer(r"\[([0-9][0-9,\s\u2013-]*)\]", body):
        parsed_spans.append(m.span())
        for part in re.split(r"[,\s]+", m.group(1).strip()):
            if not part:
                continue
            rng = re.fullmatch(r"(\d+)[\u2013-](\d+)", part)
            if rng:
                cited.update(range(int(rng.group(1)), int(rng.group(2)) + 1))
            else:
                assert part.isdigit(), f"unparsed citation fragment {part!r}"
                cited.add(int(part))
    # Anything bracketed that holds a digit and was not parsed above must be a
    # markdown label, not a citation. Surface it rather than assume.
    for m in re.finditer(r"\[[^\]]*\d[^\]]*\]", body):
        if m.span() in parsed_spans:
            continue
        label = m.group(0)
        assert body[m.end():m.end() + 1] == "(", (
            f"{label} is bracketed, holds a digit, is not a link label and was "
            "not read as a citation, so this audit cannot say whether it is one")
    return cited


def test_every_bracketed_citation_resolves_to_a_reference():
    """A citation to a reference that does not exist is a broken claim."""
    refs = _reference_list()
    assert refs, "the paper has no parsable reference list"
    dangling = sorted(_body_citations() - set(refs))
    assert not dangling, (
        f"the body cites {dangling} and the reference list has "
        f"{sorted(refs)}, so those claims point at nothing")


def test_every_reference_is_cited_at_least_once():
    """An uncited reference is either a missing citation or padding."""
    refs = _reference_list()
    uncited = sorted(set(refs) - _body_citations())
    assert not uncited, (
        f"references {uncited} are listed and never cited in the body, so "
        "either a claim is missing its marker or the entry should be cut")


def test_the_reference_numbering_is_contiguous_from_one():
    refs = _reference_list()
    assert sorted(refs) == list(range(1, len(refs) + 1)), (
        f"the reference numbers are {sorted(refs)}, which is not 1..n")


def test_the_catalytic_assignment_carries_a_citation():
    """Which residue is nucleophile is not readable from coordinates.

    The identities of residues 477, 515 and 588 are read from the committed
    structure. Their mechanistic roles are not, and this paper objects
    elsewhere to residue claims it cannot trace, so the assignment has to name
    its source and say what kind of source it is.
    """
    text = PAPER.read_text()
    i = text.index("Asp477 as nucleophile")
    window = text[i : i + 900]
    assert re.search(r"\[\d+\]", window), (
        "the catalytic triad assignment no longer cites a source within the "
        "passage that makes it")
    assert "QM/MM" in window and "computational" in window, (
        "the assignment no longer says that its source is computational rather "
        "than experimental, which overstates how well established it is")


def test_the_qmmm_reference_is_the_one_that_was_verified():
    """Pinned to the DOI, so swapping the source cannot pass silently."""
    refs = _reference_list()
    match = [r for r in refs.values() if "10.1039/C8OB00284C" in r]
    assert len(match) == 1, (
        "the QM/MM source verified for the catalytic assignment "
        "(doi 10.1039/C8OB00284C) is not in the reference list exactly once")
    entry = match[0]
    for part in ("Jaña GA", "Org Biomol Chem", "2018", "16(14):2438-47"):
        assert part in entry, f"the verified citation no longer states {part!r}: {entry}"


# --------------------------------------------------------------------------
# Provenance of the published numbers, as opposed to their correctness.
#
# Every defect found in this repository until now was a real computation
# wearing a label it had not earned. The angle-11 drift figure was something
# else: a number with no computation behind it at all, which no label check and
# no wording check would ever have caught. These two tests are for that class.
# --------------------------------------------------------------------------

BENCHMARK_SHEET = REPO_ROOT / "FALSIFICATION_BENCHMARK.md"


def _script_summary():
    """Run the provenance script and read its own closing counts.

    Thirteen of its checks read the project's commit history, so the totals are
    lower in a copy without .git: 378 against 391 when this was written. A
    checkout without history cannot check the paper's figures, so this skips
    there rather than reporting a defect that is really a missing clone.
    """
    import subprocess
    if not (REPO_ROOT / ".git").exists():
        pytest.skip("no git history in this checkout, so the script's totals "
                    "are lower than the paper's by the history-derived checks")
    proc = subprocess.run(
        ["python3", str(REPO_ROOT / "scripts" / "compute_paper_numbers.py")],
        capture_output=True, text=True, timeout=600, cwd=str(REPO_ROOT))
    assert proc.returncode == 0, (
        "scripts/compute_paper_numbers.py does not run, so the paper's claim "
        f"about it cannot be checked: {proc.stderr[-800:]}")
    out = proc.stdout
    tail = re.search(r"(\d+) numbers emitted: (\d+) recomputed here, "
                     r"(\d+) read from a recorded run", out)
    disc = re.search(r"DISCREPANCIES \((\d+)\)", out)
    assert tail and disc, (
        "the script's closing summary no longer has the shape this test reads, "
        "so the paper's figures cannot be checked against it")
    return tuple(int(g) for g in tail.groups()) + (int(disc.group(1)),)


def test_the_paper_states_the_counts_the_script_actually_emits():
    """The paper quotes four figures about its own provenance script.

    They were 389, 348, 41 and 25 while the script emitted 391, 346, 45 and 26.
    Those four were not invented: the script emits them at ce8242c, the commit
    the draft that quoted them accompanied, and they went stale at f82bc3a,
    which added checks. The script is deterministic, verified by three
    byte-identical consecutive runs at one commit, so a stale self-description
    is the easiest claim in a paper to check and the easiest to leave rotting.
    An earlier version of this test only asked whether each live figure
    appeared somewhere in the paper, which a paper stating both the current and
    the stale counts would have passed. It now reads the four numbers out of
    the sentence that describes the current version.
    """
    emitted, recomputed, recorded, discrepancies = _script_summary()
    paper = PAPER.read_text(encoding="utf-8")
    sentence = re.search(
        r"it emits (\d+) values, of which (\d+) are recomputed from primary data "
        r"and (\d+) are read from recorded runs and labelled as such, alongside "
        r"(\d+) failed cross-checks",
        paper)
    assert sentence, (
        "the paper no longer carries the sentence that states the provenance "
        "script's four counts in the expected form, so this guard cannot read "
        "them; restore the wording or rewrite this test deliberately")
    stated = tuple(int(g) for g in sentence.groups())
    assert stated == (emitted, recomputed, recorded, discrepancies), (
        f"the paper states {stated} as the script's current output and a live "
        f"run gives {(emitted, recomputed, recorded, discrepancies)}, so its "
        f"account of its own provenance is stale")


def _benchmark_rows():
    """The candidate table in the benchmark sheet, as cells."""
    rows = []
    for line in BENCHMARK_SHEET.read_text(encoding="utf-8").splitlines():
        if not line.startswith("|") or "`" not in line:
            continue
        cells = [c.strip().strip("*").strip() for c in line.strip("|").split("|")]
        seq = [c.strip("`") for c in cells if c.startswith("`")]
        if seq:
            rows.append((seq[0], cells))
    return rows


def _benchmark_json_value(seq, field):
    """The benchmark json's value for one sequence and field, or None.

    The file is a handful of named blocks rather than a table, so this walks it
    and returns the field from whichever block names this sequence.
    """
    path = REPO_ROOT / "demo_data" / "results" / "falsification_benchmark.json"
    if not path.exists():
        return None
    doc = json.loads(path.read_text(encoding="utf-8"))

    def walk(node):
        if isinstance(node, dict):
            if node.get("sequence") == seq and field in node:
                yield node[field]
            for value in node.values():
                yield from walk(value)
        elif isinstance(node, list):
            for value in node:
                yield from walk(value)

    for found in walk(doc):
        if isinstance(found, (int, float)):
            return float(found)
    return None


def test_the_benchmark_table_does_not_disagree_with_the_ledger():
    """Two published cells in that table came from nowhere.

    The linear control is not in the ledger at all, and its amylase score and
    selectivity difference were in no artefact. The naive polyarginine row
    carried an amylase score of -6.500 and a difference of -1.09 while the
    ledger recorded -7.580 and -0.006 for the same sequence. A table that
    disagrees with the run record is the 0.92 class: a number a reader cannot
    trace.
    """
    with LEDGER.open(newline="", encoding="utf-8") as fh:
        ledger = {r["sequence"]: r for r in csv.DictReader(fh)}
    columns = {2: "vina_3aic_kcal", 3: "vina_1smd_kcal", 4: "selectivity_ddg"}
    checked = 0
    for seq, cells in _benchmark_rows():
        row = ledger.get(seq)
        for idx, col in columns.items():
            if idx >= len(cells):
                continue
            cell = cells[idx].replace("$", "").replace("\\mathbf{", "").strip("}")
            value = re.search(r"-?\d+\.?\d*", cell)
            if row is None:
                # Not in the ledger. The benchmark json may still record the
                # field for this sequence, so check there before calling the
                # cell unsourced.
                recorded_elsewhere = _benchmark_json_value(seq, col)
                if recorded_elsewhere is not None:
                    assert value is not None and \
                        abs(float(value.group(0)) - recorded_elsewhere) <= 0.011, (
                            f"the benchmark table states {cells[idx]!r} for "
                            f"{seq} {col} and the benchmark json records "
                            f"{recorded_elsewhere}")
                    continue
                assert value is None or "not recorded" in cells[idx], (
                    f"the benchmark table publishes {cells[idx]!r} for {seq} "
                    f"in {col}, and no artefact in the repository records it: "
                    f"{seq} is not in the ledger and the benchmark json has no "
                    f"{col} for it")
                continue
            if value is None:
                continue
            stated, recorded = float(value.group(0)), float(row[col])
            assert abs(stated - recorded) <= 0.011, (
                f"the benchmark table states {stated} for {seq} {col} and the "
                f"ledger records {recorded}")
            checked += 1
    assert checked >= 12, (
        f"only {checked} cells were compared against the ledger, so this test "
        "is not covering the table it is meant to guard")


def test_the_documents_partition_the_ledger_the_way_the_ledger_does():
    """260 rows are 220 designs, one naive control and 39 scrambled controls.

    The paper called all 221 non-scrambled rows designed, which counts the
    naive control `naive_baseline_CRRRWWRRRRC` among the designs, while the
    README states the partition correctly. Both documents describe the same
    file, so they cannot disagree about what is in it.
    """
    import csv
    with (REPO_ROOT / "data" / "contenders_ledger.csv").open(encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    scrambled = [r for r in rows if r["falsification_status"] == "SCRAMBLED_NULL_CONTROL"]
    naive = [r for r in rows if r["contender_id"].startswith("naive_baseline")]
    designed = len(rows) - len(scrambled) - len(naive)
    assert (len(rows), len(scrambled), len(naive), designed) == (260, 39, 1, 220), (
        f"the ledger now holds {len(rows)} rows, {len(scrambled)} scrambled "
        f"controls, {len(naive)} naive controls and {designed} designs, so the "
        f"partition quoted in the paper and the README needs restating")
    paper = PAPER.read_text(encoding="utf-8")
    assert "221 designed" not in paper, (
        "the paper calls 221 rows designed, which counts the naive control "
        "among the designs; the non-scrambled set is 220 designs plus that "
        "control")
    readme = (REPO_ROOT / "README.md").read_text(encoding="utf-8")
    assert "220 design" in readme, (
        "the README no longer states the 220 design partition that the ledger "
        "holds")
