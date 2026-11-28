"""Tests on the data figures, so a figure cannot disagree with the manuscript.

A wrong figure is the thing a reader believes fastest, and the figures that
were here before were PNGs built inside a notebook that this environment cannot
run, so nothing checked them against the data at all.

Each figure now embeds a machine-readable provenance block holding every value
it plots. These tests parse that block and recompute each value independently,
from the ledger and the recorded runs, so the figure, the manuscript and the
data have to agree or the suite fails.
"""

import csv
import hashlib
import json
import math
import re
import statistics
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
RENDERER = REPO / "scripts" / "render_figures_svg.py"
FIGURES = REPO / "manuscript" / "figures"
LEDGER = REPO / "data" / "contenders_ledger.csv"
DISPROOF = REPO / "data" / "exhaustive_disproof_report.json"
POSES = REPO / "data" / "poses"

LEAD = "CAHWRPWHC"
PARENT = "CPHWRPWHC"
NULL_STATUS = "SCRAMBLED_NULL_CONTROL"

EXPECTED = ("fig_null_vs_lead.svg", "fig_alanine_scan.svg",
            "fig_affinity_vs_selectivity.svg")


def _rows():
    with LEDGER.open(newline="") as fh:
        return list(csv.DictReader(fh))


def _provenance(name):
    text = (FIGURES / name).read_text()
    match = re.search(r"<!-- figure-provenance\n(.*?)\n-->", text, re.S)
    assert match, f"{name} carries no provenance block"
    return json.loads(match.group(1))


def _digest(path):
    digest = hashlib.md5()
    for line in Path(path).read_bytes().splitlines(keepends=True):
        if line.startswith((b"ATOM", b"HETATM")):
            digest.update(line[30:54])
    return digest.hexdigest()


def _independent_repeats():
    angle = next(a for a in json.loads(DISPROOF.read_text())["angles"]
                 if a.get("vina_scores"))
    campaign = _digest(POSES / "contender_mutant_P2A_CAHWRPWHC_3aic.pdbqt")
    keep = []
    for seed, score in zip([42, 101, 202, 303, 404], angle["vina_scores"]):
        pose = POSES / f"disproof_seed_{seed}_CAHWRPWHC_3AIC_docked.pdbqt"
        if _digest(pose) != campaign:
            keep.append(score)
    return keep


# ------------------------------------------------------------- the files exist

@pytest.mark.parametrize("name", EXPECTED)
def test_each_figure_is_committed_and_is_svg(name):
    path = FIGURES / name
    assert path.exists(), f"{name} has not been rendered into manuscript/figures"
    text = path.read_text()
    assert text.startswith("<svg"), f"{name} is not an SVG document"
    assert text.rstrip().endswith("</svg>"), f"{name} is truncated"
    assert "<title>" in text, f"{name} has no accessible title"


def test_the_committed_figures_are_what_the_data_now_renders():
    """The renderer is deterministic, so a stale figure is a failure.

    This is the test that makes the figures trustworthy: if the ledger changes
    and nobody re-renders, the committed file no longer matches and the suite
    says so.
    """
    result = subprocess.run([sys.executable, str(RENDERER), "--check"],
                            capture_output=True, text=True, cwd=str(REPO))
    assert result.returncode == 0, (
        f"the committed figures no longer match the data they are drawn from. "
        f"Re-run scripts/render_figures_svg.py.\n{result.stdout}\n{result.stderr}")


@pytest.mark.parametrize("name", EXPECTED)
def test_each_figure_names_its_source_and_row_count(name):
    stamp = _provenance(name)
    assert stamp["source"] == "data/contenders_ledger.csv"
    assert stamp["source_rows"] == len(_rows()), (
        f"{name} was drawn from a different number of ledger rows than the "
        f"ledger now holds")
    assert stamp["source_sha256"] == hashlib.sha256(
        LEDGER.read_bytes()).hexdigest()[:16], (
        f"{name} was drawn from a different version of the ledger")


# ------------------------------------------------ figure 1, the null histogram

def test_the_null_histogram_plots_the_null_the_paper_reports():
    stamp = _provenance("fig_null_vs_lead.svg")
    rows = _rows()
    nulls = [float(r["vina_3aic_kcal"]) for r in rows
             if r["falsification_status"] == NULL_STATUS]
    lead = float(next(r for r in rows if r["sequence"] == LEAD)["vina_3aic_kcal"])
    parent = float(next(r for r in rows if r["sequence"] == PARENT)["vina_3aic_kcal"])

    assert stamp["null_n"] == len(nulls)
    assert stamp["lead"] == lead
    assert stamp["parent"] == parent
    assert stamp["null_mean"] == pytest.approx(statistics.fmean(nulls), abs=5e-6)
    assert stamp["null_sample_sd"] == pytest.approx(statistics.stdev(nulls), abs=5e-6)
    assert stamp["null_best"] == pytest.approx(min(nulls), abs=5e-6)
    assert stamp["nulls_at_least_as_good_as_lead"] == sum(
        1 for v in nulls if v <= lead)


def test_the_histogram_bars_account_for_every_control():
    stamp = _provenance("fig_null_vs_lead.svg")
    assert sum(stamp["bin_counts"]) == stamp["null_n"], (
        "the drawn bars do not add up to the number of controls, so the figure "
        "is dropping or double counting members of the null")
    assert len(stamp["bin_edges"]) == len(stamp["bin_counts"]) + 1
    nulls = [float(r["vina_3aic_kcal"]) for r in _rows()
             if r["falsification_status"] == NULL_STATUS]
    assert stamp["bin_edges"][0] <= min(nulls)
    assert stamp["bin_edges"][-1] >= max(nulls)


def test_the_histogram_shows_the_lead_inside_the_null_not_beyond_it():
    """The figure's whole point. If this ever flips, the prose must change too."""
    stamp = _provenance("fig_null_vs_lead.svg")
    assert stamp["lead"] > stamp["null_best"], (
        "the lead now scores better than every control, which would be a "
        "different result from the one the manuscript reports")
    assert stamp["nulls_at_least_as_good_as_lead"] > 0, (
        "no control now matches the lead; section 3.3b says 23 of 39 do")


# ------------------------------------------------- figure 2, the alanine scan

def test_the_alanine_figure_plots_the_deltas_the_ledger_gives():
    stamp = _provenance("fig_alanine_scan.svg")
    rows = _rows()
    parent = float(next(r for r in rows if r["sequence"] == PARENT)["vina_3aic_kcal"])
    muts = sorted((r for r in rows if "_mutant_" in r["contender_id"]),
                  key=lambda r: float(r["vina_3aic_kcal"]) - parent)
    labels = [r["contender_id"].split("_mutant_")[1].split("_")[0] for r in muts]
    deltas = [float(r["vina_3aic_kcal"]) - parent for r in muts]

    assert stamp["parent"] == parent
    assert stamp["positions"] == labels
    assert stamp["deltas"] == pytest.approx(deltas, abs=5e-6)
    assert stamp["rows_plotted"] == len(muts)


def test_the_alanine_figure_agrees_with_the_recorded_scan():
    """The ledger and the recorded run have to say the same thing."""
    stamp = _provenance("fig_alanine_scan.svg")
    angle = next(a for a in json.loads(DISPROOF.read_text())["angles"]
                 if a.get("mutants"))
    recorded = {k: v["delta_ddg"] for k, v in angle["mutants"].items()}
    drawn = dict(zip(stamp["positions"], stamp["deltas"]))
    assert set(drawn) == set(recorded), (
        f"the figure draws {sorted(drawn)} but the run recorded {sorted(recorded)}")
    for position, value in recorded.items():
        assert drawn[position] == pytest.approx(value, abs=1e-3), (
            f"{position}: figure {drawn[position]:.4f}, recorded run {value}")
    assert stamp["hotspot_gate"] == pytest.approx(angle["hotspot_cutoff_kcal"])
    assert stamp["largest_loss"] == pytest.approx(
        angle["largest_affinity_loss_kcal"], abs=1e-3)


def test_the_alanine_figure_reports_no_hotspot():
    stamp = _provenance("fig_alanine_scan.svg")
    assert stamp["positions_reaching_gate"] == 0, (
        "a position now reaches the hotspot threshold, which contradicts the "
        "manuscript's finding and the figure's own headline")
    assert max(stamp["deltas"]) < stamp["hotspot_gate"]


def test_the_alanine_noise_band_is_the_four_run_spread():
    stamp = _provenance("fig_alanine_scan.svg")
    repeats = _independent_repeats()
    assert stamp["repeat_runs_used"] == len(repeats)
    assert stamp["run_to_run_sample_sd"] == pytest.approx(
        statistics.stdev(repeats), abs=5e-6)
    assert stamp["lead_is_position"] == "P2A", (
        "the lead is the P2A substitution of the parent, and the figure should "
        "say so: a reader otherwise takes them for two separate molecules")


# ------------------------------------- figure 3, affinity against selectivity

def test_the_scatter_plots_every_row_once():
    stamp = _provenance("fig_affinity_vs_selectivity.svg")
    assert stamp["rows_plotted"] == len(_rows())


def test_the_dominance_claim_on_the_scatter_is_computed_not_asserted():
    stamp = _provenance("fig_affinity_vs_selectivity.svg")
    rows = _rows()
    pts = [(float(r["vina_3aic_kcal"]), float(r["selectivity_ddg"]), r) for r in rows]
    frontier = [r for x, y, r in pts
                if not any(a <= x and b <= y and (a < x or b < y)
                           for a, b, _ in pts)]
    assert stamp["non_dominated_count"] == len(frontier)
    assert stamp["non_dominated_sequences"] == sorted(r["sequence"] for r in frontier)
    assert stamp["lead_is_non_dominated"] == any(
        r["sequence"] == LEAD for r in frontier)

    lead_x, lead_y = next((x, y) for x, y, r in pts if r["sequence"] == LEAD)
    assert stamp["candidates_better_on_both_axes"] == sum(
        1 for x, y, _ in pts if x < lead_x and y < lead_y)
    assert stamp["candidates_better_on_affinity_alone"] == sum(
        1 for x, _, _ in pts if x < lead_x)


def test_the_pareto_rank_column_is_not_a_dominance_rank():
    """Why the figure computes dominance instead of reading the column.

    The column is written as 1 for a CORROBORATED row and 2 otherwise, so it
    restates the status. Pinning the disagreement keeps the figure's caption
    honest, and if the column is ever fixed this test says so.
    """
    stamp = _provenance("fig_affinity_vs_selectivity.svg")
    assert stamp["rows_carrying_pareto_rank_1"] != stamp["non_dominated_count"], (
        "pareto_rank now agrees with the dominance relation; the figure's note "
        "saying it does not should be removed")


# --------------------------------------------- the figure that was taken out

def test_the_acid_switch_figure_is_gone_and_its_removal_is_explained():
    """Both neutral arms of that result reuse a file of unestablished origin.

    It is removed rather than redrawn, and the README has to say so, because a
    reader who remembers the figure needs to know why it went.
    """
    assert not (FIGURES / "fig4_stephan_acid_switch.png").exists(), (
        "the acid-switch figure is back; both of its neutral arms reuse the "
        "duplicated pose file, so it cannot be drawn from this data")
    readme = (REPO / "README.md").read_text()
    assert "fig4_stephan_acid_switch" in readme, (
        "the README does not mention the removed acid-switch figure, so its "
        "disappearance is unexplained")


def test_the_live_figure_directory_holds_no_handmade_png():
    """Named for what it checks.

    This was called test_no_superseded_png_figure_is_still_committed, which
    claims more than it does: the superseded PNGs are still committed, on
    purpose, under manuscript/superseded/, and the README says so. What this
    checks is narrower and is the thing worth holding: the live figure
    directory has exactly one source, the ledger.
    """
    stale = sorted(p.name for p in FIGURES.glob("*.png"))
    assert not stale, (
        f"these figures are still committed but are not generated from the "
        f"data: {stale}. Either render them from the ledger or remove them.")
