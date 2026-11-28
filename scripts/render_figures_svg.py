#!/usr/bin/env python3
"""Draw the repository's data figures as SVG, from the data, with no plotting library.

Why this exists. The previous figures were PNGs produced inside a notebook that
needs matplotlib, seaborn and pandas, so they could not be regenerated in the
build and had quietly drifted from the data behind them. These are drawn from
the committed ledger every time, in the standard library alone, so a figure
cannot say something the data does not.

Three rules this file keeps.

1. Nothing is typed in. Every plotted value is read from data/contenders_ledger.csv
   at render time. The only literals are the two stated gates, which are design
   decisions rather than measurements, and they are labelled as such.
2. Every figure embeds its own provenance: the source file, its row count, the
   rows selected, and the values plotted, as a machine-readable comment. The
   tests in tests/test_figures.py read that block and compare it with what
   scripts/compute_paper_numbers.py recomputes, so a stale figure fails the suite.
3. No gradients, no shadows, no decorative colour. Light background, dark text.

Usage:  python3 scripts/render_figures_svg.py [--out DIR] [--check]
        --check renders in memory and reports whether the committed files match.
"""

import argparse
import csv
import hashlib
import html
import json
import math
import statistics
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
LEDGER = REPO / "data" / "contenders_ledger.csv"
DISPROOF = REPO / "data" / "exhaustive_disproof_report.json"
POSES = REPO / "data" / "poses"
CAMPAIGN_POSE = POSES / "contender_mutant_P2A_CAHWRPWHC_3aic.pdbqt"
DEFAULT_OUT = REPO / "manuscript" / "figures"

LEAD = "CAHWRPWHC"
PARENT = "CPHWRPWHC"
NULL_STATUS = "SCRAMBLED_NULL_CONTROL"

# The affinity gate is a design decision rather than a measurement, and it is
# written as a literal in the campaign scripts, so it is quoted from there. The
# hotspot threshold is read out of the recorded run instead, since the run
# records it. Nothing else here is a constant.
AFFINITY_GATE = -7.5              # scripts/run_campaign.py, "affinity_ok"

INK = "#1b1b1b"
MUTED = "#6b6b6b"
GRID = "#d8d8d8"
MARK = "#20506e"
ACCENT = "#9c2f2f"


# ----------------------------------------------------------------- SVG canvas

class Canvas:
    """A sheet of SVG with a few primitives. Coordinates are SVG user units."""

    def __init__(self, width, height, title, desc=""):
        self.w, self.h = width, height
        self.title, self.desc = title, desc
        self.body = []
        self.provenance = {}

    def text(self, x, y, s, size=11, anchor="start", weight="normal",
             fill=INK, rotate=None):
        tr = f' transform="rotate({rotate} {x:.1f} {y:.1f})"' if rotate else ""
        self.body.append(
            f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" '
            f'text-anchor="{anchor}" font-weight="{weight}" fill="{fill}" '
            f'font-family="DejaVu Sans, Helvetica, Arial, sans-serif"{tr}>'
            f'{html.escape(str(s))}</text>')

    def label(self, x, y, s, size=10, anchor="start", fill=INK, pad=2.0):
        """Text on a white backing box, so an annotation stays readable over marks.

        The box is sized from an average glyph width rather than real metrics,
        which is why it is slightly transparent: a box a little wide or narrow
        then reads as a halo instead of a misplaced panel.
        """
        s = str(s)
        w = 0.55 * size * len(s)
        x0 = x if anchor == "start" else (x - w if anchor == "end" else x - w / 2.0)
        self.body.append(
            f'<rect x="{x0 - pad:.1f}" y="{y - size * 0.82 - pad:.1f}" '
            f'width="{w + 2 * pad:.1f}" height="{size + 2 * pad:.1f}" '
            f'fill="#ffffff" fill-opacity="0.85" stroke="none"/>')
        self.text(x, y, s, size=size, anchor=anchor, fill=fill)

    def line(self, x1, y1, x2, y2, stroke=INK, width=1.0, dash=None):
        d = f' stroke-dasharray="{dash}"' if dash else ""
        self.body.append(
            f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
            f'stroke="{stroke}" stroke-width="{width}"{d}/>')

    def rect(self, x, y, w, h, fill="none", stroke=INK, width=1.0):
        self.body.append(
            f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" '
            f'fill="{fill}" stroke="{stroke}" stroke-width="{width}"/>')

    def circle(self, cx, cy, r, fill="none", stroke=INK, width=0.8):
        self.body.append(
            f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="{r:.2f}" fill="{fill}" '
            f'stroke="{stroke}" stroke-width="{width}"/>')

    def render(self):
        out = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{self.w}" '
               f'height="{self.h}" viewBox="0 0 {self.w} {self.h}" role="img" '
               f'aria-label="{html.escape(self.title)}">',
               f'<title>{html.escape(self.title)}</title>']
        if self.desc:
            out.append(f'<desc>{html.escape(self.desc)}</desc>')
        block = json.dumps(self.provenance, indent=2, sort_keys=True)
        assert "--" not in block, "an SVG comment cannot contain a double hyphen"
        out.append("<!-- figure-provenance\n" + block + "\n-->")
        out.append(f'<rect width="{self.w}" height="{self.h}" fill="#ffffff"/>')
        out += self.body
        out.append("</svg>")
        return "\n".join(out) + "\n"


def auto_bins(values):
    """Bin count by numpy's auto rule: the larger of Freedman-Diaconis and Sturges.

    A bin count picked by eye is a free parameter, and at 39 points a fine
    binning invents structure: at 14 bins this histogram shows two peaks that
    are not in the data. The rule is recorded in the figure's provenance so a
    reader can check it rather than take the shape on trust.
    """
    xs = sorted(values)
    n = len(xs)
    if n < 2:
        return 1
    span = xs[-1] - xs[0]
    if span <= 0:
        return 1

    def quantile(p):
        h = (n - 1) * p
        lo = int(math.floor(h))
        return xs[lo] + (h - lo) * (xs[min(lo + 1, n - 1)] - xs[lo])

    iqr = quantile(0.75) - quantile(0.25)
    fd = math.ceil(span / (2 * iqr * n ** (-1 / 3))) if iqr > 0 else 1
    sturges = math.ceil(math.log2(n)) + 1
    return max(1, fd, sturges)


def count_ticks(top, target=5):
    """Whole-number ticks from 0 to at least `top`, at most about `target` of them."""
    step = max(1, int(math.ceil(top / float(target))))
    for candidate in (1, 2, 5, 10, 20, 25, 50, 100):
        if candidate >= step:
            step = candidate
            break
    return list(range(0, int(top) + step, step))


def nice_ticks(lo, hi, target=6):
    """Round tick positions spanning [lo, hi] inclusive."""
    if hi <= lo:
        hi = lo + 1.0
    raw = (hi - lo) / float(target)
    mag = 10.0 ** math.floor(math.log10(raw))
    step = mag * 10
    for mult in (1, 2, 2.5, 5):
        if raw <= mag * mult:
            step = mag * mult
            break
    ticks, v = [], math.floor(lo / step) * step
    while v < hi - step * 1e-9:
        ticks.append(round(v, 10))
        v += step
    ticks.append(round(v, 10))
    if ticks[-1] < hi - step * 1e-9:
        ticks.append(round(ticks[-1] + step, 10))
    return ticks


class Axes:
    """A plot frame with linear scales and labelled ticks."""

    def __init__(self, canvas, left, top, width, height,
                 xlim, ylim, xlabel, ylabel, xticks=None, yticks=None):
        self.c = canvas
        self.left, self.top, self.width, self.height = left, top, width, height
        self.x0, self.x1 = xlim
        self.y0, self.y1 = ylim

        canvas.rect(left, top, width, height, stroke=GRID, width=1.0)
        for xv in (xticks if xticks is not None else nice_ticks(*xlim)):
            if not self.x0 - 1e-9 <= xv <= self.x1 + 1e-9:
                continue
            px = self.px(xv)
            canvas.line(px, top, px, top + height, stroke=GRID, width=0.6)
            canvas.line(px, top + height, px, top + height + 4, stroke=MUTED)
            canvas.text(px, top + height + 16, _fmt(xv), size=10,
                        anchor="middle", fill=MUTED)
        for yv in (yticks if yticks is not None else nice_ticks(*ylim)):
            if not self.y0 - 1e-9 <= yv <= self.y1 + 1e-9:
                continue
            py = self.py(yv)
            canvas.line(left, py, left + width, py, stroke=GRID, width=0.6)
            canvas.line(left - 4, py, left, py, stroke=MUTED)
            canvas.text(left - 8, py + 3.5, _fmt(yv), size=10,
                        anchor="end", fill=MUTED)
        canvas.text(left + width / 2.0, top + height + 36, xlabel, size=11,
                    anchor="middle")
        canvas.text(left - 44, top + height / 2.0, ylabel, size=11,
                    anchor="middle", rotate=-90)

    def px(self, xv):
        return self.left + (xv - self.x0) / (self.x1 - self.x0) * self.width

    def py(self, yv):
        return self.top + self.height - (yv - self.y0) / (self.y1 - self.y0) * self.height


def _display(path):
    """Repository-relative where possible, absolute otherwise. Never raises."""
    resolved = Path(path).resolve()
    try:
        return str(resolved.relative_to(REPO))
    except ValueError:
        return str(resolved)


def _fmt(v):
    if abs(v - round(v)) < 1e-9:
        return str(int(round(v)))
    s = f"{v:.2f}".rstrip("0").rstrip(".")
    return s if s not in ("-0", "") else "0"


# ---------------------------------------------------------------------- data

def read_ledger():
    with LEDGER.open(newline="") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        raise SystemExit(f"{LEDGER} has no rows")
    return rows


def _f(row, key):
    return float(row[key])


def null_scores(rows):
    return [_f(r, "vina_3aic_kcal") for r in rows
            if r["falsification_status"] == NULL_STATUS]


def one(rows, sequence):
    hits = [r for r in rows if r["sequence"] == sequence]
    if len(hits) != 1:
        raise SystemExit(f"expected one row for {sequence}, found {len(hits)}")
    return hits[0]


def coordinate_digest(path):
    """md5 over the coordinate columns alone, so REMARK lines cannot change it."""
    digest = hashlib.md5()
    for line in Path(path).read_bytes().splitlines(keepends=True):
        if line.startswith((b"ATOM", b"HETATM")):
            digest.update(line[30:54])
    return digest.hexdigest()


def disproof_angles():
    return json.loads(DISPROOF.read_text())["angles"]


def hotspot_gate():
    """The threshold the alanine scan itself recorded."""
    for angle in disproof_angles():
        if "hotspot_cutoff_kcal" in angle:
            return float(angle["hotspot_cutoff_kcal"])
    raise SystemExit("no recorded hotspot threshold in the disproof report")


def independent_repeat_scores():
    """Repeat runs of one ligand, less the run that reports the campaign pose.

    One of the five nominal seed runs carries a docked pose byte-identical to
    the campaign pose, so it replicates nothing and is dropped here. The test
    for that is the pose, not the score, so this stays correct if the run is
    ever redone.
    """
    seeds = [42, 101, 202, 303, 404]
    for angle in disproof_angles():
        if not angle.get("vina_scores"):
            continue
        scores = angle["vina_scores"]
        if len(scores) != len(seeds):
            raise SystemExit("the recorded sweep no longer has one score per seed")
        campaign = coordinate_digest(CAMPAIGN_POSE)
        keep = []
        for seed, score in zip(seeds, scores):
            pose = POSES / f"disproof_seed_{seed}_CAHWRPWHC_3AIC_docked.pdbqt"
            if coordinate_digest(pose) != campaign:
                keep.append(score)
        return keep
    raise SystemExit("no recorded seed sweep in the disproof report")


def non_dominated(rows):
    """Rows no other row beats on affinity and selectivity at once.

    The ledger carries a pareto_rank column, but it is written as 1 for a
    CORROBORATED row and 2 otherwise, so it restates the status rather than
    ranking anything. This computes the dominance relation the name implies.
    """
    pts = [(_f(r, "vina_3aic_kcal"), _f(r, "selectivity_ddg"), r) for r in rows]
    return [r for x, y, r in pts
            if not any(a <= x and b <= y and (a < x or b < y) for a, b, _ in pts)]


def source_stamp(rows, selected, note):
    return {"source": str(LEDGER.relative_to(REPO)),
            "source_sha256": hashlib.sha256(LEDGER.read_bytes()).hexdigest()[:16],
            "source_rows": len(rows),
            "rows_plotted": selected,
            "generated_by": "scripts/render_figures_svg.py",
            "note": note}


# ------------------------------------------------------------------- figures

def fig_null_vs_lead(rows):
    """The central figure: where the lead falls inside its own null distribution."""
    nulls = null_scores(rows)
    lead = _f(one(rows, LEAD), "vina_3aic_kcal")
    parent = _f(one(rows, PARENT), "vina_3aic_kcal")
    mean = statistics.fmean(nulls)
    sd = statistics.stdev(nulls)
    at_least_as_good = sum(1 for v in nulls if v <= lead)

    lo, hi = math.floor(min(nulls) * 4) / 4.0, math.ceil(max(nulls) * 4) / 4.0
    nbins = auto_bins(nulls)
    edges = [lo + (hi - lo) * i / nbins for i in range(nbins + 1)]
    counts = [0] * nbins
    for v in nulls:
        idx = min(int((v - lo) / (hi - lo) * nbins), nbins - 1)
        counts[idx] += 1

    c = Canvas(760, 430,
               "Scrambled null distribution with the lead marked",
               "Histogram of 39 composition-matched scrambled controls against "
               "GTF-SI, with the lead and its parent marked in place.")
    c.text(40, 30, "The lead sits inside its own null distribution", size=15,
           weight="bold")
    c.text(40, 50, f"39 scrambled controls, mean {mean:.3f} kcal/mol, "
                   f"sample SD {sd:.3f}. {at_least_as_good} of 39 score at least "
                   f"as well as the lead.", size=11, fill=MUTED)

    ax = Axes(c, 92, 74, 620, 268, (lo, hi), (0, max(counts) + 1),
              "AutoDock Vina affinity against GTF-SI (kcal/mol), lower is better",
              "Scrambled controls",
              yticks=count_ticks(max(counts) + 1))
    bw = (ax.px(edges[1]) - ax.px(edges[0]))
    for i, n in enumerate(counts):
        if not n:
            continue
        x = ax.px(edges[i])
        y = ax.py(n)
        c.rect(x + 1, y, bw - 2, ax.py(0) - y, fill="#c9d8e4", stroke=MARK, width=0.8)

    for value, label, colour in ((lead, f"lead {LEAD} {lead:.3f}", ACCENT),
                                 (parent, f"parent {PARENT} {parent:.3f}", "#4a4a4a")):
        x = ax.px(value)
        c.line(x, ax.top, x, ax.py(0), stroke=colour, width=1.6)
    c.line(ax.px(mean), ax.top, ax.px(mean), ax.py(0), stroke=MUTED, width=1.2,
           dash="4 3")

    # Labels placed below the frame, so no line of text can cover a bar.
    ly = ax.top + ax.height + 52
    c.line(92, ly - 4, 106, ly - 4, stroke=ACCENT, width=1.6)
    c.text(112, ly, f"lead {LEAD} at {lead:.3f}", size=10)
    c.line(300, ly - 4, 314, ly - 4, stroke="#4a4a4a", width=1.6)
    c.text(320, ly, f"parent {PARENT} at {parent:.3f}", size=10)
    c.line(520, ly - 4, 534, ly - 4, stroke=MUTED, width=1.2, dash="4 3")
    c.text(540, ly, f"null mean {mean:.3f}", size=10)
    c.text(40, 412, f"Drawn from data/contenders_ledger.csv, {len(rows)} rows, "
                    f"{len(nulls)} of them scrambled controls.", size=9, fill=MUTED)

    c.provenance = source_stamp(rows, len(nulls),
                                "one-sided permutation comparison of the lead "
                                "against its scrambled controls")
    c.provenance.update({
        "null_n": len(nulls), "null_mean": round(mean, 6),
        "null_sample_sd": round(sd, 6), "null_best": round(min(nulls), 6),
        "lead": lead, "parent": parent,
        "nulls_at_least_as_good_as_lead": at_least_as_good,
        "bin_edges": [round(e, 6) for e in edges], "bin_counts": counts,
        "binning": "the larger of the Freedman-Diaconis and Sturges bin counts, "
                   "which is numpy's auto rule; not a number chosen by eye"})
    return "fig_null_vs_lead.svg", c


def fig_alanine_scan(rows):
    """Every single-position alanine substitution, against the gate it had to clear."""
    parent_row = one(rows, PARENT)
    parent = _f(parent_row, "vina_3aic_kcal")
    muts = sorted((r for r in rows if "_mutant_" in r["contender_id"]),
                  key=lambda r: _f(r, "vina_3aic_kcal") - parent)
    labels = [r["contender_id"].split("_mutant_")[1].split("_")[0] for r in muts]
    deltas = [_f(r, "vina_3aic_kcal") - parent for r in muts]
    gate = hotspot_gate()
    repeats = independent_repeat_scores()
    noise_sd = statistics.stdev(repeats)
    lead_position = next(
        (lab for lab, r in zip(labels, muts) if r["sequence"] == LEAD), None)

    lo = min(min(deltas), -noise_sd) - 0.15
    hi = max(max(deltas), gate) + 0.15
    c = Canvas(760, 430, "Single-position alanine scan against its own gate",
               "Change in Vina affinity for each alanine substitution of the "
               "parent peptide, against the stated hotspot threshold.")
    c.text(40, 30, "No position reaches the hotspot threshold the scan set",
           size=15, weight="bold")
    subtitle = (f"Change against the parent {PARENT} at {parent:.3f} kcal/mol. "
                f"Positive is a loss of affinity, which is what a hotspot means.")
    if lead_position:
        subtitle += f" {lead_position} is the lead, {LEAD}."
    c.text(40, 50, subtitle, size=11, fill=MUTED)

    ax = Axes(c, 92, 74, 620, 268, (lo, hi), (-0.5, len(muts) - 0.5),
              "Change in Vina affinity against the parent (kcal/mol)", "",
              yticks=[])
    band_lo, band_hi = ax.px(-noise_sd), ax.px(noise_sd)
    c.rect(band_lo, ax.top, band_hi - band_lo, ax.height, fill="#f0f0f0",
           stroke="none", width=0)
    c.line(ax.px(0), ax.top, ax.px(0), ax.top + ax.height, stroke=MUTED, width=1.0)
    gx = ax.px(gate)
    c.line(gx, ax.top, gx, ax.top + ax.height, stroke=ACCENT, width=1.4, dash="5 3")
    c.text(gx - 6, ax.top + 14, f"hotspot gate +{gate:.2f}", size=10,
           anchor="end", fill=ACCENT)

    bh = ax.height / len(muts) * 0.55
    for i, (lab, d) in enumerate(zip(labels, deltas)):
        cy = ax.py(len(muts) - 1 - i)
        is_lead = lab == lead_position
        x = min(ax.px(0), ax.px(d))
        w = abs(ax.px(d) - ax.px(0))
        c.rect(x, cy - bh / 2.0, max(w, 0.8), bh,
               fill="#e3cccc" if is_lead else "#c9d8e4",
               stroke=ACCENT if is_lead else MARK, width=0.8)
        c.text(ax.left - 8, cy + 3.5, lab, size=10, anchor="end",
               fill=ACCENT if is_lead else INK)
        tx = ax.px(d) + (6 if d >= 0 else -6)
        c.text(tx, cy + 3.5, f"{d:+.3f}", size=9,
               anchor="start" if d >= 0 else "end", fill=MUTED)

    c.text(40, 388, "The grey band is plus or minus one run-to-run standard "
                    f"deviation of a single ligand ({noise_sd:.3f} kcal/mol over "
                    f"{len(repeats)} repeat runs), so every bar inside it is "
                    "smaller than the scatter of one molecule redocked.",
           size=9, fill=MUTED)
    c.text(40, 412, f"Drawn from data/contenders_ledger.csv, {len(rows)} rows, "
                    f"{len(muts)} alanine substitutions.", size=9, fill=MUTED)

    c.provenance = source_stamp(rows, len(muts),
                                "alanine scan deltas against the parent peptide")
    c.provenance.update({
        "parent": parent, "positions": labels,
        "deltas": [round(d, 6) for d in deltas],
        "hotspot_gate": gate,
        "hotspot_gate_source": "data/exhaustive_disproof_report.json",
        "positions_reaching_gate": sum(1 for d in deltas if d >= gate),
        "largest_loss": round(max(deltas), 6),
        "substitutions_that_improved": sum(1 for d in deltas if d < 0),
        "lead_is_position": lead_position,
        "repeat_runs_used": len(repeats),
        "run_to_run_sample_sd": round(noise_sd, 6),
        "deltas_inside_the_run_to_run_band": sum(
            1 for d in deltas if abs(d) <= noise_sd)})
    return "fig_alanine_scan.svg", c


def fig_affinity_vs_selectivity(rows):
    """The leaderboard's two axes, with the lead's place in it shown honestly."""
    pts = [(_f(r, "vina_3aic_kcal"), _f(r, "selectivity_ddg"), r) for r in rows]
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    lead_row, parent_row = one(rows, LEAD), one(rows, PARENT)
    lead_x, lead_y = _f(lead_row, "vina_3aic_kcal"), _f(lead_row, "selectivity_ddg")
    better_both = [r for x, y, r in pts if x < lead_x and y < lead_y]
    frontier = non_dominated(rows)
    lead_on_frontier = any(r["sequence"] == LEAD for r in frontier)

    c = Canvas(760, 470,
               "Affinity against selectivity for every evaluated candidate",
               "Each candidate's GTF-SI affinity against its GTF-SI minus "
               "alpha-amylase selectivity, with the lead marked.")
    headline = ("The lead is not on the frontier of its own leaderboard"
                if not lead_on_frontier else
                "The lead is one of the non-dominated candidates")
    c.text(40, 30, headline, size=15, weight="bold")
    c.text(40, 50, f"{len(better_both)} of {len(rows)} candidates score better than "
                   f"the lead on both axes at once. {len(frontier)} are dominated by "
                   f"nothing in the set.", size=11, fill=MUTED)

    ax = Axes(c, 92, 74, 620, 300,
              (math.floor(min(xs) * 2) / 2.0, math.ceil(max(xs) * 2) / 2.0),
              (math.floor(min(ys) * 2) / 2.0, math.ceil(max(ys) * 2) / 2.0),
              "AutoDock Vina affinity against GTF-SI (kcal/mol), lower is better",
              "Selectivity, GTF-SI minus alpha-amylase (kcal/mol)")
    on_frontier = {id(r) for r in frontier}
    for x, y, r in pts:
        if r["sequence"] in (LEAD, PARENT):
            continue
        is_null = r["falsification_status"] == NULL_STATUS
        if id(r) in on_frontier:
            c.circle(ax.px(x), ax.py(y), 4.2, fill="none", stroke=INK, width=1.3)
            continue
        c.circle(ax.px(x), ax.py(y), 2.6,
                 fill="none" if is_null else "#c9d8e4",
                 stroke=MUTED if is_null else MARK, width=0.7)
    gx = ax.px(AFFINITY_GATE)
    c.line(gx, ax.top, gx, ax.top + ax.height, stroke=ACCENT, width=1.2, dash="5 3")
    c.label(gx + 6, ax.top + 14, f"affinity gate {AFFINITY_GATE}", size=10,
            fill=ACCENT)
    for row, label, colour in ((lead_row, f"lead {LEAD}", ACCENT),
                               (parent_row, f"parent {PARENT}", "#4a4a4a")):
        x, y = _f(row, "vina_3aic_kcal"), _f(row, "selectivity_ddg")
        c.circle(ax.px(x), ax.py(y), 5.0, fill="none", stroke=colour, width=1.8)
        # Both sit in the densest part of the cloud, so the backing box is what
        # keeps the annotation readable without moving the mark off its value.
        c.label(ax.px(x) + 9, ax.py(y) + 3.5, label, size=10, fill=colour)

    ly = ax.top + ax.height + 52
    c.circle(98, ly - 4, 2.6, fill="#c9d8e4", stroke=MARK, width=0.7)
    c.text(108, ly, "designed candidate", size=10)
    c.circle(250, ly - 4, 2.6, fill="none", stroke=MUTED, width=0.7)
    c.text(260, ly, "scrambled null control", size=10)
    c.circle(430, ly - 4, 4.2, fill="none", stroke=INK, width=1.3)
    c.text(440, ly, "dominated by nothing in the set", size=10)
    c.text(40, 434, "Points are not coloured by the ledger's falsification status, "
                    "because 36 of the 39 rows it labels SUBTHRESHOLD_AFFINITY clear "
                    "the affinity gate and fail a different criterion.",
           size=9, fill=MUTED)
    c.text(40, 452, f"Drawn from data/contenders_ledger.csv, all {len(rows)} rows. "
                    "Dominance is computed here, not read from the pareto_rank "
                    "column, which restates the status.", size=9, fill=MUTED)

    c.provenance = source_stamp(rows, len(rows),
                                "affinity against selectivity for every row")
    c.provenance.update({
        "lead_affinity": lead_x, "lead_selectivity": lead_y,
        "parent_affinity": _f(parent_row, "vina_3aic_kcal"),
        "parent_selectivity": _f(parent_row, "selectivity_ddg"),
        "candidates_better_on_both_axes": len(better_both),
        "candidates_better_on_affinity_alone": sum(1 for x, _, _ in pts if x < lead_x),
        "non_dominated_count": len(frontier),
        "non_dominated_sequences": sorted(r["sequence"] for r in frontier),
        "lead_is_non_dominated": lead_on_frontier,
        "rows_carrying_pareto_rank_1": sum(1 for r in rows if r["pareto_rank"] == "1"),
        "affinity_gate": AFFINITY_GATE,
        "affinity_gate_source": "scripts/run_campaign.py"})
    return "fig_affinity_vs_selectivity.svg", c


FIGURES = (fig_null_vs_lead, fig_alanine_scan, fig_affinity_vs_selectivity)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default=str(DEFAULT_OUT))
    ap.add_argument("--check", action="store_true",
                    help="do not write; report whether committed files match")
    args = ap.parse_args(argv)

    rows = read_ledger()
    out = Path(args.out)
    stale = []
    for builder in FIGURES:
        name, canvas = builder(rows)
        markup = canvas.render()
        target = out / name
        if args.check:
            current = target.read_text() if target.exists() else None
            state = "matches" if current == markup else "STALE"
            if current != markup:
                stale.append(name)
            print(f"{name}: {state}")
        else:
            out.mkdir(parents=True, exist_ok=True)
            target.write_text(markup)
            print(f"wrote {_display(target)} ({len(markup)} bytes)")
    if args.check:
        if stale:
            print(f"\nResult: FAIL. {len(stale)} figure(s) differ from the data: "
                  + ", ".join(stale))
            return 1
        print("\nResult: pass. Every figure matches the data it is drawn from.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
