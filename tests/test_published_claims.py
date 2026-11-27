"""Tests on the claims the public pages make, as opposed to the paper's own prose.

These exist because the five-seed reproducibility claim was live on the GitHub
Pages site and in the README while the artefacts behind it showed four runs, not
five: the docked pose attributed to seed 42 is byte-identical to the campaign
pose across all nine output modes, so it is not an independent replicate.

Two things are pinned here.

1. The claim cannot come back. A file may refer to the old N=5 figure as
   history, which a reader of the old page needs, but it may not assert five
   independent seeds, and it may not mention N=5 without the exclusion nearby.
2. The replacement figures are recomputed from the artefacts on every run, and
   the exclusion itself is derived from the pose files rather than from a score
   typed into this file. If someone re-dockes seed 42 properly, these tests
   change their expectations on their own and fail only the published text.
"""

import hashlib
import itertools
import json
import re
import math
import statistics
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
DISPROOF = REPO / "data" / "exhaustive_disproof_report.json"
POSES = REPO / "data" / "poses"
CAMPAIGN_POSE = POSES / "contender_mutant_P2A_CAHWRPWHC_3aic.pdbqt"

# The files a reader of the project actually sees.
PUBLISHED = ("docs/index.html", "README.md", "dashboard.py",
             "scripts/build_educational_sops.py")

# Strings that assert five independent seeds. Any of them is the defect itself.
BANNED = ("N=5 independent",
          "(42, 101, 202, 303, 404)",
          "Multi-Seed Conformer & Docking (N=5)",
          "Multi-Seed Reproducibility (N=5)",
          "Multi-Seed N=5",
          "N=5 seeds")


def _coordinate_digest(path):
    """md5 over the coordinate columns of every atom record, all modes.

    Only columns 31 to 54, so provenance REMARK lines never change the digest
    and two files agree here only if their geometry agrees to the last decimal.
    """
    digest = hashlib.md5()
    for line in Path(path).read_bytes().splitlines(keepends=True):
        if line.startswith((b"ATOM", b"HETATM")):
            digest.update(line[30:54])
    return digest.hexdigest()


# Meeko glue pseudo-atoms close an opened macrocycle. They are not atoms, and
# counting them inflates every distance: it is the difference between the
# published 9.24 A and the 9.11 A this project now reports.
GLUE_TYPES = {"G", "G0", "G1", "G2", "G3"}


def _mode_one_heavy_atoms(path):
    """Real heavy atoms of the first output mode, in file order.

    Hydrogens out, glue pseudo-atoms out.
    """
    coords, started = [], False
    for line in Path(path).read_text().splitlines():
        if line.startswith("MODEL"):
            if started:
                break
            started = True
        if line.startswith(("ATOM", "HETATM")):
            autodock_type = line[76:79].strip().upper()
            element = autodock_type[:2].strip() or line[12:16].strip()[0]
            if element.upper().startswith("H"):
                continue
            if autodock_type in GLUE_TYPES:
                continue
            coords.append((float(line[30:38]), float(line[38:46]),
                           float(line[46:54])))
    return coords


def _seed_runs():
    """The nominal seed sweep, with each run's score and docked pose."""
    angles = json.loads(DISPROOF.read_text())["angles"]
    angle = next(a for a in angles if a.get("vina_scores"))
    seeds = [42, 101, 202, 303, 404]
    assert len(angle["vina_scores"]) == len(seeds), (
        "the recorded sweep no longer holds one score per nominal seed; this "
        "test derives the exclusion from the poses and needs the pairing")
    runs = []
    for seed, score in zip(seeds, angle["vina_scores"]):
        pose = POSES / f"disproof_seed_{seed}_CAHWRPWHC_3AIC_docked.pdbqt"
        runs.append({"seed": seed, "score": score, "pose": pose,
                     "digest": _coordinate_digest(pose)})
    return runs


def _independent_runs():
    """The runs whose docked pose is not the campaign pose.

    This is the whole basis of the correction, so it is computed here rather
    than asserted: a run that merely reports the campaign's geometry is not a
    replicate of anything, whatever seed it is filed under.
    """
    campaign = _coordinate_digest(CAMPAIGN_POSE)
    return [r for r in _seed_runs() if r["digest"] != campaign]


def _pairwise_rmsd(poses):
    out = []
    for a, b in itertools.combinations(poses, 2):
        pa, pb = _mode_one_heavy_atoms(a), _mode_one_heavy_atoms(b)
        assert len(pa) == len(pb), "pose files disagree on heavy-atom count"
        out.append(math.sqrt(sum((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2
                                 for (x1, y1, z1), (x2, y2, z2) in zip(pa, pb))
                             / len(pa)))
    return out


def _centroid_sd(poses):
    """The sweep's own definition: the mean of the per-axis population SDs."""
    centroids = []
    for pose in poses:
        atoms = _mode_one_heavy_atoms(pose)
        centroids.append([sum(a[i] for a in atoms) / len(atoms) for i in range(3)])
    return statistics.fmean(
        statistics.pstdev([c[i] for c in centroids]) for i in range(3))


# --------------------------------------------------------------- the exclusion

def test_exactly_one_nominal_seed_run_carries_the_campaign_geometry():
    """The premise of the correction, checked against the files themselves."""
    runs = _seed_runs()
    campaign = _coordinate_digest(CAMPAIGN_POSE)
    duplicated = [r for r in runs if r["digest"] == campaign]
    assert len(duplicated) == 1, (
        f"expected one nominal seed run to carry the campaign geometry, found "
        f"{len(duplicated)}: {[r['seed'] for r in duplicated]}. If this is now "
        f"zero the seed was re-docked and the published text should say five "
        f"again; if it is more than one, the correction understates the problem")
    assert duplicated[0]["seed"] == 42, (
        f"the duplicated run is seed {duplicated[0]['seed']}, not 42; the "
        f"published text names seed 42 and now names the wrong one")


def test_the_excluded_runs_input_conformer_is_the_campaigns():
    """Seed 42 rebuilds the campaign's conformer, and the other four do not.

    This pins the fact in both directions, because the one-directional version
    misleads. The five seeds' conformers are pairwise distinct, which is what an
    earlier version of this test checked, and that reads as though the seed 42
    arm had a build of its own. It did not. Three of this repository's four
    embedding call sites use ETKDGv3 with the seed fixed at 42, so asking for
    seed 42 reproduces that input exactly. The fourth, in
    run_falsification_campaign.py, uses the default embedding method and a lower
    MMFF cap, so it is not the same recipe and the match also says the file did
    not come from there. Either way the input says nothing about whether the
    docking was independent, and the docked pose is what rules the arm out.
    """
    seeds = {seed: _coordinate_digest(POSES / f"disproof_seed_{seed}_CAHWRPWHC.pdbqt")
             for seed in (42, 101, 202, 303, 404)}
    campaign = _coordinate_digest(POSES / "contender_mutant_P2A_CAHWRPWHC.pdbqt")

    assert seeds[42] == campaign, (
        f"seed 42's input conformer {seeds[42]} no longer matches the campaign's "
        f"{campaign}; the page explains the match by the fixed embedding seed")
    for seed in (101, 202, 303, 404):
        assert seeds[seed] != campaign, (
            f"seed {seed} rebuilt the campaign's conformer, so it is not an "
            f"independent build either")
    assert len(set(seeds.values())) == 5, (
        f"the five input conformers are no longer pairwise distinct: {seeds}")


def test_the_embedding_seed_is_fixed_at_42_in_the_builder():
    """The reason the page gives for the match has to be in the code."""
    source = (REPO / "scripts" / "run_exhaustive_disproof.py").read_text()
    assert "randomSeed = 42" in source or "randomSeed=42" in source, (
        "the builder no longer fixes RDKit's embedding seed at 42, so the page's "
        "explanation for seed 42 reproducing the campaign conformer is stale")


def test_the_campaign_scripts_build_recipe_is_not_the_seeded_one():
    """The page says the conformer match is expected for the ETKDGv3 paths and
    not for the campaign script's own. Both halves have to hold in the source.

    Three call sites build with ETKDGv3 and a fixed seed of 42. The fourth, in
    run_falsification_campaign.py, calls EmbedMolecule with a bare randomSeed
    and no parameter object, then optimises with a lower iteration cap. If that
    ever becomes the same recipe, the page's second clause is wrong.
    """
    etkdg = REPO / "scripts" / "run_exhaustive_disproof.py"
    engine = REPO / "biophysics_engine.py"
    campaign = REPO / "scripts" / "run_falsification_campaign.py"

    for path in (etkdg, engine):
        text = path.read_text()
        assert "AllChem.ETKDGv3()" in text, f"{path.name} no longer builds with ETKDGv3"

    text = campaign.read_text()
    assert "AllChem.EmbedMolecule(mol, randomSeed=42)" in text, (
        "run_falsification_campaign no longer embeds with a bare randomSeed, so the "
        "page's claim that its recipe differs needs checking again")
    assert "ETKDGv3" not in text, (
        "run_falsification_campaign now uses ETKDGv3, so it is the same recipe as the "
        "seeded paths and the page's second clause is stale")
    assert "maxIters=300" in text and "maxIters=500" in engine.read_text(), (
        "the two paths no longer differ in their optimiser cap")


def test_four_runs_remain_independent():
    runs = _independent_runs()
    assert len(runs) == 4, f"expected four independent runs, found {len(runs)}"
    assert len({r["digest"] for r in runs}) == 4, (
        "two of the remaining runs share a docked pose, so fewer than four are "
        "independent and the published count is still wrong")


# ------------------------------------------------------- the claim cannot return

@pytest.mark.parametrize("relative", PUBLISHED)
@pytest.mark.parametrize("banned", BANNED)
def test_no_published_file_asserts_five_independent_seeds(relative, banned):
    text = (REPO / relative).read_text()
    assert banned not in text, (
        f"{relative} asserts five independent seeds again: {banned!r}. One of "
        f"the five nominal runs reports the campaign pose, so the independent "
        f"count is four")


@pytest.mark.parametrize("relative", PUBLISHED)
def test_any_surviving_mention_of_the_old_figure_keeps_its_reason(relative):
    """History is allowed. A bare reappearance is not."""
    text = (REPO / relative).read_text()
    # "Across 5 seeds" slipped past an earlier version of this list, which looked
    # for the N=5 form only, and left the withdrawn five-run figures asserted in
    # the project page's verdict paragraph.
    for marker in ("N=5", "N = 5", "5 seeds", "five seeds", "five runs"):
        start = 0
        while True:
            at = text.find(marker, start)
            if at < 0:
                break
            window = text[max(0, at - 400):at + 400]
            assert "exclud" in window, (
                f"{relative} mentions {marker!r} at offset {at} with no nearby "
                f"note that a run is excluded, so a reader cannot tell the "
                f"figure is historical")
            start = at + len(marker)


# -------------------------------------------------- the replacement figures

def test_the_published_four_run_score_matches_the_artefacts():
    scores = [r["score"] for r in _independent_runs()]
    mean = statistics.fmean(scores)
    population_sd = statistics.pstdev(scores)
    sample_sd = statistics.stdev(scores)

    readme = (REPO / "README.md").read_text()
    page = (REPO / "docs" / "index.html").read_text()
    for text, name in ((readme, "README.md"), (page, "docs/index.html")):
        assert f"{mean:.2f}" in text, (
            f"{name} does not carry the recomputed four-run mean {mean:.2f}")
        assert f"{population_sd:.2f}" in text, (
            f"{name} does not carry the recomputed four-run population SD "
            f"{population_sd:.2f}")
    assert f"{sample_sd:.2f}" in readme, (
        f"README.md does not carry the four-run sample SD {sample_sd:.2f}, so a "
        f"reader cannot tell which convention the headline figure uses")


def test_the_published_four_run_pose_spread_matches_the_artefacts():
    poses = [r["pose"] for r in _independent_runs()]
    rmsds = _pairwise_rmsd(poses)
    mean_rmsd = statistics.fmean(rmsds)
    worst = max(rmsds)
    centroid_sd = _centroid_sd(poses)

    readme = (REPO / "README.md").read_text()
    page = (REPO / "docs" / "index.html").read_text()
    for text, name in ((readme, "README.md"), (page, "docs/index.html")):
        assert f"{mean_rmsd:.2f}" in text, (
            f"{name} does not carry the recomputed four-run mean pairwise RMSD "
            f"{mean_rmsd:.2f} A")
    assert f"{worst:.2f}" in readme, (
        f"README.md does not carry the recomputed worst pair {worst:.2f} A")
    assert f"{centroid_sd:.2f}" in readme, (
        f"README.md does not carry the recomputed four-run centroid SD "
        f"{centroid_sd:.2f} A")


def test_the_published_pose_figures_are_over_real_heavy_atoms_only():
    """The selection has to be stated, because two selections are in play.

    Both appear in the README on purpose: the figure this project reports and
    the figure the earlier runs reported. A reader can only tell them apart if
    the page says which atoms each counts.
    """
    poses = [r["pose"] for r in _independent_runs()]
    heavy = len(_mode_one_heavy_atoms(poses[0]))
    readme = (REPO / "README.md").read_text()
    page = (REPO / "docs" / "index.html").read_text()
    assert f"{heavy} real heavy atoms" in readme, (
        f"README.md does not say the pose distances are over the {heavy} real "
        f"heavy atoms")
    assert "glue pseudo-atoms" in readme and "glue pseudo-atoms" in page, (
        "neither page explains that the earlier distances counted Meeko's glue "
        "pseudo-atoms, so the change in the figures looks unexplained")


def test_the_correction_did_not_silently_improve_the_result():
    """The point of the fix is the claim, not the number.

    If the four-run figures ever drift far from the five-run ones, the honest
    reading changes from 'one bad member' to 'the sweep says something else',
    and that needs prose, not a swapped number.
    """
    angles = json.loads(DISPROOF.read_text())["angles"]
    all_scores = next(a for a in angles if a.get("vina_scores"))["vina_scores"]
    four = [r["score"] for r in _independent_runs()]
    shift = abs(statistics.fmean(all_scores) - statistics.fmean(four))
    assert shift < 0.1, (
        f"dropping the duplicated run moves the mean by {shift:.3f} kcal/mol, "
        f"which is no longer a rounding matter: the published prose says the "
        f"result barely moves and would need rewriting")


# ---------------------------------------------------------------------------
# The permutation p value
#
# The project published p = 0.5897 for the scrambled null. That is 23/39, the
# bare proportion of controls scoring at least as well as the lead, and it was
# labelled a one-sided empirical rank p in the code, the README, the page, the
# dashboard, the benchmark file, four notebooks and two test assertions. The
# estimator for a permutation test counts the observation among the null,
# (k+1)/(n+1) = 24/40 = 0.600, because under the null the lead is one draw from
# the same distribution as its scrambles. The conclusion does not change: both
# figures are far from significance. The label did change, and a reader who
# recomputes the number should land on the same one the text prints.
# ---------------------------------------------------------------------------

# A file may state the old figure where it says what that figure was, so the old
# number is not banned outright. What is banned is writing it in the form of a p
# value: "p = 0.5897" and its variants. A historical note can say "0.5897, the raw
# proportion 23/39" and carry the same information without asserting it as p. An
# earlier version of this test allowed any line holding an explaining word, which
# let a line both assert the figure and explain it away.
_P_AS_A_P_VALUE = re.compile(r"\$?p\s*(?:=|:)\s*\**\$?\s*0\.5897", re.IGNORECASE)

_P_SURFACES = (
    "README.md",
    "FALSIFICATION_BENCHMARK.md",
    "dashboard.py",
    "docs/index.html",
    "falsification_engine.py",
    "scripts/run_exhaustive_disproof.py",
    "scripts/build_educational_sops.py",
    "scripts/build_sop_notebook.py",
    "scripts/build_readme_visuals.py",
    "paper/negative_result_GTF-SI.md",
)


def test_the_published_p_value_is_the_permutation_estimator():
    """No surface states the raw proportion in the form of a p value."""
    offenders = []
    for rel in _P_SURFACES:
        path = REPO / rel
        if not path.exists():
            continue
        for number, line in enumerate(path.read_text().splitlines(), start=1):
            if "0.5897" in line and _P_AS_A_P_VALUE.search(line):
                offenders.append(f"{rel}:{number}: {line.strip()[:120]}")
    assert not offenders, (
        "0.5897 is written as a p value; it is the raw proportion 23/39, and the "
        "permutation p is (k+1)/(n+1) = 0.600:\n" + "\n".join(offenders))


def test_no_notebook_states_the_old_p_value():
    """The notebooks are generated, and the generators were corrected, so a hit
    here means a notebook was committed without being regenerated. That already
    happened once: the Angle 9 label was corrected in the generator and the
    committed notebook kept the old one."""
    stale = [str(p.relative_to(REPO)) for p in (REPO / "notebooks").glob("*.ipynb")
             if "0.5897" in p.read_text()]
    assert not stale, f"regenerate these notebooks from their builders: {stale}"


def test_the_battery_recomputes_p_from_its_own_counts():
    """The recorded report carries the old figure, and it is left alone as the
    record of that run. The code must not agree with it."""
    recorded = json.loads(DISPROOF.read_text())
    angles = recorded.get("angles") or recorded.get("disproof_angles") or []
    angle5 = next(a for a in angles if a.get("angle") == 5)
    k = angle5["scrambles_at_least_as_good_as_lead"]
    n = angle5["n_isomers"]
    assert angle5["p_value"] == round(k / n, 4), (
        "the recorded report no longer carries the raw proportion, so this test has "
        "lost its subject; if the report was regenerated, that is the thing to check")
    source = (REPO / "scripts" / "run_exhaustive_disproof.py").read_text()
    assert "(better_or_equal + 1) / (n_isomers + 1)" in source, (
        "angle 5 no longer computes (k+1)/(n+1)")
    assert round((k + 1) / (n + 1), 4) == 0.600
