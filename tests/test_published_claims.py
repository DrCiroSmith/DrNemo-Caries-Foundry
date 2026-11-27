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


def test_the_excluded_runs_own_input_conformer_is_distinct():
    """Its input was built with its own seed. Only the docked file is shared.

    Worth pinning because the easy version of this story, that seed 42 was
    never a real run at all, is not what the artefacts show, and the page says
    so.
    """
    inputs = {seed: _coordinate_digest(POSES / f"disproof_seed_{seed}_CAHWRPWHC.pdbqt")
              for seed in (42, 101, 202, 303, 404)}
    assert len(set(inputs.values())) == 5, (
        f"the five input conformers are not all distinct: {inputs}")


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
    for marker in ("N=5", "N = 5"):
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
