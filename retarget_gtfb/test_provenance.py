#!/usr/bin/env python3
"""
Tests for result naming and provenance.

Why this file exists
--------------------
In the repository this package corrects, two docked results appear under five
filenames: the 3AIC result at -8.191 kcal/mol under three names, and the 1SMD result
at -7.343 under two. The -8.191 is the lead's published headline number. Those files
carry no seed, no command line, no run id and no timestamp, so nothing inside them
can say whether that was one run written several times or several runs that silently
collapsed onto one.

A third duplicated content in that directory is the input ligand conformer under
three names. That one is expected, because one builder at one seed gives one file,
and it is not a defect. The distinction between the benign case and the two real
ones cannot be made from the files, which is exactly the problem.

These tests are the guarantee that the same thing cannot happen in this package. The
central one, test_real_repository_duplicates_are_caught, takes the actual duplicated
files out of the repository and shows the check catching them.

Run with: python3 -m pytest retarget_gtfb/test_provenance.py -q
"""
from __future__ import annotations

import hashlib
import importlib.util
import re
import shutil
from pathlib import Path

import pytest

HERE = Path(__file__).parent.resolve()
REPO_ROOT = HERE.parent


def _load(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, HERE / filename)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


rd = _load("rd_under_test", "run_docking.py")
an = _load("an_under_test", "analyze.py")


# --------------------------------------------------------------------------- #
# The real failure
# --------------------------------------------------------------------------- #
def test_real_repository_duplicates_are_caught(tmp_path):
    """Copy the repository's own duplicated pose files out and show the audit catch them.

    Skips once data/poses has been cleaned up, which is the intended end state. While
    the duplicates are still there, this test is the proof that the check works on the
    real failure and not only on a fixture I wrote to be caught.
    """
    poses = REPO_ROOT / "data" / "poses"
    if not poses.is_dir():
        pytest.skip("data/poses is not present in this checkout")

    by: dict[str, list[Path]] = {}
    for f in sorted(poses.glob("*.pdbqt")):
        by.setdefault(hashlib.sha256(f.read_bytes()).hexdigest(), []).append(f)
    duplicated = {h: v for h, v in by.items() if len(v) > 1}
    if not duplicated:
        pytest.skip("data/poses no longer holds one content under several names")

    staged = tmp_path / "staged"
    staged.mkdir()
    for group in duplicated.values():
        for f in group:
            shutil.copy2(f, staged / f.name)

    res = an.audit_pose_directory(staged)
    assert res["passed"] is False
    assert res["contents_under_more_than_one_name"] == len(duplicated)
    # at least one of them is the three-names-one-content case
    assert max(len(v) for v in res["detail"].values()) >= 3
    # and every staged file is accounted for
    assert res["files_hashed"] == sum(len(v) for v in duplicated.values())


def test_clean_pose_directory_passes(tmp_path):
    for i in range(4):
        (tmp_path / f"pose_{i}.pdbqt").write_text(f"ATOM  {i}\nREMARK unique {i}\n")
    res = an.audit_pose_directory(tmp_path)
    assert res["passed"] is True
    assert res["contents_under_more_than_one_name"] == 0
    assert res["files_hashed"] == 4


# --------------------------------------------------------------------------- #
# The same failure expressed as result records
# --------------------------------------------------------------------------- #
def _record(identity: str, stem: str, run_id: str, pose_sha: str, **over) -> dict:
    rec = {
        "ok": True, "identity": identity, "path_stem": stem, "run_id": run_id,
        "pose_sha256": pose_sha, "vina_command_string": "vina --receptor r --ligand l",
        "dock_seed": 1, "exhaustiveness": 16, "box_centre": [1.0, 2.0, 3.0],
        "box_size": [28.0, 28.0, 28.0], "vina_version_text": "AutoDock Vina v1.2.7",
        "receptor_md5": "a" * 32, "ligand_md5": "b" * 32,
        "started_utc": "2026-11-27T00:00:00Z", "finished_utc": "2026-11-27T00:01:00Z",
        "vina_kcal": -8.1,
    }
    rec.update(over)
    return rec


def test_three_runs_one_content_is_a_hard_failure():
    """The exact shape of the existing repository's failure, as result records."""
    same = "c" * 64
    rows = [
        _record("contender|LEAD|8FJC|p|wide_28|d1|c42|ex16|rAAA", "stem_a", "run1", same),
        _record("contender|LEAD|8FJC|p|wide_28|d2|c42|ex16|rAAA", "stem_b", "run2", same),
        _record("contender|LEAD|8FJC|p|wide_28|d3|c42|ex16|rAAA", "stem_c", "run3", same),
    ]
    res = an.provenance_audit(rows)
    assert res["passed"] is False
    assert res["hard_failures"], "a shared content across three runs must be a hard failure"
    shared = res["identical_pose_content_under_distinct_runs"]
    assert isinstance(shared, list) and len(shared) == 1
    assert shared[0]["distinct_runs_sharing_this_content"] == 3
    assert sorted(shared[0]["run_ids"]) == ["run1", "run2", "run3"]
    assert sorted(shared[0]["path_stems"]) == ["stem_a", "stem_b", "stem_c"]


def test_distinct_runs_with_distinct_content_pass():
    rows = [_record(f"id{i}", f"stem{i}", f"run{i}", f"{i:064d}") for i in range(5)]
    res = an.provenance_audit(rows)
    assert res["passed"] is True
    assert res["hard_failures"] == []
    assert res["identical_pose_content_under_distinct_runs"] == "none"
    assert res["records_missing_provenance_fields"] == 0


def test_one_path_holding_two_results_is_a_hard_failure():
    rows = [_record("id1", "same_stem", "run1", "d" * 64),
            _record("id2", "same_stem", "run2", "e" * 64)]
    res = an.provenance_audit(rows)
    assert res["passed"] is False
    assert "same_stem" in res["one_path_holding_more_than_one_result"]


def test_reused_run_id_is_a_hard_failure():
    rows = [_record("id1", "stem1", "shared_run", "d" * 64),
            _record("id2", "stem2", "shared_run", "e" * 64)]
    res = an.provenance_audit(rows)
    assert res["passed"] is False
    assert "shared_run" in res["reused_run_ids"]


def test_identity_recorded_twice_is_noted_but_not_fatal():
    rows = [_record("id1", "stem1", "run1", "d" * 64),
            _record("id1", "stem1b", "run2", "e" * 64)]
    res = an.provenance_audit(rows)
    assert res["passed"] is True
    assert "id1" in res["identities_recorded_more_than_once"]
    assert any("more than one record" in n for n in res["notes"])


def test_missing_provenance_is_reported_and_not_fatal():
    rows = [_record("id1", "stem1", "run1", "d" * 64, vina_command_string=None,
                    receptor_md5=None)]
    res = an.provenance_audit(rows)
    assert res["records_missing_provenance_fields"] == 1
    assert "vina_command_string" in res["records_missing_provenance_detail"]["id1"]
    assert res["passed"] is True, "an unattributable row is a caveat, not a stop"


def test_failed_rows_are_not_audited():
    rows = [{"ok": False, "error": "vina exit 1"}]
    res = an.provenance_audit(rows)
    assert res["records_audited"] == 0
    assert res["passed"] is True


# --------------------------------------------------------------------------- #
# Deterministic naming
# --------------------------------------------------------------------------- #
def _job(**over) -> dict:
    j = {"tier": 2, "kind": "contender", "id": "mutant_P2A_CAHWRPWHC",
         "receptor_id": "8FJC", "receptor": "/r/8FJC.pdbqt", "receptor_tag": "1111aaaa",
         "centre": [1.0, 2.0, 3.0], "box_id": "wide_28", "size": [28.0, 28.0, 28.0],
         "dock_seed": 1, "conformer_seed": 42, "exhaustiveness": 16,
         "protonation": "as_built_no_pH_model_free_termini", "prebuilt_ligand": None}
    j.update(over)
    return j


def test_path_is_deterministic():
    assert rd.path_stem(_job(), "ligtag01") == rd.path_stem(_job(), "ligtag01")


@pytest.mark.parametrize("field,value", [
    ("tier", 3),
    ("id", "wild_type_CPHWRPWHC"),
    ("receptor_id", "8FJC_CA"),
    ("protonation", "explicit_pH_5_2"),
    ("box_id", "tight_auto"),
    ("dock_seed", 2),
    ("conformer_seed", 101),
    ("exhaustiveness", 32),
    ("receptor_tag", "2222bbbb"),
])
def test_every_run_parameter_changes_the_path(field, value):
    """Each parameter he named must be visible in the path, or two runs can collide."""
    base = rd.path_stem(_job(), "ligtag01")
    changed = rd.path_stem(_job(**{field: value}), "ligtag01")
    assert base != changed, f"changing {field} left the output path unchanged"


def test_ligand_content_changes_the_path():
    assert rd.path_stem(_job(), "ligtag01") != rd.path_stem(_job(), "ligtag02")


def test_path_has_no_hand_chosen_label():
    stem = rd.path_stem(_job(), "ligtag01")
    for part in ("t2", "contender", "mutant_P2A_CAHWRPWHC", "8FJC", "wide_28",
                 "d1", "c42", "ex16", "r1111aaaa", "Lligtag01"):
        assert part in stem, f"{part} is missing from the path"


def test_path_is_filesystem_safe_and_bounded():
    stem = rd.path_stem(_job(id="X" * 400), "ligtag01")
    assert len(stem) <= 150
    assert all(ch.isalnum() or ch in "_-." for ch in stem)


def test_long_names_still_give_distinct_paths():
    a = rd.path_stem(_job(id="Y" * 400), "ligtag01")
    b = rd.path_stem(_job(id="Y" * 399 + "Z"), "ligtag01")
    assert a != b, "the hash fallback must preserve distinctness past the length cap"


def test_identity_excludes_the_tier_so_resume_works_across_tiers():
    assert rd.identity(_job(tier=2)) == rd.identity(_job(tier=6))
    assert rd.job_key(_job(tier=2)) != rd.job_key(_job(tier=6))


def test_identical_runs_requested_by_two_tiers_are_docked_once():
    jobs = [_job(tier=2), _job(tier=6), _job(tier=3, box_id="tight_auto")]
    deduped = rd.dedupe_jobs(jobs)
    assert len(deduped) == 2
    shared = [j for j in deduped if j["also_requested_by_tiers"]][0]
    assert shared["also_requested_by_tiers"] == [6]


# --------------------------------------------------------------------------- #
# Regression: the macrocycle glue-atom read
#
# An earlier version of has_glue_atoms() read the AutoDock type with a two-column
# slice, which truncates CG0 to CG, and on the strength of it I told Dr Roi that
# Meeko had left the ring rigid. It had not. The type is the last token on the line.
# --------------------------------------------------------------------------- #
GLUE_LINE = ("ATOM      9  G   UNL     1       1.000   2.000   3.000  "
             "0.00  0.00    +0.000 G0\n")
CG_LINE = ("ATOM     10  CG  UNL     1       4.000   5.000   6.000  "
           "0.00  0.00    +0.000 CG0\n")
PLAIN_LINE = ("ATOM      1  N   CYS A   1       0.000   0.000   0.000  "
              "0.00  0.00    -0.300 N\n")


def test_glue_atoms_are_detected_whatever_the_column(tmp_path):
    f = tmp_path / "open_ring.pdbqt"
    f.write_text(PLAIN_LINE + CG_LINE + GLUE_LINE)
    assert rd.has_glue_atoms(f) is True


def test_glue_atoms_detected_when_the_type_column_is_shifted(tmp_path):
    f = tmp_path / "shifted.pdbqt"
    f.write_text(PLAIN_LINE + CG_LINE.replace("    +0.000 CG0", "  +0.000 CG0"))
    assert rd.has_glue_atoms(f) is True


def test_a_rigid_ring_reads_as_rigid(tmp_path):
    f = tmp_path / "rigid.pdbqt"
    f.write_text(PLAIN_LINE * 3 + "TORSDOF 14\n")
    assert rd.has_glue_atoms(f) is False


def test_a_cg_named_atom_is_not_mistaken_for_glue(tmp_path):
    """A side-chain CG carbon is named CG and typed C. It is not a glue pseudo-atom."""
    f = tmp_path / "sidechain.pdbqt"
    f.write_text("ATOM      5  CG  LEU A   2       1.000   1.000   1.000  "
                 "0.00  0.00    +0.020 C\n")
    assert rd.has_glue_atoms(f) is False


def test_the_committed_poses_do_have_an_opened_ring():
    """The claim I got wrong, now asserted against the files themselves."""
    poses = sorted((REPO_ROOT / "data" / "poses").glob("*.pdbqt"))
    if not poses:
        pytest.skip("data/poses is not present in this checkout")
    nine_mers = [p for p in poses if rd.has_glue_atoms(p)]
    assert nine_mers, ("Meeko did open the macrocycle in the committed poses. If this "
                       "assertion ever fails, check the detector before the conclusion.")


# --------------------------------------------------------------------------- #
# Refusing to overwrite
# --------------------------------------------------------------------------- #
def test_docking_refuses_to_overwrite_an_unclaimed_result(tmp_path, monkeypatch):
    """An existing pose file at a deterministic path is not quietly replaced.

    The check runs before Vina is invoked, so this test needs no Vina.
    """
    monkeypatch.setattr(rd, "POSE_DIR", tmp_path)
    (tmp_path / "already_here.pdbqt").write_text("ATOM      1  C   UNL     1\n")
    out = rd.dock_one(ligand_pdbqt=tmp_path / "lig.pdbqt", receptor="r.pdbqt",
                      centre=(1, 2, 3), size=(28, 28, 28), dock_seed=1,
                      tag="already_here", cpu=1, exhaustiveness=16)
    assert out["ok"] is False
    assert out.get("overwrite_refused") is True
    assert "--overwrite" in out["error"]


def test_overwrite_flag_lets_it_through(tmp_path, monkeypatch):
    """With the flag set the run proceeds past the check and fails on the absent Vina.

    Reaching a Vina failure is the proof that the overwrite guard was passed, which is
    all this test is for.
    """
    monkeypatch.setattr(rd, "POSE_DIR", tmp_path)
    (tmp_path / "already_here.pdbqt").write_text("ATOM      1  C   UNL     1\n")
    (tmp_path / "lig.pdbqt").write_text("ATOM      1  C   UNL     1\n")
    (tmp_path / "r.pdbqt").write_text("ATOM      1  N   ALA A   1\n")
    try:
        out = rd.dock_one(ligand_pdbqt=tmp_path / "lig.pdbqt",
                          receptor=str(tmp_path / "r.pdbqt"),
                          centre=(1, 2, 3), size=(28, 28, 28), dock_seed=1,
                          tag="already_here", cpu=1, exhaustiveness=16,
                          allow_overwrite=True)
    except FileNotFoundError:
        return  # vina is not installed here; the guard was still passed
    assert out.get("overwrite_refused") is not True


# --------------------------------------------------------------------------- #
# The measurements that cost no docking runs
# --------------------------------------------------------------------------- #
def test_minimum_detectable_effect_is_bounded_by_the_null_size():
    null = [-8.0 - 0.1 * i for i in range(39)]  # 39 members, best is -11.8
    res = an.minimum_detectable_effect(null, alpha=0.05)
    assert res["n_null"] == 39
    assert res["smallest_attainable_p"] == round(1 / 40, 4)
    # at n = 39 and alpha = 0.05 the observation must beat every null member
    assert res["nulls_the_observation_may_not_beat"] == 0
    assert res["score_threshold_kcal"] == min(null)
    assert res["minimum_detectable_effect_kcal"] > 0


def test_a_null_too_small_to_ever_reach_significance_says_so():
    res = an.minimum_detectable_effect([-8.0, -8.1, -8.2], alpha=0.05)
    assert res["smallest_attainable_p"] == 0.25
    assert "cannot return a significant result" in res["verdict"]
    assert "minimum_detectable_effect_kcal" not in res


def test_size_slope_is_withheld_when_the_range_is_too_narrow():
    entries = [{"heavy_atoms": 84 + (i % 3), "score": -8.0 - 0.01 * i} for i in range(40)]
    res = an.size_confound(entries)
    assert res.get("slope_withheld") is True
    assert "slope_kcal_per_heavy_atom" not in res


def test_size_slope_is_reported_over_a_real_range():
    # a deliberate 0.05 kcal/mol per heavy atom relationship over 60 atoms of span
    entries = [{"heavy_atoms": 60 + i, "score": -(5.0 + 0.05 * (60 + i))} for i in range(60)]
    res = an.size_confound(entries)
    assert res["slope_kcal_per_heavy_atom"] == pytest.approx(-0.05, abs=0.001)
    assert res["r_squared"] > 0.99


def test_spearman_matches_known_cases():
    assert an.spearman([1, 2, 3, 4], [1, 2, 3, 4]) == 1.0
    assert an.spearman([1, 2, 3, 4], [4, 3, 2, 1]) == -1.0
    assert an.spearman([1, 1, 1], [1, 2, 3]) is None  # no variance to correlate
    assert an.spearman([1, 2], [1, 2]) is None  # too few points


def test_the_real_ledger_size_slope_is_small_but_not_zero():
    """His own 260-contender ledger, as the reference the narrow-range guard points to."""
    import csv
    ledger = REPO_ROOT / "data" / "contenders_ledger.csv"
    if not ledger.exists():
        pytest.skip("contenders_ledger.csv is not present in this checkout")
    rows = list(csv.DictReader(ledger.open()))
    entries = [{"heavy_atoms": int(float(r["heavy_atoms"])),
                "score": float(r["vina_3aic_kcal"])} for r in rows]
    res = an.size_confound(entries)
    assert res["n"] == len(rows)
    assert res["heavy_atom_range"][1] - res["heavy_atom_range"][0] >= 5
    assert res["slope_kcal_per_heavy_atom"] < 0, "a bigger ligand scores better in Vina"
    # the 2-heavy-atom mismatch between the lead and null A is worth this much
    assert abs(2 * res["slope_kcal_per_heavy_atom"]) < 0.15, (
        "if this ever exceeds the claimed lead-null gap, size alone explains the gap")


def test_the_duplication_count_is_pinned():
    """Pin the reconciled figures, because I got this arithmetic wrong once.

    I first reported eight filenames over five distinct results. The eight filenames
    cover three duplicated contents, of which only two are docked results, and those
    two sit under five filenames. This test fixes the numbers so neither the protocol
    nor the manuscript can drift from them.
    """
    poses = REPO_ROOT / "data" / "poses"
    if not poses.is_dir():
        pytest.skip("data/poses is not present in this checkout")

    by: dict[str, list[Path]] = {}
    for f in sorted(poses.glob("*.pdbqt")):
        by.setdefault(hashlib.sha256(f.read_bytes()).hexdigest(), []).append(f)
    dups = {h: v for h, v in by.items() if len(v) > 1}
    if not dups:
        pytest.skip("data/poses no longer holds one content under several names")

    docked, inputs = [], []
    for h, group in dups.items():
        text = group[0].read_text(errors="replace")
        (docked if "REMARK VINA RESULT" in text else inputs).append(group)

    assert len(by) == 1239, f"distinct contents changed: {len(by)}"
    assert sum(len(v) for v in by.values()) == 1244
    assert len(dups) == 3, "three contents appear under more than one name"
    assert len(docked) == 2, "only two of the three duplicated contents are docked results"
    assert len(inputs) == 1, "the third is an input conformer, which is benign"
    assert sum(len(g) for g in docked) == 5, "the two docked results sit under five names"
    # and the headline number is one of them
    affinities = set()
    for group in docked:
        for line in group[0].read_text(errors="replace").splitlines():
            if line.startswith("REMARK VINA RESULT:"):
                affinities.add(float(line.split()[3]))
                break
    assert -8.191 in affinities, "the lead's published 3AIC number is one of the duplicates"


def test_an_input_conformer_duplicate_is_told_apart_from_a_docked_one(tmp_path):
    """The benign case and the real one must be distinguishable by content, not by name."""
    conformer = ("ATOM      1  N   CYS A   1       0.0  0.0  0.0  0.00  0.00  -0.3 N\n"
                 "TORSDOF 15\n")
    pose = ("MODEL 1\nREMARK VINA RESULT:    -8.191      0.000      0.000\n"
            "ATOM      1  N   CYS A   1       0.0  0.0  0.0  0.00  0.00  -0.3 N\n"
            "ENDMDL\n")
    for name in ("a.pdbqt", "b.pdbqt", "c.pdbqt"):
        (tmp_path / name).write_text(conformer)
    for name in ("d.pdbqt", "e.pdbqt"):
        (tmp_path / name).write_text(pose)

    res = an.audit_pose_directory(tmp_path)
    assert res["contents_under_more_than_one_name"] == 2
    docked = [names for h, names in res["detail"].items()
              if "REMARK VINA RESULT" in (tmp_path / names[0]).read_text()]
    assert len(docked) == 1 and len(docked[0]) == 2


# --------------------------------------------------------------------------
# The margin at issue, and what replication it would cost.
#
# These exist because an earlier revision of PROTOCOL.md pre-registered a
# 0.15 kcal/mol "claimed advantage" that has no source anywhere in the
# repository. It was my number, not his, and it was about to be gated
# against. The sourced quantity is the lead's margin over the parent it was
# derived from, and it is pinned here so it cannot drift back to a guess.
# --------------------------------------------------------------------------

def _null_from_ledger() -> list[float]:
    """The 39 scrambled nulls from his own ledger, as every power test uses them."""
    import csv
    ledger = REPO_ROOT / "data" / "contenders_ledger.csv"
    if not ledger.exists():
        pytest.skip("contenders_ledger.csv is not present in this checkout")
    rows = list(csv.DictReader(ledger.open()))
    return [float(r["vina_3aic_kcal"]) for r in rows
            if "scrambled" in r["contender_id"]]


def test_the_margin_at_issue_is_sourced_from_the_ledger():
    import csv
    ledger = REPO_ROOT / "data" / "contenders_ledger.csv"
    if not ledger.exists():
        pytest.skip("contenders_ledger.csv is not present in this checkout")
    rows = list(csv.DictReader(ledger.open()))
    by_id = {r["contender_id"]: r for r in rows}
    lead = float(by_id["contender_mutant_P2A_CAHWRPWHC"]["vina_3aic_kcal"])
    parent = float(by_id["contender_g01_CPHWRPWHC"]["vina_3aic_kcal"])
    assert lead == -8.191
    assert parent == -8.062
    assert round(parent - lead, 4) == 0.129, (
        "the margin the protocol gates against must be this subtraction and nothing "
        "else; 0.15 kcal/mol appears nowhere in the repository")


def test_no_unsourced_margin_is_pre_registered():
    """0.15 may appear in PROTOCOL.md only inside the note withdrawing it."""
    text = (HERE / "PROTOCOL.md").read_text()
    for line in text.splitlines():
        if "0.15 " in line or "0.15 kcal" in line:
            assert "no source" in line or "withdraw" in line.lower(), (
                f"PROTOCOL.md still asserts 0.15 kcal/mol outside the withdrawal "
                f"note: {line.strip()[:120]}")


def test_seeds_needed_scales_as_the_square_root_of_the_margin():
    null = _null_from_ledger()
    a = an.seeds_needed(0.129, null)
    b = an.seeds_needed(0.258, null)
    assert a["seeds_per_control_to_reach_the_boundary"] == 36
    assert a["seeds_per_control_for_the_stated_power"] == 82
    assert a["docking_runs_to_reach_the_boundary"] == 36 * 40
    # halving the required resolution quarters the seeds, to rounding
    assert b["seeds_per_control_to_reach_the_boundary"] == 9


def test_seeds_needed_refuses_a_non_positive_margin():
    null = _null_from_ledger()
    assert "reason" in an.seeds_needed(0.0, null)
    assert "reason" in an.seeds_needed(-0.1, null)


def test_the_boundary_claim_is_about_controls_not_all_sample_sizes():
    """The protocol must not claim the effect is unreachable in principle."""
    text = (HERE / "PROTOCOL.md").read_text()
    assert "however\n   many control sequences are added" in text or \
           "however many control sequences are added" in text
    # The phrase may survive only as quoted history: inside the amendment log,
    # where it records what the document used to say, or next to the word that
    # marks it as overreach. A bare live assertion of it is the failure.
    offenders = []
    for l in text.splitlines():
        if "at any sample size" not in l:
            continue
        quoted = '"outside this design\'s resolution at any sample size"' in l
        if not (quoted or "overreach" in l or "withdraw" in l.lower()):
            offenders.append(l)
    assert not offenders, f"overreaching wording survives: {offenders[:2]}"


def test_the_protocol_carries_a_dated_amendment_log():
    """A pre-registered plan that changed silently is worth nothing."""
    text = (HERE / "PROTOCOL.md").read_text()
    assert "## 0. Amendment log" in text
    # every amendment names a commit, so the diff can be read rather than trusted
    import re
    rows = [l for l in text.splitlines()
            if l.startswith("| ") and re.search(r"\| `[0-9a-f]{7}` \|", l)]
    assert len(rows) >= 7, f"only {len(rows)} amendments carry a commit hash"
    assert "No result from this package existed at any of these" in text


def test_both_contrasts_are_named_and_point_opposite_ways():
    """0.129 and 0.152 measure different comparisons and must never be conflated."""
    text = (HERE / "PROTOCOL.md").read_text()
    assert "Contrast A, lead against parent" in text
    assert "Contrast B, lead against the scrambled null" in text
    assert "+0.129 kcal/mol" in text and "0.152 kcal/mol worse" in text
    assert "23 of the 39 controls score at least as well" in text
    # the gate must name its contrasts, not just its numbers
    rule12 = [l for l in text.splitlines() if l.startswith("| 12 |")]
    assert rule12, "decision rule 12 is missing"
    low = rule12[0].lower()
    assert "contrast a" in low and "contrast b" in low, (
        "rule 12 gives numbers without naming which comparison each measures")


def test_the_price_of_the_replication_is_stated_with_the_gate():
    text = (HERE / "PROTOCOL.md").read_text()
    rule12 = [l for l in text.splitlines() if l.startswith("| 12 |")][0]
    for token in ("36 seeds per control", "1,440", "82 seeds", "3,280"):
        assert token in rule12, f"rule 12 states the gate without its price: {token}"


# --------------------------------------------------------------------------
# The repository's N=5 seed sweep contains at most four independent runs.
#
# Seed 42's docked pose and its input conformer are both byte-identical to
# the campaign's own artefacts, so any dispersion over all five is computed
# over a set with one non-independent member. Pinned because the choice of
# four versus five flips the sign of the between-sequence variance estimate.
# --------------------------------------------------------------------------

def _mode_coord_digest(path: Path) -> str:
    coords = [l[30:54] for l in path.read_text(errors="replace").splitlines()
              if l.startswith(("ATOM", "HETATM"))]
    return hashlib.sha256("".join(coords).encode()).hexdigest()


def test_seed_42_is_not_an_independent_replicate():
    poses = REPO_ROOT / "data" / "poses"
    if not poses.exists():
        pytest.skip("data/poses is not present in this checkout")
    campaign = poses / "contender_mutant_P2A_CAHWRPWHC_3aic.pdbqt"
    seed42 = poses / "disproof_seed_42_CAHWRPWHC_3AIC_docked.pdbqt"
    others = [poses / f"disproof_seed_{s}_CAHWRPWHC_3AIC_docked.pdbqt"
              for s in (101, 202, 303, 404)]
    for f in [campaign, seed42, *others]:
        if not f.exists():
            pytest.skip(f"{f.name} is not present in this checkout")
    assert _mode_coord_digest(seed42) == _mode_coord_digest(campaign), (
        "seed 42's pose must be identical to the campaign pose; if this ever fails, "
        "the four-versus-five replicate argument has to be re-derived")
    digests = {f.name: _mode_coord_digest(f) for f in others}
    assert len(set(digests.values())) == 4, "the other four seeds must be distinct"
    assert _mode_coord_digest(campaign) not in digests.values()
    # and the input conformer duplicates in the same pattern
    lig_campaign = poses / "contender_mutant_P2A_CAHWRPWHC.pdbqt"
    lig_seed42 = poses / "disproof_seed_42_CAHWRPWHC.pdbqt"
    if lig_campaign.exists() and lig_seed42.exists():
        assert lig_campaign.read_bytes() == lig_seed42.read_bytes(), (
            "seed 42 reused the campaign's input conformer as well as its pose")



def _production_call_sites_passing_a_seed():
    """Every non-test call of the docking entry point that could pass a seed.

    The seed's positional index is read from the function's own signature rather
    than hardcoded, so reordering its parameters cannot make this pass by
    accident. A call counts as possibly seeded if it names `seed`, reaches the
    seed's position positionally, or forwards `**kwargs`. Files under tests/ are
    excluded: a test fixture that seeds the engine deliberately is not evidence
    about what produced the deposited poses.
    """
    import ast
    target = "run_docking_vina_and_vinardo"
    engine = REPO_ROOT / "biophysics_engine.py"
    seed_index = None
    engine_tree = ast.parse(engine.read_text(), filename=str(engine))
    for node in ast.walk(engine_tree):
        if isinstance(node, ast.FunctionDef) and node.name == target:
            names = [a.arg for a in node.args.args]
            if "seed" in names:
                seed_index = names.index("seed") + 1
    found = []
    for path in sorted(REPO_ROOT.rglob("*.py")):
        rel = path.relative_to(REPO_ROOT)
        if rel.parts and rel.parts[0] in ("tests", ".git"):
            continue
        try:
            tree = ast.parse(path.read_text(), filename=str(path))
        except SyntaxError as exc:
            # Never swallow this. A file that will not parse is a file whose call
            # sites were not examined, and reporting it as clean is how a seeded
            # call site would hide. Found by a negative control whose injected
            # call was itself invalid syntax, which this clause silently passed.
            found.append(f"{rel}: WOULD NOT PARSE ({exc.msg} at line {exc.lineno})")
            continue
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            name = getattr(node.func, "attr", None) or getattr(node.func, "id", None)
            if name != target:
                continue
            kw = {k.arg for k in node.keywords if k.arg}
            forwards = any(k.arg is None for k in node.keywords)
            positional = seed_index is not None and len(node.args) >= seed_index
            if "seed" in kw or forwards or positional:
                found.append(f"{rel}:{node.lineno}")
    return found

def test_the_disproof_sweep_is_seeded_even_though_the_engine_is_not():
    """Defect 2 is scoped to the production pipeline, and must stay scoped."""
    sweep = REPO_ROOT / "scripts" / "run_exhaustive_disproof.py"
    engine = REPO_ROOT / "biophysics_engine.py"
    if not (sweep.exists() and engine.exists()):
        pytest.skip("the campaign scripts are not present in this checkout")
    sweep_text = sweep.read_text()
    assert "seeds = [42, 101, 202, 303, 404]" in sweep_text
    assert "--seed" in sweep_text, "the sweep passes an explicit seed to Vina"
    assert "randomSeed" in sweep_text, "and seeds the conformer builder too"
    # The claim is NOT that the string "--seed" is absent from the engine. The
    # engine may gain an optional seed parameter without touching defect 2, and
    # on the manuscript branch it already has. The falsifiable claim is that no
    # production call site passes one, so every docking run the repository has
    # actually produced was stochastic. Asserted over the call graph, because a
    # string test would fail on a change that does not affect the finding and
    # would pass on a call site that does.
    seeded = _production_call_sites_passing_a_seed()
    assert seeded == [], (
        "a production call site now seeds Vina, so defect 2 must be revisited: "
        f"{seeded}")
    protocol = (HERE / "PROTOCOL.md").read_text()
    assert "scoped to" in protocol and "run_docking_vina_and_vinardo" in protocol, (
        "the protocol must say which code path the unseeded finding applies to")


def test_what_a_larger_null_buys_is_quantified_not_asserted():
    protocol = (HERE / "PROTOCOL.md").read_text()
    assert "A larger null buys 25 percent and then nothing" in protocol
    for token in ("2.203 standard deviations", "j = 0 at", "0.0250, 0.0450, 0.0490"):
        assert token in protocol, f"the n-dependence is not quantified: {token}"


# Meeko's macrocycle glue markers. Not atoms, so never in an r.m.s.d. The closure
# carbons typed CG0 are real carbons and are kept; note that columns 78 and 79 alone
# truncate CG0 to CG, which is the truncation behind the withdrawn rigid-ring claim.
GLUE_PSEUDO_ATOM_TYPES = {"G0", "G1", "G2", "G3"}


def _mode1_coords(path: Path, heavy_only: bool) -> list[tuple[float, float, float]]:
    """Coordinates of the first MODEL record only, in file order, no fitting."""
    out: list[tuple[float, float, float]] = []
    started = False
    for line in path.read_text().splitlines():
        if line.startswith("MODEL"):
            if started:
                break
            started = True
            continue
        if line.startswith("ENDMDL") and started:
            break
        if line.startswith(("ATOM", "HETATM")):
            atom_type = line[77:79].strip()
            if atom_type in GLUE_PSEUDO_ATOM_TYPES:
                continue
            if heavy_only and atom_type.startswith("H"):
                continue
            out.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return out


def _pairwise_rmsd(poses: dict[int, list[tuple[float, float, float]]]) -> tuple[float, float, tuple[int, int]]:
    """Mean and worst in-place r.m.s.d. over every pair, and which pair is worst."""
    import itertools
    import math

    vals = []
    for a, b in itertools.combinations(sorted(poses), 2):
        pa, pb = poses[a], poses[b]
        assert len(pa) == len(pb), (a, b, len(pa), len(pb))
        total = sum(
            (x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2
            for (x1, y1, z1), (x2, y2, z2) in zip(pa, pb)
        )
        vals.append((math.sqrt(total / len(pa)), a, b))
    worst = max(vals)
    return sum(v for v, _, _ in vals) / len(vals), worst[0], (worst[1], worst[2])


def test_dropping_the_non_independent_seed_does_not_weaken_the_pose_finding():
    """Section 4.1b says the pose conclusion survives seed 42's removal. Recompute it.

    The point this pins is directional, not a figure: if the mean pairwise spread ever
    falls when the non-independent member is dropped, then the published
    score-reproduces-pose-does-not finding would partly rest on that member, and 4.1b
    would be wrong to say the correction is only to the label and the score summary.
    """
    poses_dir = REPO_ROOT / "data" / "poses"
    if not poses_dir.exists():
        pytest.skip("data/poses is not present in this checkout")
    files = {
        s: poses_dir / f"disproof_seed_{s}_CAHWRPWHC_3AIC_docked.pdbqt"
        for s in (42, 101, 202, 303, 404)
    }
    for f in files.values():
        if not f.exists():
            pytest.skip(f"{f.name} is not present in this checkout")

    # all five docked files are genuinely different poses, so the pairwise set is real
    assert len({_mode_coord_digest(f) for f in files.values()}) == 5, (
        "4.1b asserts the pairwise comparison runs over five distinct poses")

    all_atom = {s: _mode1_coords(f, heavy_only=False) for s, f in files.items()}
    assert len(all_atom[42]) == 102, (
        "real-atom count changed; the 102 in 4.1b is stale. 104 records minus two "
        "G0 glue markers")

    mean5, worst5, _ = _pairwise_rmsd(all_atom)
    four = {s: v for s, v in all_atom.items() if s != 42}
    mean4, worst4, worst_pair = _pairwise_rmsd(four)

    assert mean4 >= mean5, (
        f"dropping seed 42 lowered the mean pairwise r.m.s.d. ({mean5:.2f} to "
        f"{mean4:.2f} A); 4.1b's claim that the pose finding does not depend on it "
        f"would need rewriting")
    assert 42 not in worst_pair, (
        "the four-seed worst pair must be between two independent seeds")
    assert abs(worst4 - worst5) < 0.02, (
        "4.1b says the worst pair is effectively unaffected by the removal")

    text = (HERE / "PROTOCOL.md").read_text()
    for figure in (f"{mean5:.2f}", f"{mean4:.2f}", f"{worst4:.2f}"):
        assert figure in text, (
            f"4.1b quotes a pairwise r.m.s.d. figure that no longer recomputes: "
            f"{figure} is absent from PROTOCOL.md")

    assert "G0" in text and "CG0" in text, (
        "4.1b must state which records are excluded; a bare atom count is not a "
        "selection a reader can reproduce")
    assert "9.24" in text and "12.19" in text, (
        "4.1b must keep the glue-inclusive figures beside the correct ones, so a "
        "reader who recomputes either can see which selection produced which")

    heavy = {s: _mode1_coords(f, heavy_only=True) for s, f in files.items()}
    assert len(heavy[42]) == 84, "heavy-atom count changed; the 84 in 4.1b is stale"
    h_mean5, _, _ = _pairwise_rmsd(heavy)
    h_mean4, h_worst4, _ = _pairwise_rmsd({s: v for s, v in heavy.items() if s != 42})
    assert h_mean4 >= h_mean5
    for figure in (f"{h_mean5:.2f}", f"{h_mean4:.2f}", f"{h_worst4:.2f}"):
        assert figure in text, f"heavy-atom figure {figure} is absent from PROTOCOL.md"


def test_the_four_run_score_summary_is_recomputed_not_asserted():
    """4.1b gives a four-run mean and spread. Derive both from the report, in order."""
    import json

    report = REPO_ROOT / "data" / "exhaustive_disproof_report.json"
    if not report.exists():
        pytest.skip("exhaustive_disproof_report.json is not present in this checkout")
    angles = json.loads(report.read_text())["angles"]
    sweep = [a for a in angles if a.get("n_seeds") == 5]
    assert len(sweep) == 1, "expected exactly one five-seed angle in the report"
    scores = sweep[0]["vina_scores"]
    assert len(scores) == 5

    # the sweep declares its seeds in this order, so the first entry is seed 42's
    src = (REPO_ROOT / "scripts" / "run_exhaustive_disproof.py").read_text()
    assert "seeds = [42, 101, 202, 303, 404]" in src, (
        "the seed order is what makes scores[0] seed 42's; it has changed")
    four = scores[1:]

    mean = sum(four) / len(four)
    pop_sd = (sum((v - mean) ** 2 for v in four) / len(four)) ** 0.5
    sample_sd = (sum((v - mean) ** 2 for v in four) / (len(four) - 1)) ** 0.5

    text = (HERE / "PROTOCOL.md").read_text()
    assert f"{mean:.3f}" in text, f"4.1b's four-run mean {mean:.3f} no longer recomputes"
    assert f"{pop_sd:.4f}" in text, f"4.1b's four-run spread {pop_sd:.4f} no longer recomputes"
    assert f"{sample_sd:.4f}" in text, (
        f"the sample standard deviation {sample_sd:.4f} is the figure the manuscript's "
        f"variance work uses; 4.1b must still state it so the two stay reconciled")
    # and the published five-run figures must still be what we are contrasting against
    assert "-8.10" in text and "0.37" in text


def test_no_file_states_a_literal_count_of_the_provenance_suite():
    """A hardcoded test count goes stale every time a test is added.

    It had gone stale twice before this guard existed, so the phrasing that invites
    one is banned outright rather than corrected again. `--check-naming` derives its
    own count from the list it runs, which is why eighteen is still spelled there.
    """
    pattern = re.compile(r"\b(?:suite|set)\s+of\s+(?:\d+|"
                         r"(?:twenty|thirty|forty|fifty|sixty|seventy|eighty|ninety)"
                         r"(?:[- ](?:one|two|three|four|five|six|seven|eight|nine))?)",
                         re.IGNORECASE)
    offenders = []
    for path in sorted(HERE.glob("*")):
        if path.is_dir() or path.suffix not in (".md", ".sh", ".py"):
            continue
        if path.name == "test_provenance.py":
            continue
        for n, line in enumerate(path.read_text().splitlines(), 1):
            if pattern.search(line):
                offenders.append(f"{path.name}:{n}: {line.strip()}")
    assert not offenders, (
        "these state a count of the provenance suite, which goes stale as soon as a "
        "test is added; say 'the full provenance suite' instead:\n" + "\n".join(offenders))


def test_the_withdrawal_count_in_the_log_is_derived_from_the_log():
    """The log's prose states how many entries withdraw a claim. Count them.

    An earlier version of that sentence said "three of those seven entries" when the
    table held nine rows. A pre-registered document whose own summary of itself is
    stale is worse than one with no summary, so the number is now checked against the
    rows rather than maintained by hand.
    """
    text = (HERE / "PROTOCOL.md").read_text()
    rows = [l for l in text.splitlines()
            if re.match(r"^\|\s*\d\d:\d\d\s*\|\s*`[0-9a-f]{7}`\s*\|", l)]
    assert len(rows) >= 9, f"the amendment log has shrunk to {len(rows)} rows"

    withdrawals = [l for l in rows if "withdraw" in l.lower()]
    words = {1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six",
             7: "seven", 8: "eight", 9: "nine", 10: "ten"}
    expected = words.get(len(withdrawals), str(len(withdrawals)))
    sentence = re.search(
        r"(\w+) of these entries withdraw a claim of mine", text, re.IGNORECASE)
    assert sentence, (
        "the log must still state how many of its entries are withdrawals; that "
        "sentence is the honest part of keeping a log at all")
    assert sentence.group(1).lower() == expected, (
        f"the log says '{sentence.group(1)}' withdrawals but {len(withdrawals)} rows "
        f"contain one; the sentence has gone stale")

    # and the log must not claim a total row count, which is the way it went stale
    assert not re.search(r"of (?:those|these) (?:\w+) entries\b(?!\s+withdraw)", text), (
        "do not restate the log's total length in prose; it goes stale on every "
        "amendment")


def test_the_glue_inclusive_selection_reproduces_the_superseded_figures():
    """The glue-inclusive selection yields 9.24 and 12.19 from these pose files.

    This began as a check on someone else's published number. That number has since
    been withdrawn and replaced, so what is pinned now is the selection itself: which
    atom filter produces which figure, from files that do not change. Keeping it means
    a future session can tell at once whether a figure it meets was computed with the
    glue included, without rediscovering the mechanism a fourth time.
    """
    import itertools
    import math

    poses_dir = REPO_ROOT / "data" / "poses"
    if not poses_dir.exists():
        pytest.skip("data/poses is not present in this checkout")
    files = {
        s: poses_dir / f"disproof_seed_{s}_CAHWRPWHC_3AIC_docked.pdbqt"
        for s in (42, 101, 202, 303, 404)
    }
    for f in files.values():
        if not f.exists():
            pytest.skip(f"{f.name} is not present in this checkout")

    def mode1_drop_h_keep_glue(path: Path):
        out, started = [], False
        for line in path.read_text().splitlines():
            if line.startswith("MODEL"):
                if started:
                    break
                started = True
                continue
            if line.startswith("ENDMDL") and started:
                break
            if line.startswith(("ATOM", "HETATM")):
                if line[77:79].strip().startswith("H"):
                    continue
                out.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
        return out

    sel = {s: mode1_drop_h_keep_glue(f) for s, f in files.items()}
    assert len(sel[42]) == 86, (
        f"the glue-inclusive heavy selection should hold 86 records, got {len(sel[42])}")

    vals = []
    for a, b in itertools.combinations(sorted(sel), 2):
        pa, pb = sel[a], sel[b]
        total = sum((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2
                    for (x1, y1, z1), (x2, y2, z2) in zip(pa, pb))
        vals.append(math.sqrt(total / len(pa)))
    mean, worst = sum(vals) / len(vals), max(vals)

    assert f"{mean:.2f}" == "9.24", f"expected the manuscript's 9.24, got {mean:.2f}"
    assert f"{worst:.2f}" == "12.19", f"expected the manuscript's 12.19, got {worst:.2f}"

    text = (HERE / "PROTOCOL.md").read_text()
    assert "9.04" in text and "11.80" in text, (
        "4.1b must state the corrected heavy-atom values beside the reproduced ones")


def _load_analyze():
    spec = importlib.util.spec_from_file_location("analyze_mod", HERE / "analyze.py")
    mod = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(mod)
    except Exception as exc:  # numpy absent, for instance
        pytest.skip(f"analyze.py will not import here: {exc}")
    return mod


def test_the_shared_reader_excludes_glue_on_every_ligand_in_the_repository():
    """One reader, checked on a glued peptide, a differently glued peptide and acarbose.

    The hazard this guards is that a selection can look correct on one ligand and be
    wrong on another: acarbose carries no glue markers at all, so a reader that admits
    them passes every acarbose check and still corrupts every peptide figure. So the
    identity records minus glue minus hydrogens is asserted per ligand rather than a
    single expected count.
    """
    mod = _load_analyze()
    poses = REPO_ROOT / "data" / "poses"
    if not poses.exists():
        pytest.skip("data/poses is not present in this checkout")

    cases = {
        "contender_mutant_P2A_CAHWRPWHC_3aic.pdbqt": 84,   # the lead, glued
        "scrambled_01_CHRPWHPWC_1smd.pdbqt": 86,           # a null, glued, different size
        "3AIC_native_acarbose_docked.pdbqt": 44,           # acarbose, no glue at all
    }
    checked = 0
    for name, expected in cases.items():
        path = poses / name
        if not path.exists():
            continue
        records = glue = hydrogens = 0
        started = False
        for line in path.read_text().splitlines():
            if line.startswith("MODEL"):
                if started:
                    break
                started = True
                continue
            if line.startswith("ENDMDL") and started:
                break
            if line.startswith(("ATOM", "HETATM")):
                records += 1
                t = mod.atom_type_of(line)
                if t in mod.GLUE_PSEUDO_ATOM_TYPES:
                    glue += 1
                if t in ("H", "HD", "HS"):
                    hydrogens += 1

        coords, elements = mod.read_pdbqt_atoms(path)
        assert len(coords) == records - glue - hydrogens, (
            f"{name}: reader selected {len(coords)} but records minus glue minus "
            f"hydrogens is {records - glue - hydrogens}")
        assert len(coords) == expected, (
            f"{name}: expected {expected} heavy atoms, reader gave {len(coords)}")
        assert "G" not in elements, f"{name}: a glue marker reached the element list"
        checked += 1

    assert checked >= 2, (
        "at least two ligands must be checked, one glued and one not; only "
        f"{checked} were present")


def test_no_module_reads_the_atom_type_with_a_fixed_short_slice():
    """Fix the reader, not the call sites. This is what stops a fourth recurrence.

    A two-character slice of the type field truncates CG0 to CG, which makes a real
    closure carbon indistinguishable from a glue marker. Three figures were reported
    wrongly from that slice. The only place allowed to touch those columns is
    analyze.py's atom_type_of, as a documented fallback for plain PDB input whose
    element field genuinely is columns 77 and 78.
    """
    # a subscript of something, so that prose describing the old bug is not flagged
    subscript = re.compile(r"\w\[\s*7[67]\s*:\s*79\s*\]")
    narrow = re.compile(r"\w\[\s*76\s*:\s*78\s*\]")
    offenders = []
    for path in sorted(HERE.glob("*.py")):
        if path.name == "test_provenance.py":
            continue
        for n, line in enumerate(path.read_text().splitlines(), 1):
            if subscript.search(line):
                offenders.append(f"{path.name}:{n}: {line.strip()}")
            if narrow.search(line) and path.name != "analyze.py":
                offenders.append(f"{path.name}:{n}: {line.strip()}")
    assert not offenders, (
        "read the AutoDock type as the last whitespace token, not by column:\n"
        + "\n".join(offenders))

    mod = _load_analyze()
    assert mod.atom_type_of(
        "ATOM      1  CA  CYS A   1       1.000   0.000   0.000  1.00  0.00    +0.152 CG0"
    ) == "CG0"
    assert mod.atom_type_of(
        "ATOM      2  G   CYS A   1       2.000   0.000   0.000  1.00  0.00    +0.000 G0"
    ) == "G0"


def test_the_pytest_free_gate_also_checks_the_atom_selection():
    """He may run without pytest, so the selection check cannot live only here."""
    text = (HERE / "analyze.py").read_text()
    assert "Atom selection" in text, (
        "analyze.py --self-test must carry the atom-selection block; RUN_ME.sh runs "
        "the self-test unconditionally but pytest only if it is installed")
    for needle in ("the glue pseudo-atom G0 is excluded",
                   "the closure carbon CG0 is kept, as carbon"):
        assert needle in text, f"the self-test no longer asserts: {needle}"


def test_the_conformer_half_of_the_seed_42_finding_stays_narrowed():
    """The input-conformer match is expected by construction and is not evidence.

    Every ETKDGv3 path that builds this ligand embeds at seed 42, so asking the sweep
    for seed 42 rebuilds the campaign conformer. 4.1b said for a while that seed 42 was
    non-independent "on both sides of the calculation", which overstated it. The
    docked pose is the evidence; the conformer is not. This test exists because that
    overstatement is easy to reintroduce from a summary.
    """
    # the document is hard-wrapped, so compare on normalised whitespace
    flat = " ".join((HERE / "PROTOCOL.md").read_text().split())
    assert "The conformer match carries no weight" in flat, (
        "4.1b must say plainly that the conformer identity is not evidence")
    assert "The docked pose match is the evidence" in flat, (
        "4.1b must say which half of the finding carries the weight")
    # the phrase may appear only where it is quoted as withdrawn history, which is how
    # the amendment log records it; it may not appear as a live claim
    live = [l for l in (HERE / "PROTOCOL.md").read_text().splitlines()
            if "on both sides of the calculation" in l
            and "withdraw" not in l.lower()]
    assert not live, (
        "that phrasing treated the conformer match as evidence. It is withdrawn, and "
        "may appear only on a line that says so:\n" + "\n".join(live))

    # and the claim rests on the builder, so the builder must still hardcode seed 42
    engine = REPO_ROOT / "biophysics_engine.py"
    if engine.exists():
        src = engine.read_text()
        assert "params.randomSeed = 42" in src, (
            "build_macrocycle_3d no longer hardcodes seed 42, so the "
            "expected-by-construction argument in 4.1b needs re-deriving")

    sweep = REPO_ROOT / "scripts" / "run_exhaustive_disproof.py"
    if sweep.exists():
        assert "params.randomSeed = s" in sweep.read_text(), (
            "the sweep no longer parameterises the conformer seed the same way")
