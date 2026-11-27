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
        "ligand_diameter_A": 21.57, "box_min_edge_A": 28.0,
        "box_clearance_per_side_A": 3.215,
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
    # A prepared ligand has coordinate columns, and the stub needs them: dock_one now
    # measures the ligand before docking, and a record with no coordinates is refused
    # as unmeasurable. Without the columns this test went green without ever reaching
    # Vina, which is not what it claims to prove.
    (tmp_path / "lig.pdbqt").write_text(
        "ATOM      1  C   UNL     1       0.000   0.000   0.000  0.00  0.00    "
        "+0.000 C \n")
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
    # Every check that sits between the overwrite guard and Vina must also have been
    # passed, or this test proves nothing about the overwrite guard.
    assert out.get("ligand_unmeasurable") is not True, (
        "the run stopped at the ligand measurement, so it never reached the overwrite "
        f"guard's downstream path: {out.get('error')}")
    assert out.get("box_too_small") is not True, (
        f"the run stopped at the box-fit check: {out.get('error')}")


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


# Third-party modules whose absence is an environment fact, not a defect. Any
# other import error, and every non-import error, is a broken analyser and must
# fail rather than skip: a skip here would report the atom-selection guard as
# satisfied while nothing had been executed.
_ABSENT_DEPENDENCY_IS_ALLOWED = {"numpy", "scipy", "pandas", "rdkit", "meeko"}


def _load_analyze():
    spec = importlib.util.spec_from_file_location("analyze_mod", HERE / "analyze.py")
    mod = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(mod)
    except ImportError as exc:
        missing = (getattr(exc, "name", "") or "").split(".")[0]
        if missing in _ABSENT_DEPENDENCY_IS_ALLOWED:
            pytest.skip(f"analyze.py needs {missing}, absent here")
        raise
    return mod


def test_the_evidence_fixtures_are_present_so_a_skip_is_loud():
    '''A skip is invisible in a green run, so absence is asserted, not tolerated.

    Eighteen tests in this file skip when a repository fixture is missing. In a
    complete checkout none of them should, but a renamed or moved directory would
    turn a block of evidence tests into skips and leave the run green. This test
    names the fixtures the evidence rests on and fails if one is gone, so the
    cause is reported once and loudly instead of silently subtracting checks.
    '''
    required = [
        REPO_ROOT / "data" / "poses",
        REPO_ROOT / "data" / "contenders_ledger.csv",
        REPO_ROOT / "data" / "exhaustive_disproof_report.json",
        REPO_ROOT / "biophysics_engine.py",
        REPO_ROOT / "scripts" / "run_exhaustive_disproof.py",
        HERE / "PROTOCOL.md",
        HERE / "analyze.py",
        HERE / "run_docking.py",
    ]
    missing = [str(path.relative_to(REPO_ROOT)) for path in required
               if not path.exists()]
    assert missing == [], (
        "evidence fixtures are missing, so tests in this file are skipping rather "
        f"than checking: {missing}")


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


def _build_local_ast():
    """The AST of run_docking.build_local, for property assertions."""
    import ast
    tree = ast.parse((HERE / "run_docking.py").read_text())
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name == "build_local":
            return node
    raise AssertionError("run_docking.build_local no longer exists")


def test_a_failed_embedding_is_raised_and_not_discarded():
    """Upstream assigns embed_res twice and never reads it, so this reads it.

    A molecule that failed to embed has no conformer. Upstream it goes on to
    Meeko regardless. The property asserted is not a phrase in the source but
    that the value is loaded somewhere in the function after assignment, and
    that the function can raise.
    """
    import ast
    fn = _build_local_ast()
    # A bare "is embed_res ever read" assertion would pass on the upstream code,
    # because its retry is itself guarded by `if embed_res != 0`. The property is
    # narrower: some `if` testing embed_res must raise, so a failed embedding
    # cannot reach Meeko.
    guarded_raise = False
    for node in ast.walk(fn):
        if not isinstance(node, ast.If):
            continue
        names = {n.id for n in ast.walk(node.test)
                 if isinstance(n, ast.Name)}
        if "embed_res" not in names:
            continue
        if any(isinstance(st, ast.Raise) for st in ast.walk(node)):
            guarded_raise = True
    assert guarded_raise, (
        "no branch on embed_res raises, so a molecule that failed to embed is "
        "still handed to Meeko, which is the upstream behaviour")


def test_the_minimisation_outcome_is_recorded_rather_than_assumed():
    """Both force fields may fail silently upstream; the outcome is now recorded.

    The fallback itself is kept, because a UFF conformer is still usable. What
    must not survive is a silent `pass`: a sequence that minimised with neither
    field is docked from raw ETKDG output while the README claims MMFF94, and
    no run in the repository records which happened.
    """
    import ast
    fn = _build_local_ast()
    src = ast.get_source_segment((HERE / "run_docking.py").read_text(), fn) or ""
    assert "MMFFGetMoleculeProperties" in src, (
        "whether MMFF could type the molecule at all is not recorded")
    for handler in [n for n in ast.walk(fn) if isinstance(n, ast.ExceptHandler)]:
        bodies = [type(st).__name__ for st in handler.body]
        assert bodies != ["Pass"], (
            "an exception handler in build_local still does nothing at all")
    assert '"minimisation"' in src and '"mmff_typed"' in src, (
        "the build report does not carry the minimisation outcome")


def test_the_seed_42_path_records_its_minimisation_as_unknown():
    """The imported builder reports nothing, so seed 42 must not claim success."""
    src = (HERE / "run_docking.py").read_text()
    i = src.index("def ligand_for(")
    body = src[i:i + 3000]
    assert '"minimisation": "unrecorded"' in body, (
        "the seed 42 path must record its minimisation as unknown, because the "
        "imported builder swallows both force-field failures")
    assert "minimisation" in body and "build_report.get" in body, (
        "the ligand metadata does not carry the minimisation outcome")


def test_no_skip_in_this_file_can_hide_a_broken_analyser():
    """A dependency may be absent; a broken analyze.py may not be tolerated."""
    import ast
    src = (HERE / "test_provenance.py").read_text()
    tree = ast.parse(src)
    for node in ast.walk(tree):
        if not (isinstance(node, ast.FunctionDef) and node.name == "_load_analyze"):
            continue
        handlers = [n for n in ast.walk(node) if isinstance(n, ast.ExceptHandler)]
        assert handlers, "_load_analyze no longer handles an import failure"
        for h in handlers:
            caught = ast.unparse(h.type) if h.type else "bare"
            assert caught == "ImportError", (
                f"_load_analyze catches {caught}; anything broader turns a broken "
                "analyser into a skip and reports its guards as satisfied")
        assert any(isinstance(n, ast.Raise) for n in ast.walk(node)), (
            "_load_analyze never re-raises, so an unexpected failure is swallowed")
        return
    raise AssertionError("_load_analyze no longer exists")


# Measurement-like identifiers. A skip hides a measurement; a sentinel invents
# one, and the invented one is worse because it survives into a figure. This
# repository already held exactly that: scripts/run_campaign.py assigned
# openmm_kj = -950.0 in an except clause, and -950.0 sits inside the ledger's
# real range of -5985.78 to -173.64, so a failed minimisation would have entered
# the data as an ordinary measurement. Verified it never fired: 260 ledger rows,
# 260 distinct values, no blanks, no -950.0, and "-950" appears in no artefact.
# Latent, not a corrupted result. This guard is so that no such line can appear
# in this package, whose runs have not happened yet.
_MEASUREMENT_WORDS = ("openmm", "vina", "vinardo", "kcal", "_kj", "rmsd", "drift",
                      "energy", "affinity", "score", "angstrom", "margin")


def _numeric_sentinels_in_except_handlers():
    import ast
    found = []
    for path in sorted(HERE.glob("*.py")):
        tree = ast.parse(path.read_text(), filename=str(path))
        for handler in [n for n in ast.walk(tree) if isinstance(n, ast.ExceptHandler)]:
            for st in ast.walk(handler):
                if not isinstance(st, (ast.Assign, ast.AnnAssign)):
                    continue
                targets = st.targets if isinstance(st, ast.Assign) else [st.target]
                names = []
                for t in targets:
                    names += [n.id for n in ast.walk(t) if isinstance(n, ast.Name)]
                    if isinstance(t, ast.Subscript):
                        try:
                            names.append(ast.unparse(t.slice).strip("'\""))
                        except Exception:
                            pass
                value = st.value
                bare = (isinstance(value, ast.Constant)
                        and isinstance(value.value, (int, float))
                        and not isinstance(value.value, bool))
                negated = (isinstance(value, ast.UnaryOp)
                           and isinstance(value.op, ast.USub)
                           and isinstance(value.operand, ast.Constant))
                if not (bare or negated):
                    continue
                if any(w in n.lower() for n in names for w in _MEASUREMENT_WORDS):
                    found.append(f"{path.name}:{st.lineno} {names} = "
                                 f"{ast.unparse(value)}")
    return found


def test_no_except_handler_invents_a_measurement():
    """A failed measurement must propagate or be recorded absent, never guessed."""
    sentinels = _numeric_sentinels_in_except_handlers()
    assert sentinels == [], (
        "an except handler assigns a number to a measurement-like name, which is "
        f"how a failure becomes a datum: {sentinels}")


def test_the_force_field_return_code_is_read_not_discarded():
    """maxIters is a cap, so running MMFF is not the same as converging.

    RDKit returns 0 on convergence and 1 when the cap was reached first. A
    conformer that stopped at the cap is partially minimised, and recording it
    as minimised would be a claim one step stronger than the evidence.
    """
    import ast
    fn = _build_local_ast()
    src = ast.get_source_segment((HERE / "run_docking.py").read_text(), fn) or ""
    assert "hit_iteration_cap" in src, (
        "a conformer that stopped at the iteration cap is recorded as if converged")
    calls_whose_value_is_dropped = []
    for node in ast.walk(fn):
        if not isinstance(node, ast.Expr) or not isinstance(node.value, ast.Call):
            continue
        name = getattr(node.value.func, "attr", None) or ""
        if name in ("MMFFOptimizeMolecule", "UFFOptimizeMolecule"):
            calls_whose_value_is_dropped.append(f"line {node.lineno}: {name}")
    assert calls_whose_value_is_dropped == [], (
        "a force-field call's return value is discarded, so convergence is "
        f"unknown and would be reported as success: {calls_whose_value_is_dropped}")
    assert "minimisation_iteration_cap" in (HERE / "run_docking.py").read_text(), (
        "the cap itself is not recorded beside the outcome")


def test_the_analysis_runs_end_to_end_on_a_synthetic_result_set(tmp_path):
    """Run analyze.py as a whole on rows built by the real job planner.

    Every function in analyze.py had a unit test and main() had none, so the glue
    that will read the user's out/results.jsonl after hours of docking had never
    executed once. This builds a result set the way run_docking.py builds one, by
    calling its own build_jobs, path_stem, identity and job_key rather than
    inventing the schema, attaches real pose geometry from data/poses, and runs
    the analyser as a subprocess exactly as RUN_ME.sh does.

    It asserts only that the pipeline completes and reports its sections. The
    verdicts on synthetic input mean nothing, and the redocking gate is expected
    to fail here because the stand-in native ligand is a peptide pose rather than
    acarbose. A gate that failed on nonsense is the right direction.

    Building this found nothing wrong with the package and one thing wrong with
    the harness: the first synthetic rows omitted finished_utc, because it is set
    by mutating prov after the literal rather than inside it, and the analyser
    correctly reported every row as lacking a provenance field. That was a defect
    in the test, not in the code, which is why the field is read from
    run_docking.py here instead of being listed by hand.
    """
    import hashlib, json, shutil, subprocess, sys, time, importlib.util

    poses = REPO_ROOT / "data" / "poses"
    root = tmp_path / "repo"
    (root / "data").mkdir(parents=True)
    shutil.copytree(HERE, root / "retarget_gtfb")
    shutil.copy(REPO_ROOT / "biophysics_engine.py", root / "biophysics_engine.py")
    shutil.copy(REPO_ROOT / "data" / "contenders_ledger.csv",
                root / "data" / "contenders_ledger.csv")
    pkg = root / "retarget_gtfb"
    shutil.rmtree(pkg / "out", ignore_errors=True)
    out = pkg / "out"
    out.mkdir()

    pool = sorted(poses.glob("*.pdbqt"))[:40]
    assert pool, "data/poses holds no pose files to build a synthetic run from"
    stand_in = out / "receptor.pdbqt"
    stand_in.write_bytes(pool[0].read_bytes())
    native = out / "native_ligand.pdbqt"
    native.write_bytes(pool[0].read_bytes())
    sha = hashlib.sha256(stand_in.read_bytes()).hexdigest()
    (out / "receptor_validation.json").write_text(json.dumps(
        {"sha256": sha, "title": "SYNTHETIC, NOT A RECEPTOR"}))
    (out / "boxes.json").write_text(json.dumps({
        "source_sha256": sha, "centre": [1.0, 2.0, 3.0],
        "receptor_pdbqt": str(stand_in), "native_ligand_pdbqt": str(native),
        "boxes": {"wide_28": {"size": [28, 28, 28]},
                  "tight_auto": {"size": [18, 18, 18]}}}))

    spec = importlib.util.spec_from_file_location("rd_dry", pkg / "run_docking.py")
    rd = importlib.util.module_from_spec(spec)
    sys.modules["rd_dry"] = rd
    try:
        spec.loader.exec_module(rd)
    except ImportError as exc:
        missing = (getattr(exc, "name", "") or "").split(".")[0]
        if missing in _ABSENT_DEPENDENCY_IS_ALLOWED:
            pytest.skip(f"run_docking needs {missing}, absent here")
        raise

    # The provenance fields the analyser requires, read from the analyser, and the
    # time fields read from the writer, so neither list is maintained by hand here.
    required = set(_load_analyze().REQUIRED_PROVENANCE)

    jobs = rd.dedupe_jobs(rd.build_jobs(json.loads((out / "boxes.json").read_text()), [1]))
    assert jobs, "the real job planner produced no tier 1 jobs"
    rd.POSE_DIR.mkdir(parents=True, exist_ok=True)
    rows = []
    for i, j in enumerate(jobs):
        tag = f"{i:08x}"
        stem = rd.path_stem(j, tag)
        pose = rd.POSE_DIR / f"{stem}.pdbqt"
        pose.write_bytes(pool[i % len(pool)].read_bytes())
        base = -8.0 - 0.05 * i
        modes = [{"mode": m + 1, "affinity_kcal": round(base + 0.37 * m, 3)}
                 for m in range(9)]
        stamp = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        digest = hashlib.sha256(pose.read_bytes()).hexdigest()
        rows.append({
            "key": rd.job_key(j), "identity": rd.identity(j), "path_stem": stem,
            "ligand_tag": tag, **j, "ok": True, "run_id": f"dry-{i}",
            "ligand_meta": {"heavy_atoms": 84, "minimisation": "mmff94_converged"},
            "recorded_utc": stamp, "started_utc": stamp, "finished_utc": stamp,
            "vina_command_string": f"vina --seed {j['dock_seed']}",
            "vina_version_text": "AutoDock Vina 1.2.7 (synthetic)",
            "exhaustiveness": j["exhaustiveness"], "dock_seed": j["dock_seed"],
            "box_centre": list(j["centre"]), "box_size": list(j["size"]),
            "receptor_md5": hashlib.md5(stand_in.read_bytes()).hexdigest(),
            "ligand_md5": hashlib.md5(pose.read_bytes()).hexdigest(),
            "ligand_diameter_A": 21.57, "box_min_edge_A": float(min(j["size"])),
            "box_clearance_per_side_A": round((float(min(j["size"])) - 21.57) / 2.0, 3),
            "seconds": 1.0, "vina_kcal": modes[0]["affinity_kcal"],
            "vina_modes": modes, "n_modes": len(modes),
            "mode_gap_kcal": round(modes[1]["affinity_kcal"] - modes[0]["affinity_kcal"], 3),
            "vinardo_kcal": round(modes[0]["affinity_kcal"] - 0.6, 3),
            "vinardo_scoring_context": "fixed_box", "ss_bond_angstrom": 2.05,
            "pose_pdbqt": str(pose), "pose_sha256": digest,
        })

    for field in sorted(required):
        assert field in rows[0], (
            f"the synthetic row omits {field}, which the analyser requires; fix the "
            "row rather than the requirement")

    (out / "results.jsonl").write_text(
        "".join(json.dumps(r) + "\n" for r in rows))

    proc = subprocess.run([sys.executable, "analyze.py"], cwd=pkg,
                          capture_output=True, text=True, timeout=900)
    assert proc.returncode == 0, (
        f"analyze.py exited {proc.returncode} on a well-formed result set.\n"
        f"stdout tail:\n{proc.stdout[-2000:]}\nstderr tail:\n{proc.stderr[-2000:]}")
    report_path = out / "analysis.json"
    assert report_path.exists(), "analyze.py exited 0 and wrote no analysis.json"
    report = json.loads(report_path.read_text())
    for section in ("provenance_audit", "redocking_calibration", "null_tests",
                    "seed_sensitivity", "noise_floor", "rows_usable", "caveats"):
        assert section in report, f"the report has no {section} section"
    assert report["rows_usable"] == len(rows)
    audit = report["provenance_audit"]
    assert audit["records_missing_provenance_fields"] == 0, (
        "well-formed rows were reported as lacking provenance: "
        f"{audit.get('records_missing_provenance_detail')}")


def test_only_the_dry_run_may_relax_the_tool_checks():
    """require_tools defaults to True and is relaxed at exactly one call site.

    Relaxing the Vina, version and package checks is safe for --dry-run and
    nowhere else. A later edit that passed require_tools=False unconditionally,
    or flipped the default, would let a real campaign start with no Vina version
    check, and every score it produced would be incomparable with the published
    ones for a reason no row would record. So this reads the default off the
    FunctionDef and the argument off the call, rather than trusting either.
    """
    import ast

    tree = ast.parse((HERE / "run_docking.py").read_text())

    default = None
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name == "preflight":
            names = [a.arg for a in node.args.args]
            assert "require_tools" in names, (
                "preflight no longer takes require_tools; this guard is stale and the "
                "dry run's promise is unpinned")
            offset = len(names) - len(node.args.defaults)
            default = node.args.defaults[names.index("require_tools") - offset]
    assert default is not None, "preflight was not found in run_docking.py"
    assert isinstance(default, ast.Constant) and default.value is True, (
        "require_tools must default to True so that anything which docks checks "
        f"its tools; the default is now {ast.dump(default)}")

    calls = [n for n in ast.walk(tree)
             if isinstance(n, ast.Call) and isinstance(n.func, ast.Name)
             and n.func.id == "preflight"]
    assert len(calls) == 1, f"expected one preflight call site, found {len(calls)}"
    kw = {k.arg: k.value for k in calls[0].keywords}
    assert "require_tools" in kw, (
        "the preflight call no longer says which checks it wants; if it relies on "
        "the default it cannot dry-run, and if it omits the keyword the reader "
        "cannot see which mode a real run uses")
    v = kw["require_tools"]
    assert isinstance(v, ast.UnaryOp) and isinstance(v.op, ast.Not), (
        "require_tools must be passed as the negation of the dry-run flag, not as a "
        f"literal; it is {ast.dump(v)}")
    inner = ast.dump(v.operand)
    assert "dry_run" in inner, (
        f"require_tools is negating something other than the dry-run flag: {inner}")


def test_the_dry_run_prints_the_plan_and_names_what_it_skipped(tmp_path):
    """--dry-run is advertised as listing the work, so run it and read the output.

    It was advertised in README_RUN_THIS.md as a way to list the jobs, and it
    could not run at all without Vina on PATH, because preflight died before the
    branch was reached. That is the one machine where a person would want to read
    the plan first. Fixing it created a second risk, that a clean dry run reads as
    a green light, so the output must also name every check it skipped.
    """
    import json, shutil, subprocess, sys, hashlib

    pkg = tmp_path / "retarget_gtfb"
    shutil.copytree(HERE, pkg)
    shutil.rmtree(pkg / "out", ignore_errors=True)
    out = pkg / "out"
    out.mkdir()
    stand_in = out / "receptor.pdbqt"
    stand_in.write_text("REMARK synthetic stand-in, not a receptor\n")
    native = out / "native_ligand.pdbqt"
    native.write_text("REMARK synthetic stand-in, not a ligand\n")
    sha = hashlib.sha256(stand_in.read_bytes()).hexdigest()
    (out / "receptor_validation.json").write_text(json.dumps(
        {"sha256": sha, "title": "SYNTHETIC, NOT A RECEPTOR"}))
    (out / "boxes.json").write_text(json.dumps({
        "source_sha256": sha, "centre": [1.0, 2.0, 3.0],
        "receptor_pdbqt": str(stand_in), "native_ligand_pdbqt": str(native),
        "boxes": {"wide_28": {"size": [28, 28, 28]},
                  "tight_auto": {"size": [18, 18, 18]}}}))

    proc = subprocess.run([sys.executable, "run_docking.py", "--dry-run", "--tiers", "1"],
                          cwd=pkg, capture_output=True, text=True, timeout=300)
    assert proc.returncode == 0, (
        f"--dry-run exited {proc.returncode}; it is advertised as listing the work "
        f"and stopping.\nstdout:\n{proc.stdout[-3000:]}\nstderr:\n{proc.stderr[-3000:]}")
    assert "docking runs in tiers" in proc.stdout, (
        f"--dry-run printed no plan:\n{proc.stdout[-3000:]}")

    if shutil.which("vina") is None:
        assert "NOT be" in proc.stdout or "NOT verified" in proc.stdout, (
            "Vina is absent here, so the dry run skipped the tool and version checks "
            "and must say so; a silent clean dry run reads as a green light to dock."
            f"\n{proc.stdout[-3000:]}")
        assert "vina" in proc.stdout.lower()
        assert "not evidence that the run can start" in proc.stdout.lower(), (
            f"the dry run did not say what it is not evidence of:\n{proc.stdout[-3000:]}")
    else:
        assert "NOT CHECKED" not in proc.stdout, (
            "Vina is on PATH here, so nothing should have been skipped")


# --------------------------------------------------------------------------- #
# Whether the ligand fits the box it was docked in
# --------------------------------------------------------------------------- #
def _single_model_pose(pattern: str) -> Path:
    """A committed single-conformer pose file, for measuring real geometry."""
    hits = sorted((REPO_ROOT / "data" / "poses").glob(pattern))
    assert hits, f"no committed pose matches {pattern}"
    return hits[0]


def test_the_heavy_atom_selection_reproduces_the_counts_rdkit_reports():
    """84 for the lead, 86 for the parent, 44 for acarbose, on the committed files.

    These are the counts RDKit gives before Meeko writes the pdbqt, so reproducing
    them from the written file is what shows the selection is right: hydrogens and
    the G pseudo-atoms out, the CG closure carbons in. A two-character column slice
    truncates CG0 to CG and gets this wrong in both directions, which is the defect
    that produced two withdrawn r.m.s.d. figures.
    """
    for pattern, expected in (("*CAHWRPWHC*split_*.pdbqt", 84),
                              ("*CPHWRPWHC*split_*.pdbqt", 86),
                              ("*acarbose*split_*.pdbqt", 44)):
        hits = sorted((REPO_ROOT / "data" / "poses").glob(pattern))
        assert hits, f"no committed single-model pose matches {pattern}"
        counts = {len(rd.heavy_atom_coords(h)) for h in hits}
        assert counts == {expected}, (
            f"{pattern}: heavy-atom counts {sorted(counts)}, expected {expected}")


def test_a_diameter_is_never_measured_across_superimposed_modes():
    """A multi-model pose file is refused rather than measured.

    A nine-mode acarbose file reads as 396 atoms spanning 21.22 A; its own single
    models are 44 atoms spanning 14.76 A. Measuring across modes gives the spread of
    the ensemble, not the size of a molecule, and it would overstate the ligand and
    refuse boxes that do in fact fit. I made this mistake by hand on these very
    files, which is why it is refused in code rather than remembered.
    """
    multi = sorted((REPO_ROOT / "data" / "poses").glob("*acarbose*docked.pdbqt"))
    multi = [m for m in multi if "split" not in m.name]
    assert multi, "no multi-model acarbose pose is committed"
    with pytest.raises(ValueError, match="models"):
        rd.ligand_diameter_A(multi[0])


def test_a_box_shorter_than_the_ligand_is_recorded_as_impossible_not_scored(tmp_path,
                                                                            monkeypatch):
    """The run is refused, and the record says by how much it fell short.

    Vina keeps the whole ligand inside the box, so when the ligand's diameter exceeds
    the shortest edge, the orientations along that edge cannot be sampled. The score
    that comes back is then not comparable with a box where they could be, and
    nothing else in the output would say so. The protocol tests box size as a factor,
    so this is the difference between measuring a box and measuring a truncation.
    """
    monkeypatch.setattr(rd, "POSE_DIR", tmp_path)
    (tmp_path / "r.pdbqt").write_text("ATOM      1  N   ALA A   1\n")
    pose = _single_model_pose("*CAHWRPWHC*split_*.pdbqt")
    diameter = rd.ligand_diameter_A(pose)
    assert diameter > 14.0, f"this test needs a ligand larger than its box; got {diameter}"

    out = rd.dock_one(ligand_pdbqt=pose, receptor=str(tmp_path / "r.pdbqt"),
                      centre=(1, 2, 3), size=(14.0, 28.0, 28.0), dock_seed=1,
                      tag="too_small", cpu=1, exhaustiveness=16)
    assert out["ok"] is False
    assert out["box_too_small"] is True
    assert out["box_min_edge_A"] == 14.0
    assert abs(out["ligand_diameter_A"] - diameter) < 0.01
    assert out["box_clearance_per_side_A"] < 0
    assert not (tmp_path / "too_small.pdbqt").exists(), (
        "a pose file was written for a run that was never docked")

    # And the same ligand in the campaign's own box is not refused. Reaching Vina,
    # which is absent here, is the proof that the fit check passed it.
    try:
        ok_box = rd.dock_one(ligand_pdbqt=pose, receptor=str(tmp_path / "r.pdbqt"),
                             centre=(1, 2, 3), size=(28.0, 28.0, 28.0), dock_seed=1,
                             tag="fits", cpu=1, exhaustiveness=16)
    except FileNotFoundError:
        return  # vina is not installed here; the fit check was passed
    assert ok_box.get("box_too_small") is not True, (
        "the check refuses a ligand that does fit the 28 A box the campaign used")


def test_an_unmeasurable_ligand_fails_its_own_run_and_not_the_campaign(tmp_path,
                                                                       monkeypatch):
    """One corrupt ligand file must not abort a run of hundreds.

    It must also not be measured from the atoms that did parse: a diameter over part
    of a molecule understates it, and an understated diameter is exactly what would
    pass a box that does not fit.
    """
    monkeypatch.setattr(rd, "POSE_DIR", tmp_path)
    truncated = tmp_path / "truncated_ligand.pdbqt"
    truncated.write_text("ATOM      1  C   UNL     1\nATOM      2  C   UNL     1\n")
    out = rd.dock_one(ligand_pdbqt=truncated, receptor=str(tmp_path / "r.pdbqt"),
                      centre=(1, 2, 3), size=(28.0, 28.0, 28.0), dock_seed=1,
                      tag="bad", cpu=1, exhaustiveness=16)
    assert out["ok"] is False
    assert out["ligand_unmeasurable"] is True
    assert "could not be" in out["error"]


def test_the_analyser_requires_the_fit_each_run_was_made_at():
    """ligand_diameter_A and box_min_edge_A are required provenance.

    The protocol compares two box sizes. A row that cannot say whether its ligand fit
    the box it was docked in cannot take part in that comparison, so the requirement
    belongs with the other fields that make a result attributable.
    """
    required = _load_analyze().REQUIRED_PROVENANCE
    for field in ("ligand_diameter_A", "box_min_edge_A"):
        assert field in required, f"{field} is not required provenance"


def test_a_box_arm_with_no_usable_runs_is_named_rather_than_counted():
    """An empty experimental arm is the finding, not a number in rows_failed.

    If every tight-box run were refused, a bare count beside a full report would not
    say that one side of the box comparison has no data.
    """
    an_mod = _load_analyze()
    rows = [
        {"ok": True, "box_id": "wide_28", "identity": "a"},
        {"ok": False, "box_id": "tight_auto", "identity": "b", "box_too_small": True,
         "ligand_diameter_A": 21.57, "box_min_edge_A": 14.0, "error": "too small"},
        {"ok": False, "box_id": "tight_auto", "identity": "c", "box_too_small": True,
         "ligand_diameter_A": 22.58, "box_min_edge_A": 14.0, "error": "too small"},
    ]
    res = an_mod.failure_breakdown(rows)
    assert res["rows_failed"] == 2
    assert res["per_box_arm"]["tight_auto"] == {"planned": 2, "usable": 0}
    assert res["per_box_arm"]["wide_28"] == {"planned": 1, "usable": 1}
    assert any("tight_auto" in n and "none usable" in n for n in res["notes"]), res["notes"]
    assert len(res["not_docked_because_the_box_was_shorter_than_the_ligand"]) == 2

    # And an arm that has data produces no such note.
    clean = an_mod.failure_breakdown([{"ok": True, "box_id": "wide_28", "identity": "a"}])
    assert clean["notes"] == []


def test_the_reference_diameter_in_config_is_what_the_committed_poses_measure():
    """config.json's observed_ligand_diameter_max_A is recomputed, not trusted.

    Defect 4 rests on this number, and a figure published in a config file is exactly
    the kind that goes stale silently: the 0.92 A drift figure in the manuscript's
    README had no source in any file, which is why every number here has to be
    recomputable from one. Tolerance is 0.01 A, the rounding in the config.
    """
    import json

    cfg = json.loads((HERE / "config.json").read_text())
    stated = cfg["docking"]["observed_ligand_diameter_max_A"]

    worst, worst_file = 0.0, None
    for f in sorted((REPO_ROOT / "data" / "poses").glob("*split_*.pdbqt")):
        try:
            d = rd.ligand_diameter_A(f)
        except ValueError:
            continue
        if d > worst:
            worst, worst_file = d, f.name
    assert worst > 0, "no committed single-model pose could be measured"
    assert abs(worst - stated) < 0.01, (
        f"config says the largest observed ligand diameter is {stated} A; the committed "
        f"poses measure {worst:.3f} A, in {worst_file}. One of them is stale, and "
        f"defect 4 and the plan-time warning both read this figure.")


def test_a_plan_that_uses_a_box_shorter_than_the_ligands_says_so_before_docking(tmp_path):
    """Requesting tier 3 warns at plan time, naming the arm and the tier.

    The per-run check refuses the runs, but it does so one at a time inside a run that
    may have been started for nine hours. The warning costs nothing and arrives before
    any docking.
    """
    import json, shutil, subprocess, sys, hashlib

    pkg = tmp_path / "retarget_gtfb"
    shutil.copytree(HERE, pkg)
    shutil.rmtree(pkg / "out", ignore_errors=True)
    out = pkg / "out"
    out.mkdir()
    stand_in = out / "receptor.pdbqt"
    stand_in.write_text("REMARK synthetic stand-in, not a receptor\n")
    native = out / "native_ligand.pdbqt"
    native.write_text("REMARK synthetic stand-in, not a ligand\n")
    sha = hashlib.sha256(stand_in.read_bytes()).hexdigest()
    (out / "receptor_validation.json").write_text(json.dumps(
        {"sha256": sha, "title": "SYNTHETIC, NOT A RECEPTOR"}))
    # A tight box as short as the one acarbose's bounding box would give here.
    (out / "boxes.json").write_text(json.dumps({
        "source_sha256": sha, "centre": [1.0, 2.0, 3.0],
        "receptor_pdbqt": str(stand_in), "native_ligand_pdbqt": str(native),
        "boxes": {"wide_28": {"size": [28, 28, 28]},
                  "tight_auto": {"size": [14.0, 20.5, 22.5]}}}))

    def plan(tiers):
        return subprocess.run([sys.executable, "run_docking.py", "--dry-run",
                               "--tiers", tiers], cwd=pkg, capture_output=True,
                              text=True, timeout=300)

    tight = plan("3")
    assert tight.returncode == 0, tight.stderr[-2000:]
    assert "WARNING" in tight.stdout, (
        f"tier 3 uses a 14.0 A edge and the plan said nothing:\n{tight.stdout[-2000:]}")
    assert "tight_auto" in tight.stdout
    assert "4.1c" in tight.stdout, "the warning does not say where the arithmetic is"

    wide = plan("1")
    assert wide.returncode == 0, wide.stderr[-2000:]
    assert "WARNING" not in wide.stdout, (
        f"tier 1 docks in the 28 A box and should draw no warning:\n{wide.stdout[-2000:]}")

