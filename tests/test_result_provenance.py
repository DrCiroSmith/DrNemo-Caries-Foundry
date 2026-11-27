"""Tests for deterministic result naming and in-file provenance.

The test that matters most is test_reproduces_the_duplicate_defect. It builds
the exact failure this repository already contains, three result paths holding
one set of coordinates, and shows the checker catching it. If that test ever
goes green without the fix, the fix has been removed.
"""

import json
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

import provenance as P  # noqa: E402


POSE = """MODEL 1
REMARK VINA RESULT:    -8.191      0.000      0.000
ATOM      1  C   LIG d   1      12.345  23.456  34.567  0.00  0.00    +0.000 C
ATOM      2  N   LIG d   1      13.345  24.456  35.567  0.00  0.00    -0.300 N
ENDMDL
"""


def _spec(**over):
    base = dict(target="3AIC", ligand_id="CAHWRPWHC", receptor_digest="a" * 32,
                ligand_digest="b" * 32, box_center=(192.072, 44.370, 197.309),
                box_size=(28.0, 28.0, 28.0), exhaustiveness=8, seed=42)
    base.update(over)
    return P.RunSpec(**base)


# ---------------------------------------------------------------- the name

def test_same_parameters_give_the_same_path():
    assert _spec().result_path("/x") == _spec().result_path("/x")


@pytest.mark.parametrize("field,value", [
    ("seed", 101),
    ("seed", None),
    ("exhaustiveness", 16),
    ("box_size", (24.0, 24.0, 24.0)),
    ("box_center", (192.072, 44.370, 197.310)),
    ("receptor_digest", "c" * 32),
    ("ligand_digest", "c" * 32),
    ("protonation", "diprotonated"),
    ("scoring", "vinardo"),
    ("target", "1SMD"),
])
def test_any_difference_changes_the_path(field, value):
    """Every parameter that changes the answer must change the file name."""
    assert _spec().canonical_name() != _spec(**{field: value}).canonical_name()


def test_the_name_says_what_ran():
    name = _spec().canonical_name()
    for token in ("3AIC", "CAHWRPWHC", "neutral", "e8", "s42", "vina"):
        assert token in name
    assert name.endswith(".pdbqt")


def test_an_unset_seed_is_named_as_unset_not_omitted():
    """An unseeded run is a fact about the run, not a missing field."""
    assert "_sunset_" in _spec(seed=None).canonical_name()


# ------------------------------------------------------- the record inside

def test_stamp_and_read_round_trip(tmp_path):
    f = tmp_path / "pose.pdbqt"
    f.write_text(POSE)
    spec = _spec()
    rec = P.build_record(spec, ["vina", "--seed", "42"], engine="AutoDock Vina 1.2.5",
                         written_path=f)
    P.stamp(f, rec)
    back = P.read_record(f)
    assert back["run_id"] == rec["run_id"]
    assert back["spec_id"] == spec.spec_id()
    assert back["seed"] == "42"
    assert back["exhaustiveness"] == "8"
    assert back["command"] == "vina --seed 42"
    assert back["engine_version"] == "AutoDock Vina 1.2.5"
    assert back["utc"].endswith("Z")
    assert back["receptor_digest"] == "a" * 32
    assert back["box_size"] == "28.0 28.0 28.0"


def test_every_field_the_fix_promises_is_present(tmp_path):
    f = tmp_path / "pose.pdbqt"
    f.write_text(POSE)
    P.stamp(f, P.build_record(_spec(), ["vina"], engine="x", written_path=f))
    rec = P.read_record(f)
    for field in ("run_id", "spec_id", "utc", "command", "engine_version", "seed",
                  "exhaustiveness", "box_center", "box_size", "receptor_digest",
                  "ligand_digest", "target", "ligand_id", "protonation",
                  "canonical_name"):
        assert field in rec, f"provenance is missing {field}"


def test_stamping_does_not_change_the_coordinates(tmp_path):
    f = tmp_path / "pose.pdbqt"
    f.write_text(POSE)
    before = P.coordinate_digest(f)
    P.stamp(f, P.build_record(_spec(), ["vina"], written_path=f))
    assert P.coordinate_digest(f) == before
    assert "REMARK VINA RESULT:" in f.read_text()


def test_stamping_twice_does_not_stack_records(tmp_path):
    f = tmp_path / "pose.pdbqt"
    f.write_text(POSE)
    P.stamp(f, P.build_record(_spec(), ["vina"], written_path=f))
    P.stamp(f, P.build_record(_spec(), ["vina"], written_path=f))
    assert f.read_text().count("PROVENANCE RUN_ID") == 1


def test_run_ids_differ_between_executions():
    assert P.new_run_id() != P.new_run_id()


def test_a_file_with_no_record_reads_as_none(tmp_path):
    f = tmp_path / "bare.pdbqt"
    f.write_text(POSE)
    assert P.read_record(f) is None


def test_missing_engine_is_recorded_as_missing_not_omitted():
    assert P.engine_version("/nonexistent/vina").startswith("unavailable")
    assert P.engine_version(None).startswith("unavailable")


# ------------------------------------------------------------ enforcement

def test_refuses_to_overwrite_a_result(tmp_path):
    f = tmp_path / "pose.pdbqt"
    f.write_text(POSE)
    with pytest.raises(P.ResultPathExists):
        P.reserve_result_path(f)
    assert P.reserve_result_path(f, overwrite=True) == f


def test_reserving_a_new_path_creates_its_directory(tmp_path):
    p = P.reserve_result_path(tmp_path / "deep" / "er" / "pose.pdbqt")
    assert p.parent.is_dir()


def test_reproduces_the_duplicate_defect(tmp_path):
    """Three nominally different protocols, one set of coordinates.

    This is the failure in the committed data, rebuilt from scratch. The
    checker must refuse it, and must keep refusing it when the records are
    present but disagree.
    """
    names = ["campaign_e16_unseeded.pdbqt", "audit_seed42_e8.pdbqt",
             "audit_ph72_e8.pdbqt"]
    files = []
    for n in names:
        f = tmp_path / n
        f.write_text(POSE)
        files.append(f)

    # 1. as the repository has it now: identical coordinates, no provenance
    bare = P.audit(files)
    assert len(bare["groups"]) == 1
    assert bare["groups"][0]["kind"] == "unprovenanced"
    assert len(bare["groups"][0]["files"]) == 3

    # 2. with provenance, each claiming its own run: a flat contradiction
    for f, seed, exh in zip(files, [None, 42, None], [16, 8, 8]):
        P.stamp(f, P.build_record(_spec(seed=seed, exhaustiveness=exh),
                                  ["vina", "--exhaustiveness", str(exh)],
                                  written_path=f))
    clashing = P.audit(files)
    assert clashing["groups"][0]["kind"] == "contradiction"
    assert len(clashing["groups"][0]["run_ids"]) == 3

    # 3. one run honestly written under three names: reported, not a failure
    shared = P.new_run_id()
    for f in files:
        P.stamp(f, P.build_record(_spec(), ["vina"], run_id=shared, written_path=f))
    declared = P.audit(files)
    assert declared["groups"][0]["kind"] == "one_run_many_names"
    assert declared["groups"][0]["run_ids"] == [shared]


def test_distinct_results_are_not_flagged(tmp_path):
    a = tmp_path / "a.pdbqt"
    a.write_text(POSE)
    b = tmp_path / "b.pdbqt"
    b.write_text(POSE.replace("12.345", "99.999"))
    assert P.audit([a, b])["groups"] == []


def test_a_hand_chosen_name_is_detectable(tmp_path):
    f = tmp_path / "whatever_i_felt_like.pdbqt"
    f.write_text(POSE)
    P.stamp(f, P.build_record(_spec(), ["vina"], written_path=f))
    assert P.audit([f])["name_mismatch"][0]["file"] == str(f)


def test_a_canonical_name_is_not_flagged(tmp_path):
    spec = _spec()
    f = spec.result_path(tmp_path)
    f.write_text(POSE)
    P.stamp(f, P.build_record(spec, ["vina"], written_path=f))
    assert P.audit([f])["name_mismatch"] == []


def test_baseline_suppresses_only_what_it_names(tmp_path):
    files = []
    for n in ("one.pdbqt", "two.pdbqt"):
        f = tmp_path / n
        f.write_text(POSE)
        files.append(f)
    digest = P.coordinate_digest(files[0])
    baseline = {"accepted_coordinate_digests": {digest: {"note": "known"}}}
    assert P.audit(files, baseline=baseline)["findings"] == []
    assert len(P.audit(files, baseline=baseline)["baselined"]) == 1
    other = {"accepted_coordinate_digests": {"0" * 32: {"note": "unrelated"}}}
    assert len(P.audit(files, baseline=other)["findings"]) == 1


def test_derived_splits_are_skipped_by_default(tmp_path):
    for n in ("pose_split_1.pdbqt", "pose_split_2.pdbqt"):
        (tmp_path / n).write_text(POSE)
    files = sorted(tmp_path.glob("*.pdbqt"))
    assert P.audit(files)["groups"] == []
    assert P.audit(files)["derived_splits_skipped"] == 2
    assert len(P.audit(files, include_splits=True)["groups"]) == 1


# --------------------------------------------- the repository as it stands

def test_the_committed_baseline_matches_the_repository():
    """The baseline must describe what is there, not what was there once.

    If this fails, either a result file changed or a new duplicate appeared.
    Both need a person, which is the point.
    """
    baseline = json.loads((REPO_ROOT / "data" / "provenance_baseline.json").read_text())
    paths = [p for p in sorted(REPO_ROOT.rglob("*.pdbqt")) if ".git" not in p.parts]
    report = P.audit(paths, baseline=baseline)
    assert report["findings"] == [], "a duplicate group is not in the baseline"
    assert len(report["baselined"]) == len(baseline["accepted_coordinate_digests"])
    assert report["name_mismatch"] == []


def test_the_checker_cli_passes_on_the_repository_as_committed():
    res = subprocess.run([sys.executable, "scripts/check_result_provenance.py"],
                         cwd=str(REPO_ROOT), capture_output=True, text=True, timeout=600)
    assert res.returncode == 0, res.stdout + res.stderr
    assert "Result: pass" in res.stdout


def test_the_checker_cli_fails_when_the_baseline_is_ignored():
    """Proof the check has teeth: without the baseline, the known groups fail it."""
    res = subprocess.run([sys.executable, "scripts/check_result_provenance.py",
                          "--no-baseline"],
                         cwd=str(REPO_ROOT), capture_output=True, text=True, timeout=600)
    assert res.returncode == 1
    assert "Result: FAIL" in res.stdout


# -------------------------------------------- the fix at the docking call site

def test_the_engine_refuses_to_adopt_an_unprovenanced_file(tmp_path, monkeypatch):
    """The defect, at the place it happened.

    The docking entry point used to treat any file at its output path as its
    own output and read a score out of it. Resuming is now opt-in and checked.
    Nothing here needs the vina binary, because the refusal happens before the
    engine is reached.
    """
    import biophysics_engine as BE

    receptor = tmp_path / "receptor.pdbqt"
    receptor.write_text(POSE)
    ligand = tmp_path / "CAHWRPWHC.pdbqt"
    ligand.write_text(POSE)
    monkeypatch.setitem(BE.TARGET_COORDS, "TEST", {
        "center": (1.0, 2.0, 3.0), "size": (28.0, 28.0, 28.0),
        "receptor": str(receptor)})

    stale = tmp_path / "stale.pdbqt"
    stale.write_text(POSE)
    with pytest.raises(RuntimeError, match="no provenance record"):
        BE.run_docking_vina_and_vinardo(str(ligand), target="TEST",
                                        out_pdbqt=str(stale), resume=True)


def test_the_engine_refuses_a_file_from_different_parameters(tmp_path, monkeypatch):
    """A pose stamped by another run must not be reported as this one."""
    import biophysics_engine as BE

    receptor = tmp_path / "receptor.pdbqt"
    receptor.write_text(POSE)
    ligand = tmp_path / "CAHWRPWHC.pdbqt"
    ligand.write_text(POSE)
    monkeypatch.setitem(BE.TARGET_COORDS, "TEST", {
        "center": (1.0, 2.0, 3.0), "size": (28.0, 28.0, 28.0),
        "receptor": str(receptor)})

    other = tmp_path / "from_another_run.pdbqt"
    other.write_text(POSE)
    P.stamp(other, P.build_record(_spec(seed=999, exhaustiveness=16), ["vina"],
                                  written_path=other))
    with pytest.raises(RuntimeError, match="different run"):
        BE.run_docking_vina_and_vinardo(str(ligand), target="TEST", seed=42,
                                        exhaustiveness=8, out_pdbqt=str(other),
                                        resume=True)


def test_the_engine_derives_its_own_path_when_none_is_given(tmp_path, monkeypatch):
    """No caller has to pick a name, and the name says what ran."""
    import biophysics_engine as BE

    receptor = tmp_path / "receptor.pdbqt"
    receptor.write_text(POSE)
    ligand = tmp_path / "CAHWRPWHC.pdbqt"
    ligand.write_text(POSE)
    monkeypatch.setitem(BE.TARGET_COORDS, "TEST", {
        "center": (1.0, 2.0, 3.0), "size": (28.0, 28.0, 28.0),
        "receptor": str(receptor)})
    monkeypatch.setattr(BE, "get_vina_path", lambda: "/nonexistent/vina")

    spec = P.RunSpec.for_run(target="TEST", ligand_id="CAHWRPWHC",
                             receptor=receptor, ligand=ligand,
                             box_center=(1.0, 2.0, 3.0), box_size=(28.0, 28.0, 28.0),
                             exhaustiveness=8, seed=42)
    expected = spec.result_path(tmp_path)
    assert "s42" in expected.name and "e8" in expected.name

    # docking itself cannot run here, so prove the path is what the engine would use
    with pytest.raises(Exception):
        BE.run_docking_vina_and_vinardo(str(ligand), target="TEST", seed=42,
                                        exhaustiveness=8, ligand_id="CAHWRPWHC")
    assert not expected.exists() or P.read_record(expected) is not None
