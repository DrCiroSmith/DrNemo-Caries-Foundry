"""The twelve angle labels, and the counts every surface quotes for them.

A label here is a verdict a reader acts on, and three defects in this
repository were of one kind: a label that stated more than its evidence, and
stated it in several files at once. Angle 12 is arithmetic on the sequence with
no structure in it, and Angle 11 integrates no dynamics, yet both were written
into the code as the literal string PASSED_IN_SILICO while the gate booleans
they had just computed were discarded. Angle 9's margin is 0.53 of the scoring
error Angle 4 is judged against, and it was rendered as a pass on the page and
in the dashboard while the README called it inconclusive.

These tests do not check that any particular label is right. They check that
every surface says the same thing, that no surface claims a pass the table does
not carry, and that a classification in the code depends on something.
"""

import ast
import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
README = REPO / "README.md"
PAGE = REPO / "docs" / "index.html"
DASHBOARD = REPO / "dashboard.py"
BATTERY = REPO / "scripts" / "run_exhaustive_disproof.py"

# The prose each label is written as in a verdict sentence. The counts are not
# listed here on purpose: they are counted off the README table, so this file
# cannot drift from it.
PROSE = {
    "INCONCLUSIVE_WITHIN_NOISE": "inconclusive within noise",
    "AT_RISK": "at risk",
    "CALIBRATION_FAILED": "calibration failed",
    "NO_HOTSPOT_IDENTIFIED": "no hotspot identified",
    "SCORE_REPRODUCIBLE_POSE_NOT": "score-reproducible-but-pose-not",
    "COMPOSITION_ARITHMETIC_ONLY": "composition arithmetic only",
    "NOT_COMPUTED": "not computed",
}

ROW = re.compile(r"^\|\s*\*\*(\d+)\*\*\s*\|.*\|\s*`([A-Z_]+)`\s*\|?\s*(?:here)?\s*\|?\s*$")


def readme_labels():
    """The label in the last cell of each of the twelve angle rows."""
    found = {}
    for line in README.read_text().splitlines():
        m = ROW.match(line.strip())
        if m:
            found[int(m.group(1))] = m.group(2)
    return found


def test_the_readme_table_still_has_twelve_angles():
    labels = readme_labels()
    assert sorted(labels) == list(range(1, 13)), (
        f"expected angles 1 to 12, parsed {sorted(labels)}. If the table was "
        "reformatted, fix the ROW pattern rather than deleting this test.")


def test_no_angle_in_the_readme_table_is_labelled_a_pass():
    """Both figures published as passes are withdrawn, with the reason in the row.

    Angle 9's -1.52 kcal/mol margin is 0.53 of the 2.85 kcal/mol scoring error
    Angle 4 is judged against. Angle 12 is Henderson-Hasselbalch arithmetic on
    the sequence. If a later run earns a pass, this test should be changed by
    whoever can point at the measurement behind it.
    """
    passes = {n: l for n, l in readme_labels().items() if l == "PASSED_IN_SILICO"}
    assert not passes, f"PASSED_IN_SILICO is back in the README table at {sorted(passes)}"


def test_the_readme_verdict_counts_match_its_own_table():
    """The verdict sentence used to say two passes and five inconclusive while the
    table carried one pass and six. A reader takes the count, not the table."""
    labels = readme_labels()
    text = README.read_text()
    verdict = next(l for l in text.splitlines() if l.startswith("> **Verdict.**"))
    counted = {}
    for label in labels.values():
        counted[label] = counted.get(label, 0) + 1
    assert sum(counted.values()) == 12
    for label, n in sorted(counted.items()):
        phrase = f"{n} {PROSE[label]}"
        assert phrase in verdict, (
            f"the table carries {n} x {label} but the verdict sentence does not "
            f"say {phrase!r}")
    assert "no passes" in verdict, (
        "the table carries no pass, so the verdict has to say so outright")


def test_the_page_carries_no_pass_tag():
    """The page rendered Angles 9 and 12 as PASSED while the README did not."""
    tags = re.findall(r'class="status-tag (\w+)">([^<]*)<', PAGE.read_text())
    passes = [t for t in tags if t[0] == "pass" or t[1].strip().upper() == "PASSED"]
    assert not passes, f"pass tags are back on the page: {passes}"
    assert len(tags) == 12, f"expected twelve angle tags on the page, found {len(tags)}"


def test_the_dashboard_claims_no_pass_in_its_verdict():
    text = DASHBOARD.read_text()
    verdict = next(l for l in text.splitlines() if "Honest Verdict" in l)
    assert "no passes" in verdict, f"dashboard verdict still claims a pass: {verdict.strip()[:160]}"
    assert "1 Passed in Silico" not in verdict


def _classification_nodes():
    """Every `"classification": ...` value in an angle function, with its angle."""
    tree = ast.parse(BATTERY.read_text())
    out = []
    for fn in [n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef)]:
        m = re.match(r"evaluate_angle(\d+)_", fn.name)
        if not m:
            continue
        for node in ast.walk(fn):
            if not isinstance(node, ast.Dict):
                continue
            for key, value in zip(node.keys, node.values):
                if isinstance(key, ast.Constant) and key.value == "classification":
                    out.append((int(m.group(1)), value))
    return out


def test_no_angle_function_writes_an_unconditional_pass():
    """A literal label is not read off the gate it names.

    Angles 11 and 12 computed their gate booleans and then wrote
    "PASSED_IN_SILICO" as a constant, so the label would have read as a pass had
    both criteria failed. A literal is still allowed where it states the angle's
    own finding, such as AT_RISK on the protease check, but not where it claims
    a pass.
    """
    offenders = []
    for angle, value in _classification_nodes():
        if isinstance(value, ast.Constant) and value.value == "PASSED_IN_SILICO":
            offenders.append(angle)
    assert not offenders, (
        f"angles {sorted(set(offenders))} write PASSED_IN_SILICO as a literal, so "
        "the label does not depend on the gate it names")


def test_angle_12_reads_the_gates_it_computes():
    """donnan_ok and switch_ok were computed and discarded. They decide the label now."""
    values = [v for a, v in _classification_nodes() if a == 12]
    assert len(values) == 1, f"expected one classification in angle 12, found {len(values)}"
    src = ast.unparse(values[0])
    assert "donnan_ok" in src and "switch_ok" in src, (
        f"angle 12's label does not consult its own gate booleans: {src}")
    assert "COMPOSITION_ARITHMETIC_ONLY" in src


@pytest.mark.parametrize("angle,substring", [(11, "MINIMISATION_ONLY_NO_DYNAMICS"),
                                             (12, "COMPOSITION_ARITHMETIC_ONLY")])
def test_the_withdrawn_passes_say_what_was_actually_done(angle, substring):
    values = [v for a, v in _classification_nodes() if a == angle]
    assert any(substring in ast.unparse(v) for v in values), (
        f"angle {angle} no longer carries {substring}")
