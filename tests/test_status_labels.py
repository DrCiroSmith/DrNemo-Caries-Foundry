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


# ---------------------------------------------------------------------------
# The seed claim
#
# The paper and README say no published score can be regenerated because the
# production pipeline leaves Vina's seed unset. That claim was written when
# biophysics_engine.py had no seed argument at all. It has one now, added by
# this review's provenance commit, so the claim had to narrow from "the engine
# passes no seed" to "every call site leaves it unset". This test pins the
# narrowed version, including the count the text prints.
# ---------------------------------------------------------------------------

ENGINE = REPO / "biophysics_engine.py"
PAPER = REPO / "paper" / "negative_result_GTF-SI.md"


def _seed_parameter_index():
    """The seed's position in the docking function's own signature.

    Read rather than assumed, so that reordering the parameter list cannot make
    a positional caller look seedless.
    """
    for path in sorted(REPO.rglob("*.py")):
        if path.relative_to(REPO).parts[0] == ".git":
            continue
        for node in ast.walk(_parse(path)):
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and \
                    node.name == "run_docking_vina_and_vinardo":
                names = [a.arg for a in node.args.posonlyargs + node.args.args]
                if "seed" in names:
                    return names.index("seed")
    raise AssertionError(
        "run_docking_vina_and_vinardo has no seed parameter, so the paper's "
        "claim about call sites leaving it unset no longer describes the code")


def _parse(path):
    """Parse a file, reporting a parse failure instead of skipping the file.

    An earlier version of this helper caught SyntaxError and continued, which
    meant a file that would not parse was counted as clean. A sibling audit
    found the same hole in its own walker by injecting a seeded call with
    invalid syntax and watching the test stay green. A silent skip in an audit
    is indistinguishable from a passing file, so it fails here instead.
    """
    try:
        return ast.parse(path.read_text(encoding="utf-8", errors="replace"))
    except SyntaxError as exc:
        raise AssertionError(
            f"{path.relative_to(REPO)} does not parse, so this audit cannot "
            f"say whether it passes a seed: {exc}") from exc


def _docking_call_sites(scope="production"):
    """Every call of run_docking_vina_and_vinardo, with whether it passes a seed.

    scope "production" excludes tests/ and the separate retarget_gtfb protocol;
    scope "all" is what an auditor's grep of the repository finds, which is the
    number the paper has to survive.
    """
    assert scope in ("production", "all")
    seed_index = _seed_parameter_index()
    sites = []
    for path in sorted(REPO.rglob("*.py")):
        parts = path.relative_to(REPO).parts
        if parts[0] == ".git":
            continue
        if scope == "production" and parts[0] in ("tests", "retarget_gtfb"):
            continue
        for node in ast.walk(_parse(path)):
            if not isinstance(node, ast.Call):
                continue
            name = getattr(node.func, "id", None) or getattr(node.func, "attr", None)
            if name != "run_docking_vina_and_vinardo":
                continue
            passes = any(k.arg == "seed" for k in node.keywords if k.arg)
            if len(node.args) > seed_index:
                passes = True
            sites.append((str(path.relative_to(REPO)), node.lineno, passes))
    return sites


def test_no_production_call_site_passes_a_seed():
    seeded = [f"{f}:{n}" for f, n, passes in _docking_call_sites() if passes]
    assert not seeded, (
        "a production call site now passes a seed, so the claim that no published "
        f"score can be regenerated needs rewriting: {seeded}")


WORDS = {1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six",
         7: "seven", 8: "eight", 9: "nine", 10: "ten", 11: "eleven",
         12: "twelve", 13: "thirteen", 14: "fourteen", 15: "fifteen",
         16: "sixteen"}


def test_the_paper_prints_the_right_number_of_call_sites():
    """The text names the count, so the count has to be checked and not recalled.

    This guard formerly compared the paper against the production count alone,
    which is the narrower thing it measures, while the sentence said "all
    twelve of its call sites". The repository holds fifteen, two of which pass
    a seed, so the sentence was false to anyone who grepped for it even though
    this test was green. The scope of the claim is now checked, not just its
    digit.
    """
    paper = PAPER.read_text()
    production = len(_docking_call_sites("production"))
    assert production in WORDS, f"{production} sites, extend the map"
    assert f"all {WORDS[production]} production call sites leave it unset" in paper, (
        f"there are {production} production call sites and the paper does not "
        f"say {WORDS[production]} production")


def test_the_paper_does_not_claim_more_than_production_of_the_unset_seed():
    """A sentence about every call site is false while test fixtures pass one."""
    paper = PAPER.read_text()
    for overclaim in ("all twelve of its call sites",
                      "every call site leaves it unset",
                      "no call site passes a seed"):
        assert overclaim not in paper, (
            f"the paper says {overclaim!r}, which the test fixtures at "
            f"{[f'{f}:{n}' for f, n, p in _docking_call_sites('all') if p]} "
            "make false")


def test_the_paper_discloses_the_seeded_test_fixtures():
    """What a grep of the repository finds has to be in the paper already."""
    every = _docking_call_sites("all")
    seeded = [(f, n) for f, n, passes in every if passes]
    assert seeded, (
        "no call site passes a seed any more, so the paper's disclosure of the "
        "seeded fixtures now describes code that is gone and should be cut")
    paper = PAPER.read_text()
    assert len(every) in WORDS, f"{len(every)} sites, extend the map"
    files = {f for f, _ in seeded}
    assert len(files) == 1, f"seeded sites now span {files}, so the text needs rewriting"
    only = files.pop()
    assert only in paper, (
        f"{len(seeded)} call sites in {only} pass a seed and the paper does not "
        "name the file, so a reader who greps finds seed= and reads the "
        "surrounding claim as false")
    outside = len(every) - len(_docking_call_sites("production"))
    assert outside in WORDS, f"{outside} sites outside production, extend the map"
    assert f"{WORDS[outside]} further call sites" in paper.lower(), (
        "the paper does not account for the call sites outside production")


def test_the_engine_still_only_seeds_when_asked():
    """The pass-through is conditional. If it ever becomes unconditional the
    pipeline is deterministic and the paper's section 3.2 is wrong."""
    source = ENGINE.read_text()
    assert "if seed is not None:" in source, (
        "biophysics_engine no longer guards its --seed pass-through, so either the "
        "pipeline is now seeded or the seed is passed as None")


ENGINE = REPO / "biophysics_engine.py"
BENCHMARK = REPO / "FALSIFICATION_BENCHMARK.md"


def _mmff_call_uses_its_return_value():
    """Whether the builder does anything with MMFFOptimizeMolecule's result.

    RDKit documents the return as 0 converged, -1 force field could not be set
    up, 1 more iterations required, so discarding it hides both failure modes.
    A bare expression statement discards it; an assignment or a comparison does
    not.
    """
    tree = _parse(ENGINE)
    found = False
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        name = getattr(node.func, "attr", None) or getattr(node.func, "id", None)
        if name != "MMFFOptimizeMolecule":
            continue
        found = True
        for stmt in ast.walk(tree):
            if isinstance(stmt, ast.Expr) and stmt.value is node:
                break
        else:
            return True
    assert found, (
        "biophysics_engine.py no longer calls MMFFOptimizeMolecule, so the "
        "paper's description of the minimisation step is stale")
    return False


def test_the_paper_matches_whether_the_minimisation_outcome_is_checked():
    """The prose and the code have to move together.

    If the builder is fixed to read the return value, this fails so that the
    paper stops saying the outcome is unrecorded. If the builder still discards
    it, this fails unless the paper says so. Today's lesson was a guard that
    enforced what I meant rather than what the sentence claimed, so this one is
    pinned to the code's actual shape.
    """
    paper = PAPER.read_text()
    present = "the builder discards the optimiser's return value"
    past = "discarded the optimiser's return value"
    if _mmff_call_uses_its_return_value():
        assert present not in paper, (
            "the builder now reads MMFFOptimizeMolecule's return value, so the "
            "paper must not say in the present tense that it is discarded")
        assert past in paper, (
            "the builder now reads the return value, but the conformers in the "
            "paper were built by code that did not. The paper has to say that in "
            "the past tense, or it implies the deposited structures carry an "
            "outcome that nothing records")
    else:
        assert present in paper, (
            "the builder still discards MMFFOptimizeMolecule's return value and "
            "the paper no longer says so, which overstates the minimisation step")
    assert "UFF" in paper, (
        "the paper no longer mentions the UFF fallback that the code's "
        "except clause contains")


def _numeric(node):
    """Whether an AST node is a plain number, including a negated one."""
    if isinstance(node, ast.Constant):
        return isinstance(node.value, (int, float)) and not isinstance(node.value, bool)
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.USub, ast.UAdd)):
        return _numeric(node.operand)
    return False


def test_no_failure_handler_substitutes_a_number_for_a_failed_measurement():
    """A failure has to be recorded as absent, never as a plausible value.

    run_campaign.py wrote openmm_kj = -950.0 in its except clause, and -950.0 is
    inside the range of real relaxed energies in the ledger, so a failed
    minimisation would have entered the ledger as an ordinary measurement. No
    committed row carries it, but the path existed. This fails on any except
    clause anywhere that assigns a bare number to a name holding a measured
    quantity.
    """
    measured = ("openmm", "vina", "vinardo", "kcal", "_kj", "rmsd", "drift",
                "energy", "affinity")
    offenders = []
    for path in sorted(REPO.rglob("*.py")):
        if path.relative_to(REPO).parts[0] in (".git", "retarget_gtfb"):
            continue
        for node in ast.walk(_parse(path)):
            if not isinstance(node, ast.ExceptHandler):
                continue
            for stmt in ast.walk(node):
                if not isinstance(stmt, ast.Assign):
                    continue
                for target in stmt.targets:
                    name = getattr(target, "id", None)
                    if name is None or not any(m in name.lower() for m in measured):
                        continue
                    if _numeric(stmt.value):
                        offenders.append(
                            f"{path.relative_to(REPO)}:{stmt.lineno} sets "
                            f"{name} = {ast.unparse(stmt.value)} on failure")
    assert not offenders, (
        "a failed measurement is being written as a number, which is "
        "indistinguishable from a real one downstream: " + "; ".join(offenders))


def test_the_builder_fails_loudly_when_embedding_fails():
    """A conformer-less molecule must not reach the writers.

    EmbedMolecule returns -1 on failure and adds no conformer. The retry's
    return was assigned and never read, so a double failure fell through to the
    PDB writer and to Meeko with nothing embedded. RDKit is not installed here,
    so what those two do with such a molecule is untested; the fix is to not
    reach them.
    """
    raises_on_embed_failure = False
    for node in ast.walk(_parse(ENGINE)):
        if not isinstance(node, ast.If):
            continue
        names = {n.id for n in ast.walk(node.test) if isinstance(n, ast.Name)}
        if "embed_res" not in names:
            continue
        if any(isinstance(s, ast.Raise) for s in ast.walk(node)):
            raises_on_embed_failure = True
    assert raises_on_embed_failure, (
        "biophysics_engine.py no longer raises on a failed embedding, so a "
        "molecule with no conformer can reach the PDB writer and Meeko")


def test_no_document_states_the_minimisation_as_a_step_count():
    """maxIterations is a cap, and six surfaces called it a step count.

    OpenMM's Simulation.minimizeEnergy documents maxIterations as "The maximum
    number of iterations to perform", documents no return value, and exposes no
    convergence status afterwards. The minimiser also halts early on its force
    tolerance. So "200-step relaxed" states a step count that nothing can
    report, for every ligand in the table.
    """
    banned = ("200-step", "200 steps", "200 step ")
    for path in (PAPER, README, PAGE, DASHBOARD, REPO / "biophysics_engine.py",
                 REPO / "gtfb_orchestrator.py", BENCHMARK,
                 REPO / "DOCUMENTATION.md"):
        if not path.exists():
            continue
        text = path.read_text(encoding="utf-8", errors="replace")
        for phrase in banned:
            if phrase in text and "not describe the output as 200-step" in text:
                continue
            assert phrase not in text, (
                f"{path.relative_to(REPO)} calls the minimisation {phrase!r}, "
                "but maxIterations=200 is an upper bound and OpenMM reports "
                "neither the iterations taken nor whether it converged")


def test_no_document_quotes_a_drift_figure_the_report_says_does_not_exist():
    """Angle 11's 0.92 A drift was in the README and in no artefact.

    The committed report's own words for that angle are that no relaxation was
    run and no drift figure exists. A number quoted beside that is unsourced.
    """
    import json
    report = json.loads((REPO / "data" / "exhaustive_disproof_report.json")
                        .read_text(encoding="utf-8"))
    angle11 = [a for a in report["angles"] if a.get("angle") == 11]
    assert angle11, "angle 11 is no longer in the committed report"
    recorded = angle11[0].get("heavy_atom_drift_rmsd_angstrom")
    if recorded is not None:
        return
    for path in (README, PAGE, PAPER, DASHBOARD):
        if not path.exists():
            continue
        for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
            if "Relaxation Drift" not in line and "relaxation drift" not in line:
                continue
            assert not re.search(r"\d\.\d+\s*(A|&Aring;|\u00c5)", line), (
                f"{path.relative_to(REPO)} quotes a drift figure for angle 11, "
                "which the committed report states does not exist: " + line.strip())


def test_no_document_states_mmff_minimisation_as_completed():
    """Attempted is not done, and four surfaces used to say done."""
    banned = ("minimised with MMFF94", "minimized with MMFF94",
              "MMFF94-minimised", "MMFF94-minimized",
              "ETKDG/MMFF94 conformer")
    for path in (PAPER, REPO / "README.md", REPO / "docs" / "index.html",
                 BENCHMARK, REPO / "dashboard.py"):
        if not path.exists():
            continue
        text = path.read_text(encoding="utf-8", errors="replace")
        for phrase in banned:
            assert phrase not in text, (
                f"{path.relative_to(REPO)} states MMFF94 minimisation as "
                f"completed ({phrase!r}), but the builder discards the "
                "optimiser's return value and nothing records which force "
                "field ran")


def test_the_benchmark_sheet_does_not_advertise_mmff94_as_a_component():
    """It listed RDKit MMFF94 as a pipeline component that ran."""
    text = BENCHMARK.read_text()
    assert "**RDKit MMFF94:**" not in text, (
        "FALSIFICATION_BENCHMARK.md again lists RDKit MMFF94 as a component "
        "without saying the outcome is unrecorded")
    assert "MMFF94 attempted" in text, (
        "FALSIFICATION_BENCHMARK.md no longer says MMFF94 is only attempted")
