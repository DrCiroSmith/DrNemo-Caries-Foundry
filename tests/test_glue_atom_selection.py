"""Every distance routine in this repository must drop Meeko's glue pseudo-atoms.

Meeko opens a macrocycle so its ring torsions can be sampled and closes it with
glue pseudo-atoms, typed G, G0, G1 and so on. They are bookkeeping and not
atoms. 1183 of the 1244 committed pose files carry a pair of them, so a
selection that keeps them inflates every distance computed from a macrocyclic
ligand: it is the whole difference between the earlier 9.24 A pose RMSD and the
9.11 A this project now reports, and between 5.58 A and 5.63 A for the mode
gap.

One test per number would not have caught this, because the defect is in the
atom selection rather than in any one figure. So this file tests the selection
itself, two ways:

  - behaviourally, by feeding each reader a pose that contains a glue atom and
    checking the reader drops it. This is the test that matters.
  - statically, as a backstop, so a new module that reads coordinates and
    computes a distance cannot be added without a glue exclusion.
"""

import ast
import textwrap
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]

GLUE_TYPES = ("G", "G0", "G1", "G2", "G3")

# A three-atom ligand: one carbon, one polar hydrogen, one glue pseudo-atom.
# A correct reader keeps the carbon alone. Columns match the PDBQT the pose
# files use, with the AutoDock type in 77 to 79.
SYNTHETIC_POSE = """MODEL 1
REMARK VINA RESULT:      -8.191      0.000      0.000
ATOM      1  C   CYS A   1       1.000   2.000   3.000  1.00  0.00    +0.152 C
ATOM      2  HD  CYS A   1      50.000  50.000  50.000  1.00  0.00    +0.100 HD
ATOM      3  G   CYS A   1     900.000 900.000 900.000  1.00  0.00    +0.000 G0
ENDMDL
MODEL 2
ATOM      1  C   CYS A   1       4.000   5.000   6.000  1.00  0.00    +0.152 C
ATOM      2  HD  CYS A   1      50.000  50.000  50.000  1.00  0.00    +0.100 HD
ATOM      3  G   CYS A   1     900.000 900.000 900.000  1.00  0.00    +0.000 G0
ENDMDL
"""

# Each reader: module, function, and how to get coordinates out of its result.
READERS = (
    ("scripts/run_exhaustive_disproof.py", "extract_heavy_atom_coords_from_pdbqt",
     lambda r: np.asarray(r[1], dtype=float)),
    ("scripts/run_real_campaign.py", "parse_heavy_atoms_pdbqt",
     lambda r: np.asarray([a[1] for a in r], dtype=float)),
    ("scripts/compute_redocking_rmsd.py", "read_atoms",
     lambda r: np.asarray(r[0], dtype=float)),
    ("scripts/compute_paper_numbers.py", "pdbqt_best_mode",
     lambda r: np.asarray([a[1:] for a in r], dtype=float)),
)

# Modules that read coordinates but report no distance, with the reason. A
# centroid used to place a ligand in the docking box is a build decision and
# changing it would change future runs, so it is named here rather than altered.
STATIC_EXEMPT = {
    "biophysics_engine.py":
        "its centroid translates a conformer into the search box and is not a "
        "reported distance; changing the selection would change future runs",
    "scripts/crop_and_clean_pdb.py":
        "crops a receptor and computes no RMSD or centroid distance",
    "scripts/build_educational_sops.py":
        "emits documentation text and computes no distance itself",
}


def _load_function(relative, name):
    """Exec one function out of a module, without importing the module.

    The campaign scripts import rdkit and openmm at module scope, which are not
    installed here, so the function is lifted out of the source instead. That
    keeps the test running in every environment, which matters for a test whose
    whole job is to stop a defect coming back.
    """
    source = (REPO / relative).read_text()
    tree = ast.parse(source)
    wanted = [n for n in tree.body
              if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))
              and n.name == name]
    assert wanted, f"{relative} no longer defines {name}"

    # Module-level constants the function may close over, such as GLUE_TYPES.
    namespace = {"np": np, "numpy": np, "Path": Path,
                 "Tuple": tuple, "List": list, "Dict": dict, "Any": object,
                 "Optional": object}
    for node in tree.body:
        if isinstance(node, (ast.Assign, ast.AnnAssign)):
            try:
                exec(ast.unparse(node), namespace)
            except Exception:
                pass  # constants needing absent imports are not needed here
    # Sibling helpers the reader may call, such as an element-from-type helper.
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            try:
                exec(ast.unparse(node), namespace)
            except Exception:
                pass  # functions needing absent imports are not on this path

    assert name in namespace, (
        f"{relative}:{name} could not be lifted out of its module; the test "
        f"needs updating rather than the exclusion being assumed")
    return namespace[name]


@pytest.fixture
def pose_with_a_glue_atom(tmp_path):
    path = tmp_path / "synthetic_with_glue.pdbqt"
    path.write_text(SYNTHETIC_POSE)
    return path


@pytest.mark.parametrize("relative,name,extract", READERS,
                         ids=[f"{r[0]}:{r[1]}" for r in READERS])
def test_each_reader_drops_the_glue_pseudo_atom(relative, name, extract,
                                                pose_with_a_glue_atom):
    reader = _load_function(relative, name)
    # A Path, because some readers call .read_text() on the argument directly.
    coords = extract(reader(pose_with_a_glue_atom))

    assert len(coords) == 1, (
        f"{relative}:{name} returned {len(coords)} atoms from a ligand with one "
        f"carbon, one polar hydrogen and one glue pseudo-atom. It should return "
        f"the carbon alone; a count of 2 means the glue atom is being counted as "
        f"an atom, which inflates every distance computed from a macrocycle")
    assert np.allclose(coords[0], [1.0, 2.0, 3.0]), (
        f"{relative}:{name} kept the wrong atom: {coords[0]}")
    # The two rejects sit at coordinates nothing else occupies, so a leak names
    # itself. The hydrogen is checked because widening a two-column type read to
    # catch the glue type turns HD into "HD" and stops excluding hydrogens: that
    # regression was introduced here once and this is what caught it.
    assert not np.any(np.all(np.isclose(coords, 50.0), axis=1)), (
        f"{relative}:{name} kept the polar hydrogen")
    assert not np.any(np.all(np.isclose(coords, 900.0), axis=1)), (
        f"{relative}:{name} kept the glue pseudo-atom")


def test_the_glue_atoms_are_really_there_to_be_dropped():
    """Guards the test above: if the fixture had no glue atom it would pass idly."""
    lines = [l for l in SYNTHETIC_POSE.splitlines() if l.startswith("ATOM")]
    types = [l[76:79].strip().upper() for l in lines]
    assert any(t in GLUE_TYPES for t in types), (
        "the synthetic pose contains no glue pseudo-atom, so the readers are "
        "not being tested for anything")
    assert "C" in types, "the synthetic pose contains no real heavy atom"


def _reads_coordinates(source):
    return "30:38" in source and "38:46" in source


def _reports_a_distance(source):
    lowered = source.lower()
    return ("rmsd" in lowered
            or "centroid_displacement" in lowered
            or "centroid_sd" in lowered
            or "centroid_distance" in lowered)


def _excludes_glue(source):
    return any(f'"{t}"' in source or f"'{t}'" in source for t in GLUE_TYPES)


def test_no_module_reads_coordinates_for_a_distance_without_excluding_glue():
    """The backstop. A new reader cannot be added without the exclusion."""
    offenders = []
    for path in sorted(REPO.rglob("*.py")):
        if any(part in (".git", "__pycache__", ".venv", "node_modules")
               for part in path.parts):
            continue
        relative = str(path.relative_to(REPO))
        if relative in STATIC_EXEMPT:
            continue
        source = path.read_text()
        if not (_reads_coordinates(source) and _reports_a_distance(source)):
            continue
        if not _excludes_glue(source):
            offenders.append(relative)
    assert not offenders, (
        f"these modules read PDBQT coordinates and report a distance but never "
        f"exclude Meeko's glue pseudo-atoms: {offenders}. Either exclude them, "
        f"or add the module to STATIC_EXEMPT in this file with the reason")


@pytest.mark.parametrize("relative", sorted(STATIC_EXEMPT))
def test_every_exemption_still_points_at_a_real_file(relative):
    """An exemption that outlives its file would hide the next defect."""
    assert (REPO / relative).exists(), (
        f"{relative} is exempted in this file but no longer exists; remove the "
        f"exemption so it cannot shelter a different module later")
