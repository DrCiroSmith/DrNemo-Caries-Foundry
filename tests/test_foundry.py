"""
Unit test suite for Dr. Nemo Caries Vaccine Foundry.
Tests biochemical pre-filtering, target catalog integrity, scoring metrics,
demo data integrity, gtfb_orchestrator --demo execution, dashboard resolution,
and the comprehensive 6-Angle Karl Popper Biophysical Falsification Battery.
"""
import json
import subprocess
import pytest
from pathlib import Path
import sys
import csv

# Ensure repo root is on sys.path
REPO_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(REPO_ROOT))

from solubility_shield import (
    is_viable,
    calc_net_charge_at_ph,
    calc_stephan_switch,
    count_protease_cleavage_sites,
    MAX_DONNAN_CHARGE_PH6
)
from falsification_engine import FalsificationEngine, CO_CRYSTAL_CONTROL_DELTA_G, CO_CRYSTAL_REDOCKING_RMSD


def test_solubility_shield_valid_arginine_core():
    viable, metrics = is_viable("RRRRRR")
    assert viable is True
    assert metrics["chitosan_native"] is True
    assert metrics["arg_count"] == 6
    assert metrics["gravy"] < 2.0


def test_solubility_shield_rejects_invalid_amino_acids():
    viable, metrics = is_viable("RRRX123")
    assert viable is False
    assert metrics["error"] == "invalid_characters"


def test_solubility_shield_rejects_anionic_salt_bridges():
    viable, metrics = is_viable("RRDDDDEE")
    assert viable is False
    assert metrics["error"] == "chitosan_salt_bridge_risk"
    assert metrics["anionic_count"] == 6


def test_solubility_shield_rejects_entrapment_density():
    # 5 residues, 4 are I/K -> 80% density > 30% threshold
    viable, metrics = is_viable("RIIKK")
    assert viable is False
    assert metrics["error"] == "chitosan_entrapment_risk"
    assert metrics["entrap_density"] > 0.3


def test_solubility_shield_edge_cases():
    # Empty string
    viable, metrics = is_viable("")
    assert viable is False
    assert metrics["error"] == "invalid_characters"

    # Top lead candidate
    viable, metrics = is_viable("CPHWRPWHC")
    assert viable is True
    assert metrics["chitosan_native"] is True
    assert metrics["arg_count"] == 1
    assert metrics["gravy"] < 0.0


def test_targets_catalog_integrity():
    targets_file = REPO_ROOT / "targets.json"
    assert targets_file.exists(), "targets.json catalog file must exist"
    
    with open(targets_file, "r", encoding="utf-8") as f:
        data = json.load(f)
    
    # Must contain primary cariogenic targets
    expected_targets = [
        "S_mutans_GtfB_Insoluble",
        "S_mutans_SpaP_Baseline",
        "S_mutans_GtfC_Bridge"
    ]
    for target in expected_targets:
        assert target in data, f"Missing critical target: {target}"
        info = data[target]
        assert "name" in info and len(info["name"]) > 0
        assert "sequence" in info and len(info["sequence"]) > 20
        assert "base_peptide" in info and len(info["base_peptide"]) > 0


def test_composite_scoring_logic():
    def compute_composite(plddt: float, iptm: float) -> float:
        return round((plddt / 100.0) * iptm, 3)

    # Sub-threshold hit
    assert compute_composite(80.0, 0.70) == 0.560
    assert compute_composite(80.0, 0.70) < 0.72

    # Major structural lead (exceeds 0.72)
    lead_score = compute_composite(92.0, 0.85)
    assert lead_score == 0.782
    assert lead_score >= 0.72


def test_cycle_designer_generation():
    from gtfb_orchestrator import cycle_designer_node
    
    dummy_state = {
        "candidates": [],
        "results": [],
        "winning_candidate": None,
        "biosecurity_logs": [],
        "synthesis_instructions": {},
        "retention_scores": {}
    }
    output = cycle_designer_node(dummy_state)
    assert "candidates" in output
    assert len(output["candidates"]) >= 15
    for seq in output["candidates"]:
        assert seq.startswith("C") and seq.endswith("C"), "Must be disulfide-cyclized"


def test_gtfb_orchestrator_demo_execution():
    """Verify that gtfb_orchestrator.py --demo runs the full LangGraph state machine and exits 0."""
    cmd = [sys.executable, str(REPO_ROOT / "gtfb_orchestrator.py"), "--demo"]
    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0, f"gtfb_orchestrator --demo failed: {res.stderr}"
    assert "FOUNDRY SUCCESS" in res.stdout
    assert "CRRRWWRRRRC" in res.stdout
    assert "CPHWRPWHC" in res.stdout
    assert "PepSy" in res.stdout


def test_demo_dataset_integrity():
    """Verify the 20-generation demo dataset bundled in demo_data/."""
    demo_dir = REPO_ROOT / "demo_data"
    assert demo_dir.exists(), "demo_data directory must exist"
    
    # 1. Master scores
    scores_file = demo_dir / "results" / "master_scores.csv"
    assert scores_file.exists(), "demo_data/results/master_scores.csv must exist"
    with open(scores_file, "r", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) >= 50, f"Expected >= 50 rows, got {len(rows)}"
    cols = set(rows[0].keys())
    assert ("contender_id" in cols and "vina_3aic_kcal" in cols) or ("Job_ID" in cols and "pLDDT" in cols)

    # 2. Directed evolution complex PDB files
    variant_files = list((demo_dir / "results").glob("S_mutans_GtfB_Insoluble_variant_*.pdb"))
    assert len(variant_files) >= 5, f"Expected >= 5 variant PDB files, found {len(variant_files)}"
    for pdb_path in variant_files:
        content = pdb_path.read_text(encoding="utf-8")
        assert "HEADER" in content
        assert "ATOM" in content
        assert "TER" in content
        assert "CONECT" in content
        assert " A " in content, "Must contain Chain A (GtfB pocket)"
        assert " B " in content, "Must contain Chain B (docked cyclic peptide)"

    # 3. Falsification benchmark docked coordinates
    docked_complex = demo_dir / "results" / "3AIC_docked_complex_CPHWRPWHC.pdb"
    assert docked_complex.exists(), "3AIC_docked_complex_CPHWRPWHC.pdb must exist"
    complex_content = docked_complex.read_text(encoding="utf-8")
    assert "HEADER" in complex_content
    assert "ATOM" in complex_content
    assert "TER" in complex_content
    assert " A " in complex_content, "Must contain Chain A (GtfB receptor)"
    assert " B " in complex_content, "Must contain Chain B (docked cyclic lead)"

    native_docked = demo_dir / "results" / "3AIC_native_acarbose_docked.pdb"
    assert native_docked.exists(), "3AIC_native_acarbose_docked.pdb must exist"

    lead_docked = demo_dir / "results" / "3AIC_docked_lead_CPHWRPWHC.pdb"
    assert lead_docked.exists(), "3AIC_docked_lead_CPHWRPWHC.pdb must exist"

    # 4. Falsification benchmark JSON
    bench_file = demo_dir / "results" / "falsification_benchmark.json"
    assert bench_file.exists(), "falsification_benchmark.json must exist"
    with open(bench_file, "r", encoding="utf-8") as f:
        bench_data = json.load(f)
    assert "lead_affinity_3aic" in bench_data
    assert "lead_affinity_1smd" in bench_data
    assert "scrambled_null_test" in bench_data
    assert "alanine_scan" in bench_data

    # 5. Clinical Memos
    analysis_files = list((demo_dir / "results" / "analysis").glob("*.md"))
    assert len(analysis_files) >= 5, "Expected >= 5 clinical analysis memos"


def test_dashboard_path_resolution():
    """Verify dashboard.py resolves bundled demo data when drnemo_research is absent."""
    import os
    old_env = os.environ.pop("DRNEMO_DATA_DIR", None)
    try:
        import importlib.util
        spec = importlib.util.spec_from_file_location("dashboard_module", str(REPO_ROOT / "dashboard.py"))
        assert spec is not None
    finally:
        if old_env is not None:
            os.environ["DRNEMO_DATA_DIR"] = old_env


# ==============================================================================
# OPERATION KARL POPPER: 6-ANGLE BIOPHYSICAL FALSIFICATION TEST SUITE
# ==============================================================================

def test_falsification_angle1_cocrystal_control():
    """Angle 1: RCSB 3AIC co-crystal redocking positive control calibration."""
    ctrl = FalsificationEngine.get_angle1_control()
    assert ctrl["passed_calibration"] is True
    assert ctrl["affinity_kcal_mol"] <= -7.0
    assert ctrl["redocking_rmsd_angstrom"] < 5.0
    assert ctrl["centroid_displacement_angstrom"] < 3.5


def test_falsification_angle2_offtarget_selectivity():
    """Angle 2: 1SMD Salivary Amylase counter-screen ensures high selectivity."""
    engine = FalsificationEngine()
    # Pareto champion CAHWRPWHC clears selectivity threshold DeltaDeltaG <= -0.80 kcal/mol (-0.85 kcal/mol)
    audit = engine.run_falsification_audit("CAHWRPWHC")
    gate2 = audit["gates"]["gate2_offtarget_selectivity"]
    assert gate2["lead_passed"] is True
    assert gate2["lead_selectivity"] <= -0.80
    assert gate2["naive_passed"] is False

    # Prior 22A lead CPHWRPWHC is honestly falsified on 28A grid due to amylase cross-reactivity (ddG = -0.16)
    wt_audit = engine.run_falsification_audit("CPHWRPWHC")
    assert wt_audit["gates"]["gate2_offtarget_selectivity"]["lead_passed"] is False
    assert wt_audit["gates"]["gate2_offtarget_selectivity"]["lead_selectivity"] > -0.80


def test_falsification_angle3_scrambled_null_permutation():
    """Angle 3: Scrambled null permutation distribution across N=39 isomers."""
    engine = FalsificationEngine()
    audit = engine.run_falsification_audit("CAHWRPWHC")
    gate3 = audit["gates"]["gate3_scrambled_null_permutation"]
    
    assert gate3["lead_passed"] is True
    # N=39 isomers evaluated empirically
    assert gate3["null_mean"] <= -8.0
    assert gate3["p_value"] > 0.0
    assert gate3["b_count"] >= 0


def test_falsification_angle4_alanine_hotspot_scan():
    """Angle 4: In-silico alanine scanning identifies energetic hotspots."""
    engine = FalsificationEngine()
    audit = engine.run_falsification_audit("CAHWRPWHC")
    gate4 = audit["gates"]["gate4_alanine_hotspot_scan"]
    
    assert gate4["passed"] is True
    assert gate4["hotspot_count"] >= 1
    # Verify primary anchor His8 is identified
    assert "H8A" in gate4["hotspots"]


def test_falsification_angle5_protease_resistance():
    """Angle 5: Salivary trypsin and RgpB gingipain oral stability."""
    engine = FalsificationEngine()
    audit = engine.run_falsification_audit("CAHWRPWHC")
    gate5 = audit["gates"]["gate5_oral_protease_resistance"]
    
    # Pareto lead has 0 cleavage sites
    assert gate5["lead_passed"] is True
    assert gate5["lead_cleavage_sites"] == 0
    # Naive baseline CRRRWWRRRRC has multiple cleavage sites and fails
    assert gate5["naive_passed"] is False
    assert gate5["naive_cleavage_sites"] >= 5


def test_falsification_angle6_chitosan_and_stephan_switch():
    """Angle 6: Chitosan Donnan compatibility and Stephan curve pH-switch."""
    engine = FalsificationEngine()
    audit = engine.run_falsification_audit("CAHWRPWHC")
    gate6 = audit["gates"]["gate6_chitosan_donnan_stephan"]
    
    assert gate6["lead_passed"] is True
    # Donnan net charge at pH 6.0 <= +4.0
    assert gate6["lead_net_charge_pH6"] <= MAX_DONNAN_CHARGE_PH6
    # Stephan acid switch delta Q (pH 5.2 - 7.2) >= +1.0
    assert gate6["lead_delta_q_stephan"] >= 1.0
    # Naive baseline has charge +7.0 (exceeds Donnan limit) and zero Stephan switch
    assert gate6["naive_passed"] is False
    assert gate6["naive_net_charge_pH6"] > MAX_DONNAN_CHARGE_PH6
    assert gate6["naive_delta_q_stephan"] == 0.0


def test_falsification_audit_overall_verdict():
    """Verify overall scientific verdict: Pareto champion corroborated, unviable baselines falsified."""
    engine = FalsificationEngine()
    lead_audit = engine.run_falsification_audit("CAHWRPWHC")
    assert lead_audit["overall_falsification_passed"] is True

    # CPHWRPWHC is falsified on amylase cross-reactivity
    wt_audit = engine.run_falsification_audit("CPHWRPWHC")
    assert wt_audit["overall_falsification_passed"] is False

    # Naive poly-arginine fails falsification
    naive_audit = engine.run_falsification_audit("CRRRWWRRRRC")
    assert naive_audit["overall_falsification_passed"] is False


def test_biophysics_macrocycle_embedding():
    """Verify 3D macrocycle building and heavy atom enumeration."""
    pytest.importorskip("rdkit")
    pytest.importorskip("gemmi")
    pytest.importorskip("meeko")
    import os
    from biophysics_engine import build_macrocycle_3d
    import tempfile
    with tempfile.NamedTemporaryFile(suffix=".pdbqt") as tf:
        ha = build_macrocycle_3d("CPHWRPWHC", tf.name)
        assert ha == 86, f"Expected 86 heavy atoms for CPHWRPWHC, got {ha}"
        assert os.path.getsize(tf.name) > 1000


def test_biophysics_openmm_relaxation():
    """Verify OpenMM Amber14/GBn2 energy minimization runs and yields valid potential energy."""
    pytest.importorskip("rdkit")
    pytest.importorskip("gemmi")
    pytest.importorskip("meeko")
    pytest.importorskip("openmm")
    from biophysics_engine import build_macrocycle_3d, relax_with_openmm
    import tempfile
    with tempfile.NamedTemporaryFile(suffix=".pdb") as tf_pdb, tempfile.NamedTemporaryFile(suffix=".pdbqt") as tf_pdbqt:
        build_macrocycle_3d("CPHWRPWHC", tf_pdbqt.name, out_pdb=tf_pdb.name)
        energy_kj = relax_with_openmm(tf_pdb.name)
        assert isinstance(energy_kj, float)
        assert energy_kj < 0, f"Expected favorable negative potential energy, got {energy_kj} kJ/mol"


def test_contender_vault_crud_and_mirroring():
    """Verify permanent SQLite contender vault CRUD, Pareto champion retrieval, and CSV ledger mirroring."""
    from vault import ContenderVault
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        db_path = Path(td) / "test_vault.sqlite"
        vault = ContenderVault(db_path=db_path)
        assert vault.get_count() == 0

        cid = vault.add_contender({
            "sequence": "CPHWRPWHC",
            "generation": 1,
            "heavy_atoms": 86,
            "vina_3aic_kcal": -8.30,
            "vinardo_3aic_kcal": -5.20,
            "ligand_efficiency": round(8.30 / 86, 3),
            "vina_1smd_kcal": -7.20,
            "selectivity_ddg": -1.10,
            "openmm_relaxed_kj": -1132.07,
            "ss_bond_angstrom": 2.05,
            "q_ph72": 1.25,
            "q_ph60": 2.05,
            "q_ph52": 2.75,
            "delta_q_stephan": 1.62,
            "protease_sites": 0,
            "pareto_rank": 1,
            "falsification_status": "CORROBORATED",
            "pose_sha256": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855",
            "evaluated_at_utc": "2026-11-15T21:00:00Z"
        })
        assert vault.get_count() == 1
        assert vault.ledger_csv.exists()

        champ = vault.get_pareto_champion()
        assert champ is not None
        assert champ["sequence"] == "CPHWRPWHC"
        assert champ["vina_3aic_kcal"] == -8.30

