"""
Unit test suite for Dr. Nemo Caries Vaccine Foundry.
Tests biochemical pre-filtering, target catalog integrity, scoring metrics,
and orchestrator state definitions.
"""
import json
import pytest
from pathlib import Path
import sys

# Ensure repo root is on sys.path
REPO_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(REPO_ROOT))

from solubility_shield import is_viable


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


def test_cycle_designer_generation_mock():
    # Simulates generation conforming to chitosan-native library rules
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
    assert len(output["candidates"]) == 15
    for seq in output["candidates"]:
        assert seq.startswith("C") and seq.endswith("C"), "Must be disulfide-cyclized"
