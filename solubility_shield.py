#!/usr/bin/env python3
"""
Solubility Shield v2.0 — Research-Guided Viability Pre-Filter for Dr. Nemo.

Updated with Chitosan Surface Interaction Rules (NotebookLM Research):
- REJECTS: High densities of I, K (Irreversible entrapment via ring H-bonds).
- REJECTS: Anionic residues D, E (Polyelectrolyte salt bridges with chitosan matrix).
- FAVORS: Arginine (R) for electronegative GtfB catalytic pocket binding and matrix repulsion.
"""
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# Thresholds calibrated for oral-environment peptides (saliva pH ~6.8)
MAX_GRAVY = 2.0          
MAX_INSTABILITY = 40.0   
MIN_PI = 3.0             
MAX_PI = 13.0            # Expanded for Arg-rich peptides

# Chitosan Interaction Constraints (NotebookLM 2026-04-19)
MAX_ANIONIC_COUNT = 2    # Max (D + E) before entrapment risk
MAX_ENTRAP_DENSITY = 0.3 # Max (I + K) percentage

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")

def is_viable(sequence: str) -> tuple:
    """
    Check if a peptide sequence is biochemically viable AND chitosan-native.
    """
    seq = sequence.upper().strip()

    if not seq or not all(c in VALID_AA for c in seq):
        return False, {"error": "invalid_characters"}

    # Chitosan-Specific Filter (NotebookLM Research)
    anionic_count = seq.count('D') + seq.count('E')
    entrap_count = seq.count('I') + seq.count('K')
    entrap_density = entrap_count / len(seq) if seq else 0
    arg_count = seq.count('R')

    if anionic_count > MAX_ANIONIC_COUNT:
        print(f"  🛡️ Chitosan Filter: REJECTED {seq} — Anionic (D/E) count {anionic_count} is too high (Salt Bridge Risk)")
        return False, {"error": "chitosan_salt_bridge_risk"}

    if entrap_density > MAX_ENTRAP_DENSITY:
        print(f"  🛡️ Chitosan Filter: REJECTED {seq} — (I/K) Density {entrap_density:.2f} > {MAX_ENTRAP_DENSITY} (Entrapment Risk)")
        return False, {"error": "chitosan_entrapment_risk"}

    try:
        analyzer = ProteinAnalysis(seq)
        gravy = analyzer.gravy()
        instability = analyzer.instability_index()
        pi = analyzer.isoelectric_point()
    except Exception as e:
        return False, {"error": str(e)}

    metrics = {
        "gravy": round(gravy, 3),
        "instability_index": round(instability, 3),
        "isoelectric_point": round(pi, 3),
        "arg_count": arg_count,
        "chitosan_native": True
    }

    if gravy > MAX_GRAVY:
        print(f"  🛡️ Solubility Shield: REJECTED {seq} — GRAVY {gravy:.2f} > {MAX_GRAVY}")
        return False, metrics

    if instability > MAX_INSTABILITY:
        print(f"  🛡️ Solubility Shield: REJECTED {seq} — Instability {instability:.1f} > {MAX_INSTABILITY}")
        return False, metrics

    return True, metrics

if __name__ == '__main__':
    # Test cases
    print(is_viable("RRRRRR")) # High Arg, should pass
    print(is_viable("DDDD"))   # High Anionic, should fail
    print(is_viable("IIII"))   # High Entrapment, should fail
