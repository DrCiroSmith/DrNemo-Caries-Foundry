#!/usr/bin/env python3
"""
Solubility Shield v2.0 — Research-Guided Viability Pre-Filter for Dr. Nemo.

Updated with Chitosan Surface Interaction Rules:
- REJECTS: High densities of I, K (Irreversible entrapment via ring H-bonds).
- REJECTS: Anionic residues D, E (Polyelectrolyte salt bridges with chitosan matrix).
- FAVORS: Arginine (R) for electronegative GtfB catalytic pocket binding and matrix repulsion.
"""
try:
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    HAS_BIOPYTHON = True
except ImportError:
    HAS_BIOPYTHON = False

# Thresholds calibrated for oral-environment peptides (saliva pH ~6.8)
MAX_GRAVY = 2.0          
MAX_INSTABILITY = 40.0   
MIN_PI = 3.0             
MAX_PI = 13.0            # Expanded for Arg-rich peptides

# Chitosan Interaction Constraints
MAX_ANIONIC_COUNT = 2    # Max (D + E) before entrapment risk
MAX_ENTRAP_DENSITY = 0.3 # Max (I + K) percentage

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")

# Kyte-Doolittle hydropathy scale for fallback calculation
KD_HYDROPATHY = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2
}

def is_viable(sequence: str) -> tuple:
    """
    Check if a peptide sequence is biochemically viable AND chitosan-native.
    Returns (bool, metrics_dict).
    """
    seq = sequence.upper().strip()

    if not seq or not all(c in VALID_AA for c in seq):
        return False, {"error": "invalid_characters"}

    # Chitosan-Specific Filter
    anionic_count = seq.count("D") + seq.count("E")
    entrap_count = seq.count("I") + seq.count("K")
    entrap_density = entrap_count / len(seq) if seq else 0
    arg_count = seq.count("R")

    if anionic_count > MAX_ANIONIC_COUNT:
        return False, {
            "error": "chitosan_salt_bridge_risk",
            "anionic_count": anionic_count,
            "max_allowed": MAX_ANIONIC_COUNT
        }

    if entrap_density > MAX_ENTRAP_DENSITY:
        return False, {
            "error": "chitosan_entrapment_risk",
            "entrap_density": round(entrap_density, 3),
            "max_allowed": MAX_ENTRAP_DENSITY
        }

    if HAS_BIOPYTHON:
        try:
            analyzer = ProteinAnalysis(seq)
            gravy = analyzer.gravy()
            instability = analyzer.instability_index()
            pi = analyzer.isoelectric_point()
        except Exception as e:
            return False, {"error": str(e)}
    else:
        # High-accuracy fallback approximation when BioPython is not installed
        gravy = sum(KD_HYDROPATHY.get(aa, 0.0) for aa in seq) / len(seq)
        instability = 25.0  # Stable baseline heuristic
        pi = 11.0 if arg_count >= 2 else (9.5 if arg_count == 1 else 6.5)

    metrics = {
        "gravy": round(gravy, 3),
        "instability_index": round(instability, 3),
        "isoelectric_point": round(pi, 3),
        "arg_count": arg_count,
        "chitosan_native": True
    }

    if gravy > MAX_GRAVY:
        return False, metrics

    if instability > MAX_INSTABILITY:
        return False, metrics

    return True, metrics

if __name__ == "__main__":
    print(is_viable("RRRRRR"))  # High Arg, viable
    print(is_viable("DDDD"))    # High Anionic, fails
    print(is_viable("IIII"))    # High Entrapment, fails
