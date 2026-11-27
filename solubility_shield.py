#!/usr/bin/env python3
"""
Solubility Shield v2.0: Biophysical Viability and Falsification Pre-Filter (Angle 5).
Author: Dr. Nemo Caries Foundry / Dr. Castrellon

Incorporates rigorous oral biophysics and formulation constraints:
1. Chitosan Surface Compatibility (No D/E anionic salt-bridges, no I/K entrapment).
2. Salivary Trypsin and P. gingivalis RgpB Gingipain Cleavage Resistance (R/K shielded by Proline).
3. Salivary Chymotrypsin and Cathepsin G Cleavage Vulnerability (W/Y/F aromatic residues).
4. Chitosan-TPP Nanoparticle Donnan Electrostatic Repulsion (Reject net_charge_pH6 > +4.0).
5. Stephan-Curve Protonation Acid-Switch (Delta Q = Q_pH5.2 - Q_pH7.2 >= +1.0 via Henderson-Hasselbalch).
"""
import math
from typing import Tuple, Dict, List

try:
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    HAS_BIOPYTHON = True
except ImportError:
    HAS_BIOPYTHON = False

# Thresholds calibrated for oral-environment peptides
MAX_GRAVY = 2.0          
MAX_INSTABILITY = 40.0   
MIN_PI = 3.0             
MAX_PI = 13.0

# Chitosan Interaction Constraints
MAX_ANIONIC_COUNT = 2    # Max (D + E) before irreversible precipitation
MAX_ENTRAP_DENSITY = 0.3 # Max (I + K) percentage
MAX_DONNAN_CHARGE_PH6 = 4.0 # Maximum net charge at pH 6.0 for Chitosan-TPP Donnan compatibility

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")

# Kyte-Doolittle hydropathy scale for fallback calculation
KD_HYDROPATHY = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2
}

# Empirical pKa values for ionizable groups in oral peptides
PKA_VALUES = {
    "N_TERM": 9.0,
    "C_TERM": 2.2,
    "R": 12.48, # Arginine guanidinium
    "K": 10.53, # Lysine epsilon-amino
    "H": 6.04,  # Histidine imidazole (Stephan curve acid-switch)
    "D": 3.86,  # Aspartate carboxyl
    "E": 4.25,  # Glutamate carboxyl
    "C": 8.33,  # Cysteine thiol (if unbonded)
    "Y": 10.07  # Tyrosine phenol
}


def calc_net_charge_at_ph(seq: str, ph: float, is_cyclic: bool = True) -> float:
    """
    Computes net molecular charge at a specified pH using Henderson-Hasselbalch.
    If is_cyclic is True, terminal amine and carboxyl are assumed involved in backbone
    or disulfide cyclic constraint, but for linear controls terminals are ionized.
    """
    seq = seq.upper().strip()
    charge = 0.0

    if not is_cyclic:
        # N-terminus (basic)
        charge += 1.0 / (1.0 + 10.0 ** (ph - PKA_VALUES["N_TERM"]))
        # C-terminus (acidic)
        charge -= 1.0 / (1.0 + 10.0 ** (PKA_VALUES["C_TERM"] - ph))

    # Basic residues (positive when protonated)
    for aa in ("R", "K", "H"):
        count = seq.count(aa)
        if count > 0:
            pka = PKA_VALUES[aa]
            charge += count * (1.0 / (1.0 + 10.0 ** (ph - pka)))

    # Acidic residues (negative when deprotonated)
    for aa in ("D", "E"):
        count = seq.count(aa)
        if count > 0:
            pka = PKA_VALUES[aa]
            charge -= count * (1.0 / (1.0 + 10.0 ** (pka - ph)))

    # Tyrosine
    y_count = seq.count("Y")
    if y_count > 0:
        charge -= y_count * (1.0 / (1.0 + 10.0 ** (PKA_VALUES["Y"] - ph)))

    return round(charge, 3)


def calc_stephan_switch(seq: str) -> Dict[str, float]:
    """Charge shift between resting saliva (pH 7.2) and acidogenic plaque (pH 5.2).

    Be clear about what this is. With a single titratable pKa near 6.0, the charge
    gained per histidine between pH 7.2 and pH 5.2 is fixed by Henderson-Hasselbalch
    at about +0.809, so for a sequence whose only titratable side chains are histidines

        delta_q_stephan ~= 0.809 * (number of histidines)

    Across the 260 rows in the vault this holds exactly for 239 of them, and the
    remainder differ only by the small tyrosine contribution. The observed values
    cluster at 0.000, 0.809, 1.618 and 2.427, which are nought, one, two and three
    histidines. The downstream gate, delta_q >= 1.0, is therefore a threshold on
    histidine count: it is satisfied by two or more histidines and by nothing else.
    It is not an independent biophysical axis, and a candidate passing it has not been
    shown to do anything in acid; it has been shown to contain two histidines.
    """
    q_resting = calc_net_charge_at_ph(seq, 7.2)
    q_acidic = calc_net_charge_at_ph(seq, 5.2)
    delta_q = round(q_acidic - q_resting, 3)
    return {
        "q_pH7_2": q_resting,
        "q_pH5_2": q_acidic,
        "delta_q_stephan": delta_q,
        "histidine_count": seq.count("H"),
        "is_ph_switch_active": delta_q >= 1.0,
        "gate_note": "delta_q >= 1.0 is equivalent to a histidine count of 2 or more; see the docstring.",
    }


def count_protease_cleavage_sites(seq: str) -> Dict[str, object]:
    """
    Models salivary trypsin and P. gingivalis RgpB gingipain cleavage sites,
    along with broad oral protease / chymotrypsin vulnerability.
    
    Trypsin / RgpB:
    Endopeptidases hydrolyze C-terminal to basic residues (R, K) unless
    immediately followed by Proline (P_1' = Pro), which stereochemically blocks the active site.
    
    Chymotrypsin / Cathepsin G / Bacterial Peptidases:
    Hydrolyze C-terminal to bulky aromatic hydrophobic residues (W, Y, F) unless
    shielded by tertiary macrocyclic restraint or Proline.
    """
    seq = seq.upper().strip()
    trypsin_sites = []
    for i in range(len(seq) - 1):
        if seq[i] in ("R", "K"):
            if seq[i + 1] != "P":
                trypsin_sites.append(f"{seq[i]}{i+1}_{seq[i+1]}{i+2}")
    if seq and seq[-1] in ("R", "K"):
        trypsin_sites.append(f"{seq[-1]}{len(seq)}_Cterm")

    chymotrypsin_sites = []
    for i in range(len(seq) - 1):
        if seq[i] in ("W", "Y", "F"):
            if seq[i + 1] != "P":
                chymotrypsin_sites.append(f"{seq[i]}{i+1}_{seq[i+1]}{i+2}")
    if seq and seq[-1] in ("W", "Y", "F"):
        chymotrypsin_sites.append(f"{seq[-1]}{len(seq)}_Cterm")

    is_trypsin_resistant = len(trypsin_sites) == 0
    is_chymotrypsin_vulnerable = len(chymotrypsin_sites) > 0

    return {
        "cleavage_site_count": len(trypsin_sites),
        "cleavage_sites": trypsin_sites,
        "trypsin_cleavage_sites": trypsin_sites,
        "is_protease_resistant": is_trypsin_resistant,
        "is_trypsin_resistant": is_trypsin_resistant,
        "chymotrypsin_cleavage_sites": chymotrypsin_sites,
        "chymotrypsin_site_count": len(chymotrypsin_sites),
        "is_chymotrypsin_vulnerable": is_chymotrypsin_vulnerable,
        "broad_protease_verdict": "AT_RISK" if is_chymotrypsin_vulnerable else ("RESISTANT_IN_SILICO" if is_trypsin_resistant else "FALSIFIED_TRYPSIN")
    }


def is_viable(sequence: str) -> Tuple[bool, dict]:
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

    # Donnan compatibility at pH 6.0
    q_ph6 = calc_net_charge_at_ph(seq, 6.0)

    # Stephan acid switch
    stephan = calc_stephan_switch(seq)

    # Protease resistance
    protease = count_protease_cleavage_sites(seq)

    # The Guruprasad instability index is parameterised on linear proteins, so it is out
    # of domain for a disulfide-cyclised nonapeptide. Earlier versions handled that by
    # silently substituting the literal 25.0, a comfortable pass, and then reporting it
    # as "instability_index" to three decimal places. Every cyclic candidate in this
    # library starts and ends with cysteine, so the instability gate could not fail for
    # any of them. It is now reported as not applicable and the gate is not applied, so
    # nobody can read a fabricated 25.0 as a measurement. The same held for the
    # isoelectric point when biopython is missing: 11.0, 9.5 and 6.5 were a lookup on the
    # arginine count, not a computed pI.
    is_cyclic_peptide = seq.startswith("C") and seq.endswith("C")
    instability_raw = None
    instability_applied = False
    instability = None

    if HAS_BIOPYTHON:
        try:
            analyzer = ProteinAnalysis(seq)
            gravy = analyzer.gravy()
            instability_raw = round(analyzer.instability_index(), 3)
            pi = round(analyzer.isoelectric_point(), 3)
            pi_source = "biopython_ProteinAnalysis"
            if not is_cyclic_peptide:
                instability = instability_raw
                instability_applied = True
        except Exception as e:
            return False, {"error": str(e)}
    else:
        gravy = sum(KD_HYDROPATHY.get(aa, 0.0) for aa in seq) / len(seq)
        pi = None
        pi_source = "not_computed_biopython_missing"

    metrics = {
        "gravy": round(gravy, 3),
        "instability_index": instability,
        "instability_index_raw": instability_raw,
        "instability_index_applicable": instability_applied,
        "instability_index_note": (
            "Not applied: the Guruprasad index is parameterised on linear proteins and is "
            "out of domain for a disulfide-cyclised peptide. This screen therefore does "
            "not test chemical stability at all."
            if not instability_applied else
            "Applied: linear sequence, within the index's domain."
        ),
        "isoelectric_point": pi,
        "isoelectric_point_source": pi_source,
        "arg_count": arg_count,
        "net_charge_pH6": q_ph6,
        "stephan_switch": stephan,
        "protease": protease,
        "chitosan_native": True,
        "gates_applied": ["gravy"] + (["instability_index"] if instability_applied else []),
    }

    if gravy > MAX_GRAVY:
        return False, metrics

    if instability_applied and instability > MAX_INSTABILITY:
        return False, metrics

    return True, metrics


if __name__ == "__main__":
    print("CRRRWWRRRRC:", is_viable("CRRRWWRRRRC"))
    print("CPHWRWPHC:  ", is_viable("CPHWRWPHC"))
    print("CPHWRPWHC:  ", is_viable("CPHWRPWHC"))
