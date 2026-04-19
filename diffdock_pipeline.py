#!/usr/bin/env python3
import json
import random
import os
import sys
from pathlib import Path

BASE_DIR = Path.home() / "drnemo_research"
PENDING_DIR = BASE_DIR / "jobs" / "pending"
TARGETS_FILE = BASE_DIR / "targets.json"

target_id = sys.argv[1] if len(sys.argv) > 1 else "P_gingivalis_Protease"

with open(TARGETS_FILE, 'r') as f:
    targets_data = json.load(f)

if target_id not in targets_data:
    print(f"❌ Error: Target '{target_id}' not found in targets.json.")
    sys.exit(1)

target_info = targets_data[target_id]
TARGET_PROTEIN = target_info["sequence"]

# Prototype chemical compounds (SMILES strings)
CHEMICAL_LIGANDS = [
    "CC(=O)OC1=CC=CC=C1C(=O)O", # Aspirin (Prototype Baseline)
    "CC1(C(N2C(S1)C(C2=O)NC(=O)C(C3=CC=CC=C3)N)C(=O)O)C", # Ampicillin
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", # Caffeine
    "OC1=C(C=C(C=C1)O)C2=CC(=C(C=C2)O)O", # Catechol (Generic phenol backbone)
    "C1=CC=C2C(=C1)NC=N2" # Benzimidazole (Common fragment)
]

num_ligands = len(CHEMICAL_LIGANDS)

print("\n" + "="*60)
print(f"👨‍🏫 TEACHER MODE: DIFFDOCK (CHEMICAL DRUG) PIPELINE")
print("="*60)
print(f"Project Target    : {target_id}")
print(f"Target Sequence   : {TARGET_PROTEIN[:15]}... (Length: {len(TARGET_PROTEIN)})")
print("-" * 60)
print("CONCEPT: The 'Small Molecule Inhibitor' strategy.")
print("Unlike OpenFold (which tests Amino-Acid Peptides), DiffDock")
print("tests actual chemical drugs (SMILES strings) against the target.")
print("We calculate the 'binding pose' of these synthetic chemicals in")
print("the pathogen's active lock, measuring Vina binding affinity.")
print("="*60 + "\n")

print(f"🧪 Injecting {num_ligands} chemical SMILES molecules into the docking queue...")

for i, smiles in enumerate(CHEMICAL_LIGANDS):
    job_id = f"{target_id}_docking_{i+1}"
    
    job_data = {
        "job_id": job_id,
        "protein_sequence": TARGET_PROTEIN,
        "smiles_ligand": smiles,
        "project": target_id,
        "pipeline": "diffdock",
        "notes": f"NVIDIA DiffDock Prototype Testing for {smiles}."
    }
    
    job_file = PENDING_DIR / f"{job_id}.json"
    
    with open(job_file, 'w') as f:
        json.dump(job_data, f, indent=4)
        
    print(f"✅ Created docking target: {job_id} (Ligand: {smiles})")

print("\n🚀 All synthetic targets loaded into the pending queue!")
print("👨‍🏫 Teacher Note: The orchestrator will now process these via DiffDock")
print("   to determine physical drug potential rather than peptide potential.")
