#!/usr/bin/env python3
import json
import random
import os
import sys
from pathlib import Path
import pandas as pd
from collections import Counter
from solubility_shield import is_viable

BASE_DIR = Path.home() / "drnemo_research"
PENDING_DIR = BASE_DIR / "jobs" / "pending"
RESULTS_DIR = BASE_DIR / "results"
SCORES_FILE = RESULTS_DIR / "master_scores.csv"
TARGETS_FILE = BASE_DIR / "targets.json"

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")

def is_duplicate(peptide_seq, target_id, df_scores):
    job_id = f"{target_id}_variant_{peptide_seq}"
    pending_file = PENDING_DIR / f"{job_id}.json"
    if pending_file.exists():
        return True
    if df_scores is not None and not df_scores.empty:
        if job_id in df_scores["Job_ID"].values:
            return True
    return False

def extract_peptide(job_id):
    if '_variant_' in job_id:
        return job_id.split('_variant_')[-1]
    return ""

def generate_hybrids(target_id):
    if not SCORES_FILE.exists():
        print("No master_scores.csv found!")
        sys.exit(1)
        
    df = pd.read_csv(SCORES_FILE)
    
    # Filter for project
    # Handle older naming conventions if needed (like "S_mutans_Baseline")
    # For now, we look for job_ids starting with the target_id
    project_mask = df["Job_ID"].apply(lambda jid: str(jid).startswith(target_id))
    project_df = df[project_mask].copy()
    
    if project_df.empty:
        print(f"No results found for {target_id}")
        sys.exit(0)
        
    # Sort by Composite Score
    if 'Composite_Score' in project_df.columns:
        project_df = project_df.sort_values(by="Composite_Score", ascending=False)
    else:
        project_df = project_df.sort_values(by="pLDDT", ascending=False)
        
    # Take top 5
    top_5 = project_df.head(5)
    seqs = [extract_peptide(jid) for jid in top_5["Job_ID"] if extract_peptide(jid)]
    composites = top_5["Composite_Score"].tolist() if "Composite_Score" in top_5.columns else [0] * len(seqs)
    
    if len(seqs) < 2:
        print(f"Not enough valid sequences for {target_id} to breed.")
        sys.exit(0)
        
    max_len = max(len(s) for s in seqs)
    padded = [s.ljust(max_len, '-') for s in seqs]
    
    consensus = []
    for pos_i in range(max_len):
        residues = [s[pos_i] for s in padded if s[pos_i] != '-']
        if residues:
            consensus.append(Counter(residues).most_common(1)[0][0])
        else:
            consensus.append('-')
            
    hybrids = set()
    
    # 1. Positional Best
    positional_best = []
    for pos_i in range(max_len):
        best_idx = 0
        best_score = -1
        for s_idx, seq in enumerate(padded):
            score = float(composites[s_idx]) if s_idx < len(composites) else 0
            if seq[pos_i] != '-' and score > best_score:
                best_score = score
                best_idx = s_idx
        positional_best.append(padded[best_idx][pos_i])
    hybrids.add(''.join(positional_best).replace('-', ''))
    
    # 2. Consensus motif
    hybrids.add(''.join(consensus).replace('-', ''))
    
    # 3-15. Random crossovers
    attempts = 0
    while len(hybrids) < 15 and attempts < 100:
        attempts += 1
        p1, p2 = random.sample(range(len(padded)), min(2, len(padded)))
        crossover_pt = random.randint(1, max_len - 1)
        child = padded[p1][:crossover_pt] + padded[p2][crossover_pt:]
        
        child_list = list(child)
        if random.random() < 0.2:
            mut_pos = random.randint(0, len(child_list) - 1)
            child_list[mut_pos] = random.choice(list(VALID_AA))
            
        child_seq = ''.join(child_list).replace('-', '')
        if child_seq and child_seq not in [s.replace('-', '') for s in padded]:
            hybrids.add(child_seq)
            
    # Write jobs
    with open(TARGETS_FILE, 'r') as f:
        targets_data = json.load(f)
        
    if target_id not in targets_data:
        print(f"Target {target_id} not in targets.json")
        sys.exit(1)
        
    protein_sequence = targets_data[target_id].get("sequence", "")
    PENDING_DIR.mkdir(parents=True, exist_ok=True)
    
    hybrids_list = list(hybrids)
    submitted = 0
    for h_seq in hybrids_list[:15]:
        if not is_duplicate(h_seq, target_id, project_df):
            ok, _ = is_viable(h_seq)
            if not ok:
                continue
                
            job_id = f"{target_id}_variant_{h_seq}"
            job_data = {
                "job_id": job_id,
                "protein_sequence": protein_sequence,
                "peptide_sequence": h_seq,
                "project": target_id,
                "notes": "Automated Combinatorial Breeding (Deep Exploitation)"
            }
            out_path = PENDING_DIR / f"{job_id}.json"
            with open(out_path, 'w') as f:
                json.dump(job_data, f, indent=4)
            print(f"✅ Queued High-Confidence Breeder Target: {job_id}")
            submitted += 1
            
    print(f"🧬 Automated Breeding generated {submitted} novel jobs for {target_id}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python automated_breeder.py <target_id>")
        sys.exit(1)
    target_id = sys.argv[1]
    print(f"\n🧬 Starting Automated Combinatorial Breeder for {target_id}...")
    generate_hybrids(target_id)
