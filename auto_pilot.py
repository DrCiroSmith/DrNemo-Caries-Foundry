#!/usr/bin/env python3
"""
Dr. Nemo — Autonomous Genetic Auto-Pilot (Phase 8)
Headless recursive evolutionary loop that breeds elite peptides,
filters through solubility shield, and queues the next generation.

Usage:
    python3 auto_pilot.py                   # Run full recursive loop
    python3 auto_pilot.py --single-generation  # Run one generation only (for orchestrator)
"""
import json
import random
import os
import sys
import time
import csv
from pathlib import Path
from datetime import datetime
import pandas as pd
from solubility_shield import is_viable

BASE_DIR = Path.home() / "drnemo_research"
PENDING_DIR = BASE_DIR / "jobs" / "pending"
RESULTS_DIR = BASE_DIR / "results"
COMPLETED_DIR = BASE_DIR / "jobs" / "completed"
SCORES_FILE = RESULTS_DIR / "master_scores.csv"
TARGETS_FILE = BASE_DIR / "targets.json"
EVOLUTION_LOG = BASE_DIR / "evolution_log.csv"

VALID_AA = list("ACDEFGHIKLMNPQRSTVWY")
MAX_GENERATIONS = 50
PLATEAU_THRESHOLD = 5   # stop after N generations without improvement
ELITE_PERCENTILE = 0.05 # top 5%
MIN_ELITES = 2
CHILDREN_PER_PROJECT = 15
MUTATION_RATE = 0.10     # 10% chance per position

def init_evolution_log():
    if not EVOLUTION_LOG.exists():
        with open(EVOLUTION_LOG, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(["Timestamp", "Generation", "Project", "Num_Elites",
                             "Elite_Avg_Score", "Elite_Best_Score", "Children_Generated",
                             "Children_Passed_Filter", "Children_Queued"])

def log_generation(gen, project, num_elites, avg_score, best_score,
                   children_gen, children_passed, children_queued):
    with open(EVOLUTION_LOG, 'a', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([datetime.now().isoformat(), gen, project, num_elites,
                         f"{avg_score:.6f}", f"{best_score:.6f}",
                         children_gen, children_passed, children_queued])

def extract_peptide(job_id):
    if "variant_" in str(job_id):
        return str(job_id).split("variant_")[-1]
    return ""

def is_duplicate(peptide_seq, target_id, df_scores):
    job_id = f"{target_id}_variant_{peptide_seq}"
    pending_file = PENDING_DIR / f"{job_id}.json"
    if pending_file.exists():
        return True
    if df_scores is not None and not df_scores.empty:
        if job_id in df_scores["Job_ID"].values:
            return True
    return False

def crossover(parent1, parent2):
    """Single-point crossover between two peptide sequences."""
    min_len = min(len(parent1), len(parent2))
    if min_len < 2:
        return parent1
    point = random.randint(1, min_len - 1)
    child = parent1[:point] + parent2[point:]
    return child

def mutate(seq, rate=MUTATION_RATE):
    """Point mutation at each position with given probability."""
    seq_list = list(seq)
    for i in range(len(seq_list)):
        if random.random() < rate:
            seq_list[i] = random.choice(VALID_AA)
    return ''.join(seq_list)

def breed_generation(project_id, df_scores, targets_data, generation_num):
    """Select elites, breed children, filter, and queue for one project."""
    # Filter to project
    project_mask = df_scores["Job_ID"].apply(lambda jid: str(jid).startswith(project_id))
    project_df = df_scores[project_mask & (df_scores["Status"] == "SUCCESS")].copy()
    
    if project_df.empty:
        return 0, 0, 0, 0, 0.0, 0.0

    if "Composite_Score" not in project_df.columns:
        return 0, 0, 0, 0, 0.0, 0.0

    project_df["Composite_Score"] = pd.to_numeric(project_df["Composite_Score"], errors="coerce")
    project_df = project_df.dropna(subset=["Composite_Score"])
    project_df = project_df.sort_values(by="Composite_Score", ascending=False)

    # Select elites: top 5% (minimum 2)
    n_elites = max(MIN_ELITES, int(len(project_df) * ELITE_PERCENTILE))
    elites = project_df.head(n_elites)
    
    elite_seqs = [extract_peptide(jid) for jid in elites["Job_ID"] if extract_peptide(jid)]
    elite_scores = elites["Composite_Score"].tolist()
    
    if len(elite_seqs) < 2:
        return len(elite_seqs), 0, 0, 0, 0.0, 0.0

    avg_score = sum(elite_scores[:len(elite_seqs)]) / len(elite_seqs)
    best_score = elite_scores[0] if elite_scores else 0.0

    protein_sequence = targets_data[project_id].get("sequence", "")
    
    # Breed children via crossover + mutation
    children = set()
    attempts = 0
    while len(children) < CHILDREN_PER_PROJECT * 3 and attempts < 300:
        attempts += 1
        p1, p2 = random.sample(elite_seqs, min(2, len(elite_seqs)))
        child = crossover(p1, p2)
        child = mutate(child)
        if child and child not in elite_seqs:
            children.add(child)

    # Filter through solubility shield
    viable_children = []
    for child_seq in children:
        ok, _ = is_viable(child_seq)
        if ok:
            viable_children.append(child_seq)

    # Queue non-duplicate viable children
    queued = 0
    PENDING_DIR.mkdir(parents=True, exist_ok=True)
    for child_seq in viable_children[:CHILDREN_PER_PROJECT]:
        if is_duplicate(child_seq, project_id, df_scores):
            continue
        job_id = f"{project_id}_variant_{child_seq}"
        job_data = {
            "job_id": job_id,
            "protein_sequence": protein_sequence,
            "peptide_sequence": child_seq,
            "project": project_id,
            "notes": f"Auto-Pilot Gen {generation_num} — Elite Breeding"
        }
        out_path = PENDING_DIR / f"{job_id}.json"
        with open(out_path, 'w') as f:
            json.dump(job_data, f, indent=4)
        queued += 1

    return len(elite_seqs), len(children), len(viable_children), queued, avg_score, best_score


def wait_for_pending_jobs(timeout_hours=6):
    """Block until all pending jobs are processed (or timeout)."""
    start = time.time()
    max_wait = timeout_hours * 3600
    while time.time() - start < max_wait:
        pending = list(PENDING_DIR.glob("*.json"))
        if not pending:
            return True
        print(f"   ⏳ Waiting... {len(pending)} jobs still pending ({(time.time()-start)/60:.1f}m elapsed)")
        time.sleep(30)
    print("⏰ Timeout waiting for pending jobs.")
    return False


def run_auto_pilot(single_generation=False):
    """Main auto-pilot loop."""
    init_evolution_log()
    
    with open(TARGETS_FILE, 'r') as f:
        targets_data = json.load(f)
    
    project_ids = list(targets_data.keys())
    
    print("=" * 60)
    print("🧬 DR. NEMO — AUTONOMOUS GENETIC AUTO-PILOT (Phase 8)")
    print("=" * 60)
    print(f"📋 Projects: {len(project_ids)}")
    print(f"🧪 Elite Percentile: Top {ELITE_PERCENTILE*100:.0f}%")
    print(f"🔄 Max Generations: {MAX_GENERATIONS}")
    print(f"📉 Plateau Threshold: {PLATEAU_THRESHOLD} generations")
    print(f"👶 Children per Project: {CHILDREN_PER_PROJECT}")
    print("=" * 60)
    
    best_scores_history = {}  # project -> list of best scores per generation
    
    max_gen = 1 if single_generation else MAX_GENERATIONS
    
    for gen in range(1, max_gen + 1):
        print(f"\n{'─'*60}")
        print(f"  🧬 GENERATION {gen}")
        print(f"{'─'*60}")
        
        if not SCORES_FILE.exists():
            print("No master_scores.csv found. Nothing to evolve from.")
            break
            
        df = pd.read_csv(SCORES_FILE)
        total_queued = 0
        plateau_all = True
        
        for proj in project_ids:
            n_elites, n_children, n_passed, n_queued, avg_sc, best_sc = \
                breed_generation(proj, df, targets_data, gen)
            
            if n_elites < 2:
                continue
                
            log_generation(gen, proj, n_elites, avg_sc, best_sc,
                          n_children, n_passed, n_queued)
            
            total_queued += n_queued
            
            print(f"  {proj}: {n_elites} elites (best={best_sc:.4f}, avg={avg_sc:.4f}) "
                  f"→ {n_children} bred → {n_passed} viable → {n_queued} queued")
            
            # Plateau detection
            if proj not in best_scores_history:
                best_scores_history[proj] = []
            best_scores_history[proj].append(best_sc)
            
            if len(best_scores_history[proj]) >= PLATEAU_THRESHOLD:
                recent = best_scores_history[proj][-PLATEAU_THRESHOLD:]
                if max(recent) - min(recent) < 0.001:
                    print(f"  📉 {proj}: Plateau detected ({PLATEAU_THRESHOLD} gens, Δ < 0.001)")
                else:
                    plateau_all = False
            else:
                plateau_all = False
        
        if total_queued == 0:
            print("\n🏁 No new children queued across any project. Evolution complete.")
            break
            
        print(f"\n  📊 Generation {gen} Summary: {total_queued} total jobs queued")
        
        if single_generation:
            print("  ℹ️  Single-generation mode. Exiting.")
            break
        
        if plateau_all and gen >= PLATEAU_THRESHOLD:
            print(f"\n🏁 All projects have plateaued for {PLATEAU_THRESHOLD} generations. Stopping.")
            break
        
        # Wait for the pipeline to process all pending jobs
        print(f"\n  ⏳ Waiting for generation {gen} to finish folding...")
        if not wait_for_pending_jobs():
            print("  ⚠️  Timeout. Proceeding to next generation anyway.")
    
    print(f"\n{'='*60}")
    print("🏁 Auto-Pilot complete.")
    print(f"{'='*60}")


if __name__ == "__main__":
    single = "--single-generation" in sys.argv
    run_auto_pilot(single_generation=single)
