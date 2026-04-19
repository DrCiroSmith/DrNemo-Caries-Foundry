#!/usr/bin/env python3
import os
import csv
import time
import subprocess
import json
import sys
from pathlib import Path
import pandas as pd

BASE_DIR = Path.home() / "drnemo_research"
VENV_PYTHON = BASE_DIR / "venv" / "bin" / "python3"
PENDING_DIR = BASE_DIR / "jobs" / "pending"
RESULTS_DIR = BASE_DIR / "results"
SCORES_FILE = RESULTS_DIR / "master_scores.csv"
TARGETS_FILE = BASE_DIR / "targets.json"

STATS_FILE = BASE_DIR / "expansion_stats.json"

BREEDING_THRESHOLD = 100
SAMPLES_PER_PROJECT = int(sys.argv[1]) if len(sys.argv) > 1 else 1625

def get_project_count(t_id):
    if not SCORES_FILE.exists():
        return 0
    try:
        df = pd.read_csv(SCORES_FILE, usecols=['Job_ID'])
        return len(df[df['Job_ID'].str.startswith(t_id)])
    except:
        return 0

def get_current_count():
    if not SCORES_FILE.exists():
        return 0
    with open(SCORES_FILE, 'r') as f:
        return sum(1 for line in f) - 1

def get_avg_score():
    if not SCORES_FILE.exists():
        return 0.0
    try:
        df = pd.read_csv(SCORES_FILE)
        if 'Composite_Score' in df.columns:
            return df['Composite_Score'].mean()
        return 0.0
    except:
        return 0.0

def run_command(cmd_list):
    print(f"🚀 Running: {' '.join(cmd_list)}")
    subprocess.run(cmd_list)

def update_stats(curr, target, start_count):
    stats = {
        "current_samples": curr,
        "target_samples": target,
        "start_samples": start_count,
        "avg_score": get_avg_score(),
        "last_update": time.ctime(),
        "status": "RUNNING"
    }
    with open(STATS_FILE, 'w') as f:
        json.dump(stats, f, indent=4)


def main():
    initial_count = get_current_count()
    
    with open(TARGETS_FILE, 'r') as f:
        targets_list = list(json.load(f).keys())
    
    target_count = SAMPLES_PER_PROJECT * len(targets_list)

    if initial_count == 0:
        print("🌱 Seeding initial population...")
        run_command([str(VENV_PYTHON), str(BASE_DIR / "mutation_scanner.py")])
        initial_count = get_current_count()
    
    print(f"🧬 Starting Autonomous Lab Expansion...")
    print(f"📊 Current Samples: {initial_count}")
    print(f"🎯 Target Samples: {target_count}")
    
    update_stats(initial_count, target_count, initial_count)
    
    while True:
        curr = get_current_count()
        update_stats(curr, target_count, initial_count)
        
        if all(get_project_count(t) >= SAMPLES_PER_PROJECT for t in targets_list):
            print("🎉 All projects reached the 1500 sample target. Expansion Complete!")
            break
            
        print(f"\n--- Progress: {curr}/{target_count} ({ (curr/target_count)*100 :.1f}%) ---")
        
        # 1. Check if pending is empty
        pending_count = len(list(PENDING_DIR.glob("*.json")))
        if pending_count == 0:
            print("🧪 Evolution Engine: Generating next-gen variants...")
            generation_triggered = False
            for t_id in targets_list:
                p_count = get_project_count(t_id)
                
                if p_count >= SAMPLES_PER_PROJECT:
                    print(f"✅ Project {t_id} reached target ({p_count}/{SAMPLES_PER_PROJECT}). Skipping.")
                    continue
                
                generation_triggered = True
                if p_count >= BREEDING_THRESHOLD:
                    print(f"🧬 Threshold met ({p_count} >= {BREEDING_THRESHOLD}). Running Automated Breeder for {t_id}")
                    run_command([str(VENV_PYTHON), str(BASE_DIR / "automated_breeder.py"), t_id])
                else:
                    print(f"🔬 Exploring ({p_count} < {BREEDING_THRESHOLD}). Running Directed Evolution for {t_id}")
                    run_command([str(VENV_PYTHON), str(BASE_DIR / "directed_evolution.py"), t_id])
            
            if not generation_triggered:
                print("🏁 No more variants to generate. Waiting for pipeline to finish...")
        
        # 2. Run Structural Pipeline
        print("🧊 Folding Engine: Predicting 3D structures...")
        run_command([str(VENV_PYTHON), str(BASE_DIR / "pipeline.py")])
        
        # 3. Run Clinical Analysis
        print("🧠 Diagnostician: Writing clinical reports...")
        run_command([str(VENV_PYTHON), str(BASE_DIR / "analyzer.py")])
        
        # Short sleep to prevent CPU thrashing
        time.sleep(5)

if __name__ == "__main__":
    main()
