#!/usr/bin/env python3
"""
Dr. Nemo — Orchestrator (Phase 8)
Monitors pending queue, runs the async pipeline, then triggers auto_pilot.
"""
import time
import subprocess
from pathlib import Path
import os
import sys

BASE_DIR = Path.home() / "drnemo_research"
VENV_PYTHON = BASE_DIR / "venv" / "bin" / "python3"
PENDING_DIR = BASE_DIR / "jobs" / "pending"
PIPELINE_SCRIPT = BASE_DIR / "pipeline.py"
ANALYZER_SCRIPT = BASE_DIR / "analyzer.py"
AUTOPILOT_SCRIPT = BASE_DIR / "auto_pilot.py"
ENV_FILE = BASE_DIR / ".env"

def create_env_template():
    if not ENV_FILE.exists():
        with open(ENV_FILE, 'w') as f:
            f.write("NGC_API_KEY=YOUR_OPENFOLD_API_KEY_HERE\n")
            f.write("# HF_TOKEN=YOUR_HUGGINGFACE_TOKEN_HERE (needed for MedGemma)\n")
        print(f"⚠️ Created missing .env file at {ENV_FILE}. Please edit it with your actual API keys!")

def main_loop():
    create_env_template()
    print("=" * 60)
    print("🤖 Dr. Nemo 24/7 Orchestrator (Phase 8 — Async Mode)")
    print("=" * 60)
    print(f"👀 Monitoring for new jobs in: {PENDING_DIR}")
    print(f"⚡ Pipeline: Async aiohttp batch dispatch")
    print(f"🧬 Auto-Pilot: {'ENABLED' if AUTOPILOT_SCRIPT.exists() else 'DISABLED'}")
    print("=" * 60)
    
    while True:
        try:
            job_files = list(PENDING_DIR.glob("*.json"))
            if job_files:
                print(f"\n[{time.strftime('%Y-%m-%d %H:%M:%S')}] 🚨 Detected {len(job_files)} pending job(s)!")
                
                # 1. Run async structural prediction pipeline
                print("⚡ Running Async Structural Pipeline...")
                subprocess.run([str(VENV_PYTHON), str(PIPELINE_SCRIPT)])
                
                # 2. Run LLM Clinical Analysis on the new results
                if ANALYZER_SCRIPT.exists():
                    print("🧠 Running Clinical Analysis (LLM Synthesizer)...")
                    subprocess.run([str(VENV_PYTHON), str(ANALYZER_SCRIPT)])
                else:
                    print("⚠️ Analyzer script not found. Skipping clinical synthesis step.")

                # 3. Auto-trigger genetic auto-pilot if available
                if AUTOPILOT_SCRIPT.exists():
                    remaining = list(PENDING_DIR.glob("*.json"))
                    if not remaining:
                        print("🧬 Triggering Autonomous Genetic Auto-Pilot...")
                        subprocess.run([str(VENV_PYTHON), str(AUTOPILOT_SCRIPT), "--single-generation"])
                    else:
                        print(f"⏳ {len(remaining)} jobs still pending, deferring auto-pilot.")

            time.sleep(10)  # Polling interval

        except KeyboardInterrupt:
            print("\n🛑 Graceful Shutdown: Dr. Nemo Orchestrator stopped by user.")
            break
        except Exception as e:
            print(f"❌ Orchestrator Error: {e}. Backing off for 30s...")
            time.sleep(30)

if __name__ == "__main__":
    main_loop()
