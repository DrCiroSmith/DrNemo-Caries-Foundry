#!/usr/bin/env python3
"""
Dr. Nemo — High-Octane Async Structural Pipeline (Phase 8)
Uses aiohttp to batch 10 concurrent API requests per API key.
Each API key gets its own thread; inside each thread, asyncio dispatches batches.
"""
import asyncio
import aiohttp
import os
import json
import shutil
import csv
import threading
import time
from pathlib import Path
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor
from dotenv import load_dotenv

BASE_DIR = Path.home() / "drnemo_research"
PENDING_DIR = BASE_DIR / "jobs" / "pending"
COMPLETED_DIR = BASE_DIR / "jobs" / "completed"
FAILED_DIR = BASE_DIR / "jobs" / "failed"
RESULTS_DIR = BASE_DIR / "results"
SCORES_FILE = RESULTS_DIR / "master_scores.csv"
SCORES_LOCK = threading.Lock()

API_URL = "https://health.api.nvidia.com/v1/biology/openfold/openfold3/predict"
BATCH_SIZE = 1          # concurrent requests per key per wave
RATE_LIMIT_SLEEP = 1.5  # seconds between waves (Math: 60/1.5 = 40 RPM Absolute Safety Cap)
MAX_BACKOFF = 60         # max exponential backoff seconds

API_KEYS = []

def init_env():
    global API_KEYS
    load_dotenv(BASE_DIR / ".env")
    for d in [PENDING_DIR, COMPLETED_DIR, FAILED_DIR, RESULTS_DIR]:
        d.mkdir(parents=True, exist_ok=True)
    if not SCORES_FILE.exists():
        with open(SCORES_FILE, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(["Timestamp", "Job_ID", "Project", "Peptide_Sequence",
                             "pLDDT", "iPTM", "Composite_Score", "Confidence",
                             "Status", "PDB_File", "Error_Msg"])
                             
    RAW_KEYS = os.getenv("NGC_API_KEYS", os.getenv("NGC_API_KEY", ""))
    API_KEYS = [k.strip() for k in RAW_KEYS.split(",") if k.strip()]

def log_result(job_id, status, project="", peptide="", plddt="", iptm="",
               composite="", conf="", pdb_file="", error=""):
    with SCORES_LOCK:
        with open(SCORES_FILE, 'a', newline='') as f:
            writer = csv.writer(f)
            writer.writerow([datetime.now().isoformat(), job_id, project, peptide,
                             plddt, iptm, composite, conf, status, pdb_file, error])


def build_payload(job_data, request_id):
    """Build the API payload from a job JSON."""
    molecules = job_data.get("molecules", [])
    if not molecules:
        p_seq = job_data.get("protein_sequence")
        if p_seq:
            molecules.append({
                "type": "protein",
                "sequence": p_seq,
                "msa": {"dummy_db": {"a3m": {"alignment": f">seq\n{p_seq}", "format": "a3m"}}}
            })
        pept_seq = job_data.get("peptide_sequence")
        if pept_seq:
            molecules.append({
                "type": "protein",
                "sequence": pept_seq,
                "msa": {"dummy_db": {"a3m": {"alignment": f">seq\n{pept_seq}", "format": "a3m"}}}
            })
    if not molecules:
        raise ValueError("No protein or peptide sequences found in job file.")
    return {
        "inputs": [{
            "input_id": request_id,
            "molecules": molecules,
            "output_format": "pdb"
        }]
    }


async def process_job_async(session, job_file, api_key, key_index):
    """Process a single job asynchronously with exponential back-off."""
    try:
        with open(job_file, 'r') as f:
            job_data = json.load(f)
    except Exception as e:
        log_result(job_file.stem, "FAILED_EXCEPTION", error=str(e))
        return job_file, False

    request_id = job_data.get("request_id", job_data.get("job_id", job_file.stem))
    project_id = job_data.get("project", "")
    peptide_seq = job_data.get("peptide_sequence", "")

    try:
        payload = build_payload(job_data, request_id)
    except ValueError as e:
        log_result(request_id, "FAILED_EXCEPTION", project=project_id,
                   peptide=peptide_seq, error=str(e))
        return job_file, False

    # Compliance logging
    try:
        from compliance_logger import log_openfold_execution
        log_openfold_execution(request_id, payload["inputs"][0])
    except ImportError:
        pass

    headers = {"Authorization": f"Bearer {api_key}", "Accept": "application/json"}
    key_num = key_index + 1
    max_retries = 20  # increased to allow proper exponential backoff over all keys

    for attempt in range(max_retries):
        try:
            print(f"🚀 [{request_id}] → NVIDIA Cloud (Key #{key_index + 1})")
            async with session.post(API_URL, json=payload, headers=headers, timeout=aiohttp.ClientTimeout(total=120)) as resp:
                if resp.status == 429:
                    b_txt = await resp.text()
                    backoff = min(2 ** (attempt + 1), MAX_BACKOFF)
                    print(f"⚠️ Key #{key_index + 1} rate limited (429). Body: {b_txt[:100]} | Backoff {backoff}s...")
                    await asyncio.sleep(backoff)
                    continue

                if resp.status != 200:
                    body = await resp.text()
                    error_msg = f"API Error {resp.status}: {body[:200]}"
                    print(f"❌ {error_msg}")
                    log_result(request_id, "FAILED_API", project=project_id,
                               peptide=peptide_seq, error=error_msg)
                    return job_file, False

                result = await resp.json()
                struct_data = result['outputs'][0]['structures_with_scores'][0]
                pdb_content = struct_data['structure']
                plddt = struct_data.get('complex_plddt_score', 0)
                iptm = struct_data.get('iptm_score', 0)
                conf = struct_data.get('confidence_score', 0)
                composite_score = (plddt / 100.0) * iptm

                pdb_filename = RESULTS_DIR / f"{request_id}.pdb"
                with open(pdb_filename, "w", encoding="utf-8") as f:
                    f.write(pdb_content)

                print(f"✅ SUCCESS: {pdb_filename.name} (pLDDT: {plddt:.2f}, iPTM: {iptm:.2f}, Composite: {composite_score:.4f})")
                log_result(request_id, "SUCCESS", project=project_id, peptide=peptide_seq,
                           plddt=plddt, iptm=iptm, composite=composite_score,
                           conf=conf, pdb_file=pdb_filename.name)
                return job_file, True

        except asyncio.TimeoutError:
            print(f"⏰ Timeout for [{request_id}]. Retrying...")
            continue
        except Exception as e:
            print(f"❌ Error processing {job_file.name}: {e}")
            log_result(request_id, "FAILED_EXCEPTION", project=project_id,
                       peptide=peptide_seq, error=str(e))
            return job_file, False

    print(f"🚨 All keys rate-limited for [{request_id}].")
    log_result(request_id, "FAILED_RATELIMIT", project=project_id,
               peptide=peptide_seq, error="All API keys rate-limited")
    return job_file, False


async def async_worker_loop(key_index, job_list):
    """Process jobs in batches of BATCH_SIZE concurrently using aiohttp."""
    key = API_KEYS[key_index]
    key_num = key_index + 1
    total = len(job_list)
    print(f"🔧 Worker #{key_num}: Processing {total} jobs in batches of {BATCH_SIZE}")

    async with aiohttp.ClientSession() as session:
        for batch_start in range(0, total, BATCH_SIZE):
            batch = job_list[batch_start:batch_start + BATCH_SIZE]
            batch_num = (batch_start // BATCH_SIZE) + 1
            print(f"   Worker #{key_num} | Batch {batch_num} | {len(batch)} jobs dispatching...")

            tasks = [process_job_async(session, jf, key, key_index) for jf in batch]
            results = await asyncio.gather(*tasks, return_exceptions=True)

            for res in results:
                if isinstance(res, Exception):
                    print(f"❌ Unexpected batch error: {res}")
                    continue
                job_file, success = res
                if success:
                    shutil.move(str(job_file), str(COMPLETED_DIR / job_file.name))
                else:
                    shutil.move(str(job_file), str(FAILED_DIR / job_file.name))

            # Rate-limit sleep between batches
            if batch_start + BATCH_SIZE < total:
                await asyncio.sleep(RATE_LIMIT_SLEEP)

    print(f"✅ Worker #{key_num}: All jobs complete.")


def thread_worker(key_index, job_list):
    """Thread entry point — runs the async loop for one API key."""
    asyncio.run(async_worker_loop(key_index, job_list))


def run_pipeline():
    init_env()
    job_files = sorted(list(PENDING_DIR.glob("*.json")))

    if not job_files:
        print("📭 No pending jobs found.")
        return

    num_keys = len(API_KEYS)
    if num_keys == 0:
        print("❌ No API keys available.")
        return

    print(f"\n{'='*60}")
    print(f"⚡ DR. NEMO — HIGH-OCTANE ASYNC PIPELINE (Phase 8)")
    print(f"{'='*60}")
    print(f"📥 {len(job_files)} jobs in queue")
    print(f"🔑 {num_keys} API keys × {BATCH_SIZE} concurrent = {num_keys * BATCH_SIZE} parallel requests")
    print(f"{'='*60}\n")

    # Round-robin distribute jobs to workers
    worker_jobs = [[] for _ in range(num_keys)]
    for i, jf in enumerate(job_files):
        worker_jobs[i % num_keys].append(jf)

    for i in range(num_keys):
        print(f"   Worker #{i+1}: {len(worker_jobs[i])} jobs assigned")

    start_time = time.time()

    with ThreadPoolExecutor(max_workers=num_keys) as executor:
        futures = []
        for i in range(num_keys):
            futures.append(executor.submit(thread_worker, i, worker_jobs[i]))
        for f in futures:
            f.result()

    elapsed = time.time() - start_time
    throughput = len(job_files) / (elapsed / 3600) if elapsed > 0 else 0
    print(f"\n{'='*60}")
    print(f"✅ Pipeline batch complete in {elapsed:.1f}s ({throughput:.0f} jobs/hr)")
    print(f"{'='*60}")


if __name__ == "__main__":
    run_pipeline()
