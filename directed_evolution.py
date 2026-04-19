#!/usr/bin/env python3
import json
import random
import os
import csv
from pathlib import Path
import sys
from dotenv import load_dotenv
from solubility_shield import is_viable

load_dotenv()

try:
    import ollama
    HAS_OLLAMA = True
except ImportError:
    HAS_OLLAMA = False
    print("⚠️ Ollama package not found. Will use statistical fallback.")

try:
    from openai import OpenAI
    HAS_OPENAI = True
except ImportError:
    HAS_OPENAI = False


BASE_DIR = Path.home() / "drnemo_research"
PENDING_DIR = BASE_DIR / "jobs" / "pending"
RESULTS_DIR = BASE_DIR / "results"
SCORES_FILE = RESULTS_DIR / "master_scores.csv"
TARGETS_FILE = BASE_DIR / "targets.json"

target_id = sys.argv[1] if len(sys.argv) > 1 else "S_mutans_Baseline"

with open(TARGETS_FILE, 'r') as f:
    targets_data = json.load(f)

if target_id not in targets_data:
    print(f"❌ Error: Target '{target_id}' not found in targets.json.")
    sys.exit(1)

target_info = targets_data[target_id]
TARGET_PROTEIN = target_info["sequence"]
BASE_PEPTIDE = target_info["base_peptide"]

# Amino Acids to sample from for mutations
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"

# Evolutionary Parameters
num_variants_per_parent = 5
top_parents_to_keep = 3

print("\n" + "="*60)
print(f"👨‍🏫 TEACHER MODE: DIRECTED EVOLUTION ALGORITHM")
print("="*60)
print(f"Project Target    : {target_id}")
print(f"Target Sequence   : {TARGET_PROTEIN[:15]}... (Length: {len(TARGET_PROTEIN)})")
print("-" * 60)
print("CONCEPT: The 'Evolutionary LLM Loop'.")
print("This algorithm feeds top historical performers to a local Ollama LLM.")
print("The LLM acts as a Biological Designer, deducing binding motifs and")
print("logically generating the next generation of peptide blockers.")
print("If the LLM fails, it falls back to a 1D statistical probability map.")
print("="*60 + "\n")

print(f"🧬 Booting Directed Evolution Engine... [Project: {target_id}]")

# 1. Read existing scores to find the "Fittest" candidates
try:
    with open(SCORES_FILE, 'r') as f:
        reader = csv.DictReader(f)
        results = []
        for row in reader:
            if row['Status'] == 'SUCCESS' and (row.get('Project', '') == target_id or target_id in row['Job_ID']):
                plddt = float(row.get('pLDDT', 0))
                iptm = float(row.get('iPTM', 0))
                try:
                    composite = float(row.get('Composite_Score'))
                except (TypeError, ValueError):
                    composite = (plddt / 100.0) * iptm
                
                row['pLDDT'] = plddt
                row['iPTM'] = iptm
                row['Composite_Score'] = composite
                results.append(row)
except FileNotFoundError:
    print("❌ No master_scores.csv found. Run mutation_scanner.py first to seed the population.")
    sys.exit(1)

if not results:
    print(f"❌ No successful computational folds found for project {target_id}. Seed population first.")
    sys.exit(1)

# 2. Sort by Composite Score descending
results.sort(key=lambda x: x['Composite_Score'], reverse=True)

# 3. Extract the parent sequences
parent_sequences = []
seen_sequences = set()
# Deduplication Shield: Scan pending jobs to avoid duplicates across workers
for _f in (PENDING_DIR.iterdir() if PENDING_DIR.exists() else []):
    if _f.suffix == ".json" and "variant_" in _f.name:
        _p_seq = _f.name.split("variant_")[-1].replace(".json", "")
        seen_sequences.add(_p_seq.upper())

# Deduplicate to find unique top parents
for row in results:
    job_id = row['Job_ID']
    seq = job_id.split('variant_')[-1] if 'variant_' in job_id else ""
    if not seq or seq in seen_sequences: continue
        
    parent_sequences.append({
        'sequence': seq,
        'Composite_Score': row['Composite_Score'],
        'iPTM': row['iPTM'],
        'pLDDT': row['pLDDT']
    })
    seen_sequences.add(seq)
        
    if len(parent_sequences) >= top_parents_to_keep:
        break

if not parent_sequences:
    print(f"⚠️ No valid variant sequences found in history. Seeding from default {BASE_PEPTIDE}.")
    parent_sequences = [{'sequence': BASE_PEPTIDE, 'iPTM': 0.0, 'pLDDT': 0.0, 'Composite_Score': 0.0}]

print(f"🏆 Found top {len(parent_sequences)} fittest parents using Composite Score:")
for i, p in enumerate(parent_sequences):
    label = "Neutral"
    if p['Composite_Score'] > 0.72: label = "BLUE (Excellent)"
    elif p['Composite_Score'] > 0.42: label = "GREEN (Good)"
    print(f"   {i+1}. {p['sequence']} | Score: {p['Composite_Score']:.3f} | {label}")

print(f"\n⚡ Generating next-generation mutations...")
total_created = 0
generated_targets = []

# --- CLOUD OR LOCAL IDEATION ENGINE ---
used_llm = False
github_token = os.environ.get("GITHUB_TOKEN")

if HAS_OPENAI and github_token:
    print("\n🧠 Invoking GitHub Models Cloud Ideation Engine...")
    model_id = os.environ.get("OLLAMA_MODEL", "gpt-4o")
    context_str = "\n".join([f"- Sequence: {p['sequence']} (Score: {p.get('Composite_Score', 0):.3f})" for p in parent_sequences])
    top_score = parent_sequences[0].get('Composite_Score', 0)
    
    # Narrow accuracy focus: check highest bar FIRST so it doesn't get masked by lower condition
    accuracy_instruction = "Perform high-precision exploration of binding motifs."
    if top_score > 0.72:
        accuracy_instruction = "Excellent affinity reached. Perform extremely conservative, high-precision single mutations to stabilize the fold and perfect the interaction."
    elif top_score > 0.5:
        accuracy_instruction = "We have reached high-accuracy candidates. Focus on identifying single-point mutations to refine the binding pocket further for maximum affinity."

    prompt = f"""You are a computational biologist evolving high-affinity peptide blockers.
Target Protein Sequence: {TARGET_PROTEIN}
Target Peptide Length: {len(BASE_PEPTIDE)}

Based on the physics simulation, here are the most successful peptide sequences so far:
{context_str}

Analyze the biochemical patterns, amino acid preferences, and potential binding motifs that make these successful.
STRATEGY: {accuracy_instruction}

Generate exactly 15 NEW unique peptide sequences (strictly length {len(BASE_PEPTIDE)}) that improve upon these patterns.
Return the result Strictly as a JSON object with a single key "evolved_sequences" containing a list of strings.
Example format:
{{
  "evolved_sequences": ["SEQ1", "SEQ2", "SEQ3"]
}}
"""
    messages = [{"role": "user", "content": prompt}]

    try:
        client = OpenAI(
            base_url="https://models.inference.ai.azure.com",
            api_key=github_token
        )
        response = client.chat.completions.create(timeout=15, 
            model=model_id,
            messages=messages,
            response_format={"type": "json_object"}
        )
        content = response.choices[0].message.content
        data = json.loads(content)
        llm_sequences = data.get("evolved_sequences", [])
        
        # Validate lengths and characters
        for seq in llm_sequences:
            seq = seq.upper().strip()
            if len(seq) == len(BASE_PEPTIDE) and all(c in AMINO_ACIDS for c in seq) and seq not in seen_sequences:
                ok, _ = is_viable(seq)
                if ok:
                    generated_targets.append((seq, f"Cloud Evolved ({model_id})"))
                seen_sequences.add(seq)
                
        if generated_targets:
            print(f"🧬 Cloud Model successfully designed {len(generated_targets)} viable candidates.")
            used_llm = True
        else:
            print(f"⚠️ Cloud Model generated sequences were invalid.")
            print(f"DEBUG LLM OUTPUT: {content}")
            
    except Exception as e:
        print(f"⚠️ Cloud Ideation failed ({e}).")
        try: print(f"DEBUG LLM OUTPUT: {content}")
        except: pass

elif HAS_OLLAMA:
    print("\n🧠 Invoking Ollama LLM Ideation Engine...")
    model_id = os.environ.get("OLLAMA_MODEL", "llama3.1")
    context_str = "\n".join([f"- Sequence: {p['sequence']} (Score: {p.get('Composite_Score', 0):.3f})" for p in parent_sequences])
    top_score = parent_sequences[0].get('Composite_Score', 0)
    
    # Narrow accuracy focus: check highest bar FIRST so it doesn't get masked by lower condition
    accuracy_instruction = "Perform high-precision exploration of binding motifs."
    if top_score > 0.72:
        accuracy_instruction = "Excellent affinity reached. Perform extremely conservative, high-precision single mutations to stabilize the fold and perfect the interaction."
    elif top_score > 0.5:
        accuracy_instruction = "We have reached high-accuracy candidates. Focus on identifying single-point mutations to refine the binding pocket further for maximum affinity."

    prompt = f"""You are a computational biologist evolving high-affinity peptide blockers.
Target Protein Sequence: {TARGET_PROTEIN}
Target Peptide Length: {len(BASE_PEPTIDE)}

Based on the physics simulation, here are the most successful peptide sequences so far:
{context_str}

Analyze the biochemical patterns, amino acid preferences, and potential binding motifs that make these successful.
STRATEGY: {accuracy_instruction}

Generate exactly 15 NEW unique peptide sequences (strictly length {len(BASE_PEPTIDE)}) that improve upon these patterns.
Return the result Strictly as a JSON object with a single key "evolved_sequences" containing a list of strings.
Example format:
{{
  "evolved_sequences": ["SEQ1", "SEQ2", "SEQ3"]
}}
"""
    messages = [{"role": "user", "content": prompt}]

    
    try:
        response = ollama.chat(model=model_id, messages=messages, format='json')
        content = response['message']['content']
        data = json.loads(content)
        llm_sequences = data.get("evolved_sequences", [])
        
        # Validate lengths and characters
        for seq in llm_sequences:
            seq = seq.upper().strip()
            if len(seq) == len(BASE_PEPTIDE) and all(c in AMINO_ACIDS for c in seq) and seq not in seen_sequences:
                ok, _ = is_viable(seq)
                if ok:
                    generated_targets.append((seq, "Ollama LLM Evolved"))
                seen_sequences.add(seq)
                
        if generated_targets:
            print(f"🧬 Ollama successfully designed {len(generated_targets)} viable candidates.")
            used_llm = True
        else:
            print(f"⚠️ Ollama generated sequences were invalid (wrong length or characters).")
            print(f"DEBUG LLM OUTPUT: {content}")
            
    except Exception as e:
        print(f"⚠️ Ollama Ideation failed ({e}).")
        # Try to print content if it was set
        try:
            print(f"DEBUG LLM OUTPUT: {content}")
        except: pass

# --- STATISTICAL FALLBACK ---
if not used_llm:
    print("\n🔢 Using Statistical Probability Map (Fallback)...")
    
    recognition_map = {}
    top_performers = results[:20]
    for row in top_performers:
        job_id = row['Job_ID']
        seq = job_id.split('variant_')[-1] if 'variant_' in job_id else ""
        if not seq: continue
        
        weight = row['Composite_Score']
        for i, aa in enumerate(seq):
            if i not in recognition_map: recognition_map[i] = {}
            recognition_map[i][aa] = recognition_map[i].get(aa, 0) + weight

    top_score = parent_sequences[0].get('Composite_Score', 0)
    for parent in parent_sequences:
        base_pep = parent['sequence']
        for _ in range(num_variants_per_parent):
            # Narrowing focus: smaller mutation rate for high-accuracy samples
            num_mutations = 1 if top_score > 0.5 else random.choice([1, 2])
            pep_list = list(base_pep)
            
            for _ in range(num_mutations):
                idx = random.randint(0, len(pep_list) - 1)
                if idx in recognition_map and random.random() < 0.5:
                    weights = recognition_map[idx]
                    choices = list(weights.keys())
                    choice_weights = list(weights.values())
                    pep_list[idx] = random.choices(choices, weights=choice_weights, k=1)[0]
                else:
                    pep_list[idx] = random.choice(AMINO_ACIDS)
                
            mutated_peptide = "".join(pep_list)
            if mutated_peptide not in seen_sequences:
                ok, _ = is_viable(mutated_peptide)
                if ok:
                    generated_targets.append((mutated_peptide, f"Statistical Evolved. Parent: {base_pep}"))
                seen_sequences.add(mutated_peptide)

# --- JOB CREATION ---
print("")
for seq, note in generated_targets:
    job_id = f"{target_id}_variant_{seq}"
    job_data = {
        "job_id": job_id,
        "protein_sequence": TARGET_PROTEIN,
        "peptide_sequence": seq,
        "project": target_id,
        "notes": note
    }
    
    job_file = PENDING_DIR / f"{job_id}.json"
    with open(job_file, 'w') as f:
        json.dump(job_data, f, indent=4)
        
    print(f"✅ Created computational target: {job_id}")
    total_created += 1

print(f"\n🚀 {total_created} evolved generation targets loaded into the pending queue!")
print("The 24/7 Orchestrator will now mathematically test the fitness of this new generation.")
