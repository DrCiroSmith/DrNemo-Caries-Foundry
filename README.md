# 🦷 Dr. Nemo: Caries Vaccine Foundry 🧬

**Industrialized Autonomous Peptide Discovery for Oral Health**

## 🚀 Overview

The **Dr. Nemo Caries Vaccine Foundry** is a high-precision computational biology pipeline designed to discover, validate, and synthesize next-generation vaccines for dental caries. The foundry focuses on **Arginine-rich cyclic peptides** that act as "molecular wedges," specifically inhibiting the **GtfB/C enzymes** of *Streptococcus mutans* to disrupt biofilm formation.

This project bridges the gap between in-silico structural prediction and physical manufacturing, utilizing a mucoadhesive **Chitosan/PVA film** delivery system.

---

## 🏗️ The 4-Tier Discovery Workflow

### 🥇 Tier 3: Structural Validation (The Brain)

- **Engine:** NVIDIA NIM OpenFold3 API.
- **Constraints:** Optimized T=25 diffusion steps for cyclic loop structural fidelity.
- **Goal:** Candidates must exceed a **0.72 Composite Score** (pLDDT * iPTM) for synthesis graduation.

### 🥈 Tier 2: Generative Discovery (The Engine)

- **Logic:** Arginine-enriched mutation scanning.
- **Safety:** Hardened **Chitosan-Safe Filtering** (Exclusion of I, K, D, E residues) to prevent mucosal entrapment and salt-bridge formation.
- **Library:** Seeded with the 19-mer GtfB inhibitory scaffold.

### 🥉 Tier 1: Ecological Scaffolding (The Foundation)

- **Targets:** High-affinity mapping of *S. mutans* GtfB/C/D, *S. wiggsiae*, and *C. albicans* O/N-Mannan receptors.
- **Discovery Strategy:** Multi-species biofilm disruption.

### 🏅 Tier 4: Frugal Manufacturing (The Physical)

- **Hardware:** PepSy (Arduino-based SPPS).
- **Protocol:** Automated translation of digital leads into pump-timing instructions with TPP crosslinker integration.

---

## 🛠️ Stack & Infrastructure

- **Agentic Logic:** LangGraph, LangChain.
- **Bioinformatics:** BioPython, Modeller, OpenFold3.
- **UI/UX:** Streamlit.
- **Deployment:** Zero-budget local/cloud hybrid (WSL2 + NVIDIA NIM).

---

## 🔐 Security & Privacy

This is a **PRIVATE** industrial research repository.

- Transitory simulation data (`results/`) and deep logs (`jobs/`) are git-ignored to maintain repository efficiency.
- Evolutionary trajectories are logged locally for HIPAA-compliant oversight.

---
*Created and maintained by the Dr. Nemo Project Team.*
