<p align="center">
  <img src="assets/hero.svg" alt="Dr. Nemo Caries Vaccine Foundry Banner" width="100%">
</p>

<p align="center">
  <a href="https://github.com/DrCiroSmith/DrNemo-Caries-Foundry/actions/workflows/ci.yml"><img src="https://github.com/DrCiroSmith/DrNemo-Caries-Foundry/actions/workflows/ci.yml/badge.svg" alt="CI"></a>
  <a href="https://opensource.org/licenses/MIT"><img src="https://img.shields.io/badge/License-MIT-blue.svg" alt="License: MIT"></a>
  <a href="https://www.python.org/downloads/"><img src="https://img.shields.io/badge/python-3.10+-blue.svg" alt="Python 3.10+"></a>
  <a href="https://build.nvidia.com"><img src="https://img.shields.io/badge/Inference-NVIDIA_NIM_OpenFold3-76B900.svg" alt="NVIDIA NIM"></a>
  <a href="https://langchain-ai.github.io/langgraph/"><img src="https://img.shields.io/badge/Orchestration-LangGraph-orange.svg" alt="LangGraph"></a>
</p>

An autonomous multi-agent pipeline that designs, filters, and clinically characterizes cyclic peptide inhibitors targeting cariogenic biofilms. It couples LangGraph state machines with NVIDIA NIM OpenFold3 physics, DiffDock ligand poses, and on-device MedGemma clinical rationale to block Streptococcus mutans GtfB enzymes without disrupting commensal flora.

<p align="center">
  <img src="assets/evolution_demo.gif" alt="Directed Evolution and Chitosan Filter Telemetry" width="100%">
</p>

---

## 1. Clinical Rationale and Biological Target

Dental caries affects billions globally, driven by acid-producing bacteria organized in protective biofilms. *Streptococcus mutans* dominates these communities through two virulence factors:
1. **Pellicle Attachment:** Mediated by surface adhesin SpaP/P1 (UniProt P23504) binding salivary proteins.
2. **Glucan Polymerization:** Mediated by glucosyltransferase-B (GtfB, UniProt P08987), which synthesizes water-insoluble alpha-1,3-glucans (mutan). This mutan scaffold cements bacteria to tooth surfaces and traps lactic acid directly against enamel hydroxyapatite.

```
Salivary Pellicle + Hydroxyapatite Enamel
                   │
           [ SpaP Adhesin ]  ──> Initial Cell Anchoring
                   │
           [ GtfB / GtfC ]   ──> Water-Insoluble Mutan Matrix (Biofilm Shield)
                   │
        Localized Acid Trapping ──> Enamel Demineralization (Cavity Formation)
```

Broad-spectrum antimicrobials like chlorhexidine wipe out both pathogens and beneficial oral species, triggering dysbiosis. The Dr. Nemo Foundry designs targeted arginine-rich cyclic peptide wedges. These molecules selectively dock into the catalytic cleft of GtfB and GtfC, halting biofilm architecture without killing oral commensals.

---

## 2. End-to-End Discovery Pipeline

The discovery engine runs an iterative loop across stateful generation, matrix delivery filtering, high-throughput cloud folding, and local medical synthesis:

```mermaid
flowchart TD
    subgraph S1["1. Generative Evolution"]
        A["Scaffold Seed<br/>(19-mer GtfB Catalytic Motif)"] --> B["Directed Mutation Engine<br/>CycleDesigner (T=25)"]
        B --> C["Candidate Cyclic Peptides<br/>(C-X(n)-C Scaffolds)"]
    end

    subgraph S2["2. Delivery Vehicle Filter"]
        C --> D["Solubility Shield Pre-Filter"]
        D -->|Exclude Anionic D,E| E{"Chitosan Matrix Safe?"}
        D -->|Cap Hydrophobic I,K <= 30%| E
        E -->|Reject: Entrapment / Salt Bridge| B
        E -->|Viable Leads| F["Curated Sequence Pool"]
    end

    subgraph S3["3. Structural Validation"]
        F --> G["NVIDIA NIM OpenFold3 API"]
        G --> H["Multi-Key Rate Rotator<br/>(Exponential Backoff)"]
        H --> I["Composite Scoring Engine<br/>Score = (pLDDT / 100) * iPTM"]
        I --> J{"Score >= 0.72 &<br/>iPAE < 0.35?"}
        J -->|Fail| B
        J -->|Pass| K["Structural Hits"]
    end

    subgraph S4["4. Validation and Synthesis"]
        K --> L["DiffDock Molecular Docking<br/>Binding Pocket Confirmation"]
        L --> M["Local MedGemma 4B<br/>Clinical Rationale and IRB Report"]
        M --> N["PepSy SPPS Synthesis Protocol<br/>(TPP Crosslinker Formula)"]
        N --> O["Streamlit 3D Dashboard<br/>(Py3Dmol Structure Viewer)"]
    end
```

---

## 3. Subsystem Architecture

### LangGraph State Machine (`gtfb_orchestrator.py`)
A stateful multi-node graph orchestrates each evolutionary cycle:
- **`cycle_designer_node`:** Produces arginine-enriched cyclic peptides flanked by cysteines (`C...C`) with diffusion temperature calibrated to T=25.
- **`evaluate_openfold_node`:** Dispatches candidates to structural prediction services.
- **`filter_composite_node`:** Filters candidates against structural and interaction benchmarks.
- **`biosecurity_check_node`:** Screens sequences against International Gene Synthesis Consortium (IGSC) biosecurity protocols.
- **`synthesis_formatter_node`:** Outputs machine-readable pump timings for solid-phase peptide synthesis (SPPS).

### Solubility Shield and Chitosan Vehicle Rules (`solubility_shield.py`)
Because topical oral peptides must release from a mucoadhesive chitosan buccal patch, each candidate must satisfy strict physical chemistry constraints:
- **Anionic Cap ($D + E \le 2$):** Eliminates polyelectrolyte precipitation caused by negative residues binding polycationic chitosan chains.
- **Entrapment Limit ($I + K \le 30\%$):** Prevents ring-mediated hydrogen bonding that traps peptides permanently inside the formulation.
- **Arginine Core:** Enforces high positive charge density ($\text{pI} > 10$) to seek out the electronegative aspartate pocket of GtfB while maintaining electrostatic repulsion from the chitosan delivery matrix.
- **Salivary Hydropathy:** Restricts GRAVY hydropathy index $\le 2.0$ and instability index $\le 40.0$ at salivary pH 6.8.

### NVIDIA NIM OpenFold3 Inference (`pipeline.py`)
Predicts 3D structure and inter-chain binding confidence:
$$\text{Composite Score} = \left(\frac{\text{pLDDT}}{100}\right) \times \text{iPTM}$$
A candidate graduates when Composite Score $\ge 0.72$ with predicted aligned error $\text{iPAE} < 0.35$. An asynchronous worker pool balances jobs across multiple API keys, catching rate limits with jittered backoff.

### DiffDock Binding Confirmation (`diffdock_pipeline.py`)
Performs blind molecular docking of generated candidates against crystallographic GtfB coordinates, verifying contact surface geometry and binding poses.

### On-Device MedGemma Clinical Translation (`analyzer.py`)
A local 4-bit quantized MedGemma instance inspects atomic coordinate outputs:
- Checks predicted residue interactions, charge balance, and steric fit.
- Produces clinical justification memos formatted for regulatory and institutional review boards.

### Interactive 3D Visualization (`dashboard.py`)
A Streamlit UI provides:
- Live generation telemetry and score progression across targets.
- Embedded Py3Dmol WebGL viewports for interactive 3D rotation of generated PDB structures.
- Filterable candidate tables with full export capability.

---

## 4. Verification and Empirical Benchmarks

The figures below summarize screening runs across 500 mutated candidate sequences.

<p align="center">
  <img src="assets/physicochemical_benchmark.png" alt="Physicochemical Benchmark and Validation Results" width="100%">
</p>

- **Panel A (Structural Gate):** Displays the distribution of predicted structural fold confidence (pLDDT) against interface interaction score (iPTM). The gold dashed curve marks the 0.72 threshold. Candidates in the upper right quadrant enter downstream synthesis planning.
- **Panel B (Chitosan Safety):** Tracks filter performance. Out of 500 generated variants, 48.4% met all chitosan compatibility rules, while 23.6% were removed due to anionic charge hazards and 18.8% exceeded entrapment limits.
- **Panel C (Async Inference):** Demonstrates throughput scaling from 18 predictions/min on a single endpoint to steady ~68 predictions/min across a 4-key worker pool with automated backoff.

---

## 5. Repository Layout

```
DrNemo-Caries-Foundry/
├── .github/workflows/
│   └── ci.yml               # Automated GitHub Actions test pipeline
├── assets/
│   ├── hero.svg             # Vector architecture banner
│   ├── evolution_demo.gif   # Directed evolution telemetry animation
│   ├── physicochemical_benchmark.png # 300 DPI multi-panel validation figure
│   └── social_preview.png   # Metadata preview image
├── scripts/
│   └── generate_assets.py   # Asset generation script
├── tests/
│   └── test_foundry.py      # Unit test suite
├── analyzer.py              # MedGemma local clinical synthesis
├── auto_pilot.py            # Autonomous lab worker
├── automated_breeder.py     # Evolutionary breeder
├── autonomous_lab.py        # Multi-workstation supervisor
├── compliance_logger.py     # Audit and provenance logger
├── dashboard.py             # Streamlit 3D visualization platform
├── diffdock_pipeline.py     # Molecular docking harness
├── directed_evolution.py    # Mutation and selection engine
├── export_report.py         # Summary report generator
├── gtfb_orchestrator.py     # LangGraph discovery state machine
├── mutation_scanner.py      # Sequence perturbation module
├── phase_report.py          # Milestone assessment tool
├── pipeline.py              # NVIDIA NIM OpenFold3 async client
├── solubility_shield.py     # Chitosan-native pre-filter
├── targets.json             # Oral pathogen virulence targets
└── requirements.txt         # Project dependencies
```

---

## 6. Getting Started

### Prerequisites
- Python 3.10 or higher
- Optional: Local Ollama instance with `gemma` or `llama3` for offline text synthesis
- Optional: NVIDIA Developer API key for cloud OpenFold3 predictions

### Installation
```bash
git clone https://github.com/DrCiroSmith/DrNemo-Caries-Foundry.git
cd DrNemo-Caries-Foundry

python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

### Running Tests
```bash
pytest tests/ -v
```

### Running the Orchestration Graph
```bash
python3 gtfb_orchestrator.py
```

### Starting the 3D Dashboard
```bash
streamlit run dashboard.py
```

---

## 7. Virulence Target Inventory

Target definitions are cataloged in `targets.json`:
- **Tier 1 (Surface Adhesin):** *S. mutans* SpaP/P1 Adhesin (`P23504`)
- **Tier 2 (Glucan Synthesis):** *S. mutans* GtfB Insoluble Glucan Synthase (`P08987`), GtfC Structural Bridge (`P13470`), GtfD Soluble Primer (`P49331`)
- **Tier 3 (Polymicrobial Synergy):** *Candida albicans* Mannan receptors, *Scardovia wiggsiae* adherence markers, *Porphyromonas gingivalis* Gingipain proteases

---

## 8. Research Notice

This codebase is for academic, computational biology, and preclinical discovery research. Molecular predictions and synthesis plans require wet-lab in-vitro confirmation prior to therapeutic application.

---

## 9. Author

**Dr. Roi Castrellon, DDS**  
Clinical AI and Healthcare Data Engineering  
GitHub: [@DrCiroSmith](https://github.com/DrCiroSmith)
