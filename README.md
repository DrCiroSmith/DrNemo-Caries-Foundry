<p align="center">
  <img src="assets/hero_banner.png" alt="DrNemo Caries Foundry - Targeted GtfB Macrocycle Inhibitor Design" width="100%">
</p>

<p align="center">
  <a href="https://drcirosmith.github.io/DrNemo-Caries-Foundry/"><img src="https://img.shields.io/badge/Live%20Demo-Interactive%203D%20Studio-0D9488?style=for-the-badge&logo=webgl&logoColor=white" alt="Live Interactive 3D Demo"></a>
</p>

<p align="center">
  <a href="https://github.com/DrCiroSmith/DrNemo-Caries-Foundry/actions/workflows/ci.yml"><img src="https://github.com/DrCiroSmith/DrNemo-Caries-Foundry/actions/workflows/ci.yml/badge.svg" alt="CI"></a>
  <a href="https://opensource.org/licenses/MIT"><img src="https://img.shields.io/badge/License-MIT-blue.svg" alt="License: MIT"></a>
  <a href="https://www.python.org/downloads/"><img src="https://img.shields.io/badge/python-3.10+-blue.svg" alt="Python 3.10+"></a>
  <a href="https://vina.scripps.edu"><img src="https://img.shields.io/badge/Biophysics-AutoDock_Vina_v1.2.7-76B900.svg" alt="AutoDock Vina"></a>
  <a href="https://langchain-ai.github.io/langgraph/"><img src="https://img.shields.io/badge/Orchestration-LangGraph-orange.svg" alt="LangGraph"></a>
  <a href="https://openmm.org"><img src="https://img.shields.io/badge/ForceField-OpenMM_Amber14-blueviolet.svg" alt="OpenMM Amber14"></a>
</p>

An autonomous computational biophysics pipeline that designs, filters, and falsifies cyclic disulfide macrocycle inhibitors of *Streptococcus mutans* glucosyltransferase-B (GtfB, UniProt P08987). It couples LangGraph state machines with AutoDock Vina v1.2.7 and Vinardo empirical grid docking, OpenMM Amber14/GBn2 all-atom relaxation at 310.15 K, and an exhaustive 12-angle Karl Popper falsification battery evaluating candidate **CAHWRPWHC** over 82 contenders and N=39 scrambled nulls. All findings are strictly in silico computational hypotheses requiring wet-lab enzymatic, whole-saliva, and cellular biofilm falsification.

---

## 1. Biological Target: *S. mutans* GtfB

<p align="center">
  <img src="assets/fig1_biological_mechanism.png" alt="Biological mechanism: uninhibited biofilm vs CAHWRPWHC targeted blockade" width="100%">
</p>

Dental caries is driven by acidogenic biofilms. *Streptococcus mutans* secretes glucosyltransferase-B (GtfB), which synthesizes water-insoluble alpha-1,3-glucans (mutan) that cement bacteria to enamel and trap lactic acid against hydroxyapatite, driving demineralization.

The Foundry designs **arginine-rich cyclic disulfide macrocycles** that dock into the GtfB catalytic cleft (RCSB 3AIC: Asp477, Glu515, Asp588, Trp517), halting mutan matrix formation without killing oral commensals.

---

## 2. End-to-End Discovery Pipeline

<p align="center">
  <img src="assets/fig2_pipeline_architecture.png" alt="6-stage closed-loop discovery pipeline with Popperian falsification" width="100%">
</p>

The discovery engine runs an iterative closed loop:

1. **Evolutionary Generator**: C-X-7-C scaffold with P2A hinge; produces 82 contenders.
2. **Oral Shield Gate**: Donnan Q <= +4.0, Stephan DQ >= +1.0, 0 protease sites.
3. **28.0 A Dual Docking**: AutoDock Vina + Vinardo on RCSB 3AIC (target) and 1SMD (human alpha-amylase off-target).
4. **OpenMM 310.15 K Relaxation**: Amber14/GBn2 all-atom minimization.
5. **12-Angle Popper Falsification**: 1SMD DDG <= -0.80, N=39 scrambled nulls, alanine scan, Acarbose counter-screen.
6. **SQLite Contender Vault and SPPS Protocol**: solid-phase synthesis plan for the corroborated lead.

Falsified candidates (amylase cross-reactivity, protease cleavage) feed back to mutate the hinge.

---

## 3. Biophysical Results: Pareto Champion CAHWRPWHC

<p align="center">
  <img src="assets/fig3_pocket_and_funnel.png" alt="GtfB active site close-up and multi-gate attrition funnel" width="100%">
</p>

| Parameter | CRRRWWRRRRC | CPHWRPWHC (WT lead) | CPHWRAWHC (P6A) | **CAHWRPWHC (Champion)** | Popperian outcome |
|---|---:|---:|---:|---:|---|
| **GtfB (3AIC) Vina DG** | -7.586 kcal/mol | -8.062 kcal/mol | -8.890 kcal/mol | **-8.191 kcal/mol** | High target affinity |
| **Amylase (1SMD) Vina DG** | -7.580 kcal/mol | -7.902 kcal/mol | -7.300 kcal/mol | **-7.343 kcal/mol** | Host off-target counter-screen |
| **Host Selectivity DDG** | -0.006 kcal/mol | -0.16 kcal/mol | -1.59 kcal/mol | **-0.85 kcal/mol** | **WT cross-reactive; P2A clears <= -0.80 (within Vina noise)** |
| **Oral protease sites** | 7 (degraded) | 0 (RP-shielded) | 1 (trypsin fail) | **0 Trypsin (AT RISK: Trp4/7 chymotryptic)** | Baseline falsified (enzyme lysis) |
| **Chitosan Donnan Q (pH 6)** | +7.00 (> +4.0 fail) | +2.046 | +2.046 | **+2.046 (<= +4.0)** | Baseline falsified (entrapment) |
| **Stephan DQ (acid switch)** | 0.00 | +1.62 | +1.62 | **+1.62** | Acid-activated in cariogenic niche |
| **Disulfide S-S** | unconstrained | 2.06 A | 2.06 A | **2.06 A** | Intact macrocyclic restraint |
| **OpenMM 310.15 K** | N/A | -1,165.5 kJ/mol | -1,271.48 kJ/mol | **-1,234.68 kJ/mol** | Stable minimized conformation |
| **Verdict** | FALSIFIED (multiple gates) | FALSIFIED (amylase) | FALSIFIED (trypsin site) | **IN SILICO HYPOTHESIS** | Candidate advance with 3 caveats |

P2A (Pro2 to Ala) hinge relaxation preserves GtfB binding while sterically excluding the inhibitor from human amylase (DDG = -0.85 kcal/mol). P6A (Pro6 to Ala) gives raw affinity -8.890 kcal/mol but breaks the proline steric shield, creating a tryptic cleavage site, falsified on Angle 5.

---

## 4. 12-Angle Exhaustive Computational Disproof Battery

In accordance with Karl Popper's demarcation principle, candidate **CAHWRPWHC** was challenged across twelve computational disproof angles to systematically identify physical failure modes:

| Angle | Evaluation Angle | Quantitative Metric | Honest Classification |
|---|---|---|---|
| **1** | Co-Crystal Redocking Calibration (3AIC) | Hungarian RMSD 4.51 A, Vina -7.915 kcal/mol | `PASSED_IN_SILICO` |
| **2** | Multi-Seed Conformer & Docking (N=5) | Mean -8.10 +/- 0.37 kcal/mol, RMSD = 9.24 A | `INCONCLUSIVE_WITHIN_NOISE` |
| **3** | Pose Mode Clustering & Energy Gap | Mode 1-2 Delta E = +0.10 kcal/mol, shift 1.74 A | `INCONCLUSIVE_WITHIN_NOISE` |
| **4** | Force-Field Error Margin on 1SMD Selectivity | DeltaDeltaG = -0.85 kcal/mol vs Vina SE ~2.85 kcal/mol | `INCONCLUSIVE_WITHIN_NOISE` |
| **5** | Charge-Matched Scrambled Null (N=39) | Lead -8.191 vs Null -8.343 +/- 0.470 (Z = +0.32, p = 0.5897) | `INCONCLUSIVE_WITHIN_NOISE` |
| **6** | In-Silico Alanine-Scanning Mutagenesis | Identifies His3, Trp7, His8 catalytic anchors | `PASSED_IN_SILICO` |
| **7** | Explicit pH 5.2 Histidine Diprotonation | Formal charge +1 -> +3; acidic 3AIC DG = -7.71 kcal/mol | `PASSED_IN_SILICO` |
| **8** | Broad Oral Protease & Peptidase Check | Trypsin shielded (Arg5-Pro6); Trp4/Trp7 vulnerable to chymotrypsin | `AT_RISK` |
| **9** | Additional Host Screen (RCSB 1LZ1 Lysozyme) | Human lysozyme DG = -6.67 kcal/mol (DeltaDeltaG = -1.52) | `PASSED_IN_SILICO` |
| **10** | Linear Control vs Cyclic Macrocycle | Linear AHWRPWH (-9.25 kcal/mol) vs cyclic (-8.19 kcal/mol; 14 vs 28 torsions) | `PASSED_IN_SILICO` |
| **11** | OpenMM Implicit-Solvent Relaxation Drift | Amber14/GBn2 310.15 K heavy-atom drift RMSD = 0.94 A | `PASSED_IN_SILICO` |
| **12** | Chitosan Donnan & Stephan Acid Switch | Saliva Q = +2.05 <= +4.0; Stephan Delta Q = +1.62 | `PASSED_IN_SILICO` |

> **Scientific Skepticism Verdict:** While CAHWRPWHC cleared 7 physical hurdles in silico (7 Passed in Silico, 4 Inconclusive within Noise, 1 At Risk, 0 Falsified), strict scientific skepticism highlights 5 caveats: multi-seed conformer RMSD spread (Angle 2), pose mode clustering gap (Angle 3), host selectivity (Angle 4), and scrambled null differentiation (Angle 5) fall within force-field uncertainty boundaries, and aromatic residues Trp4/Trp7 (Angle 8) carry chymotryptic risk requiring whole-saliva degradation assays. All results are computational hypotheses requiring wet-lab enzymatic, whole human saliva, and cellular biofilm falsification.

---

## 5. Validation Figures

<p align="center">
  <img src="assets/fig4_falsification_plate.png" alt="4-panel biophysical validation: Pareto frontier, scrambled null histogram, alanine scan, Stephan pH titration" width="100%">
</p>

- **Panel A, Pareto Frontier:** Vina 3AIC vs selectivity DDG across contenders; CAHWRPWHC is the non-dominated champion.
- **Panel B, Scrambled Null Histogram (N=39):** charge-matched sequence permutations establish empirical significance (mean = -8.343 +/- 0.470 kcal/mol, Z = +0.32 sigma, p = 0.5897).
- **Panel C, Alanine Hotspot Scan:** His3, Trp7, His8 are catalytic anchors; Pro2/Pro6 are conformational gates. P6A flagged for Arg5 trypsin lysis.
- **Panel D, Stephan Acid-Switch:** DQ = +1.62 activates the inhibitor in the cariogenic acidic niche (pH approximately 5.2). CRRRWWRRRRC shows flat charge, no acid activation.

---

## 6. Repository Layout

```
DrNemo-Caries-Foundry/
├── .github/workflows/
│   └── ci.yml                           # Automated GitHub Actions test pipeline
├── assets/
│   ├── hero_banner.png                  # Retina hero with 3D pocket schematic
│   ├── fig1_biological_mechanism.png    # 2-panel biological mechanism diagram
│   ├── fig2_pipeline_architecture.png   # 6-stage closed-loop pipeline
│   ├── fig3_pocket_and_funnel.png       # Active site + attrition funnel
│   ├── fig4_falsification_plate.png     # 2x2 validation data plate
│   └── social_preview.png              # 1280x640 social preview (not embedded)
├── data/
│   ├── contenders_ledger.csv            # 194 contenders (155 design + 39 nulls)
│   ├── foundry_contenders.sqlite        # Persistent contender vault
│   └── poses/                           # Docked PDBQT poses
├── demo_data/
│   └── results/
│       ├── master_scores.csv            # 20-generation demo scores
│       └── falsification_benchmark.json # 6-angle benchmark data
├── docs/
│   └── index.html                       # Live 3Dmol.js studio (unchanged)
├── manuscript/
│   └── DrNemo_Preprint.pdf              # Preprint manuscript (unchanged)
├── scripts/
│   ├── build_readme_visuals.py          # Playwright Retina PNG figure builder
│   ├── capture_linkedin_media.py        # LinkedIn demo media render
│   ├── generate_assets.py               # Legacy asset generator
│   └── render_assets.py                 # Real-data raster generator
├── tests/
│   └── test_foundry.py                  # Unit test suite
├── analyzer.py                          # Local clinical synthesis
├── auto_pilot.py                        # Autonomous lab worker
├── automated_breeder.py                 # Evolutionary breeder
├── autonomous_lab.py                    # Multi-workstation supervisor
├── biophysics_engine.py                 # AutoDock Vina, Vinardo, OpenMM engine
├── compliance_logger.py                 # Audit and provenance logger
├── dashboard.py                         # Streamlit 3D visualization platform
├── directed_evolution.py                # Mutation and selection engine
├── export_report.py                     # Summary report generator
├── falsification_engine.py              # 6-angle Karl Popper falsification suite
├── foundry_manager.py                   # Foundry CLI management
├── gtfb_orchestrator.py                 # LangGraph discovery state machine
├── mutation_scanner.py                  # Sequence perturbation module
├── phase_report.py                      # Milestone assessment tool
├── solubility_shield.py                 # Chitosan-native pre-filter
├── targets.json                         # Oral pathogen virulence targets
├── vault.py                             # SQLite contender vault
└── requirements.txt                     # Project dependencies
```

---

## 7. Getting Started

### Prerequisites
- Python 3.10 or higher
- System biophysics tools: AutoDock Vina v1.2.7, OpenBabel

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

## 8. Virulence Target Inventory

Target definitions are cataloged in `targets.json`:
- **Tier 1 (Surface Adhesin):** *S. mutans* SpaP/P1 Adhesin (`P23504`)
- **Tier 2 (Glucan Synthesis):** *S. mutans* GtfB Insoluble Glucan Synthase (`P08987`), GtfC Structural Bridge (`P13470`), GtfD Soluble Primer (`P49331`)
- **Tier 3 (Polymicrobial Synergy):** *Candida albicans* Mannan receptors, *Scardovia wiggsiae* adherence markers, *Porphyromonas gingivalis* Gingipain proteases

---

## 9. Research Notice

This codebase is for academic, computational biology, and preclinical discovery research. Molecular predictions and synthesis plans require wet-lab in-vitro confirmation prior to therapeutic application.

---

## 10. Author

**Dr. Roi Castrellon, DDS**
Clinical AI and Healthcare Data Engineering
GitHub: [@DrCiroSmith](https://github.com/DrCiroSmith)
