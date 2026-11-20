# DrNemo Non-Technical Standard Operating Procedures (SOP) Curriculum
## Streptococcus mutans GtfB Macrocyclic Inhibitor Campaign
**Principal Investigator:** Dr. Roi Castrellon, DDS  
**Affiliation:** DrNemo Oral Biophysics Research, Miami, FL, USA (`roisasha@gmail.com`)  
**Target Receptor:** RCSB PDB 3AIC (2.10 A Resolution) | Counter-Screen: Human Salivary Amylase 1SMD (1.60 A Resolution)  
**Live WebGL 3D Studio:** [Interactive WebGL Viewer](https://drcirosmith.github.io/DrNemo-Caries-Foundry/)  
**Workstation Dashboard:** `streamlit run dashboard.py` (Local Port 8501)  

---

## 🎯 Educational Curriculum Overview
This curriculum translates advanced computational biophysics, LangGraph autonomous workflows, and Karl Popper falsification protocols into clear, intuitive Standard Operating Procedures (SOPs). Written specifically for clinicians, dentists, students, recruiters, and executive leaders, each notebook pairs plain-English clinical analogies with **100% executable Python calculations** against real empirical data.

All notebooks adhere strictly to:
- Zero fabricated numbers: all calculations execute against `data/contenders_ledger.csv`, `demo_data/results/falsification_benchmark.json`, `solubility_shield.py`, and `demo_data/results/3AIC_docked_complex_CAHWRPWHC.pdb`.
- Zero em-dashes in prose strings (following canonical `SOUL.md` standards).
- Full self-contained execution with embedded publication-quality visualizations.

---

## 📚 Curriculum Table of Contents

| SOP ID | Notebook File | Notebook Title | Plain-English Summary | Status |
| :--- | :--- | :--- | :--- | :--- |
| **SOP-00** | [`SOP_00_Interactive_Page_and_Platform_Tour.ipynb`](./SOP_00_Interactive_Page_and_Platform_Tour.ipynb) | Interactive Page & Platform Tour | Step-by-step tour of the WebGL 3D Studio and local 9-tab dashboard with live 3-minute presentation script. | 🟢 Corroborated & Executed |
| **SOP-01 (Guide)** | [`SOP_01_Guide_What_Is_Karl_Popper_Falsification.ipynb`](./SOP_01_Guide_What_Is_Karl_Popper_Falsification.ipynb) | Philosophy of Science: Karl Popper Falsification | Why we break our own AI designs and how 82 contenders filter down to 2 Pareto Rank-1 survivors. | 🟢 Corroborated & Executed |
| **SOP-01 (Protocol)** | [`SOP_01_Biophysical_Falsification_Protocol.ipynb`](./SOP_01_Biophysical_Falsification_Protocol.ipynb) | Biophysical Falsification Protocol | Rigorous technical characterization reproducing publication figures from raw coordinates. | 🟢 Corroborated & Executed |
| **SOP-02** | [`SOP_02_Molecular_Docking_and_Amylase_Selectivity.ipynb`](./SOP_02_Molecular_Docking_and_Amylase_Selectivity.ipynb) | Lock-and-Key Physics & Amylase Selectivity | Understanding 28 A docking boxes, negative kcal/mol, and avoiding human saliva friendly fire. | 🟢 Corroborated & Executed |
| **SOP-03** | [`SOP_03_Scrambled_Null_Trap_and_Alanine_Scan.ipynb`](./SOP_03_Scrambled_Null_Trap_and_Alanine_Scan.ipynb) | Scrambled Null Trap & Alanine Hotspot Scan | Catching false positives with 39 scrambled letters and identifying key-teeth in the catalytic cleft. | 🟢 Corroborated & Executed |
| **SOP-04** | [`SOP_04_Saliva_Survival_Stephan_pH_Switch_and_Proteases.ipynb`](./SOP_04_Saliva_Survival_Stephan_pH_Switch_and_Proteases.ipynb) | Saliva Survival: Stephan Curve & Proteases | Plaque acid dynamics, Henderson-Hasselbalch charge surge, and Arg5-Pro6 steric protease shielding. | 🟢 Corroborated & Executed |
| **SOP-05** | [`SOP_05_Translational_Readiness_and_Clinical_Formulation.ipynb`](./SOP_05_Translational_Readiness_and_Clinical_Formulation.ipynb) | Translational Readiness & Clinical Formulation | Bridging bench-to-chairside with chitosan nanocarriers and master multi-objective Pareto radar profiles. | 🟢 Corroborated & Executed |

---

## 🛠️ How to Execute and Validate the Notebooks

All notebooks run cleanly inside the local project virtual environment:

```bash
# Activate virtual environment
source .venv/bin/activate

# Execute all SOP notebooks in-place with zero cell errors
jupyter nbconvert --to notebook --execute --inplace notebooks/SOP_00_Interactive_Page_and_Platform_Tour.ipynb
jupyter nbconvert --to notebook --execute --inplace notebooks/SOP_01_Guide_What_Is_Karl_Popper_Falsification.ipynb
jupyter nbconvert --to notebook --execute --inplace notebooks/SOP_02_Molecular_Docking_and_Amylase_Selectivity.ipynb
jupyter nbconvert --to notebook --execute --inplace notebooks/SOP_03_Scrambled_Null_Trap_and_Alanine_Scan.ipynb
jupyter nbconvert --to notebook --execute --inplace notebooks/SOP_04_Saliva_Survival_Stephan_pH_Switch_and_Proteases.ipynb
jupyter nbconvert --to notebook --execute --inplace notebooks/SOP_05_Translational_Readiness_and_Clinical_Formulation.ipynb

# Run repository test suite
pytest -q
```

---

## 🛡️ Epistemological and Clinical Guarantees
1. **Evidence Before Claims:** Every statement in these SOPs is directly corroborated by fresh executed Python output.
2. **Zero Fabricated Surrogates:** All numbers are derived from verified empirical files (`falsification_benchmark.json`, `contenders_ledger.csv`, and PDB coordinates).
3. **Non-Antibiotic Paradigm:** Targets the glucan-synthesizing GtfB enzyme specifically, preserving the healthy protective oral microbiome while halting enamel decay.
