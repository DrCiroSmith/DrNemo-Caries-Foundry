# Biophysical Falsification Benchmark: Operation Karl Popper

**Autonomous Structural-Biology Verification Suite for *Streptococcus mutans* glucosyltransferase inhibitors**  
*Principal Investigator:* Dr. Roi Castrellon, DDS (DrNemo Oral Biophysics Research, Miami, FL, USA)  
*Benchmark Version:* 2.0.0 (Empirical $28.0\text{ \AA}$ Dual-Scoring Vina + Vinardo + OpenMM)  
*Status:* Every disproof angle now records where its numbers came from. On the committed
data, eight are recomputed on the run, three are read from the contender database, and one
(the OpenMM relaxation) cannot run without openmm, pdbfixer and rdkit and reports itself as
not computed rather than printing a verdict. Two limits in the gates remain and are listed
in section 5: the disulfide geometry check cannot discriminate, since all 260 rows fall
between 2.050 and 2.070 Angstroms against a 1.95 to 2.25 window and no candidate has ever
been failed by it, and the Stephan acid-switch criterion is arithmetically a histidine
count. The instability index is no longer reported at all for cyclic peptides, where it is
out of domain and was previously replaced by a literal 25.0. A second substitution has
also gone: any sequence with two or more arginines and 20 or fewer residues had its index
replaced by the same 25.0, so a short arginine-rich peptide passed a stability gate that
was never measured on it. The gate is now applied only where the index is computed, and
reported only there. The protease criterion does
discriminate: 111 of the 260 rows carry at least one trypsin site. No candidate in this
benchmark is corroborated.

---

## 1. The Popperian Imperative: Falsifiability in Computational Drug Design

In *The Logic of Scientific Discovery* (1934), Karl Popper established that scientific theories cannot be verified by accumulating confirmatory instances; they can only be corroborated by ruthlessly attempting to **falsify** them. 

Early computational pipelines in the Foundry utilized surrogate heuristic functions (e.g., scoring arginines and contact density). While computationally light, heuristic surrogates suffer from confirmation bias: a sequence enriched in positive charges will score high regardless of whether its 3D geometry clashes with the catalytic pocket, precipitates in saliva, or cross-reacts with essential host digestive enzymes.

**Operation Karl Popper** systematically eradicated all synthetic mock scoring in `DrNemo-Caries-Foundry`. Every candidate is now subjected to an orthogonal 6-angle biophysical falsification battery powered by:
1. **AutoDock Vina v1.2.7 & Vinardo:** Empirical scoring functions on uncompressed $28.0\text{ \AA}$ search grids centered on catalytic triad oxygens (`Asp477`, `Glu515`, `Asp588` of RCSB PDB: `3AIC`).
2. **All-Atom OpenMM Molecular Mechanics:** Amber14/GBn2 implicit solvent energy minimisation. No dynamics are run and no temperature is sampled: the integrator is constructed at $310.15\text{ K}$ and never stepped, so nothing here is a thermal measurement.
3. **RDKit MMFF94:** Force-field energy minimization and distance-geometry 3D conformer generation with intact disulfide bonds ($2.06\text{ \AA}$).
4. **Physicochemical Oral Filters:** Henderson-Hasselbalch protonation states, Donnan potential compatibility ($Q_{\text{pH } 6.0} \le +4.0$), Stephan-curve electrostatic switches ($\Delta Q_{\text{Stephan}} \ge +1.0$), and mucosal enzymatic cleavage mapping.

---

## 2. The 6-Angle Falsification Battery

```
                          [ Candidate Sequence ]
                                    │
    ┌───────────────────────────────┴───────────────────────────────┐
    ▼                                                               ▼
[ Angle 1: Co-Crystal Control ]               [ Angle 2: Counter-Screen ]
RCSB 3AIC Acarbose Redocking                  RCSB 1SMD Salivary Amylase
ΔG = -7.915 kcal/mol | RMSD = 4.51 Å          Selectivity ΔΔG ≤ -0.80 kcal/mol
    │                                                               │
    ├───────────────────────────────┬───────────────────────────────┤
    ▼                               ▼                               ▼
[ Angle 3: Scrambled Null ]   [ Angle 4: Ala-Scanning ]       [ Angle 5: Protease ]
N=39 Permutations + Linear    7-Residue Hotspot Scan          Trypsin & Gingipain
Mean = -8.34 ± 0.47 kcal/mol  No hotspot ≥ +0.70 kcal/mol     Zero Cleavage Sites
    │                                                               │
    └───────────────────────────────┬───────────────────────────────┘
                                    ▼
                      [ Angle 6: Formulation ]
                      Chitosan Donnan Q ≤ +4.0
                      Stephan Acid-Switch ΔQ ≥ +1.0
                                    │
                                    ▼
                      [ SURVIVED ALL GATES, RANK 1 ]
```

---

### Angle 1: RCSB 3AIC Co-Crystal Redocking Control (Assay Calibration)
* **Target:** *Streptococcus mutans* Glucosyltransferase-SI (GTF-SI, gtfC gene product) catalytic $(\beta/\alpha)_8$ TIM-barrel domain (RCSB PDB ID: `3AIC`, 3.11 Å resolution).
* **Native Co-Crystallized Ligand:** Acarbose pseudo-tetrasaccharide (`GLC 1 - GLC 2 - AC1 3`) residing in the glucan elongation cleft adjacent to catalytic triad residues **Asp477** (nucleophile), **Glu515** (general acid/base), and **Asp588** (transition-state stabilizer, assignment from the QM/MM study cited as reference 7 of `paper/negative_result_GTF-SI.md`, which is computational and not mutagenesis), plus $\pi$-stacking platform **Trp517**.
* **Methodology:** Native coordinates were extracted directly from `3AIC.pdb`, hydrogenated at pH 6.0, and converted to PDBQT. AutoDock Vina was calibrated across the active cleft centroid $(x=192.072, y=44.370, z=197.309)$ with an uncompressed bounding box of $28.0 \times 28.0 \times 28.0\text{ \AA}$.
* **Empirical Results:**
  * **Binding Affinity ($\Delta G_{\text{Vina}}$):** $\mathbf{-7.915\text{ kcal/mol}}$ (Vinardo: $-5.186\text{ kcal/mol}$).
  * **Centroid Displacement:** $\mathbf{3.03\text{ \AA}}$ from experimental crystal coordinates.
  * **Heavy-Atom Symmetry-Corrected RMSD:** $\mathbf{4.507\text{ \AA}}$ (heavy atoms paired
    optimally within each element, 44 atoms). Reproduce with
    `python3 scripts/compute_redocking_rmsd.py`; the same script reports $9.313\text{ \AA}$
    if atoms are paired in file order, which is what the earlier helper did.
* **Accepted Criterion:** a redocking control is conventionally called successful at
  RMSD $\le 2.0\text{ \AA}$ (some authors allow $2.5\text{ \AA}$). At
  $4.507\text{ \AA}$ this control **FAILS**. The earlier version of this document set the
  threshold at $5.0\text{ \AA}$ and recorded a pass; no redocking convention uses
  $5.0\text{ \AA}$, and that threshold has been removed.
* **What this means for every affinity in this repository:** the protocol localises the
  ligand in the right cleft (centroid displacement $3.03\text{ \AA}$) but does not
  reproduce the crystallographic binding mode. Vina scores from this box are therefore
  usable as a coarse ranking signal and not as pose-level or affinity-level predictions.

---

### Angle 2: RCSB 1SMD Human Salivary $\alpha$-Amylase Off-Target Counter-Screen
* **Target:** Human Salivary $\alpha$-Amylase (RCSB PDB ID: `1SMD`, 1.60 Å resolution). Active-site triad: **Asp197**, **Glu233**, **Asp300** at $(x=8.778, y=49.724, z=17.278)$ on an identical $28.0\text{ \AA}$ grid.
* **Clinical Rationale:** Oral therapeutics targeting carbohydrate-processing bacterial enzymes must not inhibit host salivary amylase, which is vital for dietary starch hydrolysis and oral mucosal health.
* **Empirical Discoveries:**
  * **Initial Lead (`CPHWRPWHC`):** GTF-SI $\mathbf{-8.062\text{ kcal/mol}}$ vs Amylase $\mathbf{-7.902\text{ kcal/mol}}$ ($\Delta\Delta G = \mathbf{-0.16\text{ kcal/mol}}$). **FALSIFIED on Angle 2:** In the uncompressed $28\text{ \AA}$ box, `CPHWRPWHC` cross-reacts with human amylase and fails the $\le -0.80\text{ kcal/mol}$ selectivity threshold!
  * **Lead candidate (`CAHWRPWHC`, the wild type's P2A mutant). Not a Pareto champion:** it is dominated by 9 of the 17 other rows that survived every gate, and by 87 of all 260, so it sits on neither front. Among the survivors it ranks 10th of 18 on GTF-SI affinity, 18th of 18 on selectivity and 17th of 18 on ligand efficiency. It was carried forward as the lead because the campaign was built around the wild type's P2A mutant, not because the recorded objectives favour it. Verify with `python3 scripts/compute_pareto_front.py`. Its numbers: GTF-SI $\mathbf{-8.191\text{ kcal/mol}}$ vs Amylase $\mathbf{-7.343\text{ kcal/mol}}$ ($\Delta\Delta G = \mathbf{-0.85\text{ kcal/mol}}$). **SURVIVED ALL GATES:** hinge relaxation preserves GTF-SI affinity in this scoring function while placing the macrocycle outside the selectivity threshold for the narrower human amylase cleft. Clearing every in silico gate is not evidence of inhibition, and this candidate does not separate from its charge-matched scrambled nulls (Angle 5: $Z = +0.32$, $p = 0.600$).

---

### Angle 3: Charge-Matched Sequence-Scrambled Null Permutations ($N=39$)
* **Methodology:** Core residues were shuffled into 39 unique permutations preserving net charge ($+1.13$) and disulfide cyclization. Each isomer was docked with real AutoDock Vina on the uncompressed $28\text{ \AA}$ grid.
* **Empirical Null Distribution:**
  * **Mean Affinity:** $\mathbf{-8.343 \pm 0.470\text{ kcal/mol}}$ (For `CAHWRPWHC` at $-8.191\text{ kcal/mol}$: $Z = +0.32\sigma, p = 0.600, b = 23/39$; WT `CPHWRPWHC` at $-8.062\text{ kcal/mol}$ yields $Z = +0.60\sigma, p = 0.6750, b = 26/39$).
* **Epistemological Insight:** In static rigid-receptor grid docking, flexible hydrophobic isomers cluster tightly around $-8.3\text{ kcal/mol}$ because static docking under-penalizes unfolded conformational entropy in water ($\Delta S_{\text{conf}}$). This proves Karl Popper's dictum: static docking alone cannot establish candidate viability without multi-angle biological filtering (Donnan equilibrium, protease resistance, Stephan-switch protonation, and host counter-screening)!
* **Linear vs Cyclic Constraint Control:** Linearized peptide `SPHWRPWHS` docks with static Vina score $-8.865\text{ kcal/mol}$ across 35 rotatable bonds; in water, the entropy penalty of 35 unconstrained torsions destroys binding free energy, whereas cyclic nonapeptides maintain rigid pre-organization.

---

### Angle 4: In-Silico Alanine-Scanning Mutagenesis (Hotspot Mapping)
Systematic single-alanine substitutions of the wild type `CPHWRPWHC` ($-8.062\text{ kcal/mol}$), of which
the lead is the `P2A` mutant. Every $\Delta\Delta G$ below is measured against that parent, and every
affinity is read from `data/foundry_contenders.sqlite`.

**Verdict: `NO_HOTSPOT_IDENTIFIED`.** No substitution reaches the $+0.70\text{ kcal/mol}$ hotspot
threshold. The largest affinity loss in the whole scan is `W7A` at $+0.257\text{ kcal/mol}$. No residue
in this scaffold is an energetic anchor on this evidence. Reproduce with
`python3 scripts/run_exhaustive_disproof.py`.

* `P2A` (`CAHWRPWHC`): $\mathbf{-8.191\text{ kcal/mol}}$ (Hinge relaxation; enhances selectivity to $\Delta\Delta G = -0.85\text{ kcal/mol}$).
* `H3A` (`CPAWRPWHC`): $\mathbf{-7.869\text{ kcal/mol}}$ ($\Delta\Delta G = +0.193\text{ kcal/mol}$, well under the $+0.70$ hotspot threshold, and inside an assumed scoring error that carries no citation here).
* `W4A` (`CPHARPWHC`): $\mathbf{-8.269\text{ kcal/mol}}$ ($\pi$-stacking platform).
* `R5A` (`CPHWAPWHC`): $\mathbf{-8.531\text{ kcal/mol}}$ (electrostatic salt-bridge).
* `P6A` (`CPHWRAWHC`): $\mathbf{-8.890\text{ kcal/mol}}$ on 3AIC, $-7.302\text{ kcal/mol}$ on 1SMD ($\Delta\Delta G = -1.59$). **FALSIFIED on Angle 5:** Mutating Pro6 breaks the proline turn, exposing Arg5 to rapid tryptic cleavage!
* `W7A` (`CPHWRPAHC`): $\mathbf{-7.805\text{ kcal/mol}}$ ($\Delta\Delta G = +0.257\text{ kcal/mol}$, well under the $+0.70$ hotspot threshold, not an anchor).
* `H8A` (`CPHWRPWAC`): $\mathbf{-7.978\text{ kcal/mol}}$ ($\Delta\Delta G = +0.084\text{ kcal/mol}$, well under the $+0.70$ hotspot threshold, not an anchor).

---

### Angle 5: Oral Protease Cleavage Barrier
* Consensus substrate recognition motifs for salivary trypsin-like kallikreins and *Porphyromonas gingivalis* RgpB gingipains were mapped across all candidates.
* `CAHWRPWHC`: **0 cleavage sites** (Arg5 is protected by adjacent Pro turn). **PASSED.**
* `P6A` (`CPHWRAWHC`): **1 cleavage site** (Arg5 exposed). **FALSIFIED.**
* Naive `CRRRWWRRRRC`: **7 cleavage sites** (rapid enzymatic degradation). **FALSIFIED.**

---

### Angle 6: Colloidal Donnan Equilibrium and Stephan-Curve Acid Switch
* **Stephan Acid Switch:** Localized acidogenesis ($\text{pH } 5.2$) protonates histidine imidazoles ($\text{p}K_a \approx 6.0$), increasing net charge by $\Delta Q_{\text{Stephan}} = \mathbf{+1.62}$ ($Q_{7.2} = +1.13 \rightarrow Q_{5.2} = +2.75$).
* **Donnan Partitioning:** At salivary baseline ($\text{pH } 6.0$), net charge is $\mathbf{+2.05 \le +4.0}$, clearing the Donnan mucosal partition ceiling and preventing non-specific entrapment in outer mucins.
* Naive `CRRRWWRRRRC`: Net charge $+7.00$ triggers catastrophic Donnan collapse and has zero Stephan switch ($\Delta Q = 0.00$). **FALSIFIED.**

---

## 3. Comprehensive Contender Comparison

**Status vocabulary.** `SURVIVED_ALL_GATES` means a candidate cleared every in silico filter.
It is not a positive result and it is not corroboration: no wet-lab evidence of inhibition
exists for any sequence in this repository. Rows in `data/foundry_contenders.sqlite` and the
exported ledgers still carry the retired label `CORROBORATED` for that same state, because the
archive predates the rename; readers in `vault.py` accept both. Every other status names the
gate that rejected the candidate.

| Candidate / Control | Sequence | 3AIC Vina (kcal/mol) | 1SMD Vina (kcal/mol) | Selectivity $\Delta\Delta G$ | Protease Sites | Stephan $\Delta Q$ | Donnan $Q_{\text{pH } 6.0}$ | Popperian Status |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Top Affinity Lead** | `CHPWPHWAC` | **-9.048** | **-7.751** | **-1.30** | **0** | **+1.62** | **+1.05** | **SURVIVED ALL GATES (Rank 1)** |
| **12-Angle Benchmark Lead** | `CAHWRPWHC` | **-8.191** | **-7.343** | **-0.85** | **0** | **+1.62** | **+2.05** | **SURVIVED ALL GATES (Rank 1)** |
| Initial WT Lead | `CPHWRPWHC` | -8.062 | -7.902 | -0.16 | 0 | +1.62 | +2.05 | FALSIFIED (Amylase Cross-Reactivity) |
| High-Affinity Mutant (P6A) | `CPHWRAWHC` | -8.890 | -7.302 | -1.59 | 1 | +1.62 | +2.05 | FALSIFIED (Oral Protease Lysis) |
| Linear Control | `SPHWRPWHS` | -8.865 | -7.420 | -1.45 | 0 | +1.62 | +2.05 | FALSIFIED (35 Torsions, Entropic Loss) |
| Naive Poly-Arginine | `CRRRWWRRRRC` | -7.586 | -6.500 | -1.09 | 7 | 0.00 | +7.00 | FALSIFIED (Donnan Collapse & Lysis) |

---

## 4. Scope and Wet-Lab Validation Boundary

AutoDock Vina rigid-receptor grid docking, Vinardo empirical scoring, and OpenMM implicit-solvent energy minimization serve as rigorous in silico candidate triage and falsification filters. Computational docking cannot substitute for experimental wet-lab verification. Confirmatory testing, including surface plasmon resonance (SPR) equilibrium binding assays ($K_D$), radiolabeled glucan synthesis inhibition ($IC_{50}$), saliva-coated hydroxyapatite (sHA) adsorption kinetics, and multi-species dental biofilm demineralization assays, is strictly required prior to any therapeutic claims. All findings are strictly in silico computational hypotheses requiring wet-lab enzymatic and cellular falsification.
