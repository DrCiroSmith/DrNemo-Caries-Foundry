> ## Provenance notice: this memo is generated prose, not a measurement
>
> This file was written by a general-purpose language model. `analyzer.py` sends five
> summary numbers from `demo_data/jobs/completed/S_mutans_GtfB_Insoluble_variant_CRRRWWRRRRC.json`
> to a chat model (`OLLAMA_MODEL`, default `llama3`) and commits the reply. Nothing in this
> repository measured the interactions described below, and the model was not given any
> structure. The text is kept here as a record of what was published, with the measurements
> that contradict it.
>
> **Checked against `S_mutans_GtfB_Insoluble_variant_CRRRWWRRRRC.pdb`, the structure deposited
> beside this memo, and `3AIC_receptor.pdb`:**
>
> - The asserted salt bridge to **Glu515 is absent**. The closest guanidinium nitrogen to a
>   Glu515 carboxylate oxygen is **5.87 A**, and **0** nitrogen-oxygen pairs
>   lie within 4.0 A. Contacts with **Asp477 are real**: 6 pairs within 4.0 A. So the
>   sentence is right about one carboxylate and wrong about the other.
> - The asserted **stacking between Trp517 and Tyr430 does not occur**. Those two rings are
>   **12.32 A** apart, too far to sandwich a third ring, and in this pose TRP5 is 20.8 A from the Trp517 ring and 21.3 A from the Tyr430 ring; TRP6 is 28.0 A from the Trp517 ring and 27.2 A from the Tyr430 ring.
> - Net cationic charge at pH 6.8 is **7.0** by this repository's own `calc_net_charge_at_ph`,
>   not the 6.0 stated below and in the job record.
>
> **The job record this memo was written from is simulated, and internally inconsistent:**
>
> - It names target `P08987`, which is GtfB (GTF-I), while the header below
>   names GTF-SI (gtfC, P13470). One memo, two different proteins.
> - It records `generation: 20`, but the 260-row vault in
>   `data/foundry_contenders.sqlite` holds generations 0 to 7 only.
> - It records a binding affinity of -10.63 kcal/mol, a value that
>   appears nowhere in that vault, whose best score is -9.048 kcal/mol.
> - Its pLDDT is 94.2 while the memo below says 94.5. The record and its own memo disagree.
>
> **No AlphaFold, DiffDock, biosecurity screen or peptide synthesis runs anywhere in this
> repository.** The pLDDT, iPTM, DiffDock RMSD, IGSC result and synthesis yield below have no
> source here and should not be cited.

# Clinical & Biophysical Synthesis Memo: S_mutans_GtfB_Insoluble_variant_CRRRWWRRRRC

**Lead Sequence:** `CRRRWWRRRRC` (Cyclic 11-mer)  
**Target:** *Streptococcus mutans* Glucosyltransferase-SI (GTF-SI, gtfC gene product, UniProt P13470)  
**Generation:** 20 | **Structural Confidence (pLDDT):** 94.5 | **iPTM:** 0.87  
**Composite Fitness Score:** 0.829  
**DiffDock Catalytic Pocket RMSD:** 1.35 Å  

---

### 1. Mechanism of Action & Enzymatic Pocket Blocking
GtfB is the primary enzymatic virulence factor of *S. mutans*, synthesizing insoluble $\alpha$-1,3-linked glucans (mutan) that anchor the cariogenic biofilm matrix to enamel hydroxyapatite. 

Candidate `CRRRWWRRRRC` presents a conformationally constrained cyclic disulfide wedge (`C...C`). The central `W-W` aromatic pair mimics the glucosyl donor transition state, stacking between Trp517 and Tyr430 in the catalytic cleft. The flanking tetrameric Arginine clusters (`RRR`) form stable bidentate salt-bridges with catalytic carboxylates Asp477 and Glu515, preventing sucrose hydrolysis and insoluble mutan polymerization.

### 2. Mucoadhesive Chitosan Nanoparticle Compatibility
- **Exclusion of I, K, D, E:** Confirmed zero Asp/Glu (no anionic repulsion/entrapment with the polycationic chitosan polymer) and zero Lysine/Isoleucine (eliminating premature cross-linking during TPP gelation).
- **Surface Charge Density:** Net cationic charge at salivary pH 6.8 is **+6.0**, ensuring active electrostatic retention within topically applied chitosan dental varnish.
- **Hydropathy Index (GRAVY):** **-2.15**, indicating high water solubility for oral formulation.

### 3. Biosecurity & Synthesis Feasibility
- **IGSC Biosecurity Screening:** Negative for toxin or virulence homologies.
- **Solid-Phase Peptide Synthesis (PepSy):** 11 amino acid coupling cycles + on-resin air oxidation disulfide bridge. Estimated yield >88% with standard Fmoc chemistry.
