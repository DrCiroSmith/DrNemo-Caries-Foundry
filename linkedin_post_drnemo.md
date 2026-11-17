DrNemo Caries Foundry - Targeted GtfB Macrocycle Biophysics

We rebuilt the DrNemo-Caries-Foundry portfolio around real, falsifiable biophysics (our empirical biophysical campaign), replacing the deprecated v1 visuals (OpenFold3, DiffDock, pLDDT/iPTM, "Caries Vaccine" framing) with experimental-structure-guided computational design.

The target is Streptococcus mutans glucosyltransferase-B (GtfB, RCSB 3AIC), the enzyme that synthesizes the water-insoluble mutan matrix cementing acid-producing biofilms to enamel. The molecule is CAHWRPWHC, a cyclic disulfide macrocycle: Cys-Ala-His-Trp-Arg-Pro-Trp-His-Cys, disulfide S-S 2.06 angstrom, Arg5-Pro6 steric shield, zero protease cleavage sites.

Results from the foundry ledger (103 contenders: 52 ranked plus N=39 charge-matched scrambled nulls):
- AutoDock Vina into GtfB (3AIC): -8.191 kcal/mol
- Vinardo empirical cross-check: -4.990 kcal/mol
- Human salivary alpha-amylase off-target (1SMD): -7.343 kcal/mol
- Host selectivity Delta-Delta-G: -0.85 kcal/mol
- OpenMM Amber14/GBn2 relaxation at 310.15 K: -1,234.68 kJ/mol
- Donnan membrane potential Q at pH 6: +2.046 (chitosan gate <= +4.0)
- Stephan acid-switch Delta-Q: +1.62 (activated in the cariogenic niche near pH 5.2)

Popperian falsification decided the outcome. The wild-type lead CPHWRPWHC was falsified on amylase cross-reactivity (Delta-Delta-G -0.16). The P6A affinity trap (-8.890 kcal/mol raw) was falsified by a newly exposed tryptic cleavage site. Only CAHWRPWHC cleared all six gates: affinity, off-target selectivity, protease stability, chitosan compatibility, acid-switch response, and disulfide restraint.

The docked complex renders in the live 3D studio (docs/index.html) and a LinkedIn demo was exported. Source, figures, and SQLite ledger at github.com/DrCiroSmith/DrNemo-Caries-Foundry.

#ComputationalBiophysics #DrugDiscovery #Biofilm #OpenScience #StructuralBiology

Wet-lab confirmation (SPR binding, saliva-coated hydroxyapatite adsorption, multi-species biofilm demineralization) remains the next gate before any therapeutic claim.
