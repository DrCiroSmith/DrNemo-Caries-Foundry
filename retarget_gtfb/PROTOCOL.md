# Pre-registered protocol: docking the macrocycle series against an experimental GtfB structure

Written 2026-11-27, before any result from this package existed. Nothing below has
been edited after seeing an output. That is the point of writing it first: the
thresholds cannot be moved to suit the numbers once the numbers arrive.

---

## 1. The question

Does the lead cyclic peptide from the DrNemo campaign bind the active site of
*Streptococcus mutans* GtfB, as judged by grid docking, better than sequence
permutations of itself?

The existing campaign asked this question of PDB 3AIC. 3AIC is not GtfB. It is the
catalytic domain of GtfC, the *gtfC* product GTF-SI, UniProt P13470. The repository,
the code, the pose filenames and the published page all name GtfB, UniProt P08987.
That is the error this package exists to correct.

## 2. What is already established, and how firmly

### 2.1 The target identity is settled

3AIC is GtfC. Three independent lines agree: the UniProt accession for the deposited
sequence, the residue numbering of the catalytic triad used throughout the existing
code (Asp477, Glu515, Asp588, Trp517, which are GtfC numbers), and the residue range
246 to 1087 present in the prepared receptor file, which matches the published GtfC
catalytic-domain construct and not the GtfB one.

### 2.2 But the strongest form of the criticism does not survive scrutiny

It would be convenient for this package to say the old docking was run against the
wrong protein and is therefore worthless. A referee holding the primary literature
would not accept that, and neither should we.

Schormann and colleagues (Acta Crystallogr F 2023, 79(Pt 5):119-127, DOI
10.1107/S2053230X23003199, PMID 37158310) solved both catalytic domains and report
that GtfB-CD and GtfC-CD share over 90 percent sequence identity and superpose with
an r.m.s.d. of 0.33 A over 4794 of 6445 aligned atoms. The catalytic triad is
identical. Every acarbose contact maps one to one:

| GtfC (3AIC numbering) | GtfB (8FJC numbering) | role |
|---|---|---|
| Arg475 | Arg449 | subsite contact |
| Asp477 | Asp451 | nucleophile |
| Asn481 | Asn455 | subsite contact |
| Glu515 | Glu489 | acid/base catalyst |
| Trp517 | Trp491 | acceptor platform |
| Arg540 | Arg514 | subsite contact |
| His587 | His561 | transition-state stabiliser region |
| Asp588 | Asp562 | transition-state stabiliser |
| Asp909 | Asp883 | subsite contact |
| Gln960 | Gln934 | subsite contact |

So the honest statement is narrower: the old work docked into a pocket that is a
close homologue of the intended one, mislabelled throughout. Retargeting is still
warranted, on two grounds that do not depend on overstating the error. First,
labelling: a paper cannot claim GtfB on GtfC coordinates. Second, and substantively,
the same study reports that acarbose makes twelve hydrogen bonds in GtfB against
eight in GtfC at a 3.2 A cutoff, adopts different ring conformations in the two
pockets, and gains a GtfB-specific contact at Asn511. The two sites do discriminate
ligands differently, so the numbers can change even though the fold does not.

### 2.3 Five further defects in the existing pipeline, found while preparing this

These matter more than the target mix-up, because each one undermines the published
numbers independently of which protein was used.

1. **The null distribution is matched to the wrong molecule.** All 39
   `SCRAMBLED_NULL_CONTROL` entries are permutations of the wild-type sequence
   CPHWRPWHC, all with 86 heavy atoms. The lead is CAHWRPWHC, with 84. Vina scores
   correlate strongly with heavy-atom count, so comparing an 84-atom ligand with an
   86-atom null confounds sequence order with size. Recomputed from the ledger
   directly, 23 of the 39 null members score at least as well as the lead, giving an
   exact one-sided empirical p of 0.6000, against the 0.5897 the repository reports
   from a normal approximation. The conclusion (no signal) is unchanged. The stated
   test is not the test that was performed.
2. **Vina is called with no seed.** `biophysics_engine.run_docking_vina_and_vinardo`
   builds its command line without `--seed`, so the search is stochastic and the
   pipeline the README describes as deterministic is not. The spread this hides has
   never been measured.
3. **The macrocycle ring is rigid.** The committed pose files contain no CG or G
   glue pseudo-atoms, which is how Meeko marks a ring it has opened for sampling.
   Every score therefore rests on one arbitrary ETKDGv3 seed-42 backbone conformer
   per sequence, and different sequences receive conformers of differing quality.
   That is a confound easily large enough to dominate the 0.15 kcal/mol difference
   between lead and null that the campaign reports.
4. **Two mutually inconsistent linear controls exist.** The cyclization benefit is
   reported twice from two different molecules.
   `demo_data/results/falsification_benchmark.json` gives the linear control as
   SPHWRPWHS at -8.865 kcal/mol, a cyclization delta of -0.803, so cyclization helps.
   `data/exhaustive_disproof_report.json`, Angle 10, gives it as AHWRPWH at -9.25
   against the cyclic -8.191, a delta of +1.059, so the linear peptide binds better.
   Neither artefact references the other. Separately, SPHWRPWHS is a Cys-to-Ser double
   mutant, so it is not a clean cyclization control at all: serine's hydroxyl changes
   hydrogen-bond typing, which Vina's scoring function does use, so that comparison
   varies two things at once. The AHWRPWH comparison is the more honest of the two and
   it is the one that goes against the project's thesis.

   The clean control was never run. It is the lead's own sequence, CAHWRPWHC, with the
   Cys1-Cys9 disulfide reduced to free thiols: identical residues, identical heavy-atom
   count, one bond changed. This package does not carry it either, for a stated reason
   rather than an oversight. Vina's reported affinity includes a normalisation on the
   number of rotatable branches, and the reduced form has more of them than the
   macrocycle, so a score difference between the two would not be attributable to the
   constraint. Running it would produce a number that cannot answer the question it
   appears to answer. The honest position is that the cyclization benefit this project
   rests on has not been tested by any control in the repository, and grid docking is
   not the method that can test it.
5. **A gate's label is unconditional while its verdict is not.** In
   `falsification_engine.py`, `gate4_alanine_hotspot_scan` and
   `gate6_chitosan_donnan_stephan` write the literal string
   `"classification": "PASSED_IN_SILICO"` regardless of outcome. The `passed` and
   `lead_passed` booleans in those same objects are computed, and `all_passed` uses the
   booleans, so the aggregate verdict is sound. The human-readable label is not, and the
   label is what a reader sees. `data/exhaustive_disproof_report.json` shows the same
   pattern at Angle 2, where `classification` reads `PASSED_IN_SILICO` and the summary
   claims "reproducible localization in the catalytic cleft" on a 0.46 A centroid
   standard deviation, while the same object records a 9.24 A mean and 12.19 A maximum
   pairwise pose RMSD across five seeds at a score spread of only 0.371 kcal/mol. A
   conserved centroid with a scrambled pose is not reproducible localisation: the
   molecule is being rearranged inside a conserved envelope. Because no `--seed` was
   passed, those five numbers also cannot be regenerated by anyone, including the
   author.

## 3. Receptor

PDB **8FJC**: GtfB catalytic domain, residues 191 to 1051 of *S. mutans* UA159,
with acarbose bound, 2.50 A, space group P4(3)22, two chains in the asymmetric unit.
Chain A is used.

Chosen over the alternatives for stated reasons, not convenience:

- **8FJ9**, same form and resolution, is apo. No bound ligand means no native
  redocking calibration control, which is the single most informative control
  available here.
- **8FK4** is the acarbose complex in which the inhibitor fits the density better,
  but it is 3.25 A with eight copies in the asymmetric unit.
- **8FKL** is the best resolved at 1.48 A, but proteolysis removed domains B and IV,
  which flank the active-site cleft. Higher resolution of the wrong construct is
  worse, not better.
- **No predicted model is used.** A homology or AlphaFold model would make every
  pose-level statement uninterpretable, and it is unnecessary because experimental
  GtfB coordinates exist.

The coordinates are downloaded by the user, not reconstructed here. `validate_receptor.py`
then checks identity against published facts: the title must name GtfB, acarbose and
P4322 and must not name GtfC or 3AIC; the resolution is read from `REMARK 2` and
compared with the config; chain A must hold at least 700 residues inside the construct
range; and every one of the eleven catalytic and contact residues above must be the
right amino acid at the right GtfB number. That last check is the strongest: at those
numbers a GtfC structure holds different residues, so it cannot pass by accident.

## 4. Protocol

Identical for lead and for every control. Any departure for any single molecule
invalidates the comparison, so there is none.

- **Builder:** `biophysics_engine.build_macrocycle_3d`, imported unchanged from the
  repository root, so ligand construction is byte-identical to the existing campaign.
  `run_docking.py` verifies this by building CAHWRPWHC with both the imported
  function and its local copy and comparing the PDBQT hashes; a mismatch stops
  the run.
- **Receptor preparation:** chain A standard residues, then
  `obabel -xr -p 7.4 --partialcharge gasteiger`, the same flags
  `scripts/prepare_receptors.py` used for 3AIC.
- **The conserved Ca2+ is a tested condition, not a discarded heteroatom.** Stripping
  every heteroatom keeps the new numbers comparable with the old ones, but the GtfB
  catalytic domain carries a conserved structural calcium adjacent to the catalytic
  triad. Schormann et al. 2023 report it in the intact domain, give the GtfC
  coordinating residues as Glu431, Asp437, Asn481 and Asp959, and state that GtfC
  Asp437 corresponds to GtfB Asp411. Asn481, one of those residues, also hydrogen bonds
  acarbose in both enzymes, so the metal sits immediately beside the inhibitor subsites.
  The package therefore builds two receptors, one without the ion and one with it, and
  tier 6 docks the same ligands in the same box against both. Two mechanisms are at
  issue and only one of them is real: Vina's scoring function has no electrostatic term
  and does not read PDBQT partial charges, so a retained ion cannot attract a cationic
  ligand by Coulomb interaction and the common claim of an "artificial cationic sink" is
  not something a Vina calculation can produce. What a retained ion does change is that
  Vina types a metal as a hydrogen-bond donor, and what a removed ion leaves is open
  volume the steric terms can fill. Both effects are local and both are measurable, so
  they are measured. If the verdict differs between the two receptors, no conclusion
  from either stands alone.
- **Ligand termini are recorded, not altered.** `Chem.MolFromSequence`, which the
  imported builder calls, produces a free N-terminal amine and a free C-terminal
  carboxylate; the ring is closed by the Cys1-Cys9 disulfide, not head to tail.
  `solubility_shield.calc_net_charge_at_ph` models the same molecule as if both termini
  were consumed by cyclization. The molecule is left exactly as the builder makes it,
  because changing it would destroy comparability with the existing numbers, and the
  capping state is written into every per-ligand record instead so the disagreement
  cannot be inherited silently.
- **Box centre:** centroid of the heavy atoms of the bound acarbose in chain A. The
  centroid of the catalytic triad carboxylate oxygens is computed independently as a
  cross-check. If the two disagree by more than 4.0 A the run stops.
- **Boxes:** 28 A cube, for comparability with the existing numbers, and a tight box
  equal to the acarbose bounding box plus 8 A per axis. Box size is a tested factor,
  not an assumption.
- **Engine:** AutoDock Vina, version string checked for 1.2.7, exhaustiveness 16,
  explicit `--seed`, five docking seeds for the primary set and three for the nulls.
- **Scoring:** Vina as primary, Vinardo as rescore. Any Vinardo value computed after
  an autobox fallback is flagged and excluded from aggregates, because it was scored
  in a different box than the pose was docked in.
- **Controls:** acarbose redocking; null A, the 39 existing wild-type permutations;
  null B, 39 newly generated permutations of the lead's own residues with the
  cysteine termini fixed; off-target counter-screens on 1SMD, human salivary
  alpha-amylase, and 1LZ1, human lysozyme; and the calcium condition above, run as
  tier 6 on a 12-member subset of each null. Tier 6 is a sensitivity check on the
  receptor and not an independent significance test: with 12 null members the smallest
  attainable exact p is 0.077.

## 5. Decision rules, fixed in advance

Each is stated with what the outcome would mean, so no result can be reinterpreted
after it arrives.

| # | Test | Threshold declared now | If it fails |
|---|---|---|---|
| 1 | **Redocking calibration.** Symmetry-corrected RMSD of redocked acarbose against its crystallographic pose. | Under 2.0 A | Over 2.0 A means the protocol cannot reproduce a pose that is known to be correct in this box. No pose-level claim about any other ligand is then supportable, whatever the scores are. The existing campaign reports 4.51 A on 3AIC and proceeded anyway. This package will not. |
| 2 | **Primary null test.** Lead CAHWRPWHC against null B, composition matched. | Exact one-sided empirical p below 0.05 | At or above 0.05 the lead is indistinguishable from a random ordering of its own residues, and the campaign has no sequence-specific signal. That is the declared negative outcome and it is an acceptable result. |
| 3 | **Secondary null test.** Wild type CPHWRPWHC against null A. | Same | Same. Reported whatever it shows. |
| 4 | **Seed sensitivity.** Spread of Vina scores across docking seeds for a single molecule. | Any score difference claimed between molecules must exceed this spread | If seed spread exceeds the lead-minus-null difference, that difference is search noise and no ranking is reportable. |
| 5 | **Conformer sensitivity.** Spread across five ETKDGv3 conformer seeds with the ring held rigid. | Same | If conformer spread exceeds the lead-minus-null difference, the campaign measured conformer-generation luck. |
| 6 | **Pose convergence.** Heavy-atom RMSD between the top poses from different docking seeds. | Reported separately from score convergence, never conflated with it | Converged scores with divergent poses means the scoring function is reaching the same number by different arrangements, which is a known failure mode of large boxes and is not evidence of a binding mode. |
| 7 | **Mode gap.** Energy difference between pose modes 1 and 2. | Reported, no threshold | A gap under about 0.3 kcal/mol means the top pose is not meaningfully preferred over the runner-up. |
| 8 | **Box-size effect.** Score and rank change between the 28 A and tight boxes. | Reported, no threshold | A reordering of lead against controls between boxes means the original ranking was an artefact of box size. |
| 9 | **Off-target selectivity.** Delta delta G against 1SMD and 1LZ1. | Must exceed about 1 kcal/mol to be discussed at all | Vina's own standard error on its reference set is around 2.85 kcal/mol. Differences well inside that are not selectivity. |
| 10 | **Geometry integrity.** Cys1 SG to CysN SG distance in every built ligand. | 2.00 to 2.15 A | Outside that range the disulfide is not a disulfide and the molecule docked was not the molecule intended. |

### What a positive result would and would not establish

If tests 1 through 5 all pass, the correct claim is: *in this docking protocol, on
this structure, the lead scores better than permutations of its own sequence, with
a converged pose, beyond seed and conformer noise.* That is a hypothesis worth
testing at the bench. It is not evidence of inhibition. Docking scores are not
binding affinities, and this is a rigid-receptor, rigid-ring calculation on a
macrocycle, which is close to the worst case for the method. No claim about
activity, potency or selectivity follows from any outcome here.

### What a negative result would establish

Equally little about the molecule, and quite a lot about the campaign. A null result
would say the computational screen does not support its lead, which is a publishable
negative finding about the method and the screen, and is the honest endpoint if that
is what the numbers show.

## 6. Known limitations, stated now rather than discovered later

1. **Docking alone establishes neither activity nor inactivity.** Both directions
   are unsupported. A score is not an affinity.
2. **The receptor is rigid.** GH70 glucansucrases have a mobile cleft; a single
   crystal conformer cannot represent it.
3. **The macrocycle ring is rigid.** Meeko did not open it. Ring conformation is
   fixed by ETKDG before docking begins, and tier 4 measures how much that matters
   rather than fixing it. Fixing it properly would require rebuilding the ligand
   preparation with explicit macrocycle sampling, which would break comparability
   with the existing numbers and is a separate exercise.
4. **One receptor conformation, one protonation state.** pH 7.4 Gasteiger charges,
   though the relevant environment in dental plaque is acidic. The existing campaign
   has the same limitation.
5. **Crystallisation artefacts.** 8FJC was grown from 2 M ammonium sulfate. Sulfate
   may occupy the negatively charged subsites that an arginine-rich ligand would
   target. Heteroatoms are stripped for docking and any within 6 A of the box centre
   are reported, because a cationic ligand scoring well in a site a sulfate occupies
   in the crystal is a finding that needs stating.
6. **Vinardo rescoring is a rescore, not an independent method.** It shares the
   pose. Agreement between Vina and Vinardo is weaker evidence than it appears.
7. **No solvation, no entropy, no explicit water.** Standard for grid docking and
   a standing limitation of every number in this package and the previous one.
8. **Null B is one of many possible nulls.** It fixes composition and net charge and
   varies order. It does not control for, say, hydrophobic moment or the position of
   the tryptophans.
9. **n equals 39 bounds the achievable p-value.** With 39 null members the smallest
   attainable exact p is 1/40, which is 0.025. A single test can clear 0.05; it
   cannot support a strong claim of significance on its own.
10. **Grid docking is being asked to do something it was not built for.** AutoDock Vina
    was parameterised and validated on small drug-like ligands. These are 1.2 kDa
    macrocyclic nonapeptides with roughly 14 rotatable branches and, because Meeko did
    not open the ring, no ring sampling at all. The existing campaign's own data show
    what that costs: at Angle 2 in `data/exhaustive_disproof_report.json`, five seeds
    land poses 9.24 A apart on average and 12.19 A apart at worst while the scores vary
    by only 0.371 kcal/mol, which is a scoring surface too flat to discriminate between
    genuinely different placements; and at Angle 10 the linear heptapeptide AHWRPWH
    scores -9.25 against the cyclic -8.191, so on this scoring function the constraint
    the whole design rests on is a penalty. Neither observation is a reason to stop
    docking, but both are reasons that a score difference of a few tenths of a
    kilocalorie between two peptides cannot carry a claim about binding. The 2.0 A
    redocking gate in the decision rules exists for this reason and gates every
    pose-level claim here. If it fails, the honest output of this package is a
    methodological negative, not a weaker version of a positive.
11. **One number in the source paper is internally inconsistent, and this package
    follows the version that is independently checkable.** Schormann et al. 2023 state
    in their introduction that the GtfC active site consists of Asp437 as nucleophile,
    Glu431 as acid/base and Asp959 as transition-state stabiliser. Those are the same
    three residues the same paragraph gives as Ca2+-coordinating, and the paper's
    discussion, its section 3.2 and its Table 4 all give the triad instead as Asp477,
    Glu515 and Asp588. Table 4's residue-by-residue PISA comparison is the checkable
    one, so the validator uses it: GtfB Asp451, Glu489, Asp562 and Trp491, every pair
    at the same 26-residue offset. The introduction's sentence appears to reuse the
    calcium-site numbers by mistake. This is recorded because the package's residue
    checks would be wrong if the other reading were right, and because a referee may
    raise it.

## 7. Provenance

- Sequences are generated mechanically from `data/contenders_ledger.csv`, not
  transcribed. `sequences.json` records the SHA-256 of the ledger it was built from.
- Null B is drawn with a fixed seed from the enumerated set of distinct permutations
  of the lead's interior positions, with the lead itself excluded, and every member
  is asserted to match the lead's composition.
- The receptor file's SHA-256 is recorded at validation and compared on every later
  run. If it changes, the run stops rather than silently mixing two receptors into
  one result set.
- Every docking record carries its receptor id, box id, docking seed, conformer seed,
  torsion count, Vinardo scoring context and full mode list.
- No number in this package is copied from any previous artefact. `analyze.py`
  recomputes all of them from the raw records, and its statistics are checked by
  `analyze.py --self-test` against cases with known answers before any real analysis
  is allowed to run.
