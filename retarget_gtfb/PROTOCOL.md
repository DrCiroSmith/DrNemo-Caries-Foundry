# Pre-registered protocol: docking the macrocycle series against an experimental GtfB structure

Written 2026-11-27, before any result from this package existed. No docking run from
this package has been executed: `out/` is empty and no `results.jsonl` or
`analysis.json` exists anywhere in the repository. That is the point of writing it
first, and the thresholds cannot be moved to suit the numbers once the numbers arrive.

This document has been amended since its first commit. Every amendment is logged
below with its time, what it said before and what it says now. A pre-registered plan
that is quietly different from its first version is worth nothing, and the git history
does not spare anyone that, because nobody reads it.

---

## 0. Amendment log

All times UTC on 2026-11-27. **No result from this package existed at any of these
times**, which is verifiable: the first docking run has not been made, so no amendment
here could have been informed by an outcome. Each entry gives the commit, so the diff
can be read directly.

| Time | Commit | Before | After |
|---|---|---|---|
| 18:17 | `12bc883` | first version | the plan as originally registered |
| 19:08 | `0ea6857` | the structural Ca2+ was stripped with every other heteroatom, following the original campaign | Ca2+ retained and stripped are both docked, as tier 6, so the assumption is tested rather than inherited |
| 19:13 | `74f5088` | the cyclization control was left unnamed | SPHWRPWHS named as a Cys-to-Ser double mutant and therefore not a clean cyclization control, with the reason it is not run |
| 19:35 | `62f751e` | defect 3 asserted that Meeko left the macrocycle ring rigid, and tier 4 was framed as a workaround for that | that claim is **withdrawn**: 1183 of 1244 pose files carry glue pseudo-atoms, so the ring was sampled. Limitation restated as partial ring sampling, tier 4 restated as a starting-conformer sensitivity test. My detector had read the atom type with a two-character slice and truncated `CG0` to `CG` |
| 19:45 | `904f2d4` | section 4.1 counted eight duplicate filenames as five distinct results, and no power gate existed | count corrected to three duplicated contents over eight filenames, only two of them docked results; decision rule 12 and section 4.2 measurement 1 added, stating the detectable effect before any run |
| 19:50 | `03d03fd` | the origin of the duplicated results was recorded as undeterminable, and a guess attributed it to the resume guard | the repository's version history settles it as one result under three names; the **cause** is recorded as unexplained and is not attributed to anyone; the resume-guard guess is **withdrawn**, since the guard cannot write one file's bytes under another file's name |
| 20:00 | `28f8d07` | the gate was stated against a claimed advantage of 0.15 kcal/mol, and the conclusion read "outside this design's resolution at any sample size" | that 0.15 figure is **withdrawn**: it has no source anywhere in the repository and was mine, not the campaign's. Replaced by the ledger-derived 0.129 kcal/mol lead-against-parent margin, computed at run time rather than stored. The conclusion is narrowed to "however many control sequences are added", with the seed count that would move the boundary stated beside it |
| 20:06 | `cd51c2c` | the document had no amendment log, and effect sizes were given without naming which comparison they measured | section 0 added, logging every amendment; section 0.1 added, fixing contrast A (lead against parent, +0.129) and contrast B (lead against the null, 0.152 the wrong way); rule 12 names both and states the seed price |
| 20:25 | `42c5c21` | the unseeded finding could be read as covering the N=5 seed sweep, and "a larger null does not move the boundary" was true of the quantile but silent on what a larger null does relax | section 4.1b added: the sweep is explicitly seeded, defect 2 is scoped to the production engine, and seed 42 is shown non-independent on both its pose and its input conformer, so at most four of the five runs are independent. The n-dependence is quantified: the requirement relaxes from 2.203 to 1.645 null standard deviations, 1.035 to 0.773 kcal/mol, and no further |
| 20:35 | `81b5197` | 4.1b showed seed 42 non-independent but did not say what that costs the two claims built on the sweep | both quantified from the files: the score summary moves to mean -8.078 with population spread 0.4116 over the four real runs, while the pose claim strengthens slightly, the mean pairwise r.m.s.d. rising from 9.25 to 9.33 A over 102 real atoms when seed 42 is dropped. My **first two attempts at that r.m.s.d. are withdrawn**: one counted the two `G0` glue markers as atoms, the other discarded the two real `CG0` carbons, both because columns 78 and 79 truncate `CG0` to `CG`. The selection is now stated in the text and recomputed by two tests |
| 20:36 | `f8dc497` | three files stated a literal count of the provenance suite, which had gone stale twice | the count is gone from all three and a test bans the phrasing; `--check-naming` keeps its eighteen because that number is derived from the list it runs |
| 20:40 | `0b4db3e` | 4.1b guessed that its r.m.s.d. figures differed from the manuscript's "most likely in mode choice or superposition" | that **guess is withdrawn**. Sweeping the plausible selections reproduces the manuscript's 9.24 and 12.19 exactly on one of them: mode 1, hydrogens removed, the two `G0` glue markers kept. Those figures therefore count two pseudo-atoms as heavy atoms; the correct values are 9.04 and 11.80. Reported to the manuscript work and pinned by a test that recomputes both of its figures from the pose files |
| 20:48 | `c2a2bd1` | this package's own analyser carried the same defect it had just reported in the manuscript: `analyze.py` read the type as `line[76:79]`, truncating `CG0` to `CG`, and excluded only hydrogens, so the two `G0` markers entered every r.m.s.d. as element G | the type is read as the last whitespace token in one function, `atom_type_of`, and glue is excluded there rather than at three call sites. **No pre-registered threshold moves.** Rule 1 is acarbose, which carries no glue (59 records, 44 heavy atoms, verified on the file); rule 10 selects by atom name `SG`; `heavy_atoms` comes from RDKit before Meeko writes the pdbqt, 84 for the lead, matching the corrected reader exactly. Only the derived `pose_rmsd_across_seeds_A` changes, so no docking needs rerunning |

Five of these entries withdraw a claim of mine rather than refine one. That record is
kept deliberately: a protocol whose amendment log contains only improvements is a
protocol whose author was not checking. Two of the four are the same mistake twice,
a two-character slice of the atom-type field, which is why the atom selection is now
spelled out wherever it matters and checked by a test rather than described in prose.

This log does not carry a row for the commit that adds these rows, because that commit
changes nothing but the log itself. Every row that records a change to the plan names
the commit that made it.

---

## 0.1 The two contrasts, which are not interchangeable

Every effect size in this document belongs to one of two comparisons. Naming the
number without naming its comparison leaves the gate ambiguous however well sourced
the number is, so both are fixed here and referred to by name throughout.

**Contrast A, lead against parent.** CAHWRPWHC at -8.191 kcal/mol against CPHWRPWHC,
the sequence it was derived from, at -8.062. Margin **+0.129 kcal/mol** in the lead's
favour. This is the only positive margin in the existing data, and it is what the
design question is actually about: whether the P2A substitution bought anything.
Source: `data/contenders_ledger.csv`, column `vina_3aic_kcal`.

**Contrast B, lead against the scrambled null.** CAHWRPWHC at -8.191 against the mean
of its 39 scrambled controls at -8.3433. The lead is **0.152 kcal/mol worse** than the
mean of its own scrambles, and **23 of the 39 controls score at least as well**, which
is 59 percent of them. Exact one-sided permutation p = 0.6000. Source: the same ledger.

The two run in opposite directions and that is the central fact of the existing
result. Contrast A is the margin the power calculation is run against, because a
power calculation needs a positive effect to be powered for. Contrast B is what the
significance test actually reports, and it does not merely fail to reach significance:
it points the other way. Any sentence in the write-up that gives one of these figures
must name which contrast it is, and no sentence may use A's sign with B's test.

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
3. **Every score rests on a single starting conformer, and two of Vina's own
   applicability limits are exceeded.** An earlier revision of this document
   asserted here that Meeko had left the macrocycle ring rigid, on the strength of
   finding no glue pseudo-atoms in the committed poses. That was wrong, and the
   error was mine: the detector read the AutoDock atom type with a two-column slice,
   which truncates `CG0` to `CG` and matches nothing. Read correctly, 1183 of the
   1244 committed pose files carry a `CG0`/`G0` pair, which is Meeko's mark for a
   macrocycle opened at one bond for sampling. Every cyclic nine-mer pose is in that
   set. The 61 without are acarbose and the linear and eleven-residue ligands, none
   of which has a macrocycle to open. The ring was sampled. The claim is withdrawn
   and `retarget_gtfb/test_provenance.py` now asserts the opposite against the files
   themselves, so it cannot be reasserted by accident.

   What survives the correction is narrower and still worth stating. Opening a ring
   at one bond gives Vina the resulting torsions, not the ring's conformational
   freedom: bond lengths, bond angles and the closure reference geometry are all
   fixed by one ETKDGv3 embedding at seed 42 per sequence, and Vina varies none of
   them. So conformer luck is reduced rather than removed, which is what tier 4
   measures, and tier 4 is a starting-conformer sensitivity test rather than a
   workaround for a rigid ring. Separately, the eleven-residue ligands are written
   at `TORSDOF 35`, above the 0 to 32 range over which Vina's affinity normalisation
   was fitted, while the nine-mers at `TORSDOF 15` are inside it, so those two
   groups are not scored on the same footing. Neither point is large enough on its
   own to account for the 0.129 kcal/mol margin the lead holds over the parent it was
   derived from; section 4 gives the measurement that settles what is.
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

### 4.1 Naming and provenance of every result

Of the 1244 pose files in `data/poses/`, 1239 hold distinct content. Three contents
appear under more than one name, and they are not equivalent:

| Content | Names | What it is |
|---|---|---|
| `27f5671b` | 3 | A docked result at 3AIC, nine modes, mode 1 at -8.191 kcal/mol. This is the lead's published headline number. |
| `77827776` | 2 | A docked result at 1SMD, nine modes, mode 1 at -7.343 kcal/mol. |
| `56e034b2` | 3 | The input ligand conformer, with no `MODEL` records and no `REMARK VINA RESULT` lines. Not a docking result. |

So **two docked results appear under five filenames.** The third group is benign as
duplication, since one builder at one ETKDG seed gives one file, but it is not
inconsequential, for the reason in 4.1b.

### 4.1b The repository's N=5 seed sweep contains at most four independent runs

`scripts/run_exhaustive_disproof.py` angle 2 is titled "Multi-Seed Conformer and
Docking Reproducibility (N=5 seeds)" and its summary states "Across N=5 independent
conformer and docking seeds". The sweep is genuinely seeded, which matters and is
stated plainly here because defect 2 above is about a different code path: at line
127 it declares `seeds = [42, 101, 202, 303, 404]`, sets `params.randomSeed = s` on
ETKDGv3, and passes `--seed str(s)` to Vina on its own command line. Defect 2's
unseeded finding is scoped to `biophysics_engine.run_docking_vina_and_vinardo`, the
production pipeline, which has no `--seed` anywhere. The sweep's figures are not
affected by it.

Seed 42 is affected by something else. Both of its artefacts are byte-identical to the
campaign's own:

| Seed | Docked pose, coordinates over all 9 modes | Input conformer |
|---|---|---|
| 42 | identical to `contender_mutant_P2A_CAHWRPWHC_3aic.pdbqt` | identical to `contender_mutant_P2A_CAHWRPWHC.pdbqt` |
| 101, 202, 303, 404 | mutually distinct | distinct |

Verified by hashing the coordinate columns of every `MODEL` record: seed 42 and the
campaign pose both digest to the same value, and the other four do not match each
other or it. The input conformers duplicate in the same pattern, which is the
`56e034b2` group above.

The two halves of that table are not equally good evidence, and an earlier version of
this section treated them as though they were. **The conformer match carries no
weight.** Every path that builds this ligand embeds with seed 42 hardcoded or
parameterised to the same value: `biophysics_engine.build_macrocycle_3d` sets
`ETKDGv3` with `params.randomSeed = 42`, and the sweep sets `params.randomSeed = s` on
the same `ETKDGv3` with the same `MMFFOptimizeMolecule` call, so asking the sweep for
seed 42 rebuilds that conformer by construction. An identical input conformer is the
expected result, not a suspicious one. That narrowing came from the manuscript work,
which reached it independently, and it is correct.

One caveat on it, which the files cannot resolve. A third path,
`scripts/run_falsification_campaign.py` line 56, embeds with a plain
`EmbedMolecule(mol, randomSeed=42)`, which is the default method and not `ETKDGv3`,
followed by `MMFFOptimizeMolecule(maxIters=300)` rather than 500. That recipe should
not produce byte-identical coordinates to the other two. So the identity is evidence
that both files came from the `ETKDGv3` path rather than from the campaign script, and
nothing in either file records which script wrote it.

**The docked pose match is the evidence.** The campaign's docking ran through
`biophysics_engine.run_docking_vina_and_vinardo`, which passes no `--seed`, so Vina
drew its own seed. A nine-mode byte-identical match between an unseeded run and a
seed-42 run is not something determinism produces. Two readings remain and the files
distinguish neither: the file is one result standing under three names, or the whole
group came from a single seeded run of the disproof script and was relabelled, in
which case the determinism is real but all three names are still one run. In both
readings seed 42's arm is not an independent replicate. The repository's claim of five
independent conformer and docking seeds is therefore **not supported**: at most four
are independent, and any dispersion computed over all five is computed over a set
containing one non-independent member.

What that costs the repository's two angle 2 claims is not the same for each, and the
difference is worth stating because one of them survives intact.

The **score** summary moves. The published figure is mean -8.10 plus or minus 0.37
kcal/mol over five. On the four independent runs (-8.565, -7.837, -8.374, -7.535) it
is mean -8.078 with population spread 0.4116, the same convention. The sample
standard deviation of those four is 0.4753, which is the figure the manuscript's
variance work already uses, so that work is unaffected.

The **pose** claim does not move. All five docked files are mutually distinct, five
distinct coordinate digests, so the pairwise comparison is over five real poses even
though one of them was not produced by seed 42.

The atom selection has to be stated, because this ligand's files make it easy to get
wrong, and because this package's own analyser got it wrong in the same way until
`c2a2bd1`. The selection now used everywhere in this package is heavy atoms excluding
Meeko glue: the type is read as the last whitespace-separated token of the record, `G0`
through `G3` are dropped as markers rather than atoms, `CG0` through `CG3` are kept as
carbon, and `H`, `HD` and `HS` are dropped. A test bans reading that field by column at
all, and the pytest-free self-test checks the selection, because `RUN_ME.sh` runs the
self-test unconditionally and pytest only if it is installed. Mode 1 holds 104 `ATOM` records, but two of them are Meeko glue markers typed
`G0`, which are not atoms and are excluded here. Two further records are the
macrocycle closure carbons typed `CG0`, which are real carbons and are kept. Columns
78 and 79 alone truncate `CG0` to `CG`, and a selection written on that two-character
field will either count the glue markers as carbon or discard two real ones. The same
truncation is what produced the withdrawn rigid-ring claim recorded in the amendment
log, so it is called out rather than left implicit.

On that selection, mode 1, in place, with no superposition and no fitting: over all
102 real atoms the mean pairwise r.m.s.d. rises from 9.25 to 9.33 A when seed 42 is
dropped, across six pairs instead of ten, with the worst pair at 11.84 A between
seeds 101 and 404, both independent. Over the 84 heavy atoms it rises from 9.04 to
9.11 A, worst pair 11.80 A, the same two seeds. So the conclusion that the score
reproduces and the pose does not strengthens very slightly when the non-independent
member is removed, rather than depending on it.

The manuscript no longer publishes a glue-inclusive figure. As of its commit `f82bc3a`
it publishes the four-run 9.11 A mean, 11.80 A worst pair and 0.39 A centroid spread
over the 84 real heavy atoms, and names 9.34, 12.19 and 0.38 in its own text as what
the earlier figures did. Its four-run values are this package's to the digit, reached
independently on both sides.

The figures it published before that, 9.24 A mean and 12.19 A worst over five, are kept
here because the selection that produces them is worth naming and because the table
below is why this package's own reader was audited. It is mode 1 with the eighteen
hydrogens removed and the two `G0` glue markers kept, 86 records:

| Selection | n | Mean over five | Worst pair |
|---|---|---|---|
| all records | 104 | 9.41 | 12.16 |
| real atoms, glue excluded | 102 | 9.25 | 11.84 |
| **hydrogens removed, glue kept** | **86** | **9.24** | **12.19** |
| heavy atoms, glue excluded | 84 | 9.04 | 11.80 |

So the figures published before `f82bc3a` counted two pseudo-atoms as heavy atoms. The
glue markers sit at the macrocycle closure and move with it, which is why including
them inflates a distance rather than leaving it alone. The correct five-pose values are
9.04 A mean and 11.80 A worst; the four-pose values, which is what both sides now
publish, are 9.11 A and 11.80 A. The conclusion never changed, both being far above the
2.0 A convention, but the figure would have been 0.39 A too large with a pseudo-atom as
the cause. Recorded here because the finding came out of this package's own
atom-selection fix, and because the same table is why that fix was needed.

So the correction needed is to the label and to the score summary, not to the finding.

This is why the `56e034b2` group is not merely untidy. It is the evidence that the
seed sweep's first member is not what the sweep says it is. An earlier revision of this section counted all eight filenames as
covering five distinct results, which was wrong arithmetic on my part, and it also
treated the benign group as part of the defect. Corrected here, and reconciled with
the independent repository-wide count made in the manuscript work.

The defect that remains is the two docked groups. Nothing inside those files records
a seed, a command line, a run identifier or a time, so from the files alone there is
no way to tell whether each was one run saved several times or several runs that
collapsed onto one, and no way to tell the benign case from the real ones without
inspecting the contents by hand.

For these two groups the question was settled from outside the files, by the
repository's version history rather than by anything in the pose data. That history
shows the campaign path created and its content replaced on the same day, and both
audit paths created later in a single commit already holding the replaced content,
never having held anything else. So no independent docking output was ever written to
those two paths: it is one result standing under three names at 3AIC, and the same
pattern gives one result under two names at 1SMD.

The **cause** is not established and is not asserted here. The files carry no seed, no
command line, no run identifier and no timestamp, and the repository's only log file
is empty, so there is no record of how the duplication arose. The supportable
statement is that the duplication is unexplained and that its origin cannot be
attested to. Nothing here claims that anyone copied anything. An earlier guess of mine
that the resume guard accounted for it is withdrawn, because the guard cannot write
one file's bytes under a different file's name.

That the files themselves could not answer this is the whole point. Version history
happened to be available in this instance and will not be in general, which is why the
ambiguity is designed out below rather than audited for afterwards.

- **Every output path is derived, never chosen.** The stem of every file a run writes
  is `t<tier>` followed by the run's kind, ligand identifier, receptor identifier,
  protonation state, box, docking seed, conformer seed, exhaustiveness, a short
  content hash of the receptor file and a short content hash of the ligand file. No
  label is supplied by hand anywhere in the package. Two runs that differ in any one
  of those parameters cannot write to the same path, and two different receptor or
  ligand files cannot either, whatever they are called.
- **The protonation state is named rather than implied.** There is one state in this
  package and it is `as_built_no_pH_model_free_termini`, because the builder goes
  through `Chem.MolFromSequence` and `AddHs`, which applies no pH model and leaves the
  termini free. Naming it means a second state added later cannot silently collide
  with this one.
- **Resume keys on scientific identity, not on the tier label,** so a run requested by
  two tiers is docked once and reported in both, and re-running after an interruption
  cannot duplicate or skip work. The record carries both.
- **Every result record carries its own provenance:** the full Vina command line, the
  seed, exhaustiveness, box centre and size, the Vina version string, the MD5 and
  SHA-256 of both the receptor and the ligand file, a run identifier unique to that
  invocation, start and finish times in UTC, and the SHA-256 of the pose file written.
  The two files sent back are therefore self-describing; nothing has to be taken on
  trust from a filename.
- **An existing result is never overwritten silently.** A pose file sitting at a
  deterministic path that no recorded result claims is an orphan from an interrupted
  or edited run. The run refuses it and says so, and `--overwrite` is the only way
  past. Nothing in the analysis needs that flag.
- **The analysis fails loudly on the failure described above.** `provenance_audit`
  reports it and the analysis exits non-zero if one pose content appears under more
  than one distinct run, if one path was written with more than one content, or if a
  run identifier is reused. The report is written first, so a nine-hour run is never
  discarded over a check that can be read afterwards.
- **It is tested against the real failure, not only against a fixture.**
  `retarget_gtfb/test_provenance.py` takes the duplicated files out of `data/poses/`
  and shows the check catching them, asserts that each of the ten run parameters
  changes the output path, and asserts that a reused run identifier and a path holding
  two contents are both caught. Thirty-eight checks, all passing, run in under a
  second. They are a precondition of the docking rather than a tier of it.

### 4.2 Measurements made from the results, costing no docking runs

Four things are computed from data the run already produces. Each is here because a
claim in the existing campaign cannot be adjudicated without it.

1. **Minimum detectable effect.** An exact one-sided permutation test on *n* nulls
   returns p = (k+1)/(n+1), so *n* bounds the question before any score exists. On
   the existing 39-member null the smallest attainable p is 0.025, and to reach
   p < 0.05 the candidate must beat every null member, which means scoring at or
   below the best of them. On the ledger's own numbers that threshold is
   -9.378 kcal/mol, which is **1.035 kcal/mol better than the null mean** of
   -8.3433 kcal/mol. The margin actually at issue is **contrast A**, 0.129 kcal/mol,
   the lead against its parent. The floor is **8.02 times** that margin. The test this
   floor belongs to is **contrast B**, the lead against the scrambled null, where the
   lead is 0.152 kcal/mol on the wrong side and 23 of 39 controls score at least as
   well. So the gate reads: the only positive margin in the data is contrast A's
   0.129, the test that would have to detect it is contrast B's, and contrast B's
   floor is eight times the margin. **This is what gates the wording of the
   conclusion.**

   The margin is computed from the ledger at run time rather than carried here as a
   constant, because an earlier revision of this section asserted a claimed advantage
   of 0.15 kcal/mol, naming no contrast, and that figure has no source anywhere in the
   repository. It
   appears in no document, no result file and no script, and the only stated threshold
   in the project is the +0.70 kcal/mol alanine hotspot gate, which measures something
   else. The number was mine and it should not have been pre-registered as his.
   Withdrawn, and replaced by a quantity with a file behind it.
   Without it, a p value above 0.05 reads as "no effect" when what it means is "no
   power", and those are different claims. Reported for every null test, both
   pairings, both summary statistics.

   **What a larger null would not buy.** The obvious remedy is to enlarge the null,
   and it does not work here, which is worth stating before anyone spends a week on
   it. The n in (k+1)/(n+1) bounds the smallest p attainable, so a larger null does
   lift that bound. But the test itself asks whether the candidate falls in the most
   extreme 5 percent of its own permutation distribution, and that boundary is a
   property of the distribution's spread, not of its size. The existing null has a
   standard deviation of 0.470 kcal/mol, so the 5 percent boundary sits about
   **0.773 kcal/mol** better than the null mean and stays there however many
   permutations are drawn. A larger null estimates that boundary more precisely; it
   does not move it, though it does relax what is *attainable* and that is worth
   quantifying rather than waving at. At n = 39 the candidate must beat every control,
   which on this null means 2.203 standard deviations better than the mean, or
   1.035 kcal/mol. As n grows the required rank relaxes: the largest number of controls
   the candidate may lose to is the largest j with (j+1)/(n+1) < 0.05, so j = 0 at
   n = 39, 8 at n = 199 and 48 at n = 999, and the required quantile rises toward 0.05
   from below (0.0250, 0.0450, 0.0490). The requirement therefore falls from
   2.203 standard deviations toward 1.645 and stops there, which in kcal/mol is 1.035
   falling to 0.773. **A larger null buys 25 percent and then nothing.** The
   0.129 kcal/mol margin is **0.275 standard deviations** from
   the null mean, and against the null the lead falls on the wrong side entirely: its
   -8.191 kcal/mol is 0.152 kcal/mol **worse** than the null mean of -8.343, with
   59 percent of null members scoring at least as well. So the honest statement is
   not that the effect is small and the sample too small to see it. It is that an
   advantage of this size is not distinguishable from sequence-order noise **however
   many control sequences are added**.

   That phrasing is deliberately narrower than "at any sample size", which overreaches
   and a reviewer would catch. Adding controls does not move the boundary, but
   averaging *m* docking seeds per control shrinks the measured spread by about the
   square root of *m*, so measurement effort does move it. Quantified at the 0.129
   margin and the null's spread of 0.4697 kcal/mol: reaching the one-sided 5 percent
   boundary needs **36 seeds per control**, about 1,440 docking runs across 39
   controls and one candidate, and 80 percent power needs **82 seeds**, about 3,280
   runs. So the honest claim is that the effect is beyond reach at this level of
   replication, with the cost of buying it stated, rather than beyond reach in
   principle.

   Whether that cost is worth paying is a separate question and the answer is
   probably no: 0.129 kcal/mol is **4.5 percent** of the 2.85 kcal/mol standard error
   Trott and Olson report over their own benchmark set, and even the 1.035 kcal/mol
   needed to beat every control is **36.3 percent** of it. Resolving a margin of that
   size would be a claim about the scoring function rather than about the peptide,
   which is why this stays a method result. `analyze.py` prints this as
   `power_by_best.asymptotic_requirement_kcal` and
   `effect_in_null_terms.advantage_in_null_sd` on every null test, so the sentence in
   the paper is read off the output rather than reasoned about afresh.
2. **Score against ligand size.** Vina's reported affinity grows with heavy-atom
   count. Regressed over all 260 contenders in the existing ledger, which span 60 to
   118 heavy atoms, the slope is **-0.0234 kcal/mol per heavy atom**
   (Pearson r = -0.282, r-squared = 0.08). The existing null is built from
   permutations of the wild type at 86 heavy atoms while the lead it is compared
   against has 84, so the published comparison carries a 2-atom size difference worth
   about **-0.047 kcal/mol**, roughly 36 percent of contrast A's 0.129 kcal/mol. Null B
   fixes the pairing; the slope is reported so that the residual size term can be
   subtracted rather than assumed away. Where a subset spans fewer than five heavy
   atoms the slope is withheld rather than printed, because a slope fitted over a
   two-atom span looks authoritative and is not.
3. **Ligand efficiency beside the raw score.** Score divided by heavy-atom count, with
   the rank correlation between the two orderings and whether they pick the same
   winner. A ranking that reverses between them is a ranking by molecular weight.
4. **Vina against Vinardo, with what that is not.** Rank correlation and mean absolute
   difference. Vinardo is a reparameterisation of the same functional form, scored on
   the same pose by the same program, so agreement between them is one method agreeing
   with itself and is not independent confirmation. It is reported because
   disagreement would be informative, not because agreement is.

These feed a single `noise_floor` block that puts the claimed advantage beside the
docking-seed spread, the embedding-seed spread, the pose spread across seeds, Vina's
reference standard error and the minimum detectable effect, so that decision rules 4,
5 and 6 below can actually be adjudicated instead of discussed. A more negative score
is better, so an advantage is the null mean minus the candidate's score, and a
negative advantage means the candidate scored worse than its own null.

## 5. Decision rules, fixed in advance

Each is stated with what the outcome would mean, so no result can be reinterpreted
after it arrives.

| # | Test | Threshold declared now | If it fails |
|---|---|---|---|
| 1 | **Redocking calibration.** Symmetry-corrected RMSD of redocked acarbose against its crystallographic pose. | Under 2.0 A | Over 2.0 A means the protocol cannot reproduce a pose that is known to be correct in this box. No pose-level claim about any other ligand is then supportable, whatever the scores are. The existing campaign reports 4.51 A on 3AIC and proceeded anyway. This package will not. |
| 2 | **Primary null test.** Lead CAHWRPWHC against null B, composition matched. | Exact one-sided empirical p below 0.05 | At or above 0.05 the lead is indistinguishable from a random ordering of its own residues, and the campaign has no sequence-specific signal. That is the declared negative outcome and it is an acceptable result. |
| 3 | **Secondary null test.** Wild type CPHWRPWHC against null A. | Same | Same. Reported whatever it shows. |
| 4 | **Seed sensitivity.** Spread of Vina scores across docking seeds for a single molecule. | Any score difference claimed between molecules must exceed this spread | If seed spread exceeds the lead-minus-null difference, that difference is search noise and no ranking is reportable. |
| 5 | **Starting-conformer sensitivity.** Spread across five ETKDGv3 embedding seeds. Meeko opens the macrocycle at one bond, so Vina samples the resulting torsions, but bond lengths, bond angles and the ring-closure geometry stay as the embedding set them. | Same | If conformer spread exceeds the lead-minus-null difference, the campaign measured conformer-generation luck. |
| 6 | **Pose convergence.** Heavy-atom RMSD between the top poses from different docking seeds. | Reported separately from score convergence, never conflated with it | Converged scores with divergent poses means the scoring function is reaching the same number by different arrangements, which is a known failure mode of large boxes and is not evidence of a binding mode. |
| 7 | **Mode gap.** Energy difference between pose modes 1 and 2. | Reported, no threshold | A gap under about 0.3 kcal/mol means the top pose is not meaningfully preferred over the runner-up. |
| 8 | **Box-size effect.** Score and rank change between the 28 A and tight boxes. | Reported, no threshold | A reordering of lead against controls between boxes means the original ranking was an artefact of box size. |
| 9 | **Off-target selectivity.** Delta delta G against 1SMD and 1LZ1. | Must exceed about 1 kcal/mol to be discussed at all | Vina's own standard error on its reference set is around 2.85 kcal/mol. Differences well inside that are not selectivity. |
| 10 | **Geometry integrity.** Cys1 SG to CysN SG distance in every built ligand. | 2.00 to 2.15 A | Outside that range the disulfide is not a disulfide and the molecule docked was not the molecule intended. |
| 11 | **Provenance integrity.** No pose content under more than one run, no path holding more than one content, no reused run identifier. | Zero of each | Any one of them means two results cannot be told apart, so no number in the set can be attributed. This gate runs before any other is read and the analysis exits non-zero when it fails. It is the gate the existing campaign has already failed: two docked results under five filenames, one of them the lead's headline -8.191. |
| 12 | **Power.** Minimum detectable effect from the null's own size and spread, against the margin being claimed, with both named by contrast. | Contrast A's margin must exceed contrast B's minimum detectable effect before any verdict is read. A figure stated without its contrast fails this rule on its own | Where the claimed difference is below the floor, the correct word is uninformative, not negative. A test that could not have found an effect of the claimed size has not ruled one out. On the existing 39-member null the floor is 1.035 kcal/mol against a sourced margin of 0.129, so this gate decides how the result is written up whichever way the p value falls. A second, harder floor is reported beside it: the 5 percent boundary of the null's own spread, 0.773 kcal/mol at a null standard deviation of 0.470, which no sample size moves. A margin below that boundary is to be reported as beyond reach at this level of replication, in the same breath as its price: at contrast A's 0.129 kcal/mol against a control spread of 0.4697, about **36 seeds per control** to reach the one-sided 5 percent boundary, roughly **1,440 docking runs**, and **82 seeds** for 80 percent power, roughly **3,280 runs**. Not "a small effect awaiting more data", and not "beyond reach in principle", because adding controls cannot move the boundary but averaging seeds can. The margin is computed from the ledger, never carried as a constant. |

### What a positive result would and would not establish

If tests 1 through 5 all pass, the correct claim is: *in this docking protocol, on
this structure, the lead scores better than permutations of its own sequence, with
a converged pose, beyond seed and conformer noise.* That is a hypothesis worth
testing at the bench. It is not evidence of inhibition. Docking scores are not
binding affinities, and this is a rigid-receptor calculation on a macrocycle whose
ring is opened at one bond from a single embedded conformer, which is close to the
worst case for the method. No claim about
activity, potency or selectivity follows from any outcome here.

### What a negative result would establish

Equally little about the molecule, and quite a lot about the campaign. A null result
would say the computational screen does not support its lead, which is a publishable
negative finding about the method and the screen, and is the honest endpoint if that
is what the numbers show.

### How a difference below the floor is to be reported

If the primary test returns p above 0.05 while the claimed difference sits below the
minimum detectable effect, the supportable sentence is that **the screen does not
discriminate the lead from matched permutations of its own residues, and is not
powered to, so it supports neither the molecule nor its absence.** It is not "the
peptide does not bind", which docking cannot establish, and it is not "the peptide
binds weakly", which the scores do not support either. A larger null would be needed
to make the test informative, and the number required follows from the floor rather
than from preference.

## 6. Known limitations, stated now rather than discovered later

1. **Docking alone establishes neither activity nor inactivity.** Both directions
   are unsupported. A score is not an affinity.
2. **The receptor is rigid.** GH70 glucansucrases have a mobile cleft; a single
   crystal conformer cannot represent it.
3. **Ring sampling is partial, not absent.** Meeko opens the macrocycle at one
   bond, so Vina samples the torsions that opening exposes. It does not sample bond
   lengths, bond angles or the ring-closure geometry, all of which are fixed by a
   single ETKDGv3 embedding per sequence. Tier 4 measures how much the choice of
   embedding matters; it does not remove the dependence. An earlier revision of this
   document stated that the ring was rigid and that no ring sampling occurred at
   all. That was wrong and is corrected in section 2.3, defect 3.
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
