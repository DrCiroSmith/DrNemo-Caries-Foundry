---
title: "A virtual screen that cannot separate its own lead from a composition-matched null: a negative result from a cyclic peptide campaign against *Streptococcus mutans* glucansucrase GTF-SI"
author: "Roi Sasha Castrellon Avilés, DDS"
date: "2026-11-27"
status: "DRAFT, not submitted. Not for circulation."
---

# Abstract

I ran a structure-based virtual screen for disulfide-cyclised peptide inhibitors of *Streptococcus mutans* glucansucrase GTF-SI (glycoside hydrolase family 70, PDB 3AIC), evaluated 221 designed macrocycles against 39 sequence-scrambled controls, and report it as a negative result. The screen does not support the lead it produced, and the reason is the method rather than the molecule.

Three findings carry the conclusion. First, the setup fails its own redocking calibration: the crystallographic acarbose redocks to a symmetry-corrected heavy-atom RMSD of 4.507 Å against a 2.0 Å standard set in advance. Second, the score the screen ranks on is not reproducible at the scale of the differences it ranks by. Four independent rebuild-and-redock runs of the lead span 1.03 kcal/mol, with a sample standard deviation of 0.475 kcal/mol and a mean pairwise pose RMSD of 9.11 Å, against a 0.129 kcal/mol gap between the lead and its parent and a 0.47 kcal/mol standard deviation across the whole null distribution. The production pipeline passes Vina no seed, so none of the published scores can be regenerated. Third, the null comparison does not survive inspection: all 39 scrambled controls are permutations of the parent peptide and none is composition-matched to the lead; the lead scores worse than 23 of 39 of them (one-sided empirical p = 0.600); and in the only composition-matched test the data permit, the parent scores worse than 26 of its own 39 permutations (p = 0.675). The null's spread across 39 different sequences, 0.470 kcal/mol, is of the same order as the spread of one single ligand docked four times, 0.475 kcal/mol, so the null's dispersion is not demonstrably sequence-driven; and because the rejection threshold of the permutation test is set by that dispersion rather than by the number of controls, enlarging the control set cannot recover the comparison.

I report two further failures of provenance. Three docked complexes written under three different nominal docking protocols, including the lead's own published pose, hold one identical set of coordinates across all nine retained modes. How that happened is not recorded in any artefact and the author does not recall it, so the duplication is unexplained and I report it as one result carried under three names rather than three runs. The lead's reported affinity therefore cannot be attributed to a stated search protocol, the project's own reproducibility test has four real runs rather than five, and its reported pH effect of 0.485 kcal/mol is the difference between a reused file and a fresh run, which is 1.02 times the single-ligand noise.

I also report that the five deposited complexes for the arginine-rich variant series are not sterically possible in the intact enzyme. Each overlaps 20 to 31 receptor residues between positions 860 and 964 at closest approaches of 0.23 to 0.38 Å, and none of the overlapped residues is present in the 311-residue receptor slice the files contain, so the deposited coordinates cannot reveal the overlap. Structural interpretations previously attached to those five poses are therefore withdrawn.

The screen does not establish that these macrocycles bind GTF-SI, and it does not establish that they do not. It establishes that this pipeline, as built, cannot tell the two apart. I report the whole of it, including the parts that were wrong, because the discriminating failure is reusable and the lead is not.

# 1. Introduction

*Streptococcus mutans* is a principal agent of dental caries, and its glucansucrases are central to that role: they polymerise sucrose into the glucans that form the adherent matrix of cariogenic biofilm. The enzyme GTF-SI synthesises both soluble and insoluble glucan and belongs to glycoside hydrolase family 70. Ito and colleagues solved its crystal structure in the free form and in complex with acarbose and maltose, confirmed that its domain order is circularly permuted relative to the GH13 alpha-amylases, and identified the acceptor subsites: Trp517 provides the platform for glycosyl acceptor binding, while Tyr430, Asn481 and Ser589 comprise subsite +1 [1]. The authors closed by noting that the structure should be useful in designing inhibitors that prevent biofilm formation by GTF-SI, which is the premise this work set out to test.

Disulfide-cyclised peptides are an attractive scaffold for an oral-cavity target. They are conformationally constrained, they can be made to carry a net positive charge that increases as pH falls during a glycolytic acid challenge, and the cyclisation offers some protection against the proteases of saliva and plaque.

I built an automated pipeline that generates C-X(n)-C macrocycles, docks them into GTF-SI, scores them with two empirical functions, filters them through oral-stability gates, and ranks the survivors. The pipeline produced a lead. This paper reports why that lead should not be believed, and what the failure shows about the design of screens of this kind.

Negative computational results are published rarely, and in silico screens that fail their own controls are published almost never. That asymmetry is costly, because a screen whose discrimination has not been demonstrated looks identical, in a figure, to one whose has. The specific diagnostics reported here cost nothing to run and would have stopped this campaign at the outset.

# 2. Methods

## 2.1 Target preparation

The receptor was the catalytic domain of GTF-SI from PDB entry 3AIC, the structure reported by Ito and colleagues [1], prepared as a single chain of 842 residues numbered 246 to 1087 and 6,643 atoms. Residues 477, 515 and 588 were confirmed by direct inspection of the committed coordinates to be aspartate, glutamate and aspartate respectively, consistent with the GH70 catalytic triad: Asp477 as nucleophile, Glu515 as general acid/base, Asp588 as transition-state stabiliser.

Two off-target receptors were prepared for selectivity tests: human salivary alpha-amylase (PDB 1SMD, 495 residues, 3,938 atoms) and human lysozyme (PDB 1LZ1, 130 residues, 1,029 atoms).

**Correction to the target designation.** The pipeline source, its output records and its figures describe the receptor as GtfB and as UniProt P08987. That designation is wrong. The *S. mutans* glucansucrase gene products map as GtfB to GTF-I (P08987), GtfC to GTF-SI (P13470) and GtfD to GTF-S (P49331). PDB 3AIC is the GtfC product, GTF-SI, and the residue numbers used throughout the pipeline and throughout this paper (Asp477, Glu515, Asp588, Tyr430, Trp517) are GTF-SI numbering. The docking was therefore run against the structure the numbering refers to; the error is in the name attached to it, in the code, the text and the figures, and it has been corrected in the repository. I state this plainly because a misnamed target in a published figure is not a cosmetic defect: a reader comparing this work to GtfB literature would be comparing it to the wrong enzyme. A separate line of work is re-running the campaign against GtfB proper; see Section 7.

## 2.2 Ligand generation

Each sequence was built from its one-letter code with RDKit, the terminal cysteine thiols were bonded to close the macrocycle, and the structure was embedded with ETKDGv3 at a fixed random seed and minimised with MMFF94. A single conformer per sequence was carried forward. The molecule was converted to PDBQT with Meeko. The resulting ligands carry 15 active torsions for a 9-mer and 35 for an 11-mer.

Meeko's macrocycle treatment, which opens the ring and closes it with glue pseudo-atoms so that ring torsions are sampled during docking, was not applied uniformly. Counting glue pseudo-atoms across the 1,218 committed peptide pose files, all 1,172 nine-mer poses and all 7 ten-mer poses carry them, and none of the 39 eleven-mer poses does. The eleven-mer rings were therefore docked rigid, at whatever ring geometry the single ETKDG conformer happened to adopt. The disulfide is present in both cases: terminal SG to SG distances are 2.05 to 2.07 Å throughout. This asymmetry is taken up in Section 3.4.

## 2.3 Docking and scoring

Docking used AutoDock Vina [2,3] in a 28.0 Å cubic box centred on the GTF-SI active site, with nine poses retained per run. The best-scoring pose was rescored with the Vinardo function [4]. No random seed was passed to Vina; this is addressed in Section 3.2.

Exhaustiveness was not uniform across the campaign, and I did not notice this until I read the campaign script line by line. Reading the values the code actually passes: 16 for the acarbose redocking calibration, for the linearised control, for the parent peptide, for the naive eleven-mer, and for all seven alanine-scan mutants including the sequence later nominated as the lead; 8 for the 39 scrambled null controls and as the default for unlabelled candidates; and 4 for the generation 2 and generation 3 candidates. The function that docks the scrambled controls carries a docstring stating exhaustiveness 4 while passing 8, so the docstring is wrong and the value that ran is 8. Section 3.4 reports what this does to the comparisons.

## 2.4 Oral-stability gates and ranking

Candidates passed through gates for protease susceptibility (trypsin, chymotrypsin and gingipain recognition sites), disulfide geometry, and net charge at pH 7.2, 6.0 and 5.2 computed by Henderson-Hasselbalch. Survivors were ranked on Vina affinity against GTF-SI and on selectivity against alpha-amylase.

## 2.5 Null controls

Thirty-nine sequence-scrambled peptides were docked as a composition-matched null. Section 3.3 reports that the composition matching does not hold for the sequence the screen nominated as its lead, and Section 3.4 reports that the docking conditions were not identical either.

## 2.6 Provenance

Every numerical value in this paper is emitted by `scripts/compute_paper_numbers.py` in the accompanying repository, which recomputes it from the artefact that produced it rather than reading it from any summary. At the version accompanying this draft it emits 389 values, of which 348 are recomputed from primary data and 41 are read from recorded runs and labelled as such, alongside 25 failed cross-checks between artefacts; the totals move as checks are added, so the script itself is authoritative and these figures are a snapshot. It also emits, by design, a list of values that still disagree with each other; those disagreements are reported in this paper rather than reconciled silently.

# 3. Results

## 3.1 The setup fails its own redocking calibration

Before any candidate result can be read, the setup must reproduce a known answer. I redocked the crystallographic acarbose from 3AIC into its own site under the production protocol.

The redocked pose scores -7.915 kcal/mol and sits 4.507 Å from the crystallographic ligand by symmetry-corrected heavy-atom RMSD, against the 2.0 Å standard set in advance. The calibration fails.

The failure is specific, and the specificity matters. Across the nine retained modes, the pose centroid lies 1.34 to 4.10 Å from the crystallographic centroid, and the best mode places a ligand atom within 2.0 Å of 33 of the 44 crystallographic heavy atoms. The search finds the right cleft and recovers most of the crystallographic envelope. What it does not do is place the atoms in the right register. Acarbose is a near-symmetric oligosaccharide, and a shift of one sugar unit along the cleft produces exactly this signature: good envelope overlap, bad RMSD.

I therefore do not claim that this setup cannot find the GTF-SI active site. The supportable claim is narrower and still disqualifying: the setup does not reproduce a known binding mode at the resolution needed to interpret a predicted one.

## 3.2 The ranking is finer than the search that produced it

The campaign ranks candidates on a single number, the best-mode Vina affinity. I asked how far that number moves when nothing about the ligand changes.

**Across independent runs.** This is the direct test, and it is the one the conclusion rests on. The project's audit battery nominally rebuilt and redocked the lead five times, each run taking its own conformer-embedding seed and its own Vina seed at exhaustiveness 8. Four of those five are genuine independent runs. The fifth is byte-identical to an earlier pose file and is not a run at all, for the reasons set out in Section 3.10, so I exclude it and report the four.

The four independent scores are -8.565, -8.374, -7.837 and -7.535 kcal/mol: a mean of -8.078, a span of 1.03 kcal/mol, and a mean pairwise pose RMSD of 9.11 Å over the six pairs, with a worst pair at 11.80 Å and a best pair at 4.45 Å, against a 2.0 Å agreement convention. The score is loosely reproducible and the pose is not reproducible at all. Because each run reseeded both the conformer and the search, this measures the reproducibility of the whole ligand-to-score path rather than of the search alone, which is the quantity that matters for a leaderboard built that way.

Two conventions need stating, because my figures differ from the project's. I report dispersion as a sample standard deviation, 0.475 kcal/mol over the four runs; the project's record gives 0.371, which is the population formula applied to five values including the duplicate. I compute pose RMSD over the 84 real heavy atoms; the project's 9.24 Å and 12.19 Å count Meeko's two glue pseudo-atoms as atoms, and I reproduce both figures to the digit when I include those atoms and all five poses. Neither difference changes the conclusion: a nonapeptide roughly 9 Å across does not return to the same binding mode.

The comparison that matters follows immediately. The gap between the lead and its parent, the difference on which the lead was selected, is 0.129 kcal/mol. The run-to-run standard deviation of a single ligand's score is 0.475 kcal/mol, 3.7 times that gap. The standard deviation of the entire null distribution is 0.47 kcal/mol, which matches the run-to-run dispersion of this one ligand to within one per cent. Section 3.3b sets out what that coincidence does and does not establish, since the two figures are not measurements of the same quantity. A ranking whose increments are smaller than the irreproducibility of its own measurement is not a ranking.

Those four runs were seeded and are reproducible. The production pipeline is not: it invokes Vina with no seed argument, so none of the scores in the leaderboard, including the lead's, can be regenerated even in principle. The irreproducibility measured above is therefore a property of exactly the procedure that produced the published numbers, measured on the published lead.

**Within one run.** Vina retains nine poses per run. Across 107 docking runs in the committed pose set, the spread between the best and worst retained pose averages 0.704 kcal/mol (median 0.655, range 0.192 to 1.532), which exceeds the 0.129 kcal/mol selection gap in 107 of 107 runs and the 0.47 kcal/mol null standard deviation in 91 of them. For the lead and its parent, the first two modes differ by 0.101 kcal/mol while lying 5.58 Å apart in heavy-atom RMSD: two conformations the scoring function cannot separate that are not the same binding mode.

I report this second test with an objection stated, because it is weaker than it looks. Vina's retained modes are selected to be structurally diverse, so they are not repeated measurements of one quantity, and the spread among them is not a standard error on the best-mode estimate. A reader is entitled to say that a diverse pose set should span a range of scores and that this says nothing about precision. I accept that, which is why the between-run result above carries the argument. What the within-run spread does show is narrower and still relevant: at the scale of the differences this campaign ranks on, the scoring function does not distinguish binding modes that differ by 5.58 Å, so a rank order built from best-mode scores is not tracking a single consistent mode of interaction.

The ligands carry 15 active torsions as nine-mers and 35 as eleven-mers. I do not assert any threshold from the literature for the exhaustiveness such a ligand requires, because I could not find one I was willing to quote. The internal evidence is sufficient on its own: repeated runs of one ligand do not converge, so the quantity used to rank candidates is not resolved at the scale of the differences being ranked.

## 3.3 No discrimination against a composition-matched null

The published control was 39 sequence-scrambled peptides. I recomputed the comparison from the primary records.

All 39 scrambles are permutations of the parent peptide CPHWRPWHC (86 heavy atoms). The sequence the screen nominated as its lead is CAHWRPWHC (84 heavy atoms), a proline-to-alanine substitution. Zero of the 39 controls are composition-matched to the lead. The published null test therefore compares the lead to permutations of a different molecule, and a comparison of that kind cannot separate sequence order from composition, which is the only thing a scrambled null exists to do.

Two comparisons are available, and I report both, with the estimator stated for each.

For the lead against the parent-derived null: Vina affinity -8.191 kcal/mol, null mean -8.343, null standard deviation 0.47, Z = +0.324. Twenty-three of 39 controls score at least as well as the lead. As a proportion that is k/n = 0.5897; as a one-sided empirical p value with the observation counted among the null, which is the correct estimator for a permutation test, it is (k+1)/(n+1) = 0.600. The value 0.5897 appeared in earlier project documentation and should not be used.

For the only composition-matched comparison that the data permit, the parent against its own 39 permutations: Vina affinity -8.062 kcal/mol, Z = +0.599, 26 of 39 controls at least as good, one-sided empirical p = (26+1)/(39+1) = 0.675.

Both Z scores are positive. In Vina's convention a more negative score is better, so both the lead and its parent score worse than the mean of their scrambled controls. Neither comparison approaches significance, and the matched comparison is further from it than the unmatched one.

I want to be exact about what this does and does not show. Read together with Sections 3.2 and 3.4, these p values do not demonstrate that the peptides fail to bind. They demonstrate that this screen's discrimination between designed sequence and scrambled sequence has not been established in either direction. A screen that cannot separate its lead from a matched null has not yet produced evidence about any of its candidates.

## 3.3b The size of difference this design could resolve

The natural reply to a null result at n = 39 is to enlarge the null. It is worth setting out why that would not rescue this comparison, because the arithmetic is short and a reader can redo all of it from the ledger.

With one score per control and the observation counted among the null, the one-sided permutation p value is (k+1)/(n+1), where k is the number of controls scoring at least as well as the candidate. At n = 39 the smallest attainable value is 1/40 = 0.025, so p < 0.05 requires k = 0. The candidate must beat every control, which means scoring at or below the best control at -9.378 kcal/mol, a margin of 1.035 kcal/mol over the null mean.

Against that threshold there are two quantities to compare, and they measure different things, so I state both.

The first is the lead against the parent it was derived from: -8.191 kcal/mol against -8.062, so the proline-to-alanine substitution bought 0.129 kcal/mol. That is the only positive margin in the campaign, and it is the one the design decision rested on.

The second is the lead against the null, which is the contrast the screen's conclusion actually rested on: -8.191 kcal/mol against a null mean of -8.343, so the lead scores 0.152 kcal/mol worse than the mean of its own scrambled controls, with 23 of 39 controls at least as good (Section 3.3). The sign of that second quantity is the finding. There is no advantage over the null to size.

The power arithmetic below therefore runs on 0.129 kcal/mol, because it is the only positive margin available to test, and the result applies with more force to a contrast whose sign is already wrong. No document in this repository states a target advantage, and I looked; the only stated threshold is the +0.70 kcal/mol alanine hotspot gate, which measures something else. The 0.129 kcal/mol is one eighth of the margin a win would require.

Enlarging the null lowers the floor on the attainable p value and relaxes the rank a candidate must achieve, but it cannot lower the score it must reach. With n controls at one score each, a candidate may lose to at most j of them, where j is the largest integer satisfying (j+1)/(n+1) < 0.05. At n = 39 that gives j = 0, so it must beat all 39. At n = 199 it may lose to 8, and at n = 999 to 48. Expressed as a fraction of n, the rank it must reach is 1/39 = 0.026, then 0.045, then 0.049: rising toward 0.05 but below it at every n. The requirement is therefore a quantile of the null, and the 5 per cent quantile is a floor beneath it that no control count reaches.

On a normal approximation that floor lies 1.645 x 0.470 = 0.773 kcal/mol better than the null mean. At n = 39 the actual requirement is stricter, 1.035 kcal/mol or 2.20 standard deviations, because beating the best of 39 asks more than beating the 5 per cent point. So enlarging the null would move the requirement from 2.20 standard deviations down toward 1.64 and no further. A margin of 0.129 kcal/mol is 0.27 null standard deviations, which falls short of even that asymptotic floor by a factor of six. This is the sense in which more controls cannot rescue the comparison: they estimate the boundary more precisely and relax the rank, while leaving the score to be beaten essentially where it is.

That raises the question of where the dispersion comes from. The null's sample standard deviation across 39 different sequences is 0.470 kcal/mol. The sample standard deviation of four independent runs of one single ligand, from Section 3.2, is 0.475 kcal/mol, a difference of 1.2 per cent. Those four are what remains of the published five once the duplicate is set aside (Section 3.10). All five were run with an explicit Vina seed, unlike the production pipeline.

The fifth value is set aside for a reason in the files, not a preference of mine, and the reason belongs here rather than only in Section 3.10, because five runs otherwise look like more data than four. All five seed runs wrote pose files of nine ranked modes each, and their best-mode scores are exactly the five values the record lists. Four of those files hold four distinct coordinate sets. The fifth, seed 42, hashes to the same coordinates as the campaign's own pose for the lead and as the pH 7.2 arm, across all nine modes: one geometry under three names, written by three stated protocols, two of which leave the Vina seed unset. Independent searches do not agree to the last decimal of 84 atom coordinates nine times over. Whether that file is the campaign pose under another name or a run whose output no longer exists cannot be established, and in neither case is it an independent replicate, so excluding it is what the artefacts support.

A reader who declines that argument and takes the record at its word loses nothing of the conclusion, because the comparison is given both ways. Read as five runs, the single-ligand standard deviation is 0.415 kcal/mol, 0.88 times the null's, and the implied between-sequence variance is +0.049 kcal squared per mol squared. Read as four, it is 0.475, 1.01 times the null's, and the implied variance is -0.005. The sign of that quantity turns on which reading is taken, which is itself reason enough not to report either sign as a result.

What does not depend on the choice is that neither reading resolves anything. An F test for equal variances gives F = 0.98 on 38 and 3 degrees of freedom, one-sided p = 0.61, under the four-run reading, and F = 1.28 on 38 and 4 degrees of freedom, p = 0.45, under the five-run reading. The 95 per cent confidence interval on the single-ligand standard deviation runs from 0.269 to 1.772 kcal/mol on four runs and from 0.248 to 1.192 on five, and both contain the null's 0.470; the interval on the null's own standard deviation, 0.384 to 0.605 kcal/mol, contains both single-ligand estimates in turn. Four or five observations cannot estimate a variance, and these data cannot distinguish a null whose dispersion is entirely run-to-run noise from one with a real sequence-dependent component.

So I will not write that the sequence-dependent part of the null's spread is zero, and an earlier version of this section did. A negative estimate of a variance component is the routine outcome when two estimates are close and one of them rests on three degrees of freedom, and it is not evidence that the component is nil.

One further objection runs the same way, and I state it rather than wait to be told. The record's own name for the sweep is multi-seed conformer and docking reproducibility, so each run reseeded the conformer as well as the search, while every null member is one conformer docked once. The single-ligand figure therefore contains a conformer-generation term that the null does not, which makes it an upper bound on the run-to-run noise inside the null rather than an estimate of it, and biases the implied sequence component downward. Separating the two terms needs a sweep that holds the conformer fixed and varies only the docking seed, and this work does not have one.

What survives all of that is narrower and still sufficient: the spread of one ligand is of the same order as the spread across 39 different sequences, so the null's dispersion is not demonstrably sequence-driven, and a reader cannot assume that the controls in Section 3.3 are measuring sequence at all.

This narrows the claim and makes it more useful. The limit is not a property of docking in general; it is a property of this protocol's one-score-per-control design. Averaging m independent seeds per control shrinks the measured spread by roughly the square root of m and moves the boundary toward the mean. For a margin of 0.129 kcal/mol merely to reach the one-sided 5 per cent boundary takes m = 36, which over 39 controls and one candidate is about 1,440 docking runs. For 80 per cent power at that margin it takes m = 82, about 3,280 runs. Both figures assume the run-to-run noise is independent between seeds and has the dispersion estimated above.

The last step matters most for how a result of this kind should be read. Even a margin of that size resolved to statistical satisfaction would be a statement about the scoring function and not about binding. The 0.129 kcal/mol is 4.5 per cent of the 2.85 kcal/mol standard error that Vina's authors report over their own benchmark set (Section 4.2b), and even the 1.035 kcal/mol needed to beat every control is 36 per cent of it.

So the supportable claim is this. A margin of the size this screen produced lies outside what this design can resolve, however many control sequences are added, and the remedy is a different measurement rather than more of the same one. I put it that way rather than saying it is unresolvable at any sample size, because control count and measurement effort are not the same lever: more controls cannot help, while more seeds per control could, at the cost set out above.

## 3.4 Candidate and control were not searched with equal effort

Having found the composition mismatch, I read the campaign script line by line rather than trusting its summary, and found a second asymmetry between candidate and control that is independent of the first.

The lead was docked at exhaustiveness 16. Its 39 scrambled controls were docked at exhaustiveness 8. The generation 2 and generation 3 candidates, which make up most of the leaderboard, were docked at exhaustiveness 4. So the lead received twice the search effort of its own null and four times that of most of the field it was ranked above.

I first wrote that this inflates the lead relative to its null, and then that the dataset let me measure the inflation and found it negligible. Both statements were wrong and I withdraw both. The measurement rested on treating the lead's campaign score as an exhaustiveness-16 observation and the audit battery's reseeded runs as exhaustiveness-8 observations of the same ligand. Section 3.10 shows that the lead's reported score cannot be attributed to a particular search protocol at all: three files carrying three different nominal protocols hold one identical set of coordinates, so there is no exhaustiveness-16 observation to compare against. No ligand in this dataset has a usable paired comparison across settings.

What survives is narrower and still a defect. Comparisons in this campaign cross groups that were searched with different effort, which means no comparison here was made under controlled conditions, and a reader cannot tell from the published output that this was so. Raising exhaustiveness cannot make Vina's expected best score worse, so the direction of the bias favours the lead over its null and the earlier generations over the later ones. I cannot put a size on it, and I am not going to assert one in either direction to make the argument tidier.

The second asymmetry is conformational. As reported in Section 2.2, the nine-mer rings were opened by Meeko and their ring torsions sampled during docking, while the eleven-mer rings were docked rigid. The naive eleven-mer CRRRWWRRRRC, used throughout this project as the unoptimised comparator that the designed nine-mers were shown to beat, is one of those rigid-ring ligands. Its score is a function of one arbitrary ring geometry that the search was never allowed to change.

Unlike the exhaustiveness asymmetry, I cannot bound this one from the data, because no eleven-mer was docked with a flexible ring and so no paired comparison exists. A rigid ring can only restrict the search, so the direction is that the naive eleven-mer's score is at best as good as it would be with the ring sampled, which favours the designed nine-mers it was used to make look good. I state the direction and leave the size unknown.

Neither asymmetry rescues the lead and neither is needed to sink it. Section 3.3 stands on the composition mismatch and the p values, and Section 3.2 on the irreproducibility, both of which are independent of how search effort was allocated.

## 3.5 The selection of the lead was not driven by affinity

The lead was nominated from a single-position alanine scan of the parent. Recomputed from the primary records, that scan shows no hotspot: no position reaches the +0.70 kcal/mol loss threshold set in advance, the largest loss is +0.257 kcal/mol (W7A) and four of seven substitutions improve the score rather than degrading it.

Within those four, the nominated lead P2A is the smallest improvement at -0.129 kcal/mol. The largest is P6A at -0.828 kcal/mol, which was not nominated because it introduces a trypsin recognition site the parent does not have. That is a defensible engineering trade, but it means the lead was selected by a protease gate and not by affinity, and the affinity differences involved are all smaller than the within-run spread reported in Section 3.2.

The lead is also not on the Pareto front of the declared objectives. Computed over all 260 docked sequences, the 221 designed candidates and the 39 scrambled controls together, three are non-dominated and the lead is not among them.

## 3.6 Selectivity is inconclusive, and was reported inconsistently

Selectivity was assessed as the difference in Vina affinity between GTF-SI and two off-targets. For the lead these differences are -0.848 kcal/mol against human salivary alpha-amylase and -1.521 kcal/mol against human lysozyme.

The pipeline compared these against an assumed Vina scoring error of 2.85 kcal/mol. I first recorded that constant as untraceable, and I withdraw that: it is the standard error Vina's authors report for predicted against experimental free energies of binding over their 190-complex test set, which they attribute in part to ligands with many active rotatable bonds [2]. The constant is therefore sourced, though it is a benchmark-wide affinity error and not a discrimination threshold for a selectivity margin, and using it as one is an assumption the project never stated. Taken at face value, the two differences are 0.30 and 0.53 of it. On any consistent reading, both are inconclusive.

They were not read consistently. The amylase comparison was judged by whether the margin cleared the assumed error, and the lysozyme comparison by whether the absolute affinity cleared a -7.5 kcal/mol cutoff. Because the tests differed in kind and not merely in threshold, the larger margin (-1.521) was recorded as a pass while the smaller one (-0.848) was recorded as inconclusive. Applying the single margin test to both moves the project's own falsification tally from two angles passed and five inconclusive to one and six. I report selectivity here as not established.

## 3.7 Pocket engagement of the three valid poses

Three deposited complexes are genuine docked poses in the receptor frame, verified clash-free against all 6,643 committed receptor atoms at closest approaches of 2.66 to 2.76 Å. For these I compared engagement of the catalytic triad with that of the co-crystallised acarbose, counting ligand heavy atoms within 4.0 Å of a carboxylate oxygen of Asp477, Glu515 or Asp588, normalised per 100 ligand heavy atoms.

| Ligand | Heavy atoms | Contacts at 477/515/588 | Per 100 atoms | Acarbose envelope overlap |
|---|---|---|---|---|
| Acarbose (crystal) | 44 | 6 / 7 / 7 | 45.5 | reference |
| CPHWRPWHC (parent) | 86 | 0 / 1 / 2 | 3.5 | 55 / 86 |
| CAHWRPWHC (lead) | 84 | 1 / 1 / 4 | 7.1 | 64 / 84 |
| CHPWPHWAC (best score) | 80 | 3 / 0 / 1 | 5.0 | 59 / 80 |

The macrocycles occupy the cleft. Between 64.0 and 76.2 per cent of their heavy atoms lie within 4.0 Å of the crystallographic acarbose envelope, and none clashes with the receptor. What they do not do is engage the catalytic carboxylates at anything like the density of the co-crystallised inhibitor: 3.5 to 7.1 contacts per 100 heavy atoms against 45.5, a 6- to 13-fold difference.

I had expected, and was prepared to argue, that an 80-atom macrocycle would not fit the GH70 cleft at all. That argument is wrong and I have dropped it rather than qualified it: the poses fit, they overlap the inhibitor envelope substantially, and they are sterically clean.

The remaining observation stands on its own but proves less than it appears to. A competitive inhibitor of a glucansucrase need not contact the catalytic triad; occluding the cleft can be sufficient, and I have not tested occlusion. Low triad contact density is consistent with weak inhibition and also consistent with a different and still productive binding mode. Given Sections 3.2 to 3.4, I cannot distinguish these, and I report this comparison as a description of the deposited poses rather than as evidence of inactivity.

## 3.8 Five deposited complexes are not sterically possible

The five arginine-rich variant complexes were the most prominent structural output of this project. They are not valid poses, and I withdraw them.

Each of the five files contains a contiguous 311-residue slice of the receptor, residues 355 to 665. Within that slice the ligand is clean, with closest approaches of 2.62 to 2.90 Å and no atom under 2.2 Å, which is what the files' own records state. Tested against the intact 842-residue receptor, each ligand overlaps 20 to 31 receptor residues lying between positions 860 and 964, with 76 to 94 ligand heavy atoms under 2.2 Å and closest approaches of 0.23 to 0.38 Å. An approach of 0.23 Å between heavy atoms is not a close contact; it is near-coincidence.

Not one of the overlapped residues is present in the deposited slice. The clash check was run against a domain slice rather than the whole protein, the slice happens to omit precisely the region the peptide occupies, and so the deposited coordinates cannot reveal the overlap. I have not been able to determine how the slice came to be chosen, and I make no claim about that.

Every structural interpretation previously attached to these five poses is withdrawn, including the per-variant triad contact counts, which measure a placement that cannot exist in the intact enzyme.

## 3.9 Mechanistic claims that are false against their own coordinates

Five analysis files accompanying those variants state, in identical wording, that the arginine cluster forms bidentate salt bridges with Asp477 and Glu515 and that the tryptophan pair stacks between Trp517 and Tyr430. Both assertions are measurable against the coordinates they are attached to, and both fail.

For Glu515, no guanidinium nitrogen in any of the five poses lies within 4.0 Å of a Glu515 carboxylate oxygen. The closest approaches are 5.87, 4.64, 4.65, 6.26 and 4.65 Å, and the count of qualifying pairs is zero in all five. Asp477 contacts are real in these coordinates (2 to 6 pairs, closest 2.85 to 3.45 Å), though Section 3.8 applies to them.

For the stacking claim, Trp517 and Tyr430 are genuine and functionally important GTF-SI residues: Ito and colleagues identify Trp517 as the glycosyl acceptor platform and place Tyr430 in subsite +1 [1]. The residue identities are confirmed in the committed coordinates. The geometry is not. The two ring centroids lie 12.32 Å apart, which is roughly three times too far to sandwich a third aromatic ring, and the nearest ligand ring centroid in any of the five poses is 12.58 Å away. No such stack exists or could exist in this structure.

These files were generated by a script that passes five numeric fields to a general-purpose language model and commits the prose it returns. The prose is literature-grounded, in that the residues it names are the right residues from the right paper; it is simply not a measurement of the structure it describes. The files have been relabelled in the repository to state their provenance and to list these contradictions. I report this because it is the most transferable failure in the project: generated text that is correct about the literature and false about the data is very hard to catch by reading, and trivial to catch by measuring.

## 3.10 Three nominally different docking protocols report one identical set of coordinates

I checked the pose library for duplicates by hashing the coordinate fields of every docked and ligand PDBQT file, excluding the per-mode splits that are derived mechanically from their parents. Of 272 files, 267 coordinate sets are distinct. The five files that duplicate an already-counted set sit in three groups, and every group involves the lead.

The first group holds three docked complexes against GTF-SI that are identical in every byte of all 96,786, not only in their coordinates: the campaign's own pose for the lead, the audit battery's seed-42 reproducibility run, and the audit battery's pH 7.2 arm. The three scripts that write these files specify three different protocols. The campaign docks at exhaustiveness 16 with no seed argument. The seed-42 run docks at exhaustiveness 8 with an explicit seed. The pH 7.2 arm docks at exhaustiveness 8 with no seed argument. Two of the three leave the seed unset, so two of the three are stochastic searches whose starting state was never fixed.

The second group is the same pattern at the second target: the audit battery's pH 7.2 arm against alpha-amylase, at -7.343 kcal/mol, is identical in every byte to the campaign's own alpha-amylase pose for the lead. The third group is the three input conformers, and that one is benign, for the reason given below. Stated as a filename census: eight filenames carry three coordinate sets, two of those sets are docked results, at -8.191 kcal/mol under three names at GTF-SI and -7.343 kcal/mol under two at alpha-amylase, and the third is an input conformer under three names. Counting the eight filenames as separate results, or counting the conformer group as part of the defect, overstates it in opposite directions, and I have done neither.

A stochastic global search with an unset seed does not reproduce another run's output to the last decimal of 84 atom coordinates across nine ranked modes, let alone reproduce the exact text of the score and chemistry records written above them. The input conformers for these ligands are also byte-identical, which is expected and benign, because all three builders embed the same sequence with the same default embedding seed of 42. Identical input does not produce identical output from an unseeded search.

The version history settles the direction and leaves the cause open. The campaign pose was first committed on 16 November at 16:35 with one content, and that content was replaced on the same day at 18:34 by the reformatting described in the Limitations. Both audit paths were created four days later, on 20 November at 10:53, in a single commit, and the content they were created with is the 18:34 state of the campaign file. Neither audit path has ever held anything else. So the bytes now at the audit paths already existed elsewhere four days before those paths came into being, and no independent docking output was ever written to them. That is as far as the artefacts go.

An earlier version of this section offered a mechanism, and I withdraw it. All three scripts guard the docking call with a test of the form: if the output file does not exist, or is smaller than a hundred bytes, run Vina. A file already sitting at the expected path is therefore adopted silently, with no check that it was produced under the protocol the script is about to report. That is a real defect, it makes the guard a resume feature for an interrupted run and a silent provenance hole for everything else, and it is worth fixing. But it cannot be the explanation here, because the three paths are different paths. The guard only adopts content already at the path it is about to write, and it cannot put one file's bytes under another file's name. I had the right defect and the wrong inference.

Nothing in the repository records how the content reached those paths. The pose files carry Vina's score and the ligand chemistry and nothing else: no timestamp, no seed, no exhaustiveness, no command line, no run identifier. The audit battery's own report carries no timestamp and no command line either. The only log file in the repository is empty. I put the question to the author, who does not recall. The duplication is therefore unexplained, and it is unexplainable from what was kept. I will not write that one file was copied to another name, and I will not write that three independent runs failed, because the artefacts establish neither.

What follows from the duplication alone is enough without a cause, and it is four things.

First, this is one docking result carried under three names, so the audit battery's reproducibility test has four independent runs and not five, which is why Section 3.2 reports four. The battery's own record states that it used five independent seeds and lists the duplicated value first among the five, and its reported dispersion of 0.371 kcal/mol is computed over those five values. That statement does not survive the file comparison.

Second, the lead's published affinity of -8.191 kcal/mol cannot be attributed to a stated exhaustiveness, because three different stated protocols carry the identical file.

Third, the published pH effect is not a measurement. The reported shift of 0.485 kcal/mol at GTF-SI is the difference between the duplicated file and a genuine fresh run of the diprotonated ligand at -7.706 kcal/mol, and 0.485 kcal/mol is 1.02 times the run-to-run standard deviation of a single ligand established in Section 3.2. The comparison against alpha-amylase moves 0.798 kcal/mol in the opposite direction, which would say that acidification improves off-target binding, and its neutral arm is the second duplicate group, so it rests on the same unexplained file. Both arms of the acid-switch claim at both targets take their neutral reference from a file of unestablished origin, and neither shift is outside the noise.

Fourth, every downstream number that reuses either neutral arm inherits the same defect, and I flag it as such wherever it appears rather than quoting it as a result.

This is a finding about the audit rather than about the peptide, and I report it for the same reason I report the rest. An audit that silently inherits a file it believes it produced can confirm anything.

## 3.11 Disagreements I have not resolved

One quantity remains inconsistent between two records. The relaxed potential energy of the lead is recorded as -1235.93 kJ/mol in the primary results database and -1234.68 kJ/mol in the benchmark file, for what is described as the same calculation. I have not been able to determine which run produced which value, and I report the disagreement rather than choosing. Neither value supports a claim about stability at any temperature, and neither supports a claim about binding. Reading the function that produces it: it takes one coordinate file, which the callers pass as the ligand alone, builds an Amber14 system with GBn2 implicit solvent, and calls energy minimisation with a 200-iteration cap. No integrator is ever stepped and the receptor is never present. The number is therefore the relaxed potential energy of an isolated macrocycle in implicit solvent. It is not a thermal measurement, it is not molecular dynamics, and it is not a biophysical validation of a complex, because no complex was ever built. The repository's own notebooks state the first of these correctly; nothing in the project should describe it as the third.

# 4. Discussion

## 4.1 What this result is

The screen produced a lead that its own controls cannot distinguish from a scrambled sequence of the same composition, using a score that varies more within a single docking run than the difference the lead was selected on, from a setup that fails its redocking calibration. Those three facts are independent of one another and each is sufficient to withhold the lead. They are also not a sampling problem: Section 3.3b shows that the margin this screen produced sits at 0.27 of the null's own standard deviation, so no number of additional control sequences would bring it inside a 5 per cent rejection region, and the dispersion that sets that region is of the same order as the run-to-run noise of a single ligand.

The useful finding is not about these peptides. It is that a screen can run to completion, produce a ranked leaderboard, pass a dozen named falsification tests, and generate plausible mechanistic prose, while never having demonstrated that it can tell signal from composition. Every diagnostic reported here is cheap. Redocking the co-crystallised ligand takes one run. Fixing a seed takes one flag. Checking that a scrambled null is a permutation of the molecule it is a control for takes one line of code. Comparing within-run pose spread against the differences being ranked takes the output already on disk. None of them was run before results were written up, and any one of them would have stopped the campaign.

## 4.2 What this result is not

This is not evidence that cyclic peptides cannot inhibit GTF-SI, and it is not evidence that these particular sequences are inactive. A negative in silico result of this kind carries almost no information about the molecules, for a reason worth stating directly: the failure is in the measuring instrument. When an assay cannot separate a positive control from a negative control, its readings on unknowns are uninformative in both directions, and interpreting them as evidence of inactivity is the same error as interpreting them as evidence of activity, with the sign reversed.

The GH70 glucansucrases remain a reasonable target for peptide inhibition, and Ito and colleagues' own conclusion about the value of the structure for inhibitor design [1] is untouched by anything here.

## 4.2b The instrument's documented properties bear on all of this

A third-party review of this project argued that AutoDock Vina is the wrong tool for a 1.2 kDa cyclic nonapeptide, and that stripping the crystallographic heteroatoms biases the screen toward cationic peptides. I checked both against Vina's own description of its scoring function [2] and against the receptor files, and the picture is more specific than either claim.

**What the scoring function contains.** Vina's conformation-dependent score is a weighted sum of three steric terms, a hydrophobic term between hydrophobic atoms, and, where applicable, hydrogen bonding; hydrogen atoms are not scored explicitly and enter only through atom typing; all interaction functions are cut off at 8 Å; and the only entropy-like term is a division of the intermolecular score by one plus a weight times the number of active rotatable bonds [2]. Three consequences follow for this screen.

First, there is no electrostatic term. Both receptor files here contain zero heteroatom records, so every non-protein atom was removed before docking, including any metal and all ordered water, and the site that remains is formally anionic: within 8 Å of the redocked acarbose centroid there are three acidic residues (Glu515, Asp588, Asp593) and no basic ones, and within 12 Å there are ten acidic and four basic, counting heavy atoms only. So the geometric part of that claim holds. The proposed mechanism does not, because a scoring function with no Coulomb term cannot preferentially attract a cationic peptide by electrostatics. The route that remains open is the hydrogen-bonding term, since removing a metal and its waters frees carboxylate oxygens as acceptors and an arginine-rich and histidine-rich peptide carries many donors. That is a weaker and different claim, and I have not tested it. My own data are consistent with the absence of an electrostatic term: the explicitly diprotonated ligand carries two more formal charges than the neutral one, its summed partial charge moves from -0.004 to +1.994, and its score at this anionic site moves by 0.485 kcal/mol, which is 1.02 times the single-ligand noise of Section 3.2. I could not establish whether the deposited structure contains a catalytic calcium at all. Ito and colleagues' abstract does not mention a metal, the paper is not in the open full-text archive, and the structure database is not reachable from the environment this analysis ran in, so I leave that open rather than assert it.

Second, the benchmark on which Vina's binding-mode accuracy was demonstrated used ligands with 0 to 32 active rotatable bonds [2]. The nine-mers here carry 15 and sit inside that range. The eleven-mers carry 35 and sit outside it, so for the naive eleven-mer comparator there is no published basis for expecting pose accuracy at all, on top of its ring having been docked rigid (Section 3.4).

Third, the rotatable-bond penalty makes one comparison in this dataset point the wrong way. The linearised control SPHWRPWHS and the cyclic parent CPHWRPWHC have the same heavy-atom count of 86 and were docked at the same exhaustiveness, but the linear peptide carries 35 active torsions against the macrocycle's 15, so the penalty weighs on it more heavily. It nonetheless scores -8.865 against the macrocycle's -8.062, an advantage of 0.803 kcal/mol, which is 1.7 times the single-ligand noise. Cyclisation is the project's central design premise, and the one controlled test of it in this dataset goes against the premise, in the direction the scoring function's own penalty should have prevented. That is a single pair of runs and I present it as a flag and not a result.

I decline the review's broader conclusion that Vina is simply the wrong tool, because it is an interpretation and not a measurement, and because the decoupling of score from pose in Section 3.2 is equally consistent with a rugged search over a degenerate scoring landscape. What is supportable is narrower and sufficient. The instrument was benchmarked on a ligand class that these eleven-mers fall outside, it contains no term for the electrostatics the project's acid-switch narrative depends on, and its own authors report a 2.85 kcal/mol standard error on affinity that is larger than every difference this screen ranked on.

## 4.3 What would make this screen informative

In order of priority:

1. **Pass redocking before interpreting anything.** Reproduce the crystallographic acarbose pose to better than 2.0 Å symmetry-corrected RMSD, or document why the protocol cannot and treat all pose-level interpretation as unavailable until it can.
2. **Make the search converge, then show it.** Raise exhaustiveness until repeated runs of the same ligand agree in score and pose, fix the random seed, and report the within-run and between-run spread alongside every ranking. A ranking whose increments are smaller than its own noise is not a ranking.
3. **Match the null to the molecule being tested.** Generate the scrambled set from the candidate's own composition, not its parent's, and report the empirical p as (k+1)/(n+1).
4. **Check steric validity against the intact receptor**, never against the slice that will be deposited.
5. **Replace generated mechanistic prose with measured geometry.** Every interaction asserted in text should be computed from the coordinates in the same step that writes the text, and the two should be reconciled automatically.
6. **Validate wet-lab before any further computational elaboration.** The informative next experiment is an enzymatic assay on three or four synthesised peptides, not another generation of the screen.

# 5. Limitations

I state these at full strength, including those that weaken the paper's own conclusions.

**The target was misnamed throughout the pipeline.** The receptor is GTF-SI, the GtfC product, and was described in the code, output records and figures as GtfB and UniProt P08987. The docking ran against the structure the residue numbering refers to, so the computational work is internally consistent, but every artefact this project published carried the wrong enzyme name. A corrected GtfB campaign is in progress (Section 7).

**One docking program, one scoring function family.** All conclusions rest on AutoDock Vina with Vinardo rescoring. Vinardo is derived from Vina [4], so the two scores are not independent, and a second program from a different lineage would be needed to separate a scoring-function artefact from a search artefact.

**One ligand conformer per sequence, and ring flexibility applied inconsistently.** Each macrocycle entered docking as a single ETKDG/MMFF94 conformer, and the conformational dependence of the results was never tested. For the nine-mers and ten-mers, Meeko opened the ring so its torsions were sampled during docking, so the starting ring geometry mattered less. For the eleven-mers it did not, and those rings were docked rigid from that single conformer. An earlier version of this manuscript stated that the ring was opened in all cases; that was wrong, and the corrected position is in Sections 2.2 and 3.4.

**Search effort was not held constant across the comparisons.** The exhaustiveness the code specifies is 16 for the lead and the alanine scan, 8 for the scrambled nulls, and 4 for most of the leaderboard (Section 3.4). Every comparison in this paper that crosses those groups compares unequal search effort as well as unequal molecules. Worse, for the lead itself the specified value cannot be taken as the value that ran, because its pose file is one of the three duplicates in Section 3.10. The direction of the bias favours the project's original conclusions, so correcting it would not rescue the lead, but no quantitative comparison here should be read as if the conditions were matched, and I give no size for the effect because the dataset contains no clean paired comparison. An earlier version of this manuscript asserted first that the asymmetry inflated the lead and then that the inflation was measurably negligible; both are withdrawn.

**The committed coordinate files are a reformatted copy of the docking output.** A single commit in the project's history rewrote coordinates in 257 pose and receptor files. I measured every change: 510 coordinates moved, each by exactly 0.010 Å, 504 positive and 6 negative, with no change larger than 0.010 Å. That is the signature of a float round-trip through a three-decimal format, not an edit of the geometry, and 0.010 Å is below every threshold used anywhere in this paper. I record it because the files are therefore not the verbatim program output, and because it bears on Section 3.10: the duplicated audit files match the reformatted version of the campaign pose, which post-dates it.

**The duplication of Section 3.10 is unexplained and the author cannot attest to its origin.** I asked, and the answer was that he does not remember. The artefacts do not fill the gap: the pose files record no timestamp, seed, exhaustiveness, command line or run identifier, the audit report records no timestamp or command line, and the repository's only log file is empty. The honest position is that the repository contains one docking result under three names at GTF-SI and one under two names at alpha-amylase, that nothing recorded distinguishes a docking run from a file write at those paths, and that no reader of this dataset can verify which it was. Anyone repeating this work should write the seed, the exhaustiveness, the command line and a run identifier into the output file itself, so that the question cannot arise. The accompanying repository now does that: result paths are derived from the run's parameters rather than chosen by hand, every new result carries that record internally, and a check that fails when two result paths share coordinates without sharing a run identifier runs in continuous integration. It does not repair this dataset, and it is not meant to. It means the same question cannot be unanswerable next time.

**The audit battery has defects of its own.** Besides the duplicate reporting of Section 3.10, the multi-seed routine substitutes a hardcoded score of -8.191 if it cannot parse one from a pose file, so a missing or malformed file yields the lead's own expected value rather than an error. Three of the twelve angles returned values from Python literals rather than from computation before I corrected them. None of the numbers I report in this paper is taken from that battery's output; all are recomputed from the primary artefacts by the script named in Section 6.

**No wet-lab validation of any kind.** No peptide was synthesised, no enzyme was assayed, no biofilm was grown. Nothing in this paper is an experimental measurement.

**A negative in silico result does not establish inactivity.** This bears repeating as a limitation and not only as a conclusion. Nothing here licenses the statement that these peptides do not inhibit GTF-SI.

**Occlusion was not tested.** Section 3.7 compares catalytic-triad contact density, which is one mechanism of competitive inhibition among several. A ligand that occupies the cleft without contacting the triad may still block substrate access. I did not test this, so the triad comparison is weaker evidence than it looks.

**The oral-stability gates are not validated.** The protease gate counts sequence recognition motifs and does not model accessibility, so a buried motif in a constrained macrocycle is counted as a liability it may not be. The pH-dependent charge calculation uses Henderson-Hasselbalch on isolated side chains and ignores local electrostatic environment; across this set it reduces to a count of histidines, since each contributes 0.809 charge units between pH 7.2 and 5.2 and 214 of 260 sequences are described exactly by that count. Neither gate has been checked against measurement.

**One stability predictor was applied out of its domain.** The instability index [5] was trained on 12 unstable and 32 stable linear proteins to predict in vivo metabolic stability, and is length-normalised. Applied to a disulfide-cyclised 9-mer it is being used far outside that domain, and for arginine-rich sequences it returns values an order of magnitude above its own nominal ceiling. It should not have been used here and no conclusion in this paper rests on it.

**Figures have not been regenerated.** The published figures still reflect pre-correction values and the superseded falsification tally. They are flagged in the repository and are not reproduced in this paper.

**One numerical disagreement is unresolved** (Section 3.11).

**I am a single author, and I wrote the pipeline whose failure I am reporting.** No independent party has reproduced these recomputations. The provenance script is provided so that reproduction requires no cooperation from me.

# 6. Data and code availability

All primary artefacts are in the accompanying repository: receptor and ligand coordinates, all 1,244 retained pose files, the 260-row results database, the falsification battery output, and the deposited complexes discussed in Sections 3.7 to 3.9.

`scripts/compute_paper_numbers.py` recomputes every number in this paper from those artefacts and prints each with its source and a tag recording whether it was recomputed or read from a recorded run. It also prints the disagreements it finds, which is how Sections 3.3, 3.4, 3.8, 3.9 and 3.10 were written. Running it is the fastest way to check this paper against its data.

# 7. Related work in progress

> **PLACEHOLDER, DO NOT PUBLISH THIS SECTION AS WRITTEN.**
>
> A separate line of work is re-running this campaign against GtfB (GTF-I, UniProt P08987), the enzyme the original pipeline named, against the GtfB residues equivalent to GTF-SI Asp477, Glu515 and Asp588. I have deliberately not written the GtfB residue numbers here, because I have not verified that mapping against a primary source myself, and a residue number I cannot trace is exactly the kind of claim this paper exists to object to. The mapping is to be stated with its citation, or not stated.
>
> **No result from that work is available to me at the time of writing, and no number from it appears anywhere in this paper.** This section is to be completed with that result, or deleted, before this manuscript goes anywhere. If the GtfB campaign is not complete at submission, delete this section rather than describing work whose outcome is unknown.

# References

According to PubMed:

1. Ito K, Ito S, Shimamura T, Weyand S, Kawarasaki Y, Misaka T, Abe K, Kobayashi T, Cameron AD, Iwata S. Crystal structure of glucansucrase from the dental caries pathogen *Streptococcus mutans*. *J Mol Biol*. 2011;408(2):177-86. [DOI](https://doi.org/10.1016/j.jmb.2011.02.028)
2. Trott O, Olson AJ. AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, efficient optimization, and multithreading. *J Comput Chem*. 2010;31(2):455-61. [DOI](https://doi.org/10.1002/jcc.21334)
3. Eberhardt J, Santos-Martins D, Tillack AF, Forli S. AutoDock Vina 1.2.0: new docking methods, expanded force field, and Python bindings. *J Chem Inf Model*. 2021;61(8):3891-8. [DOI](https://doi.org/10.1021/acs.jcim.1c00203)
4. Quiroga R, Villarreal MA. Vinardo: a scoring function based on AutoDock Vina improves scoring, docking, and virtual screening. *PLoS One*. 2016;11(5):e0155183. [DOI](https://doi.org/10.1371/journal.pone.0155183)
5. Guruprasad K, Reddy BV, Pandit MW. Correlation between stability of a protein and its dipeptide composition: a novel approach for predicting in vivo stability of a protein from its primary sequence. *Protein Eng*. 1990;4(2):155-61. [DOI](https://doi.org/10.1093/protein/4.2.155)
6. Ito K, Ito S, Shimamura T, Kawarasaki Y, Abe K, Misaka T, Kobayashi T, Iwata S. Crystallization and preliminary X-ray analysis of a glucansucrase from the dental caries pathogen *Streptococcus mutans*. *Acta Crystallogr Sect F Struct Biol Cryst Commun*. 2010;66(Pt 9):1086-8. [DOI](https://doi.org/10.1107/S1744309110029714)
