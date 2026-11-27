# How to run this, and what it does to your machine

Plain language first, because you should not have to read the code to know what you
are about to run.

---

## What this is

The existing campaign docked the peptide series against PDB 3AIC, which the
repository calls GtfB. 3AIC is GtfC. An experimental GtfB structure does exist:
PDB 8FJC, the GtfB catalytic domain with acarbose bound at 2.50 A. This package
re-runs the campaign against it, with the controls the original run was missing,
and tries to break whatever answer comes out.

`PROTOCOL.md` is the scientific version, written before any result existed, with
every threshold fixed in advance. Read it when you want to know why a number means
what it means. This file is only about operating the thing.

## What it does to your computer

- **It writes only inside this one folder,** `retarget_gtfb/`. Specifically
  `retarget_gtfb/receptor/` for the downloaded structure and the two prepared
  receptor files, and `retarget_gtfb/out/` for everything else. Nothing is written
  to your home directory, nothing to any other part of the repository, nothing to
  any system location.
- **It changes no file that already exists in the repository.** It reads
  `biophysics_engine.py`, `data/contenders_ledger.csv` and two prepared receptors in
  `demo_data/receptors/`, and writes to none of them.
- **It makes exactly one network request,** in `fetch_receptor.sh`: it downloads
  8FJC from the RCSB PDB, or from the PDBe mirror if RCSB does not answer. Nothing
  else in the package touches the network. Nothing is uploaded anywhere.
- **It installs nothing** and modifies no Python environment.
- **It runs AutoDock Vina many times,** which will use eight CPU cores for several
  hours. That is the whole cost. Your machine will be busy and warm.
- **It can be stopped at any moment.** Press Ctrl+C. Every finished docking run is
  already written to `out/results.jsonl`, and re-running resumes from where it
  stopped instead of starting over.

## What you need installed

AutoDock Vina 1.2.7 (both `vina` and `vina_split`), Open Babel (`obabel`), and
Python with rdkit, meeko and numpy. These are the same prerequisites your existing
pipeline lists, so if `biophysics_engine.py` runs on your machine, this will too.

The run refuses to start if any of them is missing or if `vina --version` does not
report 1.2.7. That is deliberate. A silent fallback to a different version would
make the new numbers incomparable with the old ones, and you would not find out.

## The commands

```bash
cd retarget_gtfb

bash fetch_receptor.sh          # downloads 8FJC.pdb, about 1 MB, a few seconds
python3 validate_receptor.py    # checks it really is GtfB, instant
python3 prepare_receptor.py     # builds the docking receptor and the boxes, seconds
python3 analyze.py --self-test  # checks the statistics code, seconds
python3 run_docking.py          # the long one
python3 analyze.py              # writes out/analysis.json
```

Or just `bash RUN_ME.sh`, which runs those in order and stops at the first failure.

## How long the docking takes

545 docking runs in total, in six tiers that run in this order:

| Tier | What it is | Runs |
|---|---|---|
| 1 | The four named contenders plus acarbose redocking, 5 seeds each | 25 |
| 2 | Both 39-member null distributions, 3 seeds each | 234 |
| 3 | The same thing in a smaller box, to test whether box size drove the result | 87 |
| 6 | The same thing against a receptor that keeps the structural Ca2+ | 87 |
| 4 | Four more starting conformers, to test whether conformer luck drove it | 88 |
| 5 | Amylase and lysozyme counter-screens | 24 |

Tier 6 is the calcium control. Preparation strips every heteroatom from the main
receptor, which is what the original campaign did, but GtfB carries a conserved
structural calcium next to the catalytic triad and one of the residues holding it also
hydrogen bonds acarbose. So the package builds the receptor both ways and docks against
both. If your Vina build will not accept a receptor with a metal in it, or Open Babel
drops the ion during conversion, the preparation step finds that out in seconds, writes
the reason into `out/boxes.json`, and tier 6 skips itself. The run does not fail and the
remaining tiers are unaffected.

At roughly a minute per run on eight cores that is about nine hours, so it is an
overnight job. The script prints its own estimate after the third run, based on your
actual timings rather than my guess, so you will know within a few minutes whether
to leave it running.

**Tier 1 is the important one and it finishes in under half an hour.** If you want
the answer before committing a whole night:

```bash
python3 run_docking.py --tiers 1,2
python3 analyze.py
```

Tiers 1 and 2 together give the redocking calibration and both null tests, which is
the result. Tiers 3 to 6 are the attempts to break it. Send me whatever you have and
I will analyse it; partial output is genuinely useful.

If you want the calcium check in the same sitting, `--tiers 1,2,6` adds about another
hour and a half. It is the one stress test that could change the answer rather than
merely qualify it, because it changes the receptor rather than the search.

Other switches: `--cpu 16` if you have more cores, `--dry-run` to list the jobs
without running anything.

## What to send back

`out/analysis.json` and `out/results.jsonl`. Between them they hold every number I
need. If `results.jsonl` is awkwardly large, `analysis.json` alone is enough to
start on.

## What will stop the run, and why each check is there

The package is built to fail loudly rather than produce a number that looks fine.
It stops if:

- the downloaded file is under 200 kB, which means a truncated download or an error
  page rather than coordinates;
- the title does not name GtfB, acarbose and space group P4322, or does name GtfC or
  3AIC;
- the resolution in the file is not 2.50 A;
- any of the eleven catalytic and acarbose-contact residues is the wrong amino acid
  at its GtfB residue number. This is the check that cannot be fooled: at those
  numbers a GtfC structure holds different residues, so a mislabelled file fails
  here even if everything else looks right;
- the receptor file's SHA-256 differs from the one recorded on the first run, which
  would mean two different receptors mixed into one result set;
- the bound acarbose centroid and the catalytic triad centroid disagree by more than
  4 A, which would mean the active site has been located wrongly;
- fewer than 5000 protein atoms survive filtering, which would mean a partial model;
- `vina --version` does not contain 1.2.7, or `vina`, `vina_split` or `obabel` is
  missing;
- the ligand builder imported from `biophysics_engine.py` and the local copy in
  `run_docking.py` do not produce byte-identical output for the same sequence and
  seed, which would mean the conformer test was measuring a difference between two
  builders instead of between two seeds.

One thing deliberately does **not** stop the run: a missing or unusable calcium
receptor. That is recorded as a limitation and tier 6 is skipped, because an absent
control is worth knowing about and is not worth losing a night's docking over.

## One thing I could not do from here, so you are doing it

I had no network access to the RCSB or PDBe download endpoints while preparing this,
so I have never seen the 8FJC coordinate file. Two consequences, both stated plainly
rather than papered over:

1. **There is no pre-computed checksum of the correct file to compare against.** I
   could have written a plausible-looking hash into the code. It would have been
   fabricated, and it would have passed silent review forever. Instead
   `validate_receptor.py` validates the file's contents against facts published in
   Schormann et al. 2023, and records the hash of whatever you actually download so
   that every later result is traceable to those exact bytes. A checksum you cannot
   verify is not provenance, it is decoration.
2. **I did not reconstruct the coordinates from any other source.** A coordinate
   file produced by a model is not a structure. You download it from the PDB, and
   the validator tells you whether what arrived is what the package expects.

If a check fails, send me the error. Do not edit the thresholds to get past it.
