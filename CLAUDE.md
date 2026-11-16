# DrNemo-Caries-Foundry (Claude Code notes)

Public portfolio repo (GitHub DrCiroSmith). Oral-biophysics peptide design vs S. mutans GtfB. Author line: "DrNemo Oral Biophysics Research, Miami, FL".

## Rules
- Honesty first: every number in docs/paper/LinkedIn must come from `data/foundry_contenders.sqlite`, `data/contenders_ledger.csv`, `falsification_benchmark.json`. No fabricated scores, no `numpy random` in figures or notebooks.
- Pose hashes are real sha256 of `data/poses/*.pdbqt` bytes. Never edit scores by hand.
- No em dashes anywhere. LinkedIn post at `/mnt/c/dev/job_packets/linkedin_drnemo/linkedin_post_drnemo.md` <= 2,350 chars.
- Never expose TFA-related code, SOPs, or PHI. Never print `.env` or tokens.
- Token economy: wait for long jobs inside ONE blocking script; read `tail -n 25` only; never cat large CSV/JSON/HTML.

## Layout
- Engines: `biophysics_engine.py`, `falsification_engine.py`, `diffdock_pipeline.py`, `solubility_shield.py`
- Campaign scripts: `scripts/run_falsification_campaign.py`, `scripts/run_real_campaign.py`, `scripts_scrambled_runner.py`
- Data: `data/poses/*.pdbqt`, `data/foundry_contenders.sqlite`, `data/contenders_ledger.csv`
- SOP notebook: `scripts/build_sop_notebook.py` -> `notebooks/SOP_01_Biophysical_Falsification_Protocol.ipynb`
- Paper: `manuscript/` (`DrNemo_Preprint.tex/.pdf`), `paper/paper.md`, `FALSIFICATION_BENCHMARK.md`, `docs/index.html`
- Tests: `tests/test_foundry.py` (run `pytest -q`)

## Reference results (verify against ledger before quoting)
acarbose -7.915, linear -8.865, WT CPHWRPWHC -8.062 (1SMD -7.902, cross-reactive), 28A Pareto winner CAHWRPWHC -8.191 (1SMD -7.343, ddG -0.85, 0 protease sites). Target: 39/39 scrambled nulls, >= 52 docked contenders.

## Docking settings
8 workers, `--cpu 3`, exhaustiveness 4, skip existing valid poses.

## Done-check
`pytest -q`, notebook executes, `git push origin main`, CI green, `drnemo-discovery.service` active.
