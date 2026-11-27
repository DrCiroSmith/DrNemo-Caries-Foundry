#!/usr/bin/env python3
"""Fail when two result files share coordinates without sharing a run.

Run it in CI, or by hand:

    python3 scripts/check_result_provenance.py
    python3 scripts/check_result_provenance.py --json
    python3 scripts/check_result_provenance.py --include-splits

Exit status 0 when every duplicate group is either recorded in the baseline or
self-declared as one run under several names, and 1 when a group is new. A new
group is either a false record, because a stochastic search does not reproduce
another run's coordinates, or an untracked copy. Both need a person to look.

The baseline at data/provenance_baseline.json holds the groups that existed
before this convention, with the reason each is kept. It is evidence for the
manuscript and is not to be regenerated.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from provenance import audit  # noqa: E402

BASELINE = REPO / "data" / "provenance_baseline.json"


def _display(path: str) -> str:
    """Repository-relative where possible, absolute otherwise. Never raises."""
    resolved = Path(path).resolve()
    try:
        return str(resolved.relative_to(REPO))
    except ValueError:
        return str(resolved)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", default=None,
                    help="directory to scan (default: the whole repository)")
    ap.add_argument("--include-splits", action="store_true",
                    help="also check the per-mode files vina_split derives")
    ap.add_argument("--no-baseline", action="store_true",
                    help="ignore the baseline and report every group")
    ap.add_argument("--json", action="store_true", dest="as_json")
    args = ap.parse_args()

    # Resolved up front, so a relative --root behaves exactly like an absolute
    # one. Before this, "--root ." reached the display code with relative paths
    # and crashed there.
    root = Path(args.root).resolve() if args.root else REPO
    paths = [p for p in sorted(root.rglob("*.pdbqt")) if ".git" not in p.parts]
    baseline = {}
    if not args.no_baseline and BASELINE.exists():
        baseline = json.loads(BASELINE.read_text())

    report = audit(paths, baseline=baseline, include_splits=args.include_splits)
    report["root"] = str(root)

    hard = [g for g in report["findings"] if g["kind"] != "one_run_many_names"]
    soft = [g for g in report["findings"] if g["kind"] == "one_run_many_names"]
    report["new_contradictions"] = len(hard)

    if args.as_json:
        print(json.dumps(report, indent=2, sort_keys=True))
        return 1 if hard or report["name_mismatch"] else 0

    print(f"Scanned {report['files_examined']} result files under {root}, "
          f"{report['distinct_coordinate_sets']} distinct coordinate sets, "
          f"{report['derived_splits_skipped']} derived splits skipped.")
    print(f"Groups sharing coordinates: {len(report['groups'])}. "
          f"Known from the baseline: {len(report['baselined'])}. "
          f"New: {len(report['findings'])}.")

    for g in report["baselined"]:
        print(f"\n  known  [{g['kind']}] {len(g['files'])} files, "
              f"digest {g['coordinate_digest'][:12]}")
        for f in g["files"]:
            print(f"           {_display(f)}")

    for g in soft:
        print(f"\n  NOTE   one run written under {len(g['files'])} names, "
              f"run id {g['run_ids'][0][:12]}. Transparent, but pick one name.")
        for f in g["files"]:
            print(f"           {f}")

    for g in hard:
        why = ("different runs cannot produce identical coordinates, so one of "
               "these records is false" if g["kind"] == "contradiction" else
               "at least one of these files carries no provenance record, so "
               "the duplication cannot be explained from the artefacts")
        print(f"\n  FAIL   [{g['kind']}] {len(g['files'])} files, "
              f"digest {g['coordinate_digest'][:12]}: {why}")
        for f in g["files"]:
            print(f"           {f}")

    for m in report["name_mismatch"]:
        print(f"\n  FAIL   {m['file']} was written under a hand-chosen name; "
              f"its own record says it should be {m['recorded_canonical_name']}")

    if hard or report["name_mismatch"]:
        print("\nResult: FAIL. A result file is in a state no run can account for.")
        return 1
    print("\nResult: pass. Every duplicate group is accounted for.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
