#!/usr/bin/env python3
"""Compute the real Pareto front over the contender vault.

Why this script exists. The vault carries a column called pareto_rank, and the
write-ups called CAHWRPWHC the "Pareto champion" and "the non-dominated
candidate". Neither claim came from a dominance computation: pareto_rank is
assigned from the falsification status (1 for a row that survived every gate,
2 or 3 otherwise), so it is a status tier under a misleading name, and nothing
anywhere in the repository tested dominance.

This script does test it, so the claim has a producer. Run:

    python3 scripts/compute_pareto_front.py

Objectives are declared explicitly below. Both are minimised, which is the
convention the figures use: a more negative docking score is better, and a more
negative selectivity DeltaDeltaG means a stronger preference for the bacterial
target over human amylase. One candidate dominates another when it is no worse
on both and strictly better on at least one.
"""
import sqlite3
import sys
from pathlib import Path
from typing import Dict, List, Sequence

REPO_ROOT = Path(__file__).parent.parent.resolve()
DB_PATH = REPO_ROOT / "data" / "foundry_contenders.sqlite"
LEAD_SEQUENCE = "CAHWRPWHC"
SURVIVED_STATUSES = ("SURVIVED_ALL_GATES", "CORROBORATED", "NOT_CORROBORATED")

# (column, direction). "min" means a smaller number is better.
OBJECTIVES: Sequence = (
    ("vina_3aic_kcal", "min"),
    ("selectivity_ddg", "min"),
)


def load_rows() -> List[Dict]:
    conn = sqlite3.connect(f"file:{DB_PATH}?mode=ro", uri=True)
    conn.row_factory = sqlite3.Row
    try:
        columns = ", ".join(name for name, _ in OBJECTIVES)
        rows = conn.execute(
            f"SELECT sequence, {columns}, heavy_atoms, falsification_status, pareto_rank "
            f"FROM contenders"
        ).fetchall()
    finally:
        conn.close()
    return [
        dict(row)
        for row in rows
        if all(row[name] is not None for name, _ in OBJECTIVES)
    ]


def dominates(a: Dict, b: Dict) -> bool:
    """True when a is no worse than b on every objective and better on at least one."""
    no_worse = True
    strictly_better = False
    for name, direction in OBJECTIVES:
        x, y = a[name], b[name]
        if direction == "min":
            no_worse = no_worse and x <= y
            strictly_better = strictly_better or x < y
        else:
            no_worse = no_worse and x >= y
            strictly_better = strictly_better or x > y
    return no_worse and strictly_better


def pareto_front(rows: Sequence[Dict]) -> List[Dict]:
    return [
        row
        for row in rows
        if not any(
            dominates(other, row)
            for other in rows
            if other["sequence"] != row["sequence"]
        )
    ]


def report(label: str, rows: Sequence[Dict]) -> bool:
    """Print the front for one population. Returns whether the lead is on it."""
    front = pareto_front(rows)
    lead = next((r for r in rows if r["sequence"] == LEAD_SEQUENCE), None)
    print(f"\n{label}: {len(rows)} candidates, {len(front)} non-dominated")
    for row in sorted(front, key=lambda r: r["vina_3aic_kcal"]):
        marker = "  <- lead" if row["sequence"] == LEAD_SEQUENCE else ""
        print(
            f"    {row['sequence']:<12} vina {row['vina_3aic_kcal']:>8.3f}  "
            f"ddg {row['selectivity_ddg']:>7.3f}  {row['falsification_status']}{marker}"
        )
    if lead is None:
        print(f"    {LEAD_SEQUENCE} is not in this population")
        return False

    on_front = any(r["sequence"] == LEAD_SEQUENCE for r in front)
    dominating = [r for r in rows if r["sequence"] != LEAD_SEQUENCE and dominates(r, lead)]
    print(f"    {LEAD_SEQUENCE} on the front: {on_front}; dominated by {len(dominating)} candidates")
    for row in sorted(dominating, key=lambda r: r["vina_3aic_kcal"])[:10]:
        print(
            f"      dominated by {row['sequence']:<12} vina {row['vina_3aic_kcal']:>8.3f}  "
            f"ddg {row['selectivity_ddg']:>7.3f}  {row['falsification_status']}"
        )
    if len(dominating) > 10:
        print(f"      ... and {len(dominating) - 10} more")
    return on_front


def rank_lead_on_each_objective(rows: Sequence[Dict]) -> None:
    """Where the lead actually sits on each objective, among the rows that survived."""
    survivors = [r for r in rows if r["falsification_status"] in SURVIVED_STATUSES]
    lead = next((r for r in survivors if r["sequence"] == LEAD_SEQUENCE), None)
    if lead is None or not survivors:
        return

    print(f"\n{LEAD_SEQUENCE}'s rank among the {len(survivors)} rows that survived every gate:")
    for name, direction in OBJECTIVES:
        values = sorted(
            (r[name] for r in survivors if r[name] is not None),
            reverse=(direction == "max"),
        )
        position = values.index(lead[name]) + 1
        print(
            f"    {name:<20} {lead[name]:>8.3f}   rank {position} of {len(values)} "
            f"(best {values[0]:>8.3f})"
        )

    # Ligand efficiency, recomputed rather than trusted: |affinity| per heavy atom,
    # where larger is better.
    scored = [
        (abs(r["vina_3aic_kcal"]) / r["heavy_atoms"], r["sequence"])
        for r in survivors
        if r.get("heavy_atoms")
    ]
    if scored:
        scored.sort(reverse=True)
        lead_le = abs(lead["vina_3aic_kcal"]) / lead["heavy_atoms"] if lead.get("heavy_atoms") else None
        if lead_le is not None:
            position = [seq for _, seq in scored].index(LEAD_SEQUENCE) + 1
            print(
                f"    ligand_efficiency    {lead_le:>8.4f}   rank {position} of {len(scored)} "
                f"(best {scored[0][0]:.4f}, {scored[0][1]})"
            )


def main() -> int:
    if not DB_PATH.exists():
        print(f"no vault at {DB_PATH}", file=sys.stderr)
        return 2

    rows = load_rows()
    print(f"Pareto analysis of {DB_PATH.relative_to(REPO_ROOT)}")
    print("Objectives: " + ", ".join(f"{name} ({direction})" for name, direction in OBJECTIVES))

    on_front_all = report("All candidates", rows)
    survivors = [r for r in rows if r["falsification_status"] in SURVIVED_STATUSES]
    on_front_survivors = report("Candidates that survived every gate", survivors)
    rank_lead_on_each_objective(rows)

    # pareto_rank is not a dominance rank. Show that it tracks the status instead.
    mismatched = [
        r for r in rows
        if (r["pareto_rank"] == 1) != (r["falsification_status"] in SURVIVED_STATUSES)
    ]
    print(
        f"\npareto_rank == 1 and 'survived every gate' disagree on {len(mismatched)} of "
        f"{len(rows)} rows, which is what you would expect of a column assigned from the "
        f"status rather than computed from dominance."
    )

    print()
    if on_front_all or on_front_survivors:
        print(f"{LEAD_SEQUENCE} is on at least one of these fronts.")
    else:
        print(
            f"{LEAD_SEQUENCE} is on neither front. It must not be described as the "
            f"non-dominated candidate or as a Pareto champion."
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
