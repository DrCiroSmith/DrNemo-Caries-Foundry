#!/usr/bin/env python3
"""
Dr. Nemo - Phase Report
=======================
Reads master_scores.csv and prints a clear per-target summary
to help you decide which targets to focus on next.

Usage:
    python3 phase_report.py
"""
import csv
from pathlib import Path
from collections import defaultdict

BASE_DIR = Path.home() / "drnemo_research"
SCORES_FILE = BASE_DIR / "results" / "master_scores.csv"


def main():
    if not SCORES_FILE.exists():
        print("❌ No master_scores.csv found. Run the pipeline first.")
        return

    # Collect per-project stats
    stats = defaultdict(lambda: {
        "total": 0, "success": 0, "failed": 0,
        "scores": [], "best_score": 0, "best_job": "",
        "best_plddt": 0, "best_iptm": 0,
        "plddts": [], "iptms": []
    })

    with open(SCORES_FILE, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            project = row.get("Project", "unknown") or "unknown"
            s = stats[project]
            s["total"] += 1

            if row["Status"] == "SUCCESS":
                s["success"] += 1
                try:
                    composite = float(row.get("Composite_Score", 0) or 0)
                    plddt = float(row.get("pLDDT", 0) or 0)
                    iptm = float(row.get("iPTM", 0) or 0)
                    s["scores"].append(composite)
                    s["plddts"].append(plddt)
                    s["iptms"].append(iptm)
                    if composite > s["best_score"]:
                        s["best_score"] = composite
                        s["best_job"] = row.get("Job_ID", "")
                        s["best_plddt"] = plddt
                        s["best_iptm"] = iptm
                except (ValueError, TypeError):
                    pass
            else:
                s["failed"] += 1

    # Print report
    print()
    print("=" * 80)
    print("🦷 DR. NEMO — PHASE REPORT")
    print("=" * 80)
    print()

    # Sort by best composite score descending
    sorted_projects = sorted(stats.items(), key=lambda x: x[1]["best_score"], reverse=True)

    total_all = sum(s["total"] for _, s in sorted_projects)
    success_all = sum(s["success"] for _, s in sorted_projects)

    print(f"📊 GLOBAL: {success_all} successful folds out of {total_all} total jobs")
    print()

    for rank, (project, s) in enumerate(sorted_projects, 1):
        avg_score = sum(s["scores"]) / len(s["scores"]) if s["scores"] else 0
        avg_plddt = sum(s["plddts"]) / len(s["plddts"]) if s["plddts"] else 0
        avg_iptm = sum(s["iptms"]) / len(s["iptms"]) if s["iptms"] else 0

        # Classification
        if s["best_score"] >= 0.72:
            verdict = "🔴 HIT! Ready for Phase 3 optimization"
        elif s["best_score"] >= 0.5:
            verdict = "🟠 PROMISING — Continue evolving"
        elif s["best_score"] >= 0.3:
            verdict = "🟡 MODERATE — Needs more variants"
        elif s["success"] > 0:
            verdict = "🔵 EARLY — Keep exploring"
        else:
            verdict = "⚪ NO DATA yet"

        print(f"{'─' * 80}")
        print(f"  #{rank}  {project}")
        print(f"{'─' * 80}")
        print(f"  Samples: {s['success']} successful / {s['total']} total ({s['failed']} failed)")
        if s["scores"]:
            print(f"  Avg Composite:  {avg_score:.4f}    Avg pLDDT: {avg_plddt:.1f}    Avg iPTM: {avg_iptm:.4f}")
            print(f"  BEST Composite: {s['best_score']:.4f}    pLDDT: {s['best_plddt']:.1f}    iPTM: {s['best_iptm']:.4f}")
            print(f"  Best candidate: {s['best_job']}")
        print(f"  Verdict: {verdict}")
        print()

    # Recommendation
    print("=" * 80)
    print("🎯 RECOMMENDATION")
    print("=" * 80)

    hits = [(p, s) for p, s in sorted_projects if s["best_score"] >= 0.72]
    promising = [(p, s) for p, s in sorted_projects if 0.5 <= s["best_score"] < 0.72]
    needs_data = [(p, s) for p, s in sorted_projects if s["success"] < 50 and s["best_score"] < 0.5]

    if hits:
        print(f"\n  🔴 PHASE 3 READY ({len(hits)} targets): Fine-tune these winners!")
        for p, s in hits:
            print(f"     → {p} (best: {s['best_score']:.4f})")

    if promising:
        print(f"\n  🟠 KEEP EVOLVING ({len(promising)} targets): These are close!")
        for p, s in promising:
            print(f"     → {p} (best: {s['best_score']:.4f}, need > 0.72)")

    if needs_data:
        print(f"\n  🔵 NEEDS MORE DATA ({len(needs_data)} targets): Run more variants")
        for p, s in needs_data:
            print(f"     → {p} ({s['success']} samples so far)")

    print()


if __name__ == "__main__":
    main()
