#!/usr/bin/env bash
# One command for the whole package. Stops at the first failure rather than
# carrying on with a bad receptor or an unverified statistics module.
#
# Writes only inside retarget_gtfb/. See README_RUN_THIS.md.
#
#   bash RUN_ME.sh              # everything, roughly nine and a half hours
#   bash RUN_ME.sh 1,2          # tiers 1 and 2 only, the result without the stress tests
#   TIERS=1 CPU=16 bash RUN_ME.sh

set -euo pipefail

cd "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

TIERS="${1:-${TIERS:-1,2,3,4,5,6}}"
CPU="${CPU:-8}"

step() { printf '\n=== %s ===\n' "$1"; }

step "1/6 Fetching the 8FJC coordinate file"
bash fetch_receptor.sh

step "2/6 Validating that it really is GtfB"
python3 validate_receptor.py

step "3/6 Preparing the receptor and deriving the grid boxes"
python3 prepare_receptor.py

step "4/6 Self-testing the statistics code before it is trusted with real data"
python3 analyze.py --self-test

step "5/6 Docking, tiers $TIERS on $CPU cores"
echo "Ctrl+C is safe at any point. Finished runs are already saved and a re-run resumes."
python3 run_docking.py --tiers "$TIERS" --cpu "$CPU"

step "6/6 Analysis"
python3 analyze.py

printf '\nDone. Send back:\n  %s\n  %s\n' "$(pwd)/out/analysis.json" "$(pwd)/out/results.jsonl"
