#!/usr/bin/env bash
# Download the 8FJC coordinate file from the PDB.
#
# This script exists because the package was prepared in an environment with no
# network access to the structural databases. The coordinates were never in hand
# here, so you fetch them yourself and validate_receptor.py checks that what
# arrived is really GtfB before anything is docked against it.
#
# Writes exactly one file: retarget_gtfb/receptor/8FJC.pdb
# Touches nothing else.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEST_DIR="$HERE/receptor"
DEST="$DEST_DIR/8FJC.pdb"

PRIMARY="https://files.rcsb.org/download/8FJC.pdb"
MIRROR="https://www.ebi.ac.uk/pdbe/entry-files/download/pdb8fjc.ent"

mkdir -p "$DEST_DIR"

if [ -s "$DEST" ]; then
    echo "$DEST already exists ($(wc -c < "$DEST") bytes). Delete it first if you want a fresh copy."
    exit 0
fi

fetch() {
    local url="$1" out="$2"
    if command -v curl >/dev/null 2>&1; then
        curl -fsSL --retry 3 --retry-delay 2 -o "$out" "$url"
    elif command -v wget >/dev/null 2>&1; then
        wget -q -t 3 -O "$out" "$url"
    else
        echo "Neither curl nor wget is installed. Download this URL in a browser and save it as:" >&2
        echo "  $DEST" >&2
        echo "  $url" >&2
        return 1
    fi
}

echo "Fetching 8FJC from the RCSB PDB"
if ! fetch "$PRIMARY" "$DEST"; then
    echo "RCSB did not respond. Trying the PDBe mirror."
    rm -f "$DEST"
    fetch "$MIRROR" "$DEST"
fi

SIZE=$(wc -c < "$DEST")
if [ "$SIZE" -lt 200000 ]; then
    echo "STOP: the downloaded file is only $SIZE bytes. A 2.5 A two-chain structure is" >&2
    echo "around 1 MB, so this is a truncated download or an error page, not coordinates." >&2
    echo "The file has been left at $DEST so you can look at it. Delete it and try again." >&2
    exit 1
fi

echo "Saved $DEST ($SIZE bytes)"
echo
echo "Header lines, so you can see for yourself what arrived:"
grep -m 4 -E '^(HEADER|TITLE|REMARK   2 RESOLUTION)' "$DEST" || true
echo
echo "Next: python3 validate_receptor.py"
