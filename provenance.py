"""Deterministic result naming and in-file provenance for docking runs.

Why this module exists
----------------------
Three result files in this repository hold one identical set of coordinates
under three different nominal protocols, and nothing recorded says how that
happened. The paper's provenance section documents it. This module is the fix,
and it has two halves, because naming alone cannot prevent the failure.

The first half is the name. A result path is derived from the run's own
parameters: target, ligand, protonation state, seed, exhaustiveness, box, and
the content of the receptor and ligand input files. Same parameters means the
same path; any difference means a different path; no human picks a label in
between. A file sitting at an unexpected path is then detectable rather than
plausible.

The second half is the record inside the file, and it is the half that makes
the failure impossible to hide. Every result carries the command line, the
seed, the exhaustiveness, the box, the engine version, the md5 of both inputs,
a run identifier unique to the execution, and a UTC timestamp. Two files that
hold the same coordinates then either carry the same run identifier, in which
case the tooling can say plainly that one run was written under two names, or
they carry different ones, and the duplication is a contradiction the tooling
catches.

The third half is enforcement, and it lives in audit() and in the CLI at
scripts/check_result_provenance.py.

Nothing here rewrites or renames an existing result file. The files that are
already committed are the evidence for the paper, and they stay as they are.
The convention applies to output written from here on.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import shlex
import subprocess
import uuid
from dataclasses import dataclass, asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

PROVENANCE_VERSION = 1
REMARK_PREFIX = "REMARK PROVENANCE"
SPLIT_RE = re.compile(r"_split_\d+\.pdbqt$")

# How much of a hash goes into a file name. Ten hex characters is 40 bits,
# which is ample for a few thousand result files and still readable.
NAME_HASH_CHARS = 10


class ResultPathExists(FileExistsError):
    """Raised when a result would overwrite an existing file without consent."""


def file_digest(path: os.PathLike | str) -> str:
    """md5 over the raw bytes of a file. Reads in chunks, never decodes."""
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def coordinate_digest(path: os.PathLike | str) -> Optional[str]:
    """md5 over the coordinate columns of every ATOM and HETATM record.

    Columns 31 to 54 inclusive, one-based, which is x, y and z in PDB and
    PDBQT format. Deliberately ignores every REMARK, so stamping provenance
    into a file does not change its coordinate digest, and two runs that
    genuinely produced the same geometry still compare equal.

    Returns None for a file with no coordinate records at all.
    """
    parts: List[str] = []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                parts.append(line[30:54])
    if not parts:
        return None
    return hashlib.md5("".join(parts).encode()).hexdigest()


def engine_version(binary: Optional[str]) -> str:
    """The docking engine's own version string, or why it is not available.

    Never raises. A missing binary is recorded as missing rather than silently
    omitted, because an absent field reads as an oversight and this one is a
    fact about the run.
    """
    if not binary:
        return "unavailable: no engine binary resolved"
    try:
        res = subprocess.run([binary, "--version"], capture_output=True,
                             text=True, timeout=30)
    except FileNotFoundError:
        return f"unavailable: {Path(binary).name} not found"
    except Exception as exc:                                  # pragma: no cover
        return f"unavailable: {type(exc).__name__}"
    out = (res.stdout or res.stderr or "").strip().splitlines()
    return out[0].strip() if out else "unavailable: engine printed no version"


@dataclass(frozen=True)
class RunSpec:
    """Everything that determines a docking result, and nothing that does not.

    Two runs with equal RunSpec fields are the same question asked twice. They
    share a canonical path, so a second answer either matches the first or
    overwrites it with consent, and cannot quietly sit beside it under another
    name.
    """

    target: str
    ligand_id: str
    receptor_digest: str
    ligand_digest: str
    box_center: Tuple[float, float, float]
    box_size: Tuple[float, float, float]
    exhaustiveness: int
    seed: Optional[int] = None
    protonation: str = "neutral"
    scoring: str = "vina"

    def canonical_fields(self) -> Dict[str, Any]:
        d = asdict(self)
        d["box_center"] = [round(float(v), 4) for v in self.box_center]
        d["box_size"] = [round(float(v), 4) for v in self.box_size]
        return d

    def spec_id(self) -> str:
        """Deterministic identity of the question. Independent of when it ran."""
        blob = json.dumps(self.canonical_fields(), sort_keys=True,
                          separators=(",", ":"))
        return hashlib.md5(blob.encode()).hexdigest()

    def stem(self) -> str:
        seed = "unset" if self.seed is None else str(self.seed)
        return "_".join((
            self.target,
            self.ligand_id,
            self.protonation,
            f"e{self.exhaustiveness}",
            f"s{seed}",
            self.scoring,
            self.spec_id()[:NAME_HASH_CHARS],
        ))

    def canonical_name(self) -> str:
        return f"{self.stem()}.pdbqt"

    def result_path(self, root: os.PathLike | str) -> Path:
        return Path(root) / self.canonical_name()

    @classmethod
    def for_run(cls, target: str, ligand_id: str, receptor: os.PathLike | str,
                ligand: os.PathLike | str, box_center: Sequence[float],
                box_size: Sequence[float], exhaustiveness: int,
                seed: Optional[int] = None, protonation: str = "neutral",
                scoring: str = "vina") -> "RunSpec":
        """Build a spec, hashing the two input files as they are on disk now."""
        return cls(
            target=target,
            ligand_id=ligand_id,
            receptor_digest=file_digest(receptor),
            ligand_digest=file_digest(ligand),
            box_center=tuple(float(v) for v in box_center),
            box_size=tuple(float(v) for v in box_size),
            exhaustiveness=int(exhaustiveness),
            seed=None if seed is None else int(seed),
            protonation=protonation,
            scoring=scoring,
        )


def new_run_id() -> str:
    """Unique to one execution, so two files can be told apart by their run."""
    return uuid.uuid4().hex


def build_record(spec: RunSpec, command: Sequence[str], run_id: Optional[str] = None,
                 engine: Optional[str] = None,
                 written_path: Optional[os.PathLike | str] = None) -> Dict[str, Any]:
    rec: Dict[str, Any] = {
        "provenance_version": PROVENANCE_VERSION,
        "run_id": run_id or new_run_id(),
        "spec_id": spec.spec_id(),
        "utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "canonical_name": spec.canonical_name(),
        "command": shlex.join(str(c) for c in command),
        "engine_version": engine if engine is not None else "unavailable: not probed",
    }
    rec.update(spec.canonical_fields())
    if written_path is not None:
        rec["written_name"] = Path(written_path).name
    return rec


def _remark_lines(record: Dict[str, Any]) -> List[str]:
    lines = []
    for key in sorted(record):
        value = record[key]
        if isinstance(value, (list, tuple)):
            value = " ".join(str(v) for v in value)
        text = str(value).replace("\n", " ").replace("\r", " ")
        lines.append(f"{REMARK_PREFIX} {key.upper()} {text}")
    return lines


def stamp(path: os.PathLike | str, record: Dict[str, Any],
          sidecar: bool = True) -> Dict[str, Any]:
    """Write the provenance record into the result file, above everything else.

    PDBQT readers skip REMARK records, so this does not change how the file
    docks, scores or parses, and it leaves the coordinate digest untouched. A
    sidecar JSON is written too, for machines, but the file itself is the
    authority: a sidecar can be separated from what it describes.
    """
    path = Path(path)
    body = path.read_text(encoding="utf-8")
    body = "\n".join(l for l in body.splitlines()
                     if not l.startswith(REMARK_PREFIX))
    header = "\n".join(_remark_lines(record))
    path.write_text(f"{header}\n{body}\n" if body else f"{header}\n",
                    encoding="utf-8")
    if sidecar:
        path.with_suffix(path.suffix + ".provenance.json").write_text(
            json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return record


def read_record(path: os.PathLike | str) -> Optional[Dict[str, str]]:
    """Read a stamped record back, or None if the file carries none."""
    rec: Dict[str, str] = {}
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith(REMARK_PREFIX):
                    rest = line[len(REMARK_PREFIX):].strip().split(" ", 1)
                    if rest:
                        rec[rest[0].lower()] = rest[1].strip() if len(rest) > 1 else ""
                elif line.startswith(("ATOM", "HETATM", "MODEL")):
                    break
    except OSError:
        return None
    return rec or None


def reserve_result_path(path: os.PathLike | str, overwrite: bool = False) -> Path:
    """Refuse to write over an existing result unless the caller says so.

    The defect this module exists to prevent was invisible partly because the
    docking scripts treated an existing file at the output path as their own
    output. Writing is now an explicit decision.
    """
    path = Path(path)
    if path.exists() and not overwrite:
        raise ResultPathExists(
            f"{path} already exists. Pass overwrite=True to replace it, or dock "
            f"to a path derived from the parameters you actually mean to use. "
            f"Adopting a file you did not write is how a result ends up under "
            f"the wrong protocol."
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def is_derived_split(path: os.PathLike | str) -> bool:
    """True for the per-mode files vina_split writes from a parent result."""
    return bool(SPLIT_RE.search(Path(path).name))


def audit(paths: Iterable[os.PathLike | str],
          baseline: Optional[Dict[str, Any]] = None,
          include_splits: bool = False) -> Dict[str, Any]:
    """Find result files that share coordinates, and say what that means.

    Returns a dict with:
      groups        every set of two or more files sharing a coordinate digest
      findings      the groups that are a problem, each classified
      baselined     the groups a baseline file records as known and accepted
      name_mismatch files whose own name disagrees with their recorded canonical
                    name, which means the path was chosen by hand after all

    Classification of a group:
      one_run_many_names     all files carry the same run id. One execution was
                             written under several names. Transparent, and
                             still worth reporting.
      contradiction          the files carry different run ids. Identical
                             coordinates cannot come from two executions of a
                             stochastic search, so one of the records is false.
      unprovenanced          at least one file carries no record, so the
                             duplication cannot be explained from the artefacts.
    """
    baseline = baseline or {}
    known = set(baseline.get("accepted_coordinate_digests", {}))

    by_digest: Dict[str, List[str]] = {}
    splits: List[str] = []
    examined: List[str] = []
    for p in paths:
        p = Path(p)
        if is_derived_split(p) and not include_splits:
            splits.append(str(p))
            continue
        d = coordinate_digest(p)
        if d:
            by_digest.setdefault(d, []).append(str(p))
            examined.append(str(p))

    groups, findings, baselined = [], [], []
    for digest, members in sorted(by_digest.items()):
        if len(members) < 2:
            continue
        members = sorted(members)
        records = {m: read_record(m) for m in members}
        run_ids = {(r or {}).get("run_id") for r in records.values()}
        if None in run_ids:
            kind = "unprovenanced"
        elif len(run_ids) == 1:
            kind = "one_run_many_names"
        else:
            kind = "contradiction"
        group = {"coordinate_digest": digest, "files": members, "kind": kind,
                 "run_ids": sorted(x for x in run_ids if x)}
        groups.append(group)
        if digest in known:
            group["baseline_note"] = baseline["accepted_coordinate_digests"][digest]
            baselined.append(group)
        else:
            findings.append(group)

    mismatches = []
    for p in sorted(set(examined)):
        rec = read_record(p)
        if rec and rec.get("canonical_name") and Path(p).name != rec["canonical_name"]:
            mismatches.append({"file": p, "recorded_canonical_name": rec["canonical_name"]})

    return {"groups": groups, "findings": findings, "baselined": baselined,
            "name_mismatch": mismatches, "derived_splits_skipped": len(splits),
            "files_examined": len(examined),
            "distinct_coordinate_sets": len(by_digest)}
