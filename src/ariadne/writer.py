"""Write ariadne output files.

All writes are atomic: a temp file in the destination directory followed by
os.replace, so an interrupted run never leaves a half-written CSV.
"""

from __future__ import annotations

import csv
import io
import json
import math
import os
import statistics
import tempfile
from collections.abc import Iterable, Sequence
from pathlib import Path

from .parse import STATUS_OK, AtomInfo, FrameRecord, PathRecord, PathStep

__all__ = [
    "COUPLINGS_HEADER",
    "PATHS_HEADER",
    "atomic_write_text",
    "write_couplings",
    "write_paths",
    "write_run_json",
    "coupling_stats",
    "write_summary",
    "read_completed_statuses",
    "read_records",
]

COUPLINGS_HEADER = [
    "frame",
    "time_ps",
    "t_da",
    "log10_t_da",
    "n_paths",
    "n_steps",
    "status",
]

PATHS_HEADER = [
    "frame",
    "path_rank",
    "t_da",
    "step",
    "atom_index",
    "resid",
    "resname",
    "atom_name",
    "segid",
    "bond_type",
]


def atomic_write_text(path: Path, text: str) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(dir=str(path.parent), prefix=f".{path.name}.")
    try:
        with os.fdopen(fd, "w", newline="") as fh:
            fh.write(text)
        os.replace(tmp_name, path)
    except BaseException:
        if os.path.exists(tmp_name):
            os.unlink(tmp_name)
        raise


def _render_csv(header: Sequence[str], rows: Iterable[Sequence[object]]) -> str:
    buffer = io.StringIO(newline="")
    writer = csv.writer(buffer, lineterminator="\n")
    writer.writerow(header)
    writer.writerows(rows)
    return buffer.getvalue()


def _format_float(value: float | None) -> str:
    return "" if value is None else repr(value)


def write_couplings(path: Path, records: Iterable[FrameRecord]) -> None:
    rows = []
    for record in sorted(records, key=lambda r: r.frame):
        if record.t_da is not None and record.t_da > 0:
            log10 = repr(math.log10(record.t_da))
        else:
            log10 = ""
        n_steps = len(record.paths[0].steps) if record.paths else 0
        rows.append(
            [
                record.frame,
                _format_float(record.time_ps),
                _format_float(record.t_da),
                log10,
                len(record.paths),
                n_steps,
                record.status,
            ]
        )
    atomic_write_text(path, _render_csv(COUPLINGS_HEADER, rows))


def write_paths(path: Path, records: Iterable[FrameRecord]) -> None:
    rows = []
    for record in sorted(records, key=lambda r: r.frame):
        for path_record in record.paths:
            for step in path_record.steps:
                rows.append(
                    [
                        record.frame,
                        path_record.path_rank,
                        _format_float(path_record.t_da),
                        step.step,
                        step.atom.index,
                        step.atom.resid,
                        step.atom.resname,
                        step.atom.atom_name,
                        step.atom.segid,
                        step.bond_type,
                    ]
                )
    atomic_write_text(path, _render_csv(PATHS_HEADER, rows))


def write_run_json(path: Path, payload: dict) -> None:
    atomic_write_text(path, json.dumps(payload, indent=2, sort_keys=True) + "\n")


def coupling_stats(records: Iterable[FrameRecord]) -> dict:
    couplings = [
        r.t_da for r in records if r.status == STATUS_OK and r.t_da is not None
    ]
    if not couplings:
        return {
            "n_ok": 0,
            "mean": None,
            "stddev": None,
            "mean_sq": None,
            "coherence": None,
        }
    mean = statistics.fmean(couplings)
    stddev = statistics.pstdev(couplings)
    mean_sq = mean * mean + stddev * stddev
    coherence = (mean * mean / mean_sq) if mean_sq else None
    return {
        "n_ok": len(couplings),
        "mean": mean,
        "stddev": stddev,
        "mean_sq": mean_sq,
        "coherence": coherence,
    }


def write_summary(path: Path, records: Sequence[FrameRecord]) -> None:
    stats = coupling_stats(records)
    counts: dict[str, int] = {}
    for record in records:
        counts[record.status] = counts.get(record.status, 0) + 1

    lines = ["PATHWAY COUPLING STATISTICS", ""]
    if stats["n_ok"]:
        lines += [
            f"  <T_DA>:        {stats['mean']!r}",
            f"  sigma(T_DA):   {stats['stddev']!r}",
            f"  <T_DA^2>:      {stats['mean_sq']!r}",
            f"  coherence:     {stats['coherence']!r}",
        ]
    else:
        lines.append("  no successful frames")
    lines += ["", "FRAMES", ""]
    for status in sorted(counts):
        lines.append(f"  {status}: {counts[status]}")
    atomic_write_text(path, "\n".join(lines) + "\n")


def _read_float(text: str) -> float | None:
    return None if text == "" else float(text)


def read_records(couplings_path: Path, paths_path: Path) -> list[FrameRecord]:
    """Reconstruct FrameRecords from the CSVs this module writes.

    --resume needs the frames computed by earlier invocations. A completed run
    deletes its scratch directory, so couplings.csv and paths.csv are the only
    surviving record of them; without reading them back, a resume that widens
    the frame range would rewrite the CSVs from this invocation's scratch alone
    and silently drop everything computed before.

    n_atoms_pruned is not a CSV column. It exists only as a consistency check
    while collecting a worker's output and is never written, so restored
    records carry 0 and no output is affected.
    """
    couplings_path = Path(couplings_path)
    if not couplings_path.exists():
        return []

    steps_by_frame: dict[int, dict[int, tuple[float, list[PathStep]]]] = {}
    paths_path = Path(paths_path)
    if paths_path.exists():
        with paths_path.open(newline="") as fh:
            for row in csv.DictReader(fh):
                frame = int(row["frame"])
                rank = int(row["path_rank"])
                ranks = steps_by_frame.setdefault(frame, {})
                _, steps = ranks.setdefault(rank, (float(row["t_da"]), []))
                steps.append(
                    PathStep(
                        step=int(row["step"]),
                        atom=AtomInfo(
                            index=int(row["atom_index"]),
                            resid=int(row["resid"]),
                            resname=row["resname"],
                            atom_name=row["atom_name"],
                            segid=row["segid"],
                        ),
                        bond_type=row["bond_type"],
                    )
                )

    records: list[FrameRecord] = []
    with couplings_path.open(newline="") as fh:
        for row in csv.DictReader(fh):
            frame = int(row["frame"])
            ranks = steps_by_frame.get(frame, {})
            paths = tuple(
                PathRecord(
                    path_rank=rank,
                    t_da=t_da,
                    steps=tuple(sorted(steps, key=lambda s: s.step)),
                )
                for rank, (t_da, steps) in sorted(ranks.items())
            )
            records.append(
                FrameRecord(
                    frame=frame,
                    time_ps=_read_float(row["time_ps"]),
                    status=row["status"],
                    t_da=_read_float(row["t_da"]),
                    n_atoms_pruned=0,
                    paths=paths,
                )
            )
    records.sort(key=lambda r: r.frame)
    return records


def read_completed_statuses(path: Path) -> dict[int, str]:
    """Map frame number -> status from an existing couplings.csv."""
    path = Path(path)
    if not path.exists():
        return {}
    statuses: dict[int, str] = {}
    with path.open(newline="") as fh:
        for row in csv.DictReader(fh):
            statuses[int(row["frame"])] = row["status"]
    return statuses
