"""Read paths.csv back into typed rows for the viewer."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

__all__ = ["PathsError", "PathStepRow", "read_paths_csv", "group_paths"]

_REQUIRED = (
    "frame",
    "path_rank",
    "t_da",
    "step",
    "resid",
    "resname",
    "atom_name",
    "segid",
    "bond_type",
)


class PathsError(ValueError):
    """Raised for an unreadable or malformed paths.csv."""


@dataclass(frozen=True)
class PathStepRow:
    frame: int
    path_rank: int
    t_da: float
    step: int
    resid: int
    resname: str
    atom_name: str
    segid: str
    bond_type: str


def read_paths_csv(path: Path) -> list[PathStepRow]:
    path = Path(path)
    try:
        handle = path.open(newline="")
    except OSError as exc:
        raise PathsError(f"cannot read {path}: {exc}") from None

    rows: list[PathStepRow] = []
    with handle:
        reader = csv.DictReader(handle)
        missing = [c for c in _REQUIRED if c not in (reader.fieldnames or [])]
        if missing:
            raise PathsError(f"{path}: missing column(s): {', '.join(missing)}")
        for lineno, row in enumerate(reader, start=2):
            try:
                rows.append(
                    PathStepRow(
                        frame=int(row["frame"]),
                        path_rank=int(row["path_rank"]),
                        t_da=float(row["t_da"]),
                        step=int(row["step"]),
                        resid=int(row["resid"]),
                        resname=row["resname"].strip(),
                        atom_name=row["atom_name"].strip(),
                        segid=row["segid"].strip(),
                        bond_type=row["bond_type"].strip(),
                    )
                )
            except (ValueError, AttributeError):
                raise PathsError(
                    f"{path} line {lineno}: could not parse row {row!r}"
                ) from None
    return rows


def group_paths(
    rows: list[PathStepRow],
) -> dict[int, dict[int, list[PathStepRow]]]:
    """Group rows as frame -> rank -> steps, each step list ordered by step."""
    grouped: dict[int, dict[int, list[PathStepRow]]] = {}
    for row in rows:
        grouped.setdefault(row.frame, {}).setdefault(row.path_rank, []).append(row)
    for ranks in grouped.values():
        for steps in ranks.values():
            steps.sort(key=lambda s: s.step)
    return grouped
