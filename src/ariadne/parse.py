"""Turn raw Pathways/pathcore artifacts into typed records.

All index arithmetic lives here. pathcore emits 1-based indices into the
pruned, hydrogen-stripped subsystem; everything downstream of this module is
0-based.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

__all__ = [
    "ParseError",
    "AtomInfo",
    "PathStep",
    "PathRecord",
    "FrameRecord",
    "read_pruned_pdb",
    "read_bond_pairs",
    "parse_pathcore_output",
    "trim_terminals",
    "classify_bond",
    "build_frame_record",
]

STATUS_OK = "ok"
STATUS_NO_PATH = "no_path"
STATUS_ERROR = "error"


class ParseError(ValueError):
    """Raised for malformed plugin or pathcore output."""


@dataclass(frozen=True)
class AtomInfo:
    index: int  # 0-based into the pruned subsystem
    resid: int
    resname: str
    atom_name: str
    segid: str


@dataclass(frozen=True)
class PathStep:
    step: int
    atom: AtomInfo
    bond_type: str  # "" | "covalent" | "hbond" | "through_space"


@dataclass(frozen=True)
class PathRecord:
    path_rank: int
    t_da: float
    steps: tuple[PathStep, ...]


@dataclass(frozen=True)
class FrameRecord:
    frame: int
    time_ps: float | None
    status: str
    t_da: float | None
    n_atoms_pruned: int
    paths: tuple[PathRecord, ...] = ()
    message: str = ""


def read_pruned_pdb(path: Path) -> tuple[AtomInfo, ...]:
    """Read the pruned subsystem PDB the plugin writes.

    Uses fixed PDB columns, not whitespace splitting: atom and residue names
    routinely contain no separating space.
    """
    atoms: list[AtomInfo] = []
    with open(path, "r") as fh:
        for lineno, line in enumerate(fh, start=1):
            if not line.startswith(("ATOM", "HETATM")):
                continue
            try:
                resid = int(line[22:26])
            except ValueError:
                raise ParseError(
                    f"{path} line {lineno}: cannot read residue id from columns 23-26"
                ) from None
            atoms.append(
                AtomInfo(
                    index=len(atoms),
                    resid=resid,
                    resname=line[17:21].strip(),
                    atom_name=line[12:16].strip(),
                    segid=line[72:76].strip(),
                )
            )
    if not atoms:
        raise ParseError(f"{path} contains no ATOM or HETATM records")
    return tuple(atoms)


def read_bond_pairs(path: Path) -> frozenset[tuple[int, int]]:
    """Read a .cb or .hb bond list, converting 1-based pairs to sorted 0-based."""
    pairs: set[tuple[int, int]] = set()
    with open(path, "r") as fh:
        for lineno, line in enumerate(fh, start=1):
            fields = line.split()
            if not fields:
                continue
            if len(fields) != 2:
                raise ParseError(
                    f"{path} line {lineno}: expected two indices, got {line.strip()!r}"
                )
            try:
                a, b = (int(f) - 1 for f in fields)
            except ValueError:
                raise ParseError(
                    f"{path} line {lineno}: non-integer index in {line.strip()!r}"
                ) from None
            pairs.add((a, b) if a <= b else (b, a))
    return frozenset(pairs)


def parse_pathcore_output(text: str) -> list[tuple[float, list[int]]]:
    """Parse pathcore stdout into (coupling, 0-based index list) tuples."""
    results: list[tuple[float, list[int]]] = []
    for lineno, line in enumerate(text.splitlines(), start=1):
        fields = line.split()
        if not fields:
            continue
        if len(fields) < 3:
            raise ParseError(
                f"pathcore output line {lineno}: expected a coupling and at least "
                f"two atom indices, got {line.strip()!r}"
            )
        try:
            coupling = float(fields[0])
        except ValueError:
            raise ParseError(
                f"pathcore output line {lineno}: coupling {fields[0]!r} is not a number"
            ) from None
        try:
            indices = [int(f) - 1 for f in fields[1:]]
        except ValueError:
            raise ParseError(
                f"pathcore output line {lineno}: non-integer atom index"
            ) from None
        results.append((coupling, indices))
    return results


def trim_terminals(
    indices: list[int],
    donor_indices: set[int],
    acceptor_indices: set[int],
) -> list[int]:
    """Drop redundant donor and acceptor atoms, as line2path does.

    Keeps only the last donor atom before the path leaves the donor, and only
    the first acceptor atom after it enters the acceptor.
    Mirrors pathways.tcl:298-324.
    """
    if len(indices) < 2:
        return list(indices)

    start = 0
    while start + 1 < len(indices) and indices[start + 1] in donor_indices:
        start += 1

    end = len(indices) - 1
    while end - 1 > start and indices[end - 1] in acceptor_indices:
        end -= 1

    return list(indices[start : end + 1])


def classify_bond(
    a: int,
    b: int,
    covalent: frozenset[tuple[int, int]],
    hbonds: frozenset[tuple[int, int]],
) -> str:
    pair = (a, b) if a <= b else (b, a)
    if pair in covalent:
        return "covalent"
    if pair in hbonds:
        return "hbond"
    return "through_space"


def build_frame_record(
    frame: int,
    time_ps: float | None,
    out_text: str,
    atoms: tuple[AtomInfo, ...],
    covalent: frozenset[tuple[int, int]],
    hbonds: frozenset[tuple[int, int]],
    donor_indices: set[int],
    acceptor_indices: set[int],
    n_atoms_pruned: int,
) -> FrameRecord:
    raw_paths = parse_pathcore_output(out_text)
    if not raw_paths:
        return FrameRecord(
            frame=frame,
            time_ps=time_ps,
            status=STATUS_NO_PATH,
            t_da=None,
            n_atoms_pruned=n_atoms_pruned,
        )

    raw_paths.sort(key=lambda item: item[0], reverse=True)

    paths: list[PathRecord] = []
    for rank, (coupling, indices) in enumerate(raw_paths):
        for idx in indices:
            if not 0 <= idx < len(atoms):
                raise ParseError(
                    f"frame {frame}: atom index {idx + 1} is out of range for a "
                    f"pruned subsystem of {len(atoms)} atoms"
                )
        trimmed = trim_terminals(indices, donor_indices, acceptor_indices)
        steps = [PathStep(step=0, atom=atoms[trimmed[0]], bond_type="")]
        for position, idx in enumerate(trimmed[1:], start=1):
            steps.append(
                PathStep(
                    step=position,
                    atom=atoms[idx],
                    bond_type=classify_bond(
                        trimmed[position - 1], idx, covalent, hbonds
                    ),
                )
            )
        paths.append(PathRecord(path_rank=rank, t_da=coupling, steps=tuple(steps)))

    return FrameRecord(
        frame=frame,
        time_ps=time_ps,
        status=STATUS_OK,
        t_da=paths[0].t_da,
        n_atoms_pruned=n_atoms_pruned,
        paths=tuple(paths),
    )
