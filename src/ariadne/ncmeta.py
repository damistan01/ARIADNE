"""Minimal NetCDF-3 reader for AMBER trajectory metadata.

Reads only what ariadne needs: frame count, atom count and the per-frame
simulation time. Deliberately stdlib-only so the tool has no scientific-stack
dependency.
"""

from __future__ import annotations

import struct
from dataclasses import dataclass
from pathlib import Path

__all__ = ["NetCDFError", "TrajectoryMeta", "read_meta"]


class NetCDFError(ValueError):
    """Raised for unreadable, truncated or unsupported trajectory files."""


_NC_FLOAT = 5
_NC_DOUBLE = 6
_TYPE_SIZE = {1: 1, 2: 1, 3: 2, 4: 4, 5: 4, 6: 8}


@dataclass(frozen=True)
class TrajectoryMeta:
    n_frames: int
    n_atoms: int
    times: tuple[float, ...]


class _Header:
    def __init__(self, fh):
        self._fh = fh

    def int32(self) -> int:
        raw = self._fh.read(4)
        if len(raw) < 4:
            raise NetCDFError("unexpected end of file while reading header")
        return struct.unpack(">i", raw)[0]

    def offset(self, size: int) -> int:
        raw = self._fh.read(size)
        if len(raw) < size:
            raise NetCDFError("unexpected end of file while reading header")
        return struct.unpack(">q" if size == 8 else ">i", raw)[0]

    def name(self) -> str:
        length = self.int32()
        raw = self._fh.read(length)
        if len(raw) < length:
            raise NetCDFError("unexpected end of file while reading a name")
        self._fh.read((4 - length % 4) % 4)  # padding
        return raw.decode("ascii", errors="replace")

    def skip_attributes(self) -> None:
        self.int32()  # tag: 0 (absent) or 0x0C (NC_ATTRIBUTE)
        for _ in range(self.int32()):
            self.name()
            nc_type = self.int32()
            nvals = self.int32()
            if nc_type not in _TYPE_SIZE:
                raise NetCDFError(f"unsupported attribute type {nc_type}")
            nbytes = _TYPE_SIZE[nc_type] * nvals
            self._fh.read(nbytes + (4 - nbytes % 4) % 4)


def read_meta(path: Path) -> TrajectoryMeta:
    path = Path(path)
    with open(path, "rb") as fh:
        magic = fh.read(4)
        if len(magic) < 4 or magic[:3] != b"CDF" or magic[3] not in (1, 2):
            raise NetCDFError(f"{path} is not a NetCDF-3 classic or 64-bit file")
        offset_size = 8 if magic[3] == 2 else 4

        head = _Header(fh)
        n_frames = head.int32()
        if n_frames < 0:
            raise NetCDFError("streaming record count is not supported")

        head.int32()  # dimension list tag
        dims: list[tuple[str, int]] = []
        for _ in range(head.int32()):
            dims.append((head.name(), head.int32()))

        head.skip_attributes()  # global attributes

        head.int32()  # variable list tag
        recsize = 0
        time_var: tuple[int, int] | None = None  # (nc_type, begin)
        for _ in range(head.int32()):
            name = head.name()
            dimids = [head.int32() for _ in range(head.int32())]
            head.skip_attributes()
            nc_type = head.int32()
            vsize = head.int32()
            begin = head.offset(offset_size)
            is_record = bool(dimids) and dims[dimids[0]][1] == 0
            if is_record:
                recsize += vsize
                if name == "time":
                    time_var = (nc_type, begin)

        atom_dims = [length for name, length in dims if name == "atom"]
        if not atom_dims:
            raise NetCDFError(
                f"{path} has no 'atom' dimension; not an AMBER trajectory"
            )
        n_atoms = atom_dims[0]

        times = _read_times(fh, path, time_var, recsize, n_frames)

    return TrajectoryMeta(n_frames=n_frames, n_atoms=n_atoms, times=times)


def _read_times(fh, path, time_var, recsize, n_frames) -> tuple[float, ...]:
    if time_var is None:
        return ()
    nc_type, begin = time_var
    if nc_type == _NC_FLOAT:
        fmt, width = ">f", 4
    elif nc_type == _NC_DOUBLE:
        fmt, width = ">d", 8
    else:
        raise NetCDFError(f"unsupported type {nc_type} for the 'time' variable")

    times = []
    for frame in range(n_frames):
        fh.seek(begin + frame * recsize)
        raw = fh.read(width)
        if len(raw) < width:
            raise NetCDFError(
                f"{path} is truncated: header declares {n_frames} frames but "
                f"data ends during frame {frame}"
            )
        times.append(struct.unpack(fmt, raw)[0])
    return tuple(times)
