"""Decide which frames to compute and how to divide them between workers."""

from __future__ import annotations

from dataclasses import dataclass

__all__ = [
    "PlanError",
    "Chunk",
    "parse_frame_spec",
    "remaining_frames",
    "split_chunks",
]


class PlanError(ValueError):
    """Raised for an invalid frame selection or worker count."""


@dataclass(frozen=True)
class Chunk:
    index: int
    frames: tuple[int, ...]


def parse_frame_spec(spec: str | None, n_frames: int) -> list[int]:
    """Expand 'start:stop:stride' into an explicit frame list.

    stop is inclusive, matching the plugin's -first/-last convention.
    """
    if spec is None or spec.strip() == "":
        return list(range(n_frames))

    parts = spec.split(":")
    if len(parts) > 3:
        raise PlanError(f"frame range {spec!r} has too many ':'-separated fields")

    try:
        numbers = [int(p) if p.strip() != "" else None for p in parts]
    except ValueError:
        raise PlanError(
            f"frame range {spec!r} must be numbers of the form start:stop:stride"
        ) from None

    if len(numbers) == 1:
        start = numbers[0]
        stop = start
        stride = 1
    else:
        start = numbers[0] if numbers[0] is not None else 0
        stop = numbers[1] if numbers[1] is not None else n_frames - 1
        stride = numbers[2] if len(numbers) == 3 and numbers[2] is not None else 1

    if start is None:
        start = 0
    if start < 0 or stop < 0:
        raise PlanError(f"frame range {spec!r} must not contain negative frames")
    if stride < 1:
        raise PlanError(f"frame range {spec!r}: stride must be >= 1")
    if start > stop:
        raise PlanError(f"frame range {spec!r}: start comes after stop")
    if stop > n_frames - 1:
        raise PlanError(
            f"frame range {spec!r} exceeds the trajectory, which has "
            f"{n_frames} frames (0-{n_frames - 1})"
        )
    return list(range(start, stop + 1, stride))


def remaining_frames(
    selected: list[int],
    statuses: dict[int, str],
    retry_errors: bool = True,
) -> list[int]:
    """Frames still to compute, given an existing run's frame statuses."""
    settled = {"ok", "no_path"} if retry_errors else {"ok", "no_path", "error"}
    done = {frame for frame, status in statuses.items() if status in settled}
    return [frame for frame in selected if frame not in done]


def split_chunks(frames: list[int], n_workers: int) -> list[Chunk]:
    """Divide frames into at most n_workers contiguous runs of the input order.

    Contiguity matters: each worker loads one slice of the trajectory, so
    interleaving frames would force every worker to load the whole file.
    """
    if n_workers < 1:
        raise PlanError("worker count must be at least 1")
    if not frames:
        return []

    n_chunks = min(n_workers, len(frames))
    base, extra = divmod(len(frames), n_chunks)

    chunks: list[Chunk] = []
    cursor = 0
    for index in range(n_chunks):
        size = base + (1 if index < extra else 0)
        chunks.append(Chunk(index=index, frames=tuple(frames[cursor : cursor + size])))
        cursor += size
    return chunks
