"""Spawn and supervise the headless VMD workers."""

from __future__ import annotations

import os
import subprocess
import time
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from pathlib import Path

from .config import PathwaysParams
from .plan import Chunk

__all__ = [
    "RunnerError",
    "WorkerResult",
    "worker_env",
    "build_worker_command",
    "write_params_tcl",
    "read_frames_tsv",
    "run_chunks",
]

_POLL_SECONDS = 0.5


class RunnerError(RuntimeError):
    """Raised when a worker cannot be started or fails irrecoverably."""


@dataclass(frozen=True)
class WorkerResult:
    chunk: Chunk
    returncode: int
    done: bool
    log: str


def worker_env() -> dict[str, str]:
    """Environment for a worker: the current one minus DISPLAY.

    Belt and braces only. `-dispdev text` already strips DISPLAY from VMD's
    Tcl env array, and that is what actually forces the plugin's headless
    path. Removing it here keeps the worker's own guard meaningful if the
    command is ever changed.
    """
    return {key: value for key, value in os.environ.items() if key != "DISPLAY"}


def build_worker_command(worker_tcl: Path, chunkdir: Path) -> list[str]:
    return ["vmd", "-dispdev", "text", "-e", str(worker_tcl), "-args", str(chunkdir)]


def write_params_tcl(
    chunkdir: Path,
    prmtop: Path,
    traj: Path,
    chunk: Chunk,
    params: PathwaysParams,
    include_water: bool,
) -> Path:
    """Generate the Tcl file the worker sources.

    Values are wrapped in braces so Tcl performs no substitution on them;
    config.py has already rejected selections containing brace, bracket,
    backslash, dollar or quote characters.
    """
    frames = " ".join(str(f) for f in chunk.frames)
    plugin_args = " ".join(params.to_plugin_args())
    text = (
        f"set PT_PRMTOP {{{prmtop}}}\n"
        f"set PT_NC {{{traj}}}\n"
        f"set PT_CHUNKDIR {{{chunkdir}}}\n"
        f"set PT_FRAMES [list {frames}]\n"
        f"set PT_DONOR {{{params.donor}}}\n"
        f"set PT_ACCEPTOR {{{params.acceptor}}}\n"
        f"set PT_BRIDGE {{{params.effective_bridge(include_water)}}}\n"
        f"set PT_PARGS [list {plugin_args}]\n"
    )
    target = Path(chunkdir) / "params.tcl"
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(text)
    return target


def read_frames_tsv(chunkdir: Path) -> dict[int, tuple[int, str, float]]:
    """Read a worker's per-frame log: frame -> (n_atoms_pruned, status, seconds).

    Tolerates a truncated final line, which a killed worker can leave behind.
    """
    path = Path(chunkdir) / "frames.tsv"
    if not path.exists():
        return {}
    rows: dict[int, tuple[int, str, float]] = {}
    for line in path.read_text().splitlines():
        fields = line.split("\t")
        if len(fields) != 4:
            continue
        try:
            rows[int(fields[0])] = (int(fields[1]), fields[2], float(fields[3]))
        except ValueError:
            continue
    return rows


def run_chunks(
    chunks: Sequence[Chunk],
    scratch: Path,
    prmtop: Path,
    traj: Path,
    params: PathwaysParams,
    include_water: bool,
    worker_tcl: Path,
    on_chunk_done: Callable[[WorkerResult], None] | None = None,
) -> list[WorkerResult]:
    """Run every chunk concurrently and wait for all of them.

    Processes are polled rather than waited on in order, so `on_chunk_done`
    fires as soon as each chunk finishes and the caller can flush partial
    output promptly.
    """
    scratch = Path(scratch)
    env = worker_env()
    pending: list[tuple[Chunk, Path, subprocess.Popen]] = []

    for chunk in chunks:
        chunkdir = scratch / f"chunk{chunk.index:02d}"
        chunkdir.mkdir(parents=True, exist_ok=True)
        write_params_tcl(chunkdir, prmtop, traj, chunk, params, include_water)
        log = open(chunkdir / "worker.log", "w")
        try:
            proc = subprocess.Popen(
                build_worker_command(worker_tcl, chunkdir),
                stdout=log,
                stderr=subprocess.STDOUT,
                stdin=subprocess.DEVNULL,
                env=env,
            )
        except OSError as exc:
            log.close()
            for _, _, started in pending:
                started.terminate()
            raise RunnerError(f"cannot start vmd: {exc}") from None
        pending.append((chunk, chunkdir, proc))

    results: list[WorkerResult] = []
    try:
        while pending:
            still_running = []
            for chunk, chunkdir, proc in pending:
                returncode = proc.poll()
                if returncode is None:
                    still_running.append((chunk, chunkdir, proc))
                    continue
                result = WorkerResult(
                    chunk=chunk,
                    returncode=returncode,
                    done=(chunkdir / "DONE").exists(),
                    log=(chunkdir / "worker.log").read_text(errors="replace"),
                )
                results.append(result)
                if on_chunk_done is not None:
                    on_chunk_done(result)
            pending = still_running
            if pending:
                time.sleep(_POLL_SECONDS)
    except KeyboardInterrupt:
        for _, _, proc in pending:
            if proc.poll() is None:
                proc.terminate()
        for _, _, proc in pending:
            try:
                proc.wait(timeout=30)
            except subprocess.TimeoutExpired:
                proc.kill()
        raise

    results.sort(key=lambda r: r.chunk.index)
    return results
