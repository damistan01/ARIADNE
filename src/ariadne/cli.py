"""Command-line interface for ariadne."""

from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

from . import config as config_mod
from . import ncmeta, parse, plan, runner, writer
from . import view as view_module

WORKER_TCL = Path(__file__).resolve().parent / "worker.tcl"

EXIT_OK = 0
EXIT_USAGE = 1
EXIT_PREFLIGHT = 2
EXIT_WORKER_FAILED = 3


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="ariadne",
        description="Batch electron-transfer Pathways analysis over an MD trajectory.",
    )
    parser.add_argument("-i", "--input", type=Path, required=True,
                        help="input file with the Pathways parameters")
    parser.add_argument("-p", "--prmtop", type=Path, required=True,
                        help="AMBER topology (parm7)")
    parser.add_argument("-y", "--trajectory", type=Path, required=True,
                        help="AMBER NetCDF trajectory")
    parser.add_argument("-o", "--output", type=Path, required=True,
                        help="output directory (created)")
    parser.add_argument("--frames", default=None,
                        help="frame selection start:stop:stride, stop inclusive")
    parser.add_argument("-n", "--workers", type=int,
                        default=min(8, os.cpu_count() or 1),
                        help="number of parallel VMD workers")
    parser.add_argument("--water", action="store_true",
                        help="keep water in the bridge selection")
    parser.add_argument("--resume", action="store_true",
                        help="continue an interrupted run in the output directory")
    parser.add_argument("--keep-raw", action="store_true",
                        help="keep the scratch directory after completion")
    parser.add_argument("--dry-run", action="store_true",
                        help="preflight and plan only; start no workers")
    return parser


def _fail(message: str, code: int) -> int:
    print(f"ariadne: {message}", file=sys.stderr)
    return code


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _vmd_version() -> str:
    try:
        out = subprocess.run(
            ["vmd", "-h"], capture_output=True, text=True, timeout=60,
            stdin=subprocess.DEVNULL,
        )
        for line in (out.stdout + out.stderr).splitlines():
            if "VMD for" in line:
                return line.strip()
    except (OSError, subprocess.SubprocessError):
        pass
    return "unknown"


def _read_da_indices(path: Path) -> tuple[set[int], set[int]]:
    """Read ref.da: 1-based donor indices on line 1, acceptor on line 2."""
    lines = path.read_text().splitlines()
    donor = {int(tok) - 1 for tok in lines[0].split()} if len(lines) > 0 else set()
    acceptor = {int(tok) - 1 for tok in lines[1].split()} if len(lines) > 1 else set()
    return donor, acceptor


def collect_results(
    scratch: Path,
    chunks,
    meta: ncmeta.TrajectoryMeta,
    donor_indices: set[int],
    acceptor_indices: set[int],
    atoms,
    covalent,
) -> list[parse.FrameRecord]:
    """Parse every chunk's artifacts into FrameRecords."""
    records: list[parse.FrameRecord] = []
    for chunk in chunks:
        chunkdir = Path(scratch) / f"chunk{chunk.index:02d}"
        frame_rows = runner.read_frames_tsv(chunkdir)
        for frame in chunk.frames:
            if frame not in frame_rows:
                continue
            n_atoms_pruned, status, _elapsed = frame_rows[frame]
            time_ps = meta.times[frame] if frame < len(meta.times) else None
            if status != "ok":
                records.append(parse.FrameRecord(
                    frame=frame, time_ps=time_ps, status=parse.STATUS_ERROR,
                    t_da=None, n_atoms_pruned=n_atoms_pruned,
                    message="worker reported an error; see worker.log",
                ))
                continue
            if n_atoms_pruned != len(atoms):
                raise parse.ParseError(
                    f"frame {frame}: the pruned subsystem has {n_atoms_pruned} atoms "
                    f"but the reference has {len(atoms)}. The bridge selection is not "
                    f"frame-invariant, so atom indices cannot be compared across "
                    f"frames. Use a bridge selection that does not depend on distance."
                )
            out_file = chunkdir / f"f{frame:06d}1.out"
            hb_file = chunkdir / f"f{frame:06d}.hb"
            out_text = out_file.read_text() if out_file.exists() else ""
            hbonds = parse.read_bond_pairs(hb_file) if hb_file.exists() else frozenset()
            try:
                records.append(parse.build_frame_record(
                    frame=frame, time_ps=time_ps, out_text=out_text, atoms=atoms,
                    covalent=covalent, hbonds=hbonds,
                    donor_indices=donor_indices, acceptor_indices=acceptor_indices,
                    n_atoms_pruned=n_atoms_pruned,
                ))
            except parse.ParseError as exc:
                records.append(parse.FrameRecord(
                    frame=frame, time_ps=time_ps, status=parse.STATUS_ERROR,
                    t_da=None, n_atoms_pruned=n_atoms_pruned, message=str(exc),
                ))
    return records


def main(argv: list[str] | None = None) -> int:
    # Dispatch subcommands by inspecting argv rather than converting to argparse
    # subparsers, which would change the existing `ariadne -i ... -o ...` form.
    argv = sys.argv[1:] if argv is None else argv
    if argv and argv[0] == "view":
        return view_module.main(argv[1:])

    args = build_arg_parser().parse_args(argv)

    # --- configuration -----------------------------------------------------
    try:
        params = config_mod.load_config(args.input)
        # The input file supplies the default frame selection so that one file
        # can describe a whole study; --frames on the command line overrides it.
        frames_spec = args.frames
        if frames_spec is None:
            frames_spec = config_mod.parse_frames_key(args.input.read_text())
    except (config_mod.ConfigError, OSError) as exc:
        return _fail(str(exc), EXIT_USAGE)

    if not args.prmtop.exists():
        return _fail(f"topology not found: {args.prmtop}", EXIT_PREFLIGHT)

    outdir = args.output
    couplings_csv = outdir / "couplings.csv"
    run_json = outdir / "run.json"

    # The clobber and resume guards run before the remaining preflight checks:
    # "this directory already holds a run" is the message most worth surfacing
    # first, since acting on it prevents data loss rather than merely failing.
    if couplings_csv.exists() and not args.resume:
        return _fail(
            f"{outdir} already holds a run; pass --resume to continue it "
            f"or choose another output directory",
            EXIT_PREFLIGHT,
        )

    if args.resume and run_json.exists():
        previous = json.loads(run_json.read_text()).get("params", {})
        current = params.as_dict()
        differences = [
            f"{key}: previous {previous[key]!r} vs now {current[key]!r}"
            for key in sorted(set(previous) & set(current))
            if previous.get(key) != current.get(key)
        ]
        if differences:
            return _fail(
                "cannot resume: parameters differ from the existing run:\n  "
                + "\n  ".join(differences),
                EXIT_PREFLIGHT,
            )

    if shutil.which("vmd") is None:
        return _fail("vmd is not on PATH", EXIT_PREFLIGHT)
    if shutil.which("pathcore") is None:
        return _fail("pathcore is not on PATH", EXIT_PREFLIGHT)
    if not args.trajectory.exists():
        return _fail(f"trajectory not found: {args.trajectory}", EXIT_PREFLIGHT)

    try:
        meta = ncmeta.read_meta(args.trajectory)
    except ncmeta.NetCDFError as exc:
        return _fail(str(exc), EXIT_PREFLIGHT)

    try:
        selected = plan.parse_frame_spec(frames_spec, meta.n_frames)
        statuses = writer.read_completed_statuses(couplings_csv) if args.resume else {}
        todo = plan.remaining_frames(selected, statuses)
    except plan.PlanError as exc:
        return _fail(str(exc), EXIT_USAGE)

    if not todo:
        print("Nothing to do: every selected frame is already computed.")
        return EXIT_OK

    # Frames computed by earlier invocations. A completed run deletes its
    # scratch directory, so the CSVs are the only surviving record of them;
    # rebuilding the outputs from this invocation's scratch alone would drop
    # every frame it did not recompute.
    preserved: list[parse.FrameRecord] = []
    if args.resume:
        recomputing = set(todo)
        preserved = [
            record
            for record in writer.read_records(couplings_csv, outdir / "paths.csv")
            if record.frame not in recomputing
        ]

    print(f"Trajectory: {meta.n_frames} frames, {meta.n_atoms} atoms")
    print(f"Selected:   {len(selected)} frames; {len(todo)} still to compute")

    # The probe frame is the lowest-numbered frame still to compute.
    probe_frame, rest = todo[0], todo[1:]
    outdir.mkdir(parents=True, exist_ok=True)
    scratch = outdir / "scratch"

    try:
        probe_chunks = plan.split_chunks([probe_frame], 1)
        rest_chunks = plan.split_chunks(rest, args.workers)
    except plan.PlanError as exc:
        return _fail(str(exc), EXIT_USAGE)

    if args.dry_run:
        print(f"Probe frame: {probe_frame}")
        for chunk in rest_chunks:
            print(f"  chunk{chunk.index:02d}: {len(chunk.frames)} frames "
                  f"({chunk.frames[0]}-{chunk.frames[-1]})")
        return EXIT_OK

    # --- probe -------------------------------------------------------------
    print(f"Probe frame {probe_frame} ...", flush=True)
    started = time.time()
    probe_results = runner.run_chunks(
        probe_chunks, scratch / "probe", args.prmtop, args.trajectory,
        params, args.water, WORKER_TCL,
    )
    probe_elapsed = time.time() - started
    probe = probe_results[0]
    if not probe.done:
        return _fail(
            f"the probe frame failed, so the full run was not started.\n"
            f"--- worker output ---\n{probe.log[-4000:]}",
            EXIT_PREFLIGHT,
        )
    print(f"Probe frame took {probe_elapsed:.1f}s")
    if rest:
        # The probe ran alone, so it does not see the memory contention that
        # slows concurrent workers -- measured at ~2x at 4 workers and ~6x at
        # 8 on this machine. Report a floor rather than a figure that reads
        # like a prediction and is wrong by 2-4x.
        n_waves = -(-len(rest) // max(1, args.workers))
        print(f"Remaining: {len(rest)} frames in {n_waves} wave(s) across "
              f"{args.workers} workers; at least "
              f"{probe_elapsed * n_waves / 60:.1f} min "
              f"(concurrent frames run slower than the probe)")

    probe_dir = scratch / "probe" / "chunk00"
    atoms = parse.read_pruned_pdb(probe_dir / "ref.pdb")
    covalent = parse.read_bond_pairs(probe_dir / "ref.cb")
    donor_indices, acceptor_indices = _read_da_indices(probe_dir / "ref.da")

    # --- full run ----------------------------------------------------------
    results = list(probe_results)
    latest_records: list[parse.FrameRecord] = []

    def flush(_result=None) -> None:
        nonlocal latest_records
        merged = preserved + collect_results(
            scratch / "probe", probe_chunks, meta,
            donor_indices, acceptor_indices, atoms, covalent,
        ) + collect_results(
            scratch / "rest", rest_chunks, meta,
            donor_indices, acceptor_indices, atoms, covalent,
        )
        latest_records = merged
        writer.write_couplings(couplings_csv, merged)
        writer.write_paths(outdir / "paths.csv", merged)
        writer.write_summary(outdir / "summary.txt", merged)

    interrupted = False
    if rest:
        try:
            results += runner.run_chunks(
                rest_chunks, scratch / "rest", args.prmtop, args.trajectory,
                params, args.water, WORKER_TCL, on_chunk_done=flush,
            )
        except KeyboardInterrupt:
            interrupted = True
            print("\nInterrupted; writing what has completed so far.", file=sys.stderr)

    flush()

    # --- provenance --------------------------------------------------------
    traj_stat = args.trajectory.stat()
    writer.write_run_json(run_json, {
        "tool": "ariadne",
        "written_at": datetime.datetime.now().astimezone().isoformat(),
        "params": params.as_dict(),
        "input_file": args.input.read_text(),
        "include_water": args.water,
        "effective_bridge": params.effective_bridge(args.water),
        "prmtop": {"path": str(args.prmtop.resolve()), "sha256": _sha256(args.prmtop)},
        "trajectory": {
            "path": str(args.trajectory.resolve()),
            "bytes": traj_stat.st_size,
            "mtime": datetime.datetime.fromtimestamp(traj_stat.st_mtime).isoformat(),
            "n_frames": meta.n_frames,
            "n_atoms": meta.n_atoms,
            "first_time_ps": meta.times[0] if meta.times else None,
            "last_time_ps": meta.times[-1] if meta.times else None,
        },
        "vmd": _vmd_version(),
        "plugin": "pathways 1.4",
        "pathcore": shutil.which("pathcore"),
        "frames_requested": len(selected),
        "frames_computed": len(todo),
        "workers": args.workers,
        "interrupted": interrupted,
    })

    failed = [r for r in results if not r.done]
    error_frames = [r for r in latest_records if r.status == parse.STATUS_ERROR]
    ok_frames = [r for r in latest_records if r.status == parse.STATUS_OK]

    # Scratch is the only record of *why* something failed, so keep it whenever
    # anything did -- including frames that errored while their worker exited
    # cleanly, which is the common case when a plugin call fails per frame.
    if not args.keep_raw and not failed and not error_frames and not interrupted:
        shutil.rmtree(scratch, ignore_errors=True)

    print(f"Wrote {outdir}/couplings.csv, paths.csv, run.json, summary.txt")
    if error_frames:
        print(f"{len(error_frames)} frame(s) failed; worker logs kept under "
              f"{scratch}", file=sys.stderr)
    if failed:
        for result in failed:
            print(f"chunk{result.chunk.index:02d} failed "
                  f"(exit {result.returncode}); see scratch/", file=sys.stderr)
        return EXIT_WORKER_FAILED
    if not ok_frames:
        print("no frame produced a pathway; the run yielded nothing usable",
              file=sys.stderr)
        return EXIT_WORKER_FAILED
    return EXIT_OK
