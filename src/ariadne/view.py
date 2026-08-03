"""The `ariadne view` subcommand: 3D pathway viewer."""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path

from . import pathsio, viewconfig, writer
from .config import ConfigError, _check_selection

VIEW_TCL = Path(__file__).resolve().parent / "view.tcl"

EXIT_OK = 0
EXIT_USAGE = 1
# VMD's own exit always reports 0, so view.tcl signals failure with SIGKILL.
# subprocess reports that as a negative returncode (-9), which the shell would
# turn into 247; normalise every signal death to one ordinary non-zero code.
EXIT_VMD_FAILED = 3

# Statuses that mean "this frame was computed", as opposed to never attempted.
# no_path and error frames are rendered as structure with no path rather than
# reported as missing, which would be a lie.
_COMPUTED_STATUSES = ("ok", "no_path", "error")

# Written by view.tcl's export_scene alongside view.dcd; read back by
# --from-export to prove the bundle holds the frames being rendered.
EXPORT_FRAMES = "frames.txt"

__all__ = [
    "ViewError",
    "RunInfo",
    "discover_run",
    "render_viewparams",
    "read_export_frames",
    "check_export_frames",
    "main",
]


class ViewError(ValueError):
    """Raised when a run directory cannot be turned into a viewable scene."""


@dataclass(frozen=True)
class RunInfo:
    rundir: Path
    prmtop: Path
    trajectory: Path
    donor: str
    acceptor: str
    bridge: str
    computed_frames: list[int]
    input_text: str
    frame_status: dict[int, str] = field(default_factory=dict)


def discover_run(rundir: Path, require_raw: bool = True) -> RunInfo:
    """Read a run directory into a viewable scene description.

    require_raw is False for --from-export, whose whole point is that the
    topology and the 2.1 GB trajectory need not be present on this machine.
    """
    rundir = Path(rundir)
    run_json = rundir / "run.json"
    if not run_json.exists():
        raise ViewError(
            f"{rundir} has no run.json; it does not look like an ARIADNE run directory"
        )
    try:
        payload = json.loads(run_json.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise ViewError(f"cannot read {run_json}: {exc}") from None

    try:
        params = payload["params"]
        prmtop = Path(payload["prmtop"]["path"])
        trajectory = Path(payload["trajectory"]["path"])
        donor = params["donor"]
        acceptor = params["acceptor"]
        bridge = payload.get("effective_bridge", "all")
    except (KeyError, TypeError) as exc:
        raise ViewError(f"{run_json} is missing expected key: {exc}") from None

    # run.json is read from disk rather than produced by parse_config, so the
    # "no brace-quoting-breaking characters" invariant that config._check_selection
    # normally enforces on donor/acceptor/bridge cannot be assumed here: a
    # hand-edited, corrupted, or foreign-tool-generated run.json could smuggle
    # {}[]\$" straight into the brace-quoted Tcl that render_viewparams emits.
    # Reuse config._check_selection (rather than a second copy of the
    # forbidden-character set) to reject that before it reaches render_viewparams.
    for key, value in (("donor", donor), ("acceptor", acceptor), ("bridge", bridge)):
        try:
            _check_selection(key, value)
        except ConfigError as exc:
            raise ViewError(f"{run_json}: {exc}") from None

    # Only checked when the raw files are actually going to be read: with
    # --from-export the scene comes out of view/, and demanding the topology
    # and trajectory here would make the bundle useless on the machine it was
    # handed to -- which is the one place it matters.
    if require_raw:
        for label, path in (("topology", prmtop), ("trajectory", trajectory)):
            if not path.exists():
                raise ViewError(
                    f"{label} recorded in run.json no longer exists: {path}\n"
                    f"pass -p/-y to override, or use --from-export"
                )

    statuses = writer.read_completed_statuses(rundir / "couplings.csv")
    computed = sorted(
        f for f, status in statuses.items() if status in _COMPUTED_STATUSES
    )

    return RunInfo(
        rundir=rundir,
        prmtop=prmtop,
        trajectory=trajectory,
        donor=donor,
        acceptor=acceptor,
        bridge=bridge,
        computed_frames=computed,
        input_text=payload.get("input_file", ""),
        frame_status=dict(statuses),
    )


def _tcl_element(value: str) -> str:
    """Render one field of a step entry, keeping empty values as {}."""
    return value if value else "{}"


def render_viewparams(
    run: RunInfo,
    frames: list[int],
    grouped: dict[int, dict[int, list[pathsio.PathStepRow]]],
    params: viewconfig.ViewParams,
    export: bool,
    from_export: bool,
) -> str:
    lines = [
        f"set VIEW_RUNDIR {{{run.rundir.resolve()}}}",
        f"set VIEW_PRMTOP {{{run.prmtop}}}",
        f"set VIEW_TRAJ {{{run.trajectory}}}",
        f"set VIEW_DONOR {{{run.donor}}}",
        f"set VIEW_ACCEPTOR {{{run.acceptor}}}",
        f"set VIEW_BRIDGE {{{run.bridge}}}",
        f"set VIEW_FRAMES [list {' '.join(str(f) for f in frames)}]",
        f"set VIEW_RADIUS {params.view_radius}",
        f"set VIEW_RANKS {params.view_ranks}",
        # Colours resolved to VMD indices here, because draw_step passes its
        # colour argument straight to `graphics ... color`, which wants an index.
        f"set VIEW_COL_COVALENT {params.color_index('covalent')}",
        f"set VIEW_COL_HBOND {params.color_index('hbond')}",
        f"set VIEW_COL_THROUGH_SPACE {params.color_index('through_space')}",
        f"set VIEW_CONTEXT_RADIUS {params.view_context_radius}",
        f"set VIEW_CLIP_FRONT {1 if params.view_clip_front else 0}",
        f"set VIEW_EXPORT {1 if export else 0}",
        f"set VIEW_FROMEXPORT {1 if from_export else 0}",
    ]
    for frame in frames:
        lines.append(
            f"set VIEW_STATUS({frame}) {{{run.frame_status.get(frame, 'unknown')}}}"
        )
        ranks = grouped.get(frame, {})
        wanted = [0] if params.view_ranks == "0" else sorted(ranks)
        present = [r for r in wanted if r in ranks]
        if not present:
            # No VIEW_RANKLIST entry at all, so view.tcl takes its "no path"
            # branch and labels the frame instead of silently drawing nothing.
            continue
        lines.append(f"set VIEW_RANKLIST({frame}) [list {' '.join(map(str, present))}]")
        for rank in present:
            entries = " ".join(
                "{"
                + " ".join(
                    (
                        str(step.resid),
                        _tcl_element(step.resname),
                        _tcl_element(step.atom_name),
                        _tcl_element(step.segid),
                        _tcl_element(step.bond_type),
                    )
                )
                + "}"
                for step in ranks[rank]
            )
            lines.append(f"set VIEW_PATH({frame},{rank}) [list {entries}]")
    return "\n".join(lines) + "\n"


def read_export_frames(rundir: Path) -> list[int]:
    """Frame numbers recorded in the bundle, in the order they were written."""
    path = Path(rundir) / "view" / EXPORT_FRAMES
    try:
        text = path.read_text()
    except OSError:
        raise ViewError(
            f"{path} is missing, so the frames stored in view/view.dcd cannot be "
            f"identified. Re-run with --export to rebuild the bundle."
        ) from None
    try:
        return [int(tok) for tok in text.split()]
    except ValueError:
        raise ViewError(f"{path} is not a list of frame numbers") from None


def check_export_frames(rundir: Path, frames: list[int]) -> None:
    """Refuse to draw bundle frame i with another frame's pathway.

    view.tcl loads whatever is in view.dcd, in the order it was exported, while
    VIEW_FRAMES and VIEW_PATH come from the current invocation. If the two
    lists differ -- export at 0:999:10, reload at 0:999:20 -- every cylinder is
    drawn on the wrong coordinates, and nothing about the picture says so.
    """
    exported = read_export_frames(rundir)
    if exported == frames:
        return
    raise ViewError(
        "the exported bundle does not hold the frames being rendered, so the "
        "pathways would be drawn on the wrong coordinates.\n"
        f"  view/{EXPORT_FRAMES} (in view.dcd): {_summarise_frames(exported)}\n"
        f"  requested now:                      {_summarise_frames(frames)}\n"
        "Re-export with --export, or pass --frames matching the bundle."
    )


def _summarise_frames(frames: list[int]) -> str:
    shown = ", ".join(str(f) for f in frames[:12])
    more = "" if len(frames) <= 12 else f", ... (+{len(frames) - 12} more)"
    return f"{len(frames)} frame(s) [{shown}{more}]"


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="ariadne view",
        description="Show ARIADNE pathways in 3D across the analysed frames.",
    )
    parser.add_argument("rundir", type=Path, help="an ARIADNE run directory")
    parser.add_argument("-p", "--prmtop", type=Path, default=None,
                        help="override the topology recorded in run.json")
    parser.add_argument("-y", "--trajectory", type=Path, default=None,
                        help="override the trajectory recorded in run.json")
    parser.add_argument("--frames", default=None,
                        help="override view_frames for this invocation")
    parser.add_argument("--export", action="store_true",
                        help="also write a self-contained view/ bundle")
    parser.add_argument("--from-export", action="store_true",
                        help="load view/ instead of the topology and trajectory")
    parser.add_argument("--text", action="store_true",
                        help="run headless: validate and export, display nothing")
    parser.add_argument("--dry-run", action="store_true",
                        help="validate and write viewparams.tcl; start no VMD")
    return parser


def _fail(message: str) -> int:
    print(f"ariadne view: {message}", file=sys.stderr)
    return EXIT_USAGE


def main(argv: list[str]) -> int:
    args = build_arg_parser().parse_args(argv)

    try:
        run = discover_run(args.rundir, require_raw=not args.from_export)
    except ViewError as exc:
        return _fail(str(exc))

    if args.prmtop is not None:
        run = RunInfo(**{**run.__dict__, "prmtop": args.prmtop})
    if args.trajectory is not None:
        run = RunInfo(**{**run.__dict__, "trajectory": args.trajectory})

    try:
        params = viewconfig.parse_view_params(run.input_text)
        spec = args.frames if args.frames is not None else params.view_frames
        frames, missing = viewconfig.resolve_view_frames(spec, run.computed_frames)
    except ConfigError as exc:
        return _fail(str(exc))

    if missing:
        shown = ", ".join(str(f) for f in missing[:10])
        more = "" if len(missing) <= 10 else f" (+{len(missing) - 10} more)"
        print(
            f"ariadne view: {len(missing)} requested frame(s) were never computed "
            f"and will be skipped: {shown}{more}",
            file=sys.stderr,
        )

    if not frames:
        return _fail("no frames left to render")

    if args.from_export:
        try:
            check_export_frames(run.rundir, frames)
        except ViewError as exc:
            return _fail(str(exc))

    try:
        rows = pathsio.read_paths_csv(run.rundir / "paths.csv")
    except pathsio.PathsError as exc:
        return _fail(str(exc))
    grouped = pathsio.group_paths(rows)

    viewparams = run.rundir / "viewparams.tcl"
    viewparams.write_text(
        render_viewparams(run, frames, grouped, params, args.export, args.from_export)
    )
    print(f"Rendering {len(frames)} frame(s); wrote {viewparams}")

    if args.dry_run:
        return EXIT_OK

    if shutil.which("vmd") is None:
        return _fail("vmd is not on PATH")

    command = ["vmd"]
    env = dict(os.environ)
    if args.text:
        command += ["-dispdev", "text"]
    elif not env.get("DISPLAY"):
        return _fail("no DISPLAY available; use --text --export instead")
    command += ["-e", str(VIEW_TCL), "-args", str(viewparams)]

    # Only detach stdin in headless mode. In graphics mode VMD's Tcl console
    # reads stdin, so DEVNULL delivers EOF the moment the script finishes and
    # VMD exits "normally" the instant the window appears -- observed directly:
    # the window opened and closed inside two seconds. Inheriting the terminal
    # keeps the session alive and doubles as the VMD console.
    stdin = subprocess.DEVNULL if args.text else None
    completed = subprocess.run(command, env=env, stdin=stdin)
    if completed.returncode < 0:
        # Killed by a signal -- normally view.tcl's deliberate SIGKILL after a
        # failed atom resolution. Returning -9 makes the shell report 247.
        return EXIT_VMD_FAILED
    return completed.returncode
