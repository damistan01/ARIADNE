"""Presentation options for the 3D viewer.

Separate from PathwaysParams on purpose: --resume compares physics parameters
and refuses to continue when they differ, which is correct. Presentation keys
must never participate in that comparison, or nudging a cylinder radius would
block resuming a ten-hour run.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

from .config import ConfigError
from .plan import PlanError, parse_frame_spec

__all__ = [
    "ViewParams",
    "parse_view_params",
    "resolve_view_frames",
    "VMD_COLORS",
    "resolve_color",
]

_VALID_RANKS = ("0", "all")
_TRUTHY = frozenset({"yes", "true", "on", "1"})
_FALSY = frozenset({"no", "false", "off", "0"})

# VMD's built-in colour table, in index order. Verified against
# `colorinfo colors` on VMD 2.0.0; the index is what `graphics ... color`
# expects, and it is what the Pathways plugin's draw_step passes straight
# through.
VMD_COLORS = (
    "blue", "red", "gray", "orange", "yellow", "tan", "silver", "green",
    "white", "pink", "cyan", "purple", "lime", "mauve", "ochre", "iceblue",
    "black", "yellow2", "yellow3", "green2", "green3", "cyan2", "cyan3",
    "blue2", "blue3", "violet", "violet2", "magenta", "magenta2", "red2",
    "red3", "orange2", "orange3",
)

# Warm defaults on purpose. VMD colours atoms by element by default -- nitrogen
# blue, oxygen red -- so drawing the path in blue or red makes it compete with
# the licorice underneath it. Orange / light orange / yellow reads as one
# gradient and collides with nothing. orange3 is the LIGHTER orange
# (0.96 0.72 0.00); orange2 is darker (0.89 0.35 0.00).
DEFAULT_COLOR_COVALENT = "orange"
DEFAULT_COLOR_HBOND = "orange3"
DEFAULT_COLOR_THROUGH_SPACE = "yellow"


def resolve_color(key: str, value: str, lineno: int | None = None) -> str:
    """Validate a colour given as a VMD name or an index, returning the name."""
    where = "" if lineno is None else f"line {lineno}: "
    text = value.strip()
    if text.isdigit():
        index = int(text)
        if not 0 <= index < len(VMD_COLORS):
            raise ConfigError(
                f"{where}{key} index {index} is out of range; "
                f"VMD has colours 0-{len(VMD_COLORS) - 1}"
            )
        return VMD_COLORS[index]
    if text not in VMD_COLORS:
        raise ConfigError(
            f"{where}{key} must be a VMD colour name or index 0-"
            f"{len(VMD_COLORS) - 1}, got {value!r}. Valid names: "
            f"{', '.join(VMD_COLORS)}"
        )
    return text


@dataclass(frozen=True)
class ViewParams:
    view_frames: str | None = None
    view_radius: float = 0.3
    view_ranks: str = "all"
    view_color_covalent: str = DEFAULT_COLOR_COVALENT
    view_color_hbond: str = DEFAULT_COLOR_HBOND
    view_color_through_space: str = DEFAULT_COLOR_THROUGH_SPACE
    view_context_radius: float = 12.0
    view_clip_front: bool = True

    def color_index(self, bond_type: str) -> int:
        """VMD colour index for a bond type, as draw_step expects."""
        name = {
            "covalent": self.view_color_covalent,
            "hbond": self.view_color_hbond,
            "through_space": self.view_color_through_space,
        }[bond_type]
        return VMD_COLORS.index(name)

    def as_dict(self) -> dict:
        return asdict(self)


def parse_view_params(text: str) -> ViewParams:
    values: dict[str, object] = {}
    for lineno, raw in enumerate(text.splitlines(), start=1):
        line = raw.split("#", 1)[0].strip()
        if not line or "=" not in line:
            continue
        key, _, value = line.partition("=")
        key = key.strip().lower()
        value = value.strip()
        if key == "view_frames":
            values["view_frames"] = value
        elif key == "view_radius":
            try:
                radius = float(value)
            except ValueError:
                raise ConfigError(
                    f"line {lineno}: view_radius must be a number, got {value!r}"
                ) from None
            if radius <= 0:
                raise ConfigError(
                    f"line {lineno}: view_radius must be greater than 0, got {radius}"
                )
            values["view_radius"] = radius
        elif key == "view_ranks":
            if value not in _VALID_RANKS:
                raise ConfigError(
                    f"line {lineno}: view_ranks must be one of "
                    f"{', '.join(_VALID_RANKS)}, got {value!r}"
                )
            values["view_ranks"] = value
        elif key in (
            "view_color_covalent",
            "view_color_hbond",
            "view_color_through_space",
        ):
            values[key] = resolve_color(key, value, lineno)
        elif key == "view_context_radius":
            try:
                radius = float(value)
            except ValueError:
                raise ConfigError(
                    f"line {lineno}: view_context_radius must be a number "
                    f"of angstroms, got {value!r}"
                ) from None
            if radius < 0:
                raise ConfigError(
                    f"line {lineno}: view_context_radius must be >= 0 "
                    f"(0 hides the cartoon entirely), got {radius}"
                )
            values["view_context_radius"] = radius
        elif key == "view_clip_front":
            lowered = value.lower()
            if lowered not in _TRUTHY | _FALSY:
                raise ConfigError(
                    f"line {lineno}: view_clip_front must be yes or no, "
                    f"got {value!r}"
                )
            values["view_clip_front"] = lowered in _TRUTHY
    return ViewParams(**values)  # type: ignore[arg-type]


def resolve_view_frames(
    spec: str | None,
    computed: list[int],
) -> tuple[list[int], list[int]]:
    """Split a render request into (frames we can show, frames never computed).

    Requesting frames that were not computed is a warning, not an error: the
    remaining frames still render.
    """
    available = sorted(set(computed))
    if spec is None or spec.strip() == "":
        return available, []

    # Sanity check: validate spec size before materializing a list.
    # Parse arithmetically to determine frame count without allocation.
    ceiling = 10**9
    parts = spec.split(":")
    if len(parts) > 3:
        raise ConfigError(f"view_frames: {spec!r} has too many ':'-separated fields")

    # Try to parse numbers; if it fails, let parse_frame_spec handle it.
    numbers = None
    try:
        numbers = [int(p) if p.strip() != "" else None for p in parts]
    except ValueError:
        # Non-numeric value; parse_frame_spec will handle and raise ConfigError
        pass

    # An open-ended stop ("0:", "0::2", ":") means "up to the last computed
    # frame". parse_frame_spec fills an omitted stop from the n_frames it is
    # given, so handing it the 10**9 ceiling would turn `view_frames = 0:` into
    # a billion-frame request that the size guard below then rejects. Resolve
    # the limit from the computed frames instead, and only fall back to the
    # ceiling when the stop is written out explicitly.
    open_stop = len(parts) > 1 and parts[1].strip() == ""
    if open_stop:
        if not available:
            # Nothing computed: there is no "end" to run to, and the caller
            # already fails on an empty selection with a clearer message.
            return [], []
        limit = available[-1] + 1
    else:
        limit = ceiling

    # Size check happens outside any broad exception handler, so ConfigError propagates.
    if numbers is not None:
        # Determine start, stop, stride with the same logic as parse_frame_spec
        if len(numbers) == 1:
            start = numbers[0] if numbers[0] is not None else 0
            stop = start
            stride = 1
        else:
            start = numbers[0] if numbers[0] is not None else 0
            stop = numbers[1] if numbers[1] is not None else limit - 1
            stride = numbers[2] if len(numbers) == 3 and numbers[2] is not None else 1

        # Calculate frame count without allocating
        if stride > 0 and stop >= start:
            frame_count = (stop - start) // stride + 1
            # Reject if absurdly large (before attempting allocation)
            if frame_count > 10**6:
                raise ConfigError(
                    f"view_frames spec {spec!r} expands to {frame_count:,} frames, "
                    f"which exceeds the sanity limit of 1,000,000"
                )

    # Now safely call parse_frame_spec with the resolved limit
    try:
        requested = parse_frame_spec(spec, limit)
    except PlanError as exc:
        raise ConfigError(f"view_frames: {exc}") from None

    have = set(available)
    selected = [f for f in requested if f in have]
    missing = [f for f in requested if f not in have]
    return selected, missing
