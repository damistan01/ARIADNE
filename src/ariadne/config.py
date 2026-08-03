"""Parse and validate the ariadne input file."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path

__all__ = ["ConfigError", "PathwaysParams", "parse_config", "parse_frames_key", "load_config"]


class ConfigError(ValueError):
    """Raised for any malformed or invalid input file."""


# Characters that would change the meaning of the generated Tcl params file.
_TCL_FORBIDDEN = set("{}[]\\$\"")

_SELECTION_KEYS = ("donor", "acceptor", "bridge")

# key -> (python type, plugin flag). Defaults mirror pathways.tcl:600-680.
_NUMERIC_KEYS: dict[str, tuple[type, str]] = {
    "npaths": (int, "-p"),
    "epsc": (float, "-epsc"),
    "epsh": (float, "-epsh"),
    "exph": (float, "-exph"),
    "r0h": (float, "-r0h"),
    "epsts": (float, "-epsts"),
    "expts": (float, "-expts"),
    "r0ts": (float, "-r0ts"),
    "hcut": (float, "-hcut"),
    "hang": (float, "-hang"),
    "tscut": (float, "-tscut"),
    "procut": (float, "-procut"),
    "withh": (int, "-withh"),
    "cda": (int, "-cda"),
}

# Keys that live in the same input file but are parsed elsewhere.
# parse_config must tolerate them rather than reject them as unknown.
_RUN_KEYS = frozenset({"frames"})
_VIEW_KEYS = frozenset({
    "view_frames",
    "view_radius",
    "view_ranks",
    "view_color_covalent",
    "view_color_hbond",
    "view_color_through_space",
    "view_context_radius",
    "view_clip_front",
})
_FOREIGN_KEYS = _RUN_KEYS | _VIEW_KEYS


@dataclass(frozen=True)
class PathwaysParams:
    donor: str
    acceptor: str
    bridge: str = "all"
    npaths: int = 1
    epsc: float = 0.6
    epsh: float = 0.36
    exph: float = 1.7
    r0h: float = 2.8
    epsts: float = 0.6
    expts: float = 1.7
    r0ts: float = 1.4
    hcut: float = 3.0
    hang: float = 30.0
    tscut: float = 6.0
    procut: float = 7.0
    withh: int = 0
    cda: int = 1

    def effective_bridge(self, include_water: bool) -> str:
        """Bridge selection actually handed to the plugin."""
        if include_water:
            return self.bridge
        return f"({self.bridge}) and not water"

    def to_plugin_args(self) -> list[str]:
        """Plugin flags that depend only on the input file.

        Excludes -d/-a/-b/-frame/-mol/-o, which the worker supplies per frame.
        """
        args: list[str] = []
        for key, (_, flag) in _NUMERIC_KEYS.items():
            args.extend([flag, _format_number(getattr(self, key))])
        return args

    def as_dict(self) -> dict:
        return asdict(self)


def _format_number(value: float | int) -> str:
    if isinstance(value, int):
        return str(value)
    return repr(value) if value != int(value) else str(int(value))


def _check_selection(key: str, value: str) -> None:
    if not value.strip():
        raise ConfigError(f"{key} must not be empty")
    bad = sorted(_TCL_FORBIDDEN & set(value))
    if bad:
        raise ConfigError(
            f"{key} contains characters that cannot be passed safely to Tcl: "
            f"{''.join(bad)}"
        )


def parse_config(text: str) -> PathwaysParams:
    values: dict[str, object] = {}
    for lineno, raw in enumerate(text.splitlines(), start=1):
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        if "=" not in line:
            raise ConfigError(f"line {lineno}: expected 'key = value', got {raw!r}")
        key, _, value = line.partition("=")
        key = key.strip().lower()
        value = value.strip()
        if key in _SELECTION_KEYS:
            _check_selection(key, value)
            values[key] = value
        elif key in _NUMERIC_KEYS:
            caster, _ = _NUMERIC_KEYS[key]
            try:
                values[key] = caster(value)
            except ValueError:
                raise ConfigError(
                    f"line {lineno}: {key} must be a {caster.__name__}, got {value!r}"
                ) from None
        elif key in _FOREIGN_KEYS:
            continue
        else:
            raise ConfigError(f"line {lineno}: unknown key {key!r}")

    for required in ("donor", "acceptor"):
        if required not in values:
            raise ConfigError(f"missing required key: {required}")

    params = PathwaysParams(**values)  # type: ignore[arg-type]
    if params.npaths < 1:
        raise ConfigError(f"npaths must be >= 1, got {params.npaths}")
    if params.withh not in (0, 1):
        raise ConfigError(f"withh must be 0 or 1, got {params.withh}")
    if params.cda not in (0, 1):
        raise ConfigError(f"cda must be 0 or 1, got {params.cda}")
    return params


def load_config(path: Path) -> PathwaysParams:
    try:
        text = Path(path).read_text()
    except OSError as exc:
        raise ConfigError(f"cannot read input file {path}: {exc}") from None
    return parse_config(text)


def parse_frames_key(text: str) -> str | None:
    """Read the `frames` key, which selects which frames to compute.

    Kept out of PathwaysParams deliberately: --resume compares parameters, and
    extending a frame range must not be treated as a conflicting run.
    """
    for raw in text.splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line or "=" not in line:
            continue
        key, _, value = line.partition("=")
        if key.strip().lower() == "frames":
            return value.strip()
    return None
