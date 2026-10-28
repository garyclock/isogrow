"""Pipeline configuration: JSON file + ``${VAR}`` expansion + validation.

Config files are plain JSON so they are easy to diff and generate. Any string
value may reference environment variables, plus two always-defined ones:

* ``${ISOGROW_DATA}`` — dataset root (see :mod:`isogrow.paths`)
* ``${ISOGROW_HOME}`` — repository root

Unknown ``${...}`` names raise a :class:`ConfigError` instead of silently
expanding to an empty string, which is how these configs used to fail.
"""
from __future__ import annotations

import json
import os
import re
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from . import paths

_VAR = re.compile(r"\$\{([A-Za-z_][A-Za-z0-9_]*)\}")

#: Keys that hold filesystem paths and are normalised to absolute paths.
PATH_KEYS = ("mmpdb_path", "csv_path", "crem_db", "ref_ligand", "receptor", "out")

REQUIRED = ("out", "crem_db", "ref_ligand", "receptor")


class ConfigError(ValueError):
    """Raised when a pipeline config is missing, malformed or inconsistent."""


def _builtin_vars() -> dict[str, str]:
    return {
        "ISOGROW_DATA": str(paths.data_dir()),
        "ISOGROW_HOME": str(paths.repo_root()),
    }


def expand_vars(value: Any, _vars: Mapping[str, str] | None = None) -> Any:
    """Recursively expand ``${VAR}`` in strings, lists and dicts."""
    table = dict(_builtin_vars())
    table.update(os.environ)
    if _vars:
        table.update(_vars)

    def sub(text: str) -> str:
        def repl(m: re.Match) -> str:
            name = m.group(1)
            if name not in table:
                raise ConfigError(
                    f"undefined variable ${{{name}}}; define it in the environment "
                    f"or pass it explicitly"
                )
            return table[name]

        return _VAR.sub(repl, text)

    if isinstance(value, str):
        return sub(value)
    if isinstance(value, list):
        return [expand_vars(v, _vars) for v in value]
    if isinstance(value, dict):
        return {k: expand_vars(v, _vars) for k, v in value.items()}
    return value


def resolve_paths(cfg: dict, base: os.PathLike | str | None = None) -> dict:
    """Turn the known path keys into absolute paths (relative to ``base``)."""
    root = Path(base) if base else Path.cwd()
    out = dict(cfg)
    for key in PATH_KEYS:
        val = out.get(key)
        if not isinstance(val, str) or not val:
            continue
        out[key] = str((root / Path(val).expanduser()).resolve()) if not Path(val).is_absolute() \
            else str(Path(val).expanduser())
    return out


def validate(cfg: Mapping[str, Any]) -> None:
    """Check the keys the pipeline actually cannot run without."""
    missing = [k for k in REQUIRED if not cfg.get(k)]
    if missing:
        raise ConfigError(f"missing required config keys: {', '.join(missing)}")
    source = cfg.get("source", "official_csv")
    if source == "official_csv" and not cfg.get("csv_path"):
        raise ConfigError("source='official_csv' requires 'csv_path'")
    if source == "local" and not (cfg.get("query") and cfg.get("mmpdb_path")):
        raise ConfigError("source='local' requires 'query' and 'mmpdb_path'")
    if source not in ("official_csv", "local"):
        raise ConfigError(f"unknown source: {source}")


def load_config(path: os.PathLike | str, overrides: Mapping[str, Any] | None = None,
                base: os.PathLike | str | None = None) -> dict:
    """Read a JSON config, expand variables, merge overrides, resolve paths."""
    p = Path(path)
    if not p.is_file():
        raise ConfigError(f"config file not found: {p}")
    try:
        raw = json.loads(p.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:  # pragma: no cover - depends on user input
        raise ConfigError(f"invalid JSON in {p}: {exc}") from exc
    if not isinstance(raw, dict):
        raise ConfigError(f"config root must be an object, got {type(raw).__name__}")

    cfg = expand_vars(raw)
    if overrides:
        cfg.update({k: v for k, v in overrides.items() if v is not None})
    # Relative paths in a config file are interpreted relative to the working
    # directory (standard CLI behaviour); ${ISOGROW_DATA} stays absolute.
    cfg = resolve_paths(cfg, base=base)
    validate(cfg)
    return cfg
