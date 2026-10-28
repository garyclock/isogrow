"""Path and external-tool resolution.

Everything that used to be a hard-coded ``D:/...`` path goes through here, so the
package works from any checkout, container or CI runner.

Resolution order for the **data directory** (where downloaded ChEMBL files and
built mmpdb databases live):

1. ``$ISOGROW_DATA``
2. ``<repo root>/data`` (i.e. the ``data/`` folder next to ``pyproject.toml``)

Resolution order for **external executables** (``mmpdb``, ``vina``):

1. an explicit argument from the caller / config file
2. the matching environment variable (``ISOGROW_MMPDB``, ``ISOGROW_VINA``)
3. ``PATH`` lookup via :func:`shutil.which`
"""
from __future__ import annotations

import os
import shutil
from pathlib import Path

ENV_DATA = "ISOGROW_DATA"
ENV_VINA = "ISOGROW_VINA"
ENV_MMPDB = "ISOGROW_MMPDB"
ENV_KCOMBU = "ISOGROW_KCOMBU"


def repo_root() -> Path:
    """Directory that contains ``pyproject.toml`` (two levels above this file)."""
    return Path(__file__).resolve().parents[2]


def data_dir(create: bool = False) -> Path:
    """Root folder for downloaded datasets and built databases."""
    env = os.environ.get(ENV_DATA)
    path = Path(env).expanduser() if env else repo_root() / "data"
    if create:
        path.mkdir(parents=True, exist_ok=True)
    return path


def resolve_exe(name: str, *candidates: str | None, env: str | None = None,
                required_msg: str | None = None) -> str:
    """Return the first usable executable path from candidates/env/PATH.

    Raises :class:`FileNotFoundError` with an actionable message when nothing is
    found and ``required_msg`` is given; otherwise returns ``name`` unchanged so
    the caller's subprocess error surfaces naturally.
    """
    for cand in candidates:
        if cand:
            # An absolute path must exist; a bare name may still be on PATH.
            if os.path.sep in cand or (os.path.altsep and os.path.altsep in cand):
                if Path(cand).expanduser().is_file():
                    return str(Path(cand).expanduser())
            else:
                found = shutil.which(cand)
                if found:
                    return found
    if env:
        value = os.environ.get(env)
        if value:
            resolved = resolve_exe(name, value)
            if resolved != name:
                return resolved
    found = shutil.which(name)
    if found:
        return found
    if required_msg:
        raise FileNotFoundError(required_msg)
    return name


def vina_exe(explicit: str | None = None) -> str:
    """Locate the AutoDock Vina executable."""
    return resolve_exe(
        "vina",
        explicit,
        env=ENV_VINA,
        required_msg=(
            "AutoDock Vina was not found. Install it (conda install -c conda-forge "
            "vina, apt install autodock-vina, or `python scripts/data/get_vina.py`) "
            "and/or set $ISOGROW_VINA to the executable path."
        ),
    )


def mmpdb_exe(explicit: str | None = None) -> str:
    """Locate the ``mmpdb`` command-line tool (installed by the ``mmpdb`` package)."""
    return resolve_exe(
        "mmpdb",
        explicit,
        env=ENV_MMPDB,
        required_msg=(
            "The mmpdb command line tool was not found. Install it with "
            "`pip install \"mmpdb @ git+https://github.com/rdkit/mmpdb.git\"` "
            "(the PyPI package of the same name does not provide the module)."
        ),
    )


def kcombu_exe(explicit: str | None = None) -> str:
    """Locate ``fkcombu`` (KCOMBU's flexible alignment tool).

    Linux/macOS only — there is no Windows build. Install with
    ``conda install -c michellab fkcombu`` or compile from
    https://pdbj.org/kcombu/download_src.html.

    The value may also be a **command string**, which is how you drive a binary
    installed inside WSL from Windows::

        set ISOGROW_KCOMBU=wsl.exe -d Ubuntu -- /home/me/miniforge3/bin/fkcombu
    """
    if explicit:
        return explicit
    env = os.environ.get(ENV_KCOMBU)
    if env:
        return env
    found = shutil.which("fkcombu")
    if found:
        return found
    raise FileNotFoundError(
        "fkcombu (KCOMBU) was not found. It is Linux/macOS only: install with "
        "`conda install -c michellab fkcombu` or build from "
        "https://pdbj.org/kcombu/download_src.html, then set $ISOGROW_KCOMBU. "
        "On Windows, point it at WSL: ISOGROW_KCOMBU=\"wsl.exe -d Ubuntu -- "
        "/home/<user>/miniforge3/bin/fkcombu\"."
    )


def default_db_path(name: str = "chembl37.mmpdb") -> Path:
    """Conventional location of a built database inside the data directory."""
    return data_dir() / "mmpdb" / name
