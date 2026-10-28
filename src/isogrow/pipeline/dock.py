"""Docking backends with a common interface; Vina is the default, others pluggable.

A backend takes (receptor_pdbqt, ligand_pdbqt, box) and returns a dict with at least
``affinity`` (kcal/mol, more negative = better).  ``unidock`` and ``external`` are
provided as thin subprocess wrappers so new engines can be added without touching the
pipeline.
"""
from __future__ import annotations

import os
import re
import subprocess
from collections.abc import Sequence

from isogrow import paths


def _parse_affinity(text: str) -> float | None:
    # REMARK VINA RESULT:    -8.5      0.000      0.000
    m = re.search(r"REMARK VINA RESULT:\s*(-?\d+\.?\d*)", text)
    if m:
        return float(m.group(1))
    # stdout table: "   1        -8.5      0.000      0.000"
    m = re.search(r"^\s*1\s+(-?\d+\.\d+)", text, re.M)
    return float(m.group(1)) if m else None


def vina_dock(
    receptor: str,
    ligand_pdbqt: str,
    center: Sequence[float],
    size: Sequence[float],
    out_pdbqt: str | None = None,
    exhaustiveness: int = 8,
    cpu: int = 4,
    vina_exe: str | None = None,
    timeout: int = 600,
) -> dict:
    try:
        exe = paths.vina_exe(vina_exe)
    except FileNotFoundError as e:
        return {"ligand": ligand_pdbqt, "affinity": None, "error": str(e)}
    cmd = [
        exe,
        "--receptor", receptor,
        "--ligand", ligand_pdbqt,
        "--center_x", f"{center[0]:.3f}", "--center_y", f"{center[1]:.3f}", "--center_z", f"{center[2]:.3f}",
        "--size_x", f"{size[0]:.3f}", "--size_y", f"{size[1]:.3f}", "--size_z", f"{size[2]:.3f}",
        "--exhaustiveness", str(exhaustiveness), "--cpu", str(cpu),
    ]
    if out_pdbqt:
        cmd += ["--out", out_pdbqt]
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except Exception as e:
        return {"ligand": ligand_pdbqt, "affinity": None, "error": str(e)}
    text = p.stdout + "\n" + p.stderr
    if out_pdbqt and os.path.exists(out_pdbqt):
        try:
            text += "\n" + open(out_pdbqt, encoding="utf-8", errors="replace").read()
        except Exception:
            pass
    return {"ligand": ligand_pdbqt, "affinity": _parse_affinity(text), "returncode": p.returncode}


def unidock_dock(receptor, ligand_pdbqt, center, size, out_pdbqt=None,
                 unidock_exe: str = "unidock", **kw) -> dict:
    """Thin Uni-Dock wrapper (GPU engine); parsed the same way."""
    cmd = [unidock_exe, "--receptor", receptor, "--gpu_batch", ligand_pdbqt,
           "--center_x", str(center[0]), "--center_y", str(center[1]), "--center_z", str(center[2]),
           "--size_x", str(size[0]), "--size_y", str(size[1]), "--size_z", str(size[2])]
    if out_pdbqt:
        cmd += ["--dir", os.path.dirname(out_pdbqt)]
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=kw.get("timeout", 600))
    except Exception as e:
        return {"ligand": ligand_pdbqt, "affinity": None, "error": str(e)}
    return {"ligand": ligand_pdbqt, "affinity": _parse_affinity(p.stdout + p.stderr)}


BACKENDS = {"vina": vina_dock, "unidock": unidock_dock}


def get_backend(name: str):
    if name not in BACKENDS:
        raise ValueError(f"unknown docking backend: {name}")
    return BACKENDS[name]
