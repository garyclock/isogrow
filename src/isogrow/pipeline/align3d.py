"""3D alignment of generated compounds onto a reference ligand.

Two backends, selected by ``backend=`` / the ``align_backend`` config key:

``mcs`` (default)
    Rigid-body superposition on the maximum common substructure:
    ``rdFMCS`` → atom map → ``rdMolAlign.AlignMol``. Pure Python, fast, no external
    binary. Internal geometry of the target is untouched.

``kcombu``
    KCOMBU's ``fkcombu`` (BSD-3-Clause, Linux/macOS only), the tool the original
    notebook used. It first builds a 2D-MCS correspondence, then **deforms** the
    target: torsion angles are stamped from the reference (``-cfstp``) and a
    steepest-descent fit explores the remaining torsions over ``nini`` random
    initial conformations. With ``-P`` it can even avoid clashing with a receptor.
    That is strictly more than a rigid fit, at the cost of one process per molecule.

Both backends report RMSD the same way — the deviation of mapped MCS atoms after
alignment, without re-fitting — so results stay comparable.
"""
from __future__ import annotations

import os
import re
import shlex
import subprocess
import tempfile
from collections.abc import Sequence
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdFMCS, rdMolAlign

from isogrow import paths

BACKENDS = ("mcs", "kcombu")


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def load_mol(sdf_or_mol, remove_hs: bool = False):
    if isinstance(sdf_or_mol, Chem.Mol):
        return sdf_or_mol
    return Chem.MolFromMolFile(sdf_or_mol, removeHs=remove_hs)


def box_from_ligand(ref_sdf: str, padding: float = 8.0):
    """Derive an AutoDock box (center, size) from a reference ligand's coordinates."""
    mol = Chem.MolFromMolFile(ref_sdf, removeHs=False)
    if mol is None:
        raise ValueError(f"cannot read {ref_sdf}")
    conf = mol.GetConformer()
    xs = [conf.GetAtomPosition(i).x for i in range(mol.GetNumAtoms())]
    ys = [conf.GetAtomPosition(i).y for i in range(mol.GetNumAtoms())]
    zs = [conf.GetAtomPosition(i).z for i in range(mol.GetNumAtoms())]
    center = ((min(xs) + max(xs)) / 2, (min(ys) + max(ys)) / 2, (min(zs) + max(zs)) / 2)
    size = (max(xs) - min(xs) + padding, max(ys) - min(ys) + padding, max(zs) - min(zs) + padding)
    return center, size


def mcs_mapping(mol, ref, timeout: int = 10):
    """Return ``(match_in_mol, match_in_ref)`` for the MCS, or None if too small."""
    mcs = rdFMCS.FindMCS(
        [mol, ref], timeout=timeout, completeRingsOnly=True, ringMatchesRingOnly=True
    )
    if mcs.numAtoms < 3:
        return None
    core = Chem.MolFromSmarts(mcs.smartsString)
    mp = mol.GetSubstructMatch(core)
    mr = ref.GetSubstructMatch(core)
    if not mp or not mr:
        return None
    return mp, mr


def mapped_rmsd(mol, ref, mapping) -> float:
    """RMSD of the mapped atoms *as they are* (no re-fitting)."""
    mp, mr = mapping
    cm, cr = mol.GetConformer(), ref.GetConformer()
    total = 0.0
    for i, j in zip(mp, mr, strict=True):
        a, b = cm.GetAtomPosition(i), cr.GetAtomPosition(j)
        total += (a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2
    return (total / len(mp)) ** 0.5


# --------------------------------------------------------------------------- #
# mcs backend (rigid)
# --------------------------------------------------------------------------- #
def align_mcs(mol, ref, timeout: int = 10) -> tuple[Chem.Mol | None, float | None]:
    """Rigid-body MCS superposition. Returns ``(mol, rmsd)``; rmsd None if no MCS."""
    mol = load_mol(mol)
    ref = load_mol(ref)
    if mol is None or ref is None:
        return None, None
    mapping = mcs_mapping(mol, ref, timeout=timeout)
    if mapping is None:
        return mol, None
    mp, mr = mapping
    rdMolAlign.AlignMol(mol, ref, atomMap=list(zip(mp, mr, strict=True)))
    return mol, mapped_rmsd(mol, ref, mapping)


#: Backwards-compatible alias (the pipeline used to call this directly).
align_to_ref = align_mcs


# --------------------------------------------------------------------------- #
# kcombu backend (flexible)
# --------------------------------------------------------------------------- #
def _cmd_prefix(exe) -> list[str]:
    """``exe`` may be a path, a bare name, a command string, or an argv list.

    A command string is how you drive a binary that lives outside this OS, e.g.
    ``ISOGROW_KCOMBU="wsl.exe -d Ubuntu -- /home/me/miniforge3/bin/fkcombu"``.
    """
    if isinstance(exe, (list, tuple)):
        return [str(x) for x in exe]
    s = str(exe)
    if os.path.isfile(s) or " " not in s:
        return [s]
    return shlex.split(s)


def _is_wsl(prefix: Sequence[str]) -> bool:
    return bool(prefix) and Path(str(prefix[0])).name.lower().startswith("wsl")


def _to_wsl_path(path, enabled: bool) -> str:
    """Translate ``D:\\a\\b.sdf`` to ``/mnt/d/a/b.sdf`` when driving WSL."""
    p = Path(path).resolve()
    if not enabled or os.name != "nt" or not p.drive:
        return str(p)
    drive = p.drive.rstrip(":").lower()
    rest = str(p).replace("\\", "/").split(":", 1)[1]
    return f"/mnt/{drive}{rest}"


def align_kcombu(
    target_path,
    ref_path,
    out_path,
    exe=None,
    rigid: bool = False,
    nini: int = 1,
    timeout: int = 180,
    extra: Sequence[str] = (),
) -> tuple[str, float | None]:
    """Run ``fkcombu`` and write the aligned target to ``out_path`` (SDF).

    Returns ``(out_path, rmsd_match)`` where ``rmsd_match`` is KCOMBU's own
    deviation over *its* atom correspondence (2D-MCS based). It is **not**
    comparable with the MCS-backend RMSD, because the two use different atom
    correspondences — compare backends by docking score instead.

    ``rigid=True`` passes ``-S R`` (rigid fit); the default ``-S F`` is KCOMBU's
    flexible fit (torsion stamping + steepest descent). ``nini`` controls how many
    random initial conformers are explored — keep it at 1 for reproducibility.

    ``exe`` may be a list or a command string to drive a non-Windows binary, e.g.
    ``"wsl.exe -d Ubuntu -- /home/me/miniforge3/bin/fkcombu"``; paths are then
    translated to ``/mnt/<drive>/...`` automatically.
    """
    prefix = _cmd_prefix(paths.kcombu_exe(exe) if not isinstance(exe, (list, tuple)) else exe)
    wsl = _is_wsl(prefix)
    cmd = [
        *prefix,
        "-T", _to_wsl_path(target_path, wsl),
        "-R", _to_wsl_path(ref_path, wsl),
        "-osdfT", _to_wsl_path(out_path, wsl),
        "-nini", str(nini),
        "-S", "R" if rigid else "F",
        *extra,
    ]
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    if Path(out_path).exists():
        Path(out_path).unlink()
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    if not Path(out_path).exists():
        tail = (proc.stdout or "")[-400:] + (proc.stderr or "")[-400:]
        raise RuntimeError(f"fkcombu produced no output (exit {proc.returncode}): {tail.strip()}")
    m = re.search(r"RMSD_MATCH\s+([0-9.]+)", proc.stdout or "")
    return str(out_path), (float(m.group(1)) if m else None)


# --------------------------------------------------------------------------- #
# unified entry points
# --------------------------------------------------------------------------- #
def align_file(target_path, ref, out_path, backend: str = "mcs", **kwargs):
    """Align one SDF file and write the result; returns ``(mol, rmsd)``."""
    if backend not in BACKENDS:
        raise ValueError(f"unknown alignment backend {backend!r}; choose from {BACKENDS}")
    ref_mol = load_mol(ref)

    if backend == "kcombu":
        ref_path = kwargs.pop("ref_path", None) or _ref_path(ref)
        work = kwargs.pop("workdir", None) or tempfile.mkdtemp(prefix="isogrow_align_")
        tmp_out = os.path.join(work, f"{Path(target_path).stem}_kcombu.sdf")
        _path, kcombu_rmsd = align_kcombu(target_path, ref_path, tmp_out, **kwargs)
        mol = Chem.MolFromMolFile(tmp_out, removeHs=False)
        if mol is None:
            raise RuntimeError(f"fkcombu output could not be read: {tmp_out}")
        # fkcombu writes coordinates but leaves the 3D flag unset; downstream
        # (Meeko) needs it, and RDKit warns on every read otherwise.
        conf = mol.GetConformer()
        if not conf.Is3D():
            conf.Set3D(True)
        # fkcombu's SDF omits some hydrogens, leaving implicit Hs behind, which
        # Meeko rejects outright. Restore them with coordinates.
        if any(atom.GetNumImplicitHs() for atom in mol.GetAtoms()):
            mol = Chem.AddHs(mol, addCoords=True)

        # KCOMBU reports RMSD over its own 2D-MCS correspondence; fall back to a
        # local MCS mapping only if that line is missing.
        rmsd = kcombu_rmsd
        if rmsd is None:
            mapping = mcs_mapping(mol, ref_mol)
            rmsd = mapped_rmsd(mol, ref_mol, mapping) if mapping else None
        os.makedirs(Path(out_path).parent, exist_ok=True)
        with Chem.SDWriter(str(out_path)) as writer:
            writer.write(mol)
        return mol, rmsd

    mol = Chem.MolFromMolFile(target_path, removeHs=False)
    if mol is None:
        raise ValueError(f"cannot read {target_path}")
    mol, rmsd = align_mcs(mol, ref_mol, timeout=kwargs.get("timeout", 10))
    with Chem.SDWriter(str(out_path)) as w:
        w.write(mol)
    return mol, rmsd


def _ref_path(ref):
    """``ref`` is normally a path; keep a hook for Mol inputs (kcombu needs a file)."""
    if isinstance(ref, (str, os.PathLike)):
        return str(ref)
    raise TypeError("the kcombu backend needs the reference ligand as a file path")


def batch_align(sdf_paths, ref_sdf: str, outdir: str, backend: str = "mcs",
                batch_size: int = 500, **kwargs):
    """Align a list of 3D SDF files; returns the written paths (failures skipped)."""
    os.makedirs(outdir, exist_ok=True)
    ref = Chem.MolFromMolFile(ref_sdf, removeHs=False)
    out = []
    for idx, path in enumerate(sdf_paths):
        bdir = os.path.join(outdir, f"batch_{idx // batch_size:03d}")
        os.makedirs(bdir, exist_ok=True)
        opath = os.path.join(bdir, f"mol_{idx:06d}_aligned.sdf")
        try:
            if backend == "mcs":
                mol = Chem.MolFromMolFile(path, removeHs=False)
                if mol is None:
                    continue
                mol, _rms = align_mcs(mol, ref)
                with Chem.SDWriter(opath) as w:
                    w.write(mol)
            else:
                align_file(path, ref_sdf, opath, backend=backend,
                           workdir=os.path.join(bdir, "_kcombu"), ref_path=ref_sdf, **kwargs)
        except Exception:  # noqa: BLE001 - one bad molecule must not stop the batch
            continue
        out.append(opath)
    return out
