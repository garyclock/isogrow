"""Ligand PDBQT preparation with Meeko (RDKit-native, replaces the old py2 ADT script)."""
from __future__ import annotations

import os

from rdkit import Chem


def sdf_to_pdbqt(sdf_path: str, out_pdbqt: str) -> str | None:
    """Convert a 3D SDF to PDBQT using Meeko. Returns the output path or None."""
    mol = Chem.MolFromMolFile(sdf_path, removeHs=False)
    if mol is None:
        return None
    try:
        from meeko import MoleculePreparation, PDBQTWriterLegacy
        setups = MoleculePreparation().prepare(mol)
        if not setups:
            return None
        pdbqt_string, ok, _err = PDBQTWriterLegacy.write_string(setups[0])
        if not ok:
            return None
    except Exception:
        return None
    os.makedirs(os.path.dirname(out_pdbqt) or ".", exist_ok=True)
    with open(out_pdbqt, "w", encoding="utf-8") as f:
        f.write(pdbqt_string)
    return out_pdbqt


def batch_pdbqt(sdf_paths, outdir: str, batch_size: int = 500) -> list[str]:
    os.makedirs(outdir, exist_ok=True)
    out = []
    for idx, path in enumerate(sdf_paths):
        bdir = os.path.join(outdir, f"batch_{idx // batch_size:03d}")
        os.makedirs(bdir, exist_ok=True)
        opath = os.path.join(bdir, f"mol_{idx:06d}.pdbqt")
        if sdf_to_pdbqt(path, opath):
            out.append(opath)
    return out
