"""3D alignment to a reference ligand, replacing the compiled KCOMBU tool.

We use RDKit maximum-common-substructure (MCS) matching to superimpose a generated
fragment onto the corresponding substructure of the reference ligand, then report RMSD.
This keeps the whole pipeline pure-Python / cross-platform.  A ``kcombu`` backend can
be plugged in later if the binary is available.
"""
from __future__ import annotations

import os

from rdkit import Chem
from rdkit.Chem import rdFMCS, rdMolAlign


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


def align_to_ref(mol, ref, timeout: int = 10) -> tuple[Chem.Mol | None, float | None]:
    """Align ``mol`` onto ``ref`` on their MCS; mutate and return (mol, rmsd)."""
    mol = load_mol(mol)
    ref = load_mol(ref)
    if mol is None or ref is None:
        return None, None
    mcs = rdFMCS.FindMCS(
        [mol, ref], timeout=timeout, completeRingsOnly=True, ringMatchesRingOnly=True
    )
    if mcs.numAtoms < 3:
        return mol, None
    core = Chem.MolFromSmarts(mcs.smartsString)
    mp = mol.GetSubstructMatch(core)
    mr = ref.GetSubstructMatch(core)
    if not mp or not mr:
        return mol, None
    rms = rdMolAlign.AlignMol(mol, ref, atomMap=list(zip(mp, mr, strict=True)))
    return mol, rms


def batch_align(sdf_paths, ref_sdf: str, outdir: str, batch_size: int = 500):
    """Align a list of 3D SDF files to ``ref_sdf``; writes aligned SDFs."""
    os.makedirs(outdir, exist_ok=True)
    ref = Chem.MolFromMolFile(ref_sdf, removeHs=False)
    out = []
    for idx, path in enumerate(sdf_paths):
        mol = Chem.MolFromMolFile(path, removeHs=False)
        if mol is None:
            continue
        mol, rms = align_to_ref(mol, ref)
        bdir = os.path.join(outdir, f"batch_{idx // batch_size:03d}")
        os.makedirs(bdir, exist_ok=True)
        opath = os.path.join(bdir, f"mol_{idx:06d}_aligned.sdf")
        with Chem.SDWriter(opath) as w:
            w.write(mol)
        out.append(opath)
    return out
