"""3D conformer generation: ETKDG embedding + light MMFF94 minimisation."""
from __future__ import annotations

import os
from collections.abc import Iterable

from rdkit import Chem
from rdkit.Chem import AllChem


def embed_3d(smiles: str, seed: int = 42, max_its: int = 20) -> Chem.Mol | None:
    """Return a 3D molecule (H-added, MMFF-minimised) or None on failure."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, useRandomCoords=False, randomSeed=seed) != 0:
        if AllChem.EmbedMolecule(mol, useRandomCoords=True, randomSeed=seed) != 0:
            return None
    try:
        props = AllChem.MMFFGetMoleculeProperties(mol)
        ff = AllChem.MMFFGetMoleculeForceField(mol, props)
        ff.Minimize(maxIts=max_its)
    except Exception:
        pass
    return mol


def batch_embed(smiles_list: Iterable[str], outdir: str, batch_size: int = 500) -> list[str]:
    """Write 3D SDFs (batch_000/mol_000000.sdf ...); returns written paths."""
    os.makedirs(outdir, exist_ok=True)
    paths: list[str] = []
    for idx, smi in enumerate(smiles_list):
        mol = embed_3d(smi)
        if mol is None:
            continue
        bdir = os.path.join(outdir, f"batch_{idx // batch_size:03d}")
        os.makedirs(bdir, exist_ok=True)
        path = os.path.join(bdir, f"mol_{idx:06d}.sdf")
        with Chem.SDWriter(path) as w:
            w.write(mol)
        paths.append(path)
    return paths
