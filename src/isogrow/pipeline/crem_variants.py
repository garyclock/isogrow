"""CReM variant expansion, mirroring the original run_fraggen.ipynb.

Each bioisostere candidate carries one or more ``[*:n]`` attachment points.  We
replace the dummy atoms with carbons (remembering their indices as protected
sites), then run CReM ``mutate_mol`` on the resulting lead to enumerate
synthesizable variants that keep those attachment points intact.
"""
from __future__ import annotations

from collections.abc import Iterable

from crem.crem import mutate_mol
from rdkit import Chem


def replace_stars_with_carbons(smiles: str):
    """Replace ``[*:n]`` with carbon, returning (lead_smiles, protected_atom_ids)."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0:
            atom.SetAtomicNum(6)
            atom.SetAtomMapNum(9999)
    new_mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol))
    if new_mol is None:
        return None, []
    replaced = []
    for atom in new_mol.GetAtoms():
        if atom.GetAtomMapNum() == 9999:
            replaced.append(atom.GetIdx())
            atom.SetAtomMapNum(0)
    return Chem.MolToSmiles(new_mol), replaced


def crem_mutate(
    lead_smiles: str,
    protected_ids: Iterable[int],
    db_path: str,
    radius: int = 2,
    min_size: int = 0,
    max_size: int = 0,
    min_inc: int = -2,
    max_inc: int = 2,
    max_replacements: int | None = None,
    replace_cycles: bool = False,
    protect_H: bool = True,
    seed: int | None = None,
) -> list[str]:
    """Run CReM mutate; returns SMILES list with the lead first (mirrors notebook).

    CReM samples replacements, so without ``seed`` two runs of the same config give
    different variant sets. Set it (config key ``crem.seed``) for reproducibility.
    """
    mol = Chem.MolFromSmiles(lead_smiles)
    if mol is None:
        return []
    mol = Chem.AddHs(mol)
    protected = set(protected_ids)
    if protect_H:
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == "H":
                for nb in atom.GetNeighbors():
                    if nb.GetIdx() in protected:
                        protected.add(atom.GetIdx())
    extra = {"seed": seed} if seed is not None else {}
    out = list(
        mutate_mol(
            mol,
            protected_ids=protected,
            db_name=db_path,
            radius=radius,
            min_size=min_size,
            max_size=max_size,
            min_inc=min_inc,
            max_inc=max_inc,
            max_replacements=max_replacements,
            replace_cycles=replace_cycles,
            **extra,
        )
    )
    return [Chem.MolToSmiles(Chem.RemoveHs(mol))] + out


def expand_candidates(candidates: Iterable[str], db_path: str, **params) -> list[dict]:
    """Expand a list of ``[*:n]`` bioisostere candidates into variant records."""
    rows: list[dict] = []
    for cand in candidates:
        lead, prot = replace_stars_with_carbons(cand)
        if lead is None:
            continue
        for v in crem_mutate(lead, prot, db_path, **params):
            rows.append({"candidate": cand, "lead": lead, "variant": v})
    return rows
