"""Fragment property computation (RDKit), SwissBioisostere-compatible deltas.

SwissBioisostere reports ΔlogP / ΔtPSA / ΔMW of the *fragment* corresponding to a
replacement.  Values are the difference between the candidate fragment and the
query fragment, computed on the fragment level (attachment points capped).

Capping convention: each attachment point ``[*:n]`` is replaced by an explicit
hydrogen ``[H]`` (the standard way to compute a fragment's properties with the
dangling bond saturated), then properties are computed with RDKit.
"""
from __future__ import annotations

import re

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors

_DUMMY_RE = re.compile(r"\[\*:\d+\]")


def cap_dummy(smiles: str, cap: str = "[H]") -> str:
    """Replace all ``[*:n]`` attachment points with ``cap`` (default H)."""
    return _DUMMY_RE.sub(cap, smiles)


def fragment_props(smiles: str, cap: bool = True) -> dict | None:
    """Compute a property dict for a fragment SMILES (may contain [*:n]).

    Returns None if the SMILES cannot be parsed.
    """
    s = cap_dummy(smiles) if cap else smiles
    mol = Chem.MolFromSmiles(s)
    if mol is None:
        return None
    return {
        "smiles": s,
        "logp": Crippen.MolLogP(mol),
        "tpsa": rdMolDescriptors.CalcTPSA(mol),
        "mw": Descriptors.MolWt(mol),          # average MW (matches SwissBioisostere)
        "hac": mol.GetNumHeavyAtoms(),
        "rings": rdMolDescriptors.CalcNumRings(mol),
        "hbd": rdMolDescriptors.CalcNumHBD(mol),
        "hba": rdMolDescriptors.CalcNumHBA(mol),
        "rotb": rdMolDescriptors.CalcNumRotatableBonds(mol),
    }


def delta_props(query_smiles: str, cand_smiles: str) -> dict | None:
    """Δ (candidate − query) for logP / tPSA / MW plus L2 extras."""
    q = fragment_props(query_smiles)
    c = fragment_props(cand_smiles)
    if q is None or c is None:
        return None
    return {
        "dlogp": c["logp"] - q["logp"],
        "dtpsa": c["tpsa"] - q["tpsa"],
        "dmw": c["mw"] - q["mw"],
        # L2 extras
        "dhac": c["hac"] - q["hac"],
        "drings": c["rings"] - q["rings"],
        "dhbd": c["hbd"] - q["hbd"],
        "dhba": c["hba"] - q["hba"],
        "drotb": c["rotb"] - q["rotb"],
    }
