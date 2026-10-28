"""Local fragment-replacement backend built on mmpdblib (L1 SwissBioisostere replica).

Confirmed mmpdblib (v3.1.4) SQLite schema::

    rule_smiles(id, smiles, num_heavies)        # smiles like '[*:1]c1cc2ccccc2o1'
    rule(id, from_smiles_id, to_smiles_id)
    rule_environment(id, rule_id, environment_fingerprint_id, radius, num_pairs)
    pair(id, rule_environment_id, compound1_id, compound2_id, constant_id)
    rule_environment_statistics(id, rule_environment_id, property_name_id,
        count, avg, std, kurtosis, skewness, min, q1, median, q3, max, paired_t, p_value)
    property_name(id, name)

``rule_environment_statistics.avg`` is the *mean* Δ(property) over all observed
matched pairs for a rule/environment -- exactly SwissBioisostere's "mean of the
differences over all observed MMPs".  mmpdblib therefore reproduces the ΔlogP /
ΔtPSA / ΔMW columns natively when properties are supplied at index time.

Build (build with ``--symmetric`` so querying only ``from_smiles == query``
covers replacements observed in either direction)::

    mmpdb fragment data.smi -o data.fragdb -j N --delimiter whitespace
    mmpdb index data.fragdb -o data.mmpdb --properties data.props --symmetric

where ``data.smi``   = ``SMILES<TAB>ID`` per line, and
      ``data.props`` = ``ID<TAB>logP<TAB>tPSA<TAB>MW`` with a header line.
"""
from __future__ import annotations

import sqlite3
import subprocess
from collections.abc import Iterable

import pandas as pd

PROP_NAMES = ["logP", "tPSA", "MW"]


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def canonical_fragment(smiles: str) -> str | None:
    """Canonicalise a fragment SMILES the way mmpdblib stores it.

    mmpdblib writes the attachment dummy (``[*:1]``) at the *root* of the fragment
    SMILES (e.g. ``[*:1]c1cc2ccccc2o1``), which differs from RDKit's default
    canonical form (``c1ccc2oc([*:1])cc2c1``).  For single-attachment fragments we
    root the canonical SMILES at the dummy so the string matches mmpdblib's.
    Multi-attachment fragments fall back to plain canonical SMILES.
    """
    from rdkit import Chem

    if not smiles or not smiles.strip():
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    dummies = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 0]
    if len(dummies) == 1:
        return Chem.MolToSmiles(mol, rootedAtAtom=dummies[0], canonical=True)
    return Chem.MolToSmiles(mol)


def fragment_forms(smiles: str) -> set[str]:
    """Plausible stored forms of a fragment (raw, canonical, rooted at each dummy).

    mmpdblib's stored ``rule_smiles`` string depends on how it rooted the fragment;
    for single-attachment fragments that is the dummy-rooted canonical form, for
    multi-attachment (linker/scaffold) it is one of the per-dummy rooted forms.  We
    query with the whole set so 1-3 attachment points all match.
    """
    from rdkit import Chem

    forms: set[str] = set()
    if not smiles:
        return forms
    forms.add(smiles)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return forms
    forms.add(Chem.MolToSmiles(mol))
    for d in [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 0]:
        try:
            forms.add(Chem.MolToSmiles(mol, rootedAtAtom=d, canonical=True))
        except Exception:
            pass
    return forms


def write_smi_and_props(records: Iterable[tuple], smi_path: str, props_path: str) -> int:
    """Write paired .smi and .props files.

    ``records`` yields ``(smiles, compound_id, logp, tpsa, mw)``.
    Returns the number of records written.
    """
    n = 0
    with open(smi_path, "w", encoding="utf-8", newline="\n") as fs, open(
        props_path, "w", encoding="utf-8", newline="\n"
    ) as fp:
        fp.write("ID\t" + "\t".join(PROP_NAMES) + "\n")
        for smi, cid, logp, tpsa, mw in records:
            fs.write(f"{smi}\t{cid}\n")
            fp.write(f"{cid}\t{logp:.4f}\t{tpsa:.4f}\t{mw:.4f}\n")
            n += 1
    return n


# --------------------------------------------------------------------------- #
# build
# --------------------------------------------------------------------------- #
def build(
    smi_path: str,
    props_path: str,
    mmpdb_path: str,
    fragdb_path: str | None = None,
    mmpdb_exe: str = "mmpdb",
    num_jobs: int = 8,
    symmetric: bool = True,
    num_cuts: int = 3,
    num_cuts_min_heavies: int = 1,
) -> str:
    """Run ``mmpdb fragment`` then ``mmpdb index``."""
    fragdb_path = fragdb_path or (mmpdb_path + ".fragdb")
    subprocess.run(
        [
            mmpdb_exe, "fragment", smi_path, "-o", fragdb_path,
            "-j", str(num_jobs), "--delimiter", "whitespace",
            "--num-cuts", str(num_cuts),
        ],
        check=True,
    )
    cmd = [
        mmpdb_exe, "index", fragdb_path, "-o", mmpdb_path,
        "--properties", props_path,
    ]
    if symmetric:
        cmd.append("--symmetric")
    subprocess.run(cmd, check=True)
    return mmpdb_path


# --------------------------------------------------------------------------- #
# query
# --------------------------------------------------------------------------- #
_QUERY_SQL = """
SELECT rt.smiles                              AS cand,
       MAX(re.num_pairs)                      AS freq,
       AVG(CASE WHEN pn.name='logP' THEN res.avg END) AS dlogp,
       AVG(CASE WHEN pn.name='tPSA' THEN res.avg END) AS dtpsa,
       AVG(CASE WHEN pn.name='MW'   THEN res.avg END) AS dmw
FROM rule r
JOIN rule_smiles rf ON r.from_smiles_id = rf.id
JOIN rule_smiles rt ON r.to_smiles_id   = rt.id
JOIN rule_environment re ON re.rule_id  = r.id
LEFT JOIN rule_environment_statistics res ON res.rule_environment_id = re.id
LEFT JOIN property_name pn ON res.property_name_id = pn.id
WHERE rf.smiles = ?
GROUP BY r.to_smiles_id
ORDER BY freq DESC
"""


def query(mmpdb_path: str, query_fragment: str, min_pairs: int = 1) -> pd.DataFrame:
    """Return candidate replacements for ``query_fragment``.

    Columns: ``Candidate Fragments, Frequency, Δ logP, Δ tPSA, Δ MW``.
    """
    q = canonical_fragment(query_fragment) or query_fragment
    forms = sorted(fragment_forms(query_fragment) | {q})
    placeholders = ",".join("?" * len(forms))
    sql = _QUERY_SQL.replace("rf.smiles = ?", f"rf.smiles IN ({placeholders})")
    conn = sqlite3.connect(f"file:{mmpdb_path}?mode=ro", uri=True)
    try:
        df = pd.read_sql_query(sql, conn, params=tuple(forms))
    finally:
        conn.close()
    if df.empty:
        return pd.DataFrame(
            columns=["Candidate Fragments", "Frequency", "Δ logP", "Δ tPSA", "Δ MW"]
        )
    df = df[df["freq"] >= min_pairs].copy()
    df = df.rename(
        columns={
            "cand": "Candidate Fragments",
            "freq": "Frequency",
            "dlogp": "Δ logP",
            "dtpsa": "Δ tPSA",
            "dmw": "Δ MW",
        }
    )
    for c in ("Δ logP", "Δ tPSA", "Δ MW"):
        df[c] = df[c].round(2)
    return df.reset_index(drop=True)
