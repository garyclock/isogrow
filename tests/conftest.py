"""Shared fixtures: a tiny mmpdb database built on the fly.

The repository deliberately ships no ChEMBL data.  Instead, a 15-molecule
matched-molecular-pair series (``tests/data/tiny.smi``) is fragmented and
indexed into a ~2 MB database once per session.  That is enough to exercise
canonicalisation, querying, Δ-property statistics and the activity
classification, and it runs in a few seconds.
"""
from __future__ import annotations

import os
import shutil
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
DATA = Path(__file__).parent / "data"
sys.path.insert(0, str(REPO / "src"))


def _mmpdb_exe() -> str | None:
    return os.environ.get("ISOGROW_MMPDB") or shutil.which("mmpdb")


def _write_smi_props(smi_out: Path, props_out: Path, extra_props: dict[str, float] | None = None):
    """Convert ``tiny.smi`` into the .smi/.props pair mmpdb expects."""
    from rdkit import Chem
    from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors

    header = ["logP", "tPSA", "MW"] + (["pchembl"] if extra_props else [])
    rows = []
    for line in (DATA / "tiny.smi").read_text(encoding="utf-8").splitlines():
        parts = line.split()
        if len(parts) < 2:
            continue
        mol = Chem.MolFromSmiles(parts[0])
        if mol is None:
            continue
        vals = [Crippen.MolLogP(mol), rdMolDescriptors.CalcTPSA(mol), Descriptors.MolWt(mol)]
        if extra_props:
            vals.append(extra_props.get(parts[1], 5.0))
        rows.append((Chem.MolToSmiles(mol), parts[1], vals))

    with open(smi_out, "w", encoding="utf-8", newline="\n") as fs, \
            open(props_out, "w", encoding="utf-8", newline="\n") as fp:
        fp.write("ID\t" + "\t".join(header) + "\n")
        for smiles, cid, vals in rows:
            fs.write(f"{smiles}\t{cid}\n")
            fp.write(f"{cid}\t" + "\t".join(f"{v:.4f}" for v in vals) + "\n")
    return len(rows)


@pytest.fixture(scope="session")
def mmpdb_exe() -> str:
    exe = _mmpdb_exe()
    if not exe:
        pytest.skip("mmpdb command line tool not installed "
                    '(pip install "mmpdb @ git+https://github.com/rdkit/mmpdb.git")')
    return exe


@pytest.fixture(scope="session")
def tiny_mmpdb(tmp_path_factory, mmpdb_exe) -> Path:
    """Structure-only database (Frequency + ΔlogP/ΔtPSA/ΔMW)."""
    from isogrow.bioisostere import mmpdb_backend as MB

    tmp = tmp_path_factory.mktemp("tiny_db")
    smi, props, db = tmp / "tiny.smi", tmp / "tiny.props", tmp / "tiny.mmpdb"
    _write_smi_props(smi, props)
    MB.build(str(smi), str(props), str(db), mmpdb_exe=mmpdb_exe,
             num_jobs=2, symmetric=True, num_cuts=3)
    return db


@pytest.fixture(scope="session")
def tiny_activity_db(tmp_path_factory, mmpdb_exe) -> Path:
    """Activity-aware database: artificial pChEMBL values on the same series.

    ``m2`` (benzothiophene analogue) is set one log unit better than ``m1``
    (benzofuran), which gives a deterministic ΔpChEMBL = +1 for that replacement.
    """
    from isogrow.bioisostere import mmpdb_backend as MB

    scores = {"m1": 5.0, "m2": 6.0, "m3": 6.0, "m4": 4.5, "m5": 4.0,
              "m6": 5.0, "m7": 4.0, "m8": 5.5, "m9": 5.0}
    tmp = tmp_path_factory.mktemp("tiny_act_db")
    smi, props, db = tmp / "act.smi", tmp / "act.props", tmp / "act.mmpdb"
    _write_smi_props(smi, props, extra_props=scores)
    MB.build(str(smi), str(props), str(db), mmpdb_exe=mmpdb_exe,
             num_jobs=2, symmetric=True, num_cuts=3)
    return db


@pytest.fixture(scope="session")
def sb_export_csv(tiny_mmpdb, tmp_path_factory) -> Path:
    """A SwissBioisostere-format CSV generated locally (we ship none)."""
    from isogrow.bioisostere import engine
    from isogrow.bioisostere import io as sbio

    query = "[*:1]c1cc2ccccc2o1"
    df = engine.search(query, str(tiny_mmpdb))
    path = tmp_path_factory.mktemp("sb_export") / "export.csv"
    sbio.write_sb_csv(df, str(path), query)
    return path


@pytest.fixture(scope="session")
def receptor() -> Path:
    return DATA / "receptor.pdbqt"


@pytest.fixture(scope="session")
def reference_ligand() -> Path:
    return DATA / "reference_ligand.sdf"


@pytest.fixture(scope="session")
def vina():
    """Path to a working Vina executable, or skip."""
    from isogrow import paths

    try:
        return paths.vina_exe()
    except FileNotFoundError:
        pytest.skip("AutoDock Vina not installed (set ISOGROW_VINA or run scripts/data/get_vina.py)")
