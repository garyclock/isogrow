"""Turning ChEMBL structures (and activities) into an mmpdblib database.

Two builders are provided:

``build_from_chemreps``
    Structures only — produces the SwissBioisostere-equivalent structure layer
    (Frequency + ΔlogP / ΔtPSA / ΔMW).

``build_from_activity``
    Joins a per-molecule activity table (e.g. median pChEMBL from the REST API)
    so ``pchembl`` is indexed as a property too; the resulting database carries
    ΔpChEMBL and therefore supports #Better/#Similar/#Worse and a success score.

Both are thin, testable functions around :func:`isogrow.bioisostere.mmpdb_backend.build`.
"""
from __future__ import annotations

import csv
import gzip
import os
import sqlite3
from collections.abc import Iterable, Iterator
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors

from . import mmpdb_backend as MB


def _opener(path: str):
    return gzip.open if path.endswith(".gz") else open


def iter_smiles_file(path: os.PathLike | str, limit: int = 0, max_heavies: int = 60,
                     min_heavies: int = 1) -> Iterator[tuple]:
    """Yield ``(smiles, id, logP, tPSA, MW)`` from a plain ``SMILES<TAB>id`` file.

    Handy for small custom series (and for the bundled ``tests/data/tiny.smi``),
    which is how you get a database without downloading ChEMBL at all.
    """
    seen: set[str] = set()
    kept = 0
    with open(path, encoding="utf-8", errors="replace") as f:
        for i, line in enumerate(f):
            if limit and kept >= limit:
                break
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            smi = parts[0]
            cid = parts[1] if len(parts) > 1 else f"mol_{i:06d}"
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            n_heavy = mol.GetNumHeavyAtoms()
            if n_heavy > max_heavies or n_heavy < min_heavies:
                continue
            can = Chem.MolToSmiles(mol)
            if can in seen:
                continue
            seen.add(can)
            kept += 1
            yield (can, cid, Crippen.MolLogP(mol),
                   rdMolDescriptors.CalcTPSA(mol), Descriptors.MolWt(mol))


def iter_chemreps(path: str, limit: int = 0, max_heavies: int = 60,
                  min_heavies: int = 1) -> Iterator[tuple]:
    """Yield ``(smiles, chembl_id, logP, tPSA, MW)`` from a ChEMBL chemreps dump.

    Rows are de-duplicated by canonical SMILES, so the database frequency counts
    are not inflated by salt/duplicate entries.
    """
    seen: set[str] = set()
    kept = 0
    with _opener(path)(path, "rt", encoding="utf-8", errors="replace") as f:
        f.readline()  # header
        for line in f:
            if limit and kept >= limit:
                break
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            cid, smi = parts[0], parts[1]
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            n_heavy = mol.GetNumHeavyAtoms()
            if n_heavy > max_heavies or n_heavy < min_heavies:
                continue
            can = Chem.MolToSmiles(mol)
            if can in seen:
                continue
            seen.add(can)
            kept += 1
            yield (can, cid, Crippen.MolLogP(mol),
                   rdMolDescriptors.CalcTPSA(mol), Descriptors.MolWt(mol))


def _build_from_records(records, out_prefix: str, jobs: int, num_cuts: int,
                        mmpdb_exe: str | None) -> Path:
    smi, props = out_prefix + ".smi", out_prefix + ".props"
    Path(out_prefix).parent.mkdir(parents=True, exist_ok=True)
    n = MB.write_smi_and_props(records, smi, props)
    print(f"[build] wrote {n} compounds -> {smi} / {props}")
    if n == 0:
        raise ValueError("no usable structures in the input")

    print(f"[build] mmpdb fragment+index (jobs={jobs}, cuts={num_cuts})")
    MB.build(smi, props, out_prefix + ".mmpdb", mmpdb_exe=mmpdb_exe,
             num_jobs=jobs, symmetric=True, num_cuts=num_cuts)
    db = Path(out_prefix + ".mmpdb")
    print(f"[build] done -> {db} ({db.stat().st_size / 1e6:.1f} MB)")
    return db


def load_alogp(chembl_db: os.PathLike | str, chembl_ids: Iterable[str],
               chunk: int = 900) -> dict[str, float]:
    """``{chembl_id: alogp}`` from a ChEMBL SQLite dump.

    Note: since ChEMBL 24 this column is computed with **RDKit**, not ACD/Labs, so for modern
    releases it equals a locally recomputed Crippen value. It is kept because older dumps
    (before 24) store a genuine ACD logP, and because indexing a stored property avoids
    recomputation. It does *not* reproduce the SwissBioisostere ΔlogP column — see
    ``docs/benchmark.md``.
    """
    ids = list(chembl_ids)
    out: dict[str, float] = {}
    con = sqlite3.connect(f"file:{chembl_db}?mode=ro", uri=True)
    try:
        for start in range(0, len(ids), chunk):
            part = ids[start:start + chunk]
            marks = ",".join("?" * len(part))
            sql = (
                "SELECT md.chembl_id, cp.alogp FROM molecule_dictionary md "
                "JOIN compound_properties cp ON cp.molregno = md.molregno "
                f"WHERE md.chembl_id IN ({marks})"
            )
            for cid, alogp in con.execute(sql, part):
                if alogp is not None:
                    out[cid] = float(alogp)
    finally:
        con.close()
    return out


def _apply_alogp(records: list[tuple], chembl_db: str, stats: dict) -> list[tuple]:
    """Replace the Crippen logP of a batch with ChEMBL's ACD logP."""
    got = load_alogp(chembl_db, [r[1] for r in records])
    out = []
    for smiles, cid, logp, tpsa, mw in records:
        if cid in got:
            logp = got[cid]
            stats["acd"] += 1
        else:
            stats["fallback"] += 1
        out.append((smiles, cid, logp, tpsa, mw))
    return out


def with_alogp(records: Iterable[tuple], chembl_db: os.PathLike | str,
               chunk: int = 900) -> Iterator[tuple]:
    """Stream records, swapping RDKit Crippen logP for ChEMBL ACD logP in batches."""
    db = str(chembl_db)
    stats = {"acd": 0, "fallback": 0}
    batch: list[tuple] = []
    for rec in records:
        batch.append(rec)
        if len(batch) >= chunk:
            yield from _apply_alogp(batch, db, stats)
            batch = []
    if batch:
        yield from _apply_alogp(batch, db, stats)
    print(f"[build] logP source: ChEMBL alogp for {stats['acd']} compounds, "
          f"RDKit Crippen fallback for {stats['fallback']}")


def build_from_chemreps(
    chemreps: os.PathLike | str,
    out_prefix: os.PathLike | str,
    limit: int = 200_000,
    jobs: int = 8,
    num_cuts: int = 3,
    max_heavies: int = 60,
    mmpdb_exe: str | None = None,
    alogp_db: os.PathLike | str | None = None,
) -> Path:
    """Fragment + index a ChEMBL chemreps dump into ``<out_prefix>.mmpdb``.

    Pass ``alogp_db`` (the ChEMBL SQLite dump) to index ChEMBL's ACD logP instead
    of RDKit's Crippen value.
    """
    RDLogger.DisableLog("rdApp.*")
    out_prefix = str(out_prefix)
    print(f"[build] reading {chemreps} (limit={limit}, max_heavies={max_heavies}, "
          f"logP={'ChEMBL alogp' if alogp_db else 'RDKit Crippen'})")
    records = iter_chemreps(str(chemreps), limit, max_heavies)
    if alogp_db:
        records = with_alogp(records, alogp_db)
    return _build_from_records(records, out_prefix, jobs, num_cuts, mmpdb_exe)


def build_from_smiles(
    smiles_path: os.PathLike | str,
    out_prefix: os.PathLike | str,
    limit: int = 0,
    jobs: int = 4,
    num_cuts: int = 3,
    max_heavies: int = 60,
    mmpdb_exe: str | None = None,
) -> Path:
    """Build a database from a plain ``SMILES<TAB>id`` file (no downloads needed)."""
    RDLogger.DisableLog("rdApp.*")
    out_prefix = str(out_prefix)
    print(f"[build] reading {smiles_path} (limit={limit or 'all'})")
    return _build_from_records(
        iter_smiles_file(str(smiles_path), limit, max_heavies), out_prefix, jobs, num_cuts, mmpdb_exe)


def load_activity_medians(path: os.PathLike | str) -> dict[str, float]:
    """Read ``molecule_chembl_id,pchembl`` (extra columns ignored)."""
    out: dict[str, float] = {}
    with open(path, encoding="utf-8") as f:
        for row in csv.DictReader(f):
            try:
                out[row["molecule_chembl_id"]] = float(row["pchembl"])
            except (KeyError, TypeError, ValueError):
                continue
    return out


def scan_chemreps_ids(gz_path: os.PathLike | str, ids: Iterable[str]) -> dict[str, str]:
    """Return ``{chembl_id: smiles}`` for the requested ids (single pass)."""
    wanted = set(ids)
    out: dict[str, str] = {}
    with _opener(str(gz_path))(str(gz_path), "rt", encoding="utf-8", errors="replace") as f:
        f.readline()
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2 and p[0] in wanted and p[0] not in out:
                out[p[0]] = p[1]
                if len(out) == len(wanted):
                    break
    return out


def build_from_activity(
    activity_table: os.PathLike | str,
    chemreps: os.PathLike | str,
    out_prefix: os.PathLike | str,
    jobs: int = 6,
    num_cuts: int = 3,
    mmpdb_exe: str | None = None,
) -> Path:
    """Build an activity-aware database (adds ΔpChEMBL to the statistics)."""
    RDLogger.DisableLog("rdApp.*")
    out_prefix = str(out_prefix)
    Path(out_prefix).parent.mkdir(parents=True, exist_ok=True)

    act = load_activity_medians(activity_table)
    print(f"[build] molecules with activity: {len(act)}")
    smis = scan_chemreps_ids(chemreps, act)
    print(f"[build] matched structures in chemreps: {len(smis)}")

    smi_path, props_path = out_prefix + ".smi", out_prefix + ".props"
    n = 0
    with open(smi_path, "w", encoding="utf-8", newline="\n") as fs, \
            open(props_path, "w", encoding="utf-8", newline="\n") as fp:
        fp.write("ID\t" + "\t".join(MB.PROP_NAMES) + "\tpchembl\n")
        for cid, pch in act.items():
            smi = smis.get(cid)
            if not smi:
                continue
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            fs.write(f"{smi}\t{cid}\n")
            fp.write(f"{cid}\t{Crippen.MolLogP(mol):.4f}\t{rdMolDescriptors.CalcTPSA(mol):.4f}\t"
                     f"{Descriptors.MolWt(mol):.4f}\t{pch:.4f}\n")
            n += 1
    print(f"[build] wrote {n} compounds -> {smi_path}")
    if n == 0:
        raise ValueError("no structures could be joined with activities")

    MB.build(smi_path, props_path, out_prefix + ".mmpdb", mmpdb_exe=mmpdb_exe,
             num_jobs=jobs, symmetric=True, num_cuts=num_cuts)
    db = Path(out_prefix + ".mmpdb")
    print(f"[build] done -> {db} ({db.stat().st_size / 1e6:.1f} MB)")
    return db
