"""Activity-aware query: #Better/#Similar/#Worse, Δ properties and a success score.

Uses an mmpdb database built with a ``pchembl`` property (see
``scripts/data/build_mmpdb.py --activity-table ...``).  For a query fragment we
take every rule whose ``from_smiles`` equals the query, dedupe the underlying
matched pairs, and classify each pair by Δpchembl = pchembl(to) - pchembl(from):

    better  : Δpchembl >= +threshold   (candidate more potent)
    similar : |Δpchembl| <  threshold
    worse   : Δpchembl <= -threshold

threshold defaults to 0.5 log units, matching SwissBioisostere.
"""
from __future__ import annotations

import math
import sqlite3
from collections import defaultdict

import pandas as pd

from .mmpdb_backend import fragment_forms

PROP_LABELS = {"logP": "Δ logP", "tPSA": "Δ tPSA", "MW": "Δ MW", "pchembl": "Δ pchembl"}


def _prop_ids(cur) -> dict[str, int]:
    cur.execute("SELECT name, id FROM property_name")
    return {n: i for n, i in cur.fetchall()}


def _from_ids(cur, query: str) -> list[int]:
    forms = fragment_forms(query) or {query}
    placeholders = ",".join("?" * len(forms))
    cur.execute(f"SELECT id FROM rule_smiles WHERE smiles IN ({placeholders})", tuple(forms))
    return [r[0] for r in cur.fetchall()]


def query_activity(mmpdb_path: str, query: str, threshold: float = 0.5, min_pairs: int = 1) -> pd.DataFrame:
    conn = sqlite3.connect(f"file:{mmpdb_path}?mode=ro", uri=True)
    cur = conn.cursor()
    pids = _prop_ids(cur)
    if "pchembl" not in pids:
        conn.close()
        raise ValueError(
            "database has no 'pchembl' property; build it with "
            "scripts/data/build_mmpdb.py --activity-table <molecules.csv>")

    ids = _from_ids(cur, query)
    if not ids:
        conn.close()
        return pd.DataFrame(columns=["Candidate Fragments", "Frequency", "# Better", "# Similar",
                                     "# Worse", "Δ logP", "Δ tPSA", "Δ MW", "Δ pchembl", "success_score"])

    qmarks = ",".join("?" * len(ids))
    sql = f"""
    WITH up AS (
      SELECT DISTINCT re.rule_id AS rule_id, p.compound1_id AS c1, p.compound2_id AS c2
      FROM pair p JOIN rule_environment re ON p.rule_environment_id = re.id
      WHERE re.rule_id IN (SELECT id FROM rule WHERE from_smiles_id IN ({qmarks}))
    )
    SELECT rt.smiles, up.c1, up.c2
    FROM up JOIN rule r ON r.id = up.rule_id
    JOIN rule_smiles rt ON r.to_smiles_id = rt.id
    """
    rows = cur.execute(sql, ids).fetchall()

    comp = set()
    for _, c1, c2 in rows:
        comp.add(c1)
        comp.add(c2)
    prop_by_comp: dict[int, dict[int, float]] = defaultdict(dict)
    if comp:
        cm = ",".join("?" * len(comp))
        for cid, pid, val in cur.execute(
            f"SELECT compound_id, property_name_id, value FROM compound_property WHERE compound_id IN ({cm})",
            tuple(comp),
        ):
            prop_by_comp[cid][pid] = val
    conn.close()

    pid_pch = pids["pchembl"]
    agg = defaultdict(lambda: {"n": 0, "vals": defaultdict(list), "better": 0, "similar": 0, "worse": 0})
    for cand, c1, c2 in rows:
        p1, p2 = prop_by_comp.get(c1), prop_by_comp.get(c2)
        if not p1 or not p2 or pid_pch not in p1 or pid_pch not in p2:
            continue
        a = agg[cand]
        a["n"] += 1
        for name in PROP_LABELS:
            pid = pids.get(name)
            if pid in p1 and pid in p2:
                a["vals"][name].append(p2[pid] - p1[pid])
        dp = p2[pid_pch] - p1[pid_pch]
        if dp >= threshold:
            a["better"] += 1
        elif dp <= -threshold:
            a["worse"] += 1
        else:
            a["similar"] += 1

    out = []
    for cand, a in agg.items():
        n = a["n"]
        if n < min_pairs:
            continue
        rec = {"Candidate Fragments": cand, "Frequency": n,
               "# Better": a["better"], "# Similar": a["similar"], "# Worse": a["worse"]}
        for name, label in PROP_LABELS.items():
            v = a["vals"].get(name, [])
            rec[label] = round(sum(v) / len(v), 2) if v else None
        # proxy for SwissBioisostere's success score: net-improvement rate x log-count.
        # (SB also weights target/scaffold diversity, which a single-target sample cannot supply.)
        rec["success_score"] = round((a["better"] - a["worse"]) / max(1, n) * math.log10(1 + n), 3)
        out.append(rec)

    df = pd.DataFrame(out)
    if df.empty:
        return df
    return df.sort_values(["success_score", "Frequency"], ascending=False).reset_index(drop=True)


def query_pair(mmpdb_path: str, from_smiles: str, to_smiles: str, threshold: float = 0.5) -> dict | None:
    """Statistics for a *specific* replacement A -> B (SwissBioisostere's 2nd query type)."""
    conn = sqlite3.connect(f"file:{mmpdb_path}?mode=ro", uri=True)
    cur = conn.cursor()
    pids = _prop_ids(cur)
    if "pchembl" not in pids:
        conn.close()
        raise ValueError("database has no 'pchembl' property")
    fids = _from_ids(cur, from_smiles)
    tids = _from_ids(cur, to_smiles)
    if not fids or not tids:
        conn.close()
        return None
    fm = ",".join("?" * len(fids))
    tm = ",".join("?" * len(tids))
    pair_sql = f"""
    SELECT DISTINCT p.compound1_id, p.compound2_id
    FROM rule r JOIN rule_environment re ON re.rule_id = r.id
    JOIN pair p ON p.rule_environment_id = re.id
    WHERE r.from_smiles_id IN ({fm}) AND r.to_smiles_id IN ({tm})
    """
    pairs = cur.execute(pair_sql, fids + tids).fetchall()
    if not pairs:
        conn.close()
        return None
    comp = {c for pr in pairs for c in pr}
    cm = ",".join("?" * len(comp))
    prop_by_comp: dict[int, dict[int, float]] = defaultdict(dict)
    for cid, pid, val in cur.execute(
        f"SELECT compound_id, property_name_id, value FROM compound_property WHERE compound_id IN ({cm})",
        tuple(comp),
    ):
        prop_by_comp[cid][pid] = val
    conn.close()

    pid_pch = pids["pchembl"]
    vals = defaultdict(list)
    better = similar = worse = 0
    n = 0
    for c1, c2 in pairs:
        p1, p2 = prop_by_comp.get(c1), prop_by_comp.get(c2)
        if not p1 or not p2 or pid_pch not in p1 or pid_pch not in p2:
            continue
        n += 1
        for name in PROP_LABELS:
            pid = pids.get(name)
            if pid in p1 and pid in p2:
                vals[name].append(p2[pid] - p1[pid])
        dp = p2[pid_pch] - p1[pid_pch]
        if dp >= threshold:
            better += 1
        elif dp <= -threshold:
            worse += 1
        else:
            similar += 1
    res = {"from": from_smiles, "to": to_smiles, "Frequency": n,
           "# Better": better, "# Similar": similar, "# Worse": worse}
    for name, label in PROP_LABELS.items():
        v = vals.get(name, [])
        res[label] = round(sum(v) / len(v), 2) if v else None
    res["success_score"] = round((better - worse) / max(1, n) * math.log10(1 + n), 3)
    return res


def add_similarity(df: pd.DataFrame, query_smiles: str, column: str = "similarity") -> pd.DataFrame:
    """Add a Tanimoto similarity (Morgan FP on H-capped fragments) vs the query."""
    from rdkit import Chem, DataStructs, RDLogger
    from rdkit.Chem import rdFingerprintGenerator

    from .props import cap_dummy

    RDLogger.DisableLog("rdApp.*")
    qmol = Chem.MolFromSmiles(cap_dummy(query_smiles))
    if qmol is None or df.empty:
        return df
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    qfp = gen.GetFingerprint(qmol)
    sims = []
    for s in df["Candidate Fragments"]:
        m = Chem.MolFromSmiles(cap_dummy(s))
        if m is None:
            sims.append(None)
            continue
        fp = gen.GetFingerprint(m)
        sims.append(round(DataStructs.TanimotoSimilarity(qfp, fp), 3))
    df = df.copy()
    df[column] = sims
    return df
