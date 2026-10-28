"""Unified search interface producing SwissBioisostere-compatible results."""
from __future__ import annotations

from collections.abc import Sequence

import pandas as pd

from . import mmpdb_backend as MB
from . import props as P
from .io import CANON_COLS

L2_COLS = ["ΔHAC", "ΔRings", "ΔHBD", "ΔHBA", "ΔRotB"]


def search(
    query_fragment: str,
    mmpdb_path: str,
    logp_range: Sequence[float] | None = None,
    tpsa_range: Sequence[float] | None = None,
    mw_range: Sequence[float] | None = None,
    min_pairs: int = 1,
    topn: int | None = None,
    add_l2: bool = False,
) -> pd.DataFrame:
    """Search observed replacements of ``query_fragment`` in a local mmpdblib DB.

    Returns a dataframe with the SwissBioisostere column layout
    (``CANON_COLS``) and, if ``add_l2``, extra property deltas.
    """
    df = MB.query(mmpdb_path, query_fragment, min_pairs=min_pairs)
    if df.empty:
        for c in CANON_COLS:
            if c not in df.columns:
                df[c] = None
        return df[CANON_COLS]

    # drop self-replacement
    qcanon = MB.canonical_fragment(query_fragment)
    if qcanon is not None:
        df = df[df["Candidate Fragments"] != qcanon].copy()

    # property filters (only on non-null rows)
    def _filt(col, rng):
        nonlocal df
        if rng is None or col not in df.columns:
            return
        mask = df[col].notna()
        keep = (df[col] >= rng[0]) & (df[col] <= rng[1])
        df = df[~mask | keep].copy()

    _filt("Δ logP", logp_range)
    _filt("Δ tPSA", tpsa_range)
    _filt("Δ MW", mw_range)

    # bioactivity classification needs ChEMBL assay deltas (available once the
    # full SQLite is imported); left null for the structure-only L1 build.
    for c in ("# Better", "# Similar", "# Worse"):
        df[c] = None

    if add_l2:
        rows = []
        for cand in df["Candidate Fragments"]:
            d = P.delta_props(query_fragment, cand) or {}
            rows.append((d.get("dhac"), d.get("drings"), d.get("dhbd"),
                         d.get("dhba"), d.get("drotb")))
        for i, name in enumerate(L2_COLS):
            df[name] = [r[i] for r in rows]

    df = df.sort_values("Frequency", ascending=False).reset_index(drop=True)
    if topn:
        df = df.head(topn).reset_index(drop=True)

    cols = CANON_COLS + (L2_COLS if add_l2 else [])
    return df[cols]
