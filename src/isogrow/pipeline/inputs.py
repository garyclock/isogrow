"""Bioisostere input layer: either the local L1 engine or an official SwissBioisostere CSV.

Both paths return ``(query, DataFrame)`` with the SwissBioisostere column layout, so the
rest of the pipeline is agnostic to where the candidates came from.
"""
from __future__ import annotations

from collections.abc import Sequence

import pandas as pd

from isogrow.bioisostere import engine
from isogrow.bioisostere.io import read_sb_csv


def filter_bioisosteres(
    df: pd.DataFrame,
    logp_range: Sequence[float] | None = None,
    tpsa_range: Sequence[float] | None = None,
    mw_range: Sequence[float] | None = None,
) -> pd.DataFrame:
    """Apply the Δ ranges used by the original notebook (nulls kept)."""
    d = df.copy()

    def _f(col, rng):
        nonlocal d
        if rng is None or col not in d.columns:
            return
        mask = d[col].isna()
        d = d[mask | ((d[col] >= rng[0]) & (d[col] <= rng[1]))].copy()

    _f("Δ logP", logp_range)
    _f("Δ tPSA", tpsa_range)
    _f("Δ MW", mw_range)
    return d.reset_index(drop=True)


def from_official_csv(csv_path: str, **filters) -> tuple[str, pd.DataFrame]:
    """Read + filter an official SwissBioisostere export (same as the notebook's step 1)."""
    query, df = read_sb_csv(csv_path)
    return query, filter_bioisosteres(df, **filters)


def from_local_db(query: str, mmpdb_path: str, topn: int | None = None, **filters):
    """Query the local L1 mmpdblib database."""
    df = engine.search(query, mmpdb_path, topn=topn, **filters)
    return query, df


def get_bioisosteres(
    source: str = "official_csv",
    query: str | None = None,
    csv_path: str | None = None,
    mmpdb_path: str | None = None,
    topn: int | None = None,
    **filters,
) -> tuple[str, pd.DataFrame]:
    """Unified entry point.

    source="official_csv" -> read a SwissBioisostere CSV export
    source="local"        -> query the local mmpdblib DB
    """
    if source == "official_csv":
        if not csv_path:
            raise ValueError("csv_path required for source='official_csv'")
        return from_official_csv(csv_path, **filters)
    if source == "local":
        if not (query and mmpdb_path):
            raise ValueError("query and mmpdb_path required for source='local'")
        return from_local_db(query, mmpdb_path, topn=topn, **filters)
    raise ValueError(f"unknown source: {source}")
