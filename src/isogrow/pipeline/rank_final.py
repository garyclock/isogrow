"""Final scoring / ranking: collect docking affinities and rank candidates."""
from __future__ import annotations

import html
import os

import pandas as pd


def rank_results(records: list[dict], meta: dict | None = None) -> pd.DataFrame:
    """Rank by affinity ascending (more negative = better), then by variant id."""
    df = pd.DataFrame(records)
    if "affinity" not in df.columns:
        return df
    df = df.sort_values("affinity", ascending=True, na_position="last").reset_index(drop=True)
    df.insert(0, "rank", range(1, len(df) + 1))
    return df


def write_csv(df: pd.DataFrame, path: str) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    df.to_csv(path, index=False)


def write_html(df: pd.DataFrame, path: str, title: str = "IsoGrow ranking",
               top: int = 100) -> str:
    cols = [c for c in ("rank", "affinity", "variant", "candidate", "lead") if c in df.columns]
    rows = []
    for _, r in df.head(top).iterrows():
        cells = "".join(f"<td>{html.escape(str(r[c]))}</td>" for c in cols)
        rows.append(f"<tr>{cells}</tr>")
    head = "".join(f"<th>{html.escape(c)}</th>" for c in cols)
    doc = f"""<!doctype html><html lang="zh"><head><meta charset="utf-8">
<title>{html.escape(title)}</title>
<style>body{{font-family:system-ui,"Microsoft YaHei",sans-serif;max-width:1100px;margin:28px auto;padding:0 18px}}
h1{{font-size:22px}} table{{border-collapse:collapse;width:100%;font-size:13px}}
th,td{{border:1px solid #d0d7de;padding:5px 8px;text-align:left}} th{{background:#f6f8fa}}
tr:nth-child(even) td{{background:#fbfcfd}} code{{font-family:ui-monospace,Consolas,monospace}}</style>
</head><body><h1>{html.escape(title)}</h1>
<p>共 {len(df)} 条；按对接亲和力升序（越负越好）。</p>
<table><tr>{head}</tr>{''.join(rows)}</table></body></html>"""
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        f.write(doc)
    return path
