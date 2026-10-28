"""Benchmark the local engine against an official SwissBioisostere export.

Metrics
-------
* ``coverage_vs_official`` — share of official candidates the local DB reproduces
* ``precision_vs_official`` — share of local candidates that the official run also found
* ``recall_at_N`` — how much of the official top-N ranking is reproduced
* ``delta_mae`` — mean absolute error of ΔlogP / ΔtPSA / ΔMW on the overlap

Coverage is dominated by *database size* (a 50k-compound subset reproduces ~6% of
a 2.4M-compound resource); precision and recall@N measure the method itself.
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import pandas as pd

from .engine import search
from .io import read_sb_csv
from .mmpdb_backend import canonical_fragment


def _canon_map(smiles) -> dict[str, str]:
    return {(canonical_fragment(s) or s): s for s in smiles}


def benchmark(official_csv: os.PathLike | str, mmpdb: os.PathLike | str,
              out_dir: os.PathLike | str, topn: int = 50) -> dict:
    """Run the comparison and write ``report.json`` + ``per_candidate.csv``."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    query, off = read_sb_csv(str(official_csv))
    off = off.copy()
    off["canon"] = [canonical_fragment(s) or s for s in off["Candidate Fragments"]]

    local = search(query, str(mmpdb), topn=None)
    local["canon"] = [canonical_fragment(s) or s for s in local["Candidate Fragments"]]

    off_set, loc_set = set(off["canon"]), set(local["canon"])
    inter = off_set & loc_set
    n_off, n_loc = len(off_set), len(loc_set)

    coverage = len(inter) / n_off if n_off else 0.0
    precision = len(inter) / n_loc if n_loc else 0.0
    off_topn = list(off["canon"].head(topn))
    recall = sum(1 for c in off_topn if c in loc_set) / max(1, len(off_topn))

    deltas = {"dlogp": [], "dtpsa": [], "dmw": []}
    off_idx, loc_idx = off.set_index("canon"), local.set_index("canon")
    for c in inter:
        o, loc = off_idx.loc[c], loc_idx.loc[c]
        try:
            deltas["dlogp"].append(abs(float(loc["Δ logP"]) - float(o["Δ logP"])))
            deltas["dtpsa"].append(abs(float(loc["Δ tPSA"]) - float(o["Δ tPSA"])))
            deltas["dmw"].append(abs(float(loc["Δ MW"]) - float(o["Δ MW"])))
        except (TypeError, ValueError):
            continue

    def mae(a):
        return round(sum(a) / len(a), 3) if a else None

    def median(a):
        return round(sorted(a)[len(a) // 2], 3) if a else None

    def p90(a):
        if not a:
            return None
        s = sorted(a)
        return round(s[min(len(s) - 1, int(round(0.9 * (len(s) - 1))))], 3)

    report = {
        "query": query,
        "mmpdb": str(mmpdb),
        "official_candidates": n_off,
        "local_candidates": n_loc,
        "overlap": len(inter),
        "coverage_vs_official": round(coverage, 4),
        "precision_vs_official": round(precision, 4),
        f"recall_at_{topn}": round(recall, 4),
        "delta_mae": {
            "dlogp": mae(deltas["dlogp"]),
            "dtpsa": mae(deltas["dtpsa"]),
            "dmw": mae(deltas["dmw"]),
            "n_compared": len(deltas["dlogp"]),
        },
        # The mean alone is misleading: a handful of official rows carry values that
        # are inconsistent with their own fragment pair (e.g. ΔMW -78 for a 4 Da
        # change) and dominate the average. The median is the honest summary.
        "delta_median": {k: median(v) for k, v in deltas.items()},
        "delta_p90": {k: p90(v) for k, v in deltas.items()},
        "missing_count": n_off - len(inter),
        "extra_count": n_loc - len(inter),
    }

    rows = []
    for rank, (c, o) in enumerate(off_idx.iterrows(), 1):
        row = {"official_rank": rank, "official_frag": o["Candidate Fragments"],
               "official_freq": o["Frequency"], "in_local": c in loc_set}
        if c in loc_set:
            loc = loc_idx.loc[c]
            row.update({"local_frag": loc["Candidate Fragments"], "local_freq": loc["Frequency"],
                        "off_dlogp": o["Δ logP"], "loc_dlogp": loc["Δ logP"],
                        "off_dtpsa": o["Δ tPSA"], "loc_dtpsa": loc["Δ tPSA"],
                        "off_dmw": o["Δ MW"], "loc_dmw": loc["Δ MW"]})
        rows.append(row)
    pd.DataFrame(rows).to_csv(out / "per_candidate.csv", index=False)
    (out / "report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2),
                                     encoding="utf-8")
    return report
