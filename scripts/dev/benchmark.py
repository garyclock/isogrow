#!/usr/bin/env python
"""Compare the local engine with an official SwissBioisostere export.

    python scripts/dev/benchmark.py --official export.csv --db data/mmpdb/chembl.mmpdb \\
        --out runs/benchmark

Writes ``report.json`` and ``per_candidate.csv``; render a readable page with
``scripts/dev/make_report.py --run runs/benchmark``.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from isogrow.bioisostere.benchmark import benchmark  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--official", required=True, help="SwissBioisostere CSV export")
    ap.add_argument("--db", required=True, help="local .mmpdb")
    ap.add_argument("--out", required=True, help="output directory")
    ap.add_argument("--topn", type=int, default=50)
    args = ap.parse_args()

    report = benchmark(args.official, args.db, args.out, topn=args.topn)
    print(json.dumps(report, ensure_ascii=False, indent=2))
    print(f"\nwritten -> {Path(args.out) / 'report.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
