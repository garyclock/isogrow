#!/usr/bin/env python
"""Fetch a compact ChEMBL activity sample through the REST API.

    python scripts/data/fetch_activities.py --target CHEMBL203 --max-rows 5000

Outputs (``$ISOGROW_DATA/chembl_sample/``):

* ``activities_<target>.csv`` — one row per measurement
* ``molecules_<target>.csv``  — per-molecule median pChEMBL
* ``meta_<target>.json``      — target info and sample statistics

The sample is enough to exercise the activity layer (#Better/#Similar/#Worse,
ΔpChEMBL, success score) without the multi-GB SQLite dump.  ChEMBL data is
licensed CC BY-SA 4.0 — keep the attribution when redistributing.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from isogrow import fetch, paths  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", default="CHEMBL203", help="target ChEMBL id (default: EGFR)")
    ap.add_argument("--std-type", default="IC50")
    ap.add_argument("--max-rows", type=int, default=5000)
    ap.add_argument("--out-dir", default=None, help="default: $ISOGROW_DATA/chembl_sample")
    args = ap.parse_args()

    out = Path(args.out_dir) if args.out_dir else paths.data_dir(create=True) / "chembl_sample"
    fetch.save_activity_sample(args.target, args.std_type, args.max_rows, out)
    print(f"\nnext:\n  python scripts/data/build_mmpdb.py --chemreps <chemreps> "
          f"--activity-table {out / f'molecules_{args.target}.csv'} --out data/mmpdb/{args.target}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
