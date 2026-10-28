#!/usr/bin/env python
"""Build an mmpdb database from a ChEMBL chemreps dump.

Structure layer (SwissBioisostere equivalent)::

    python scripts/data/build_mmpdb.py --chemreps data/chembl/chembl_37_chemreps.txt.gz \\
        --out data/mmpdb/chembl37_sub50k --limit 50000 --jobs 6

Activity-aware layer (adds ΔpChEMBL, enables #Better/#Similar/#Worse)::

    python scripts/data/build_mmpdb.py --chemreps data/chembl/chembl_37_chemreps.txt.gz \\
        --activity-table data/chembl_sample/molecules_CHEMBL203.csv \\
        --out data/mmpdb/egfr_activity

Parameters follow SwissBioisostere: ``--num-cuts 3`` (side chain / linker /
scaffold), and fragments up to 12 heavy atoms are kept by the MP rules encoded in
the indexer.  Database size grows steeply with both ``--limit`` and
``--num-cuts``: 50k compounds at 3 cuts is already ~6.5 GB.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from isogrow import paths  # noqa: E402
from isogrow.bioisostere import chembl  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser()
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--chemreps", help="ChEMBL chemreps dump (.txt or .txt.gz)")
    src.add_argument("--smi", help="plain SMILES<TAB>id file (no download needed)")
    ap.add_argument("--out", default=None, help="output prefix (default: $ISOGROW_DATA/mmpdb/chembl)")
    ap.add_argument("--limit", type=int, default=200_000, help="0 = all compounds")
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--num-cuts", type=int, default=3)
    ap.add_argument("--max-heavies", type=int, default=60)
    ap.add_argument("--activity-table", default=None,
                    help="CSV with molecule_chembl_id,pchembl (see fetch_activities.py)")
    ap.add_argument("--chembl-db", default=None,
                    help="ChEMBL SQLite dump: index its alogp (ACD logP) instead of Crippen")
    ap.add_argument("--mmpdb-exe", default=None)
    args = ap.parse_args()

    prefix = args.out or str(paths.data_dir(create=True) / "mmpdb" / "chembl")
    if prefix.endswith(".mmpdb"):
        prefix = prefix[: -len(".mmpdb")]
    exe = paths.mmpdb_exe(args.mmpdb_exe)

    if args.activity_table:
        db = chembl.build_from_activity(args.activity_table, args.chemreps, prefix,
                                        jobs=args.jobs, num_cuts=args.num_cuts, mmpdb_exe=exe)
    elif args.smi:
        db = chembl.build_from_smiles(args.smi, prefix, limit=args.limit, jobs=args.jobs,
                                      num_cuts=args.num_cuts, max_heavies=args.max_heavies,
                                      mmpdb_exe=exe)
    else:
        db = chembl.build_from_chemreps(args.chemreps, prefix, limit=args.limit,
                                        jobs=args.jobs, num_cuts=args.num_cuts,
                                        max_heavies=args.max_heavies, mmpdb_exe=exe,
                                        alogp_db=args.chembl_db)

    print(f"\nready: {db}\n  isogrow search --db {db} --query '[*:1]c1cc2ccccc2o1'")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
