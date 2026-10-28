#!/usr/bin/env python
"""One command from "nothing downloaded" to "a database you can search".

Steps, each skipped when already done:

1. download the ChEMBL ``chemreps`` dump (~280 MB) unless a local file is given
2. fragment + index it into an mmpdb database (optionally activity-aware)
3. run a verification query and print the ``isogrow search`` command to use it

Examples
--------
Small, fast, no download (useful for a smoke test)::

    python scripts/data/build_database.py --smi tests/data/tiny.smi --out data/mmpdb/tiny

The usual iteration-sized database::

    python scripts/data/build_database.py --limit 50000 --jobs 6 --out data/mmpdb/chembl50k

Everything ChEMBL has, side-chain replacements only (see docs/data.md for the cost)::

    python scripts/data/build_database.py --limit 0 --num-cuts 1 --jobs 8 \\
        --out data/mmpdb/chembl_full_cuts1
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from isogrow import fetch, paths  # noqa: E402
from isogrow.bioisostere import chembl, engine  # noqa: E402

DEFAULT_QUERY = "[*:1]c1cc2ccccc2o1"  # benzofuran, the query used throughout the docs


def _size(path: Path) -> str:
    mb = path.stat().st_size / 1e6
    return f"{mb / 1000:.1f} GB" if mb >= 1000 else f"{mb:.1f} MB"


def _resolve_source(args) -> tuple[str, str] | None:
    """Return ``(source_path, out_prefix)`` or None after printing an error."""
    if args.smi:
        print(f"[1/3] source: SMILES file {args.smi} (no download)")
        return args.smi, args.out or str(paths.data_dir(create=True) / "mmpdb" / "smiles")

    if args.chemreps:
        if not Path(args.chemreps).is_file():
            print(f"error: --chemreps {args.chemreps} does not exist")
            return None
        print(f"[1/3] source: {args.chemreps}")
        return args.chemreps, args.out or str(
            paths.data_dir(create=True) / "mmpdb" / f"chembl{args.chembl_version}")

    chemreps = (paths.data_dir(create=True) / "chembl"
                / f"chembl_{args.chembl_version}_chemreps.txt.gz")
    if chemreps.is_file():
        print(f"[1/3] using existing {chemreps} ({_size(chemreps)})")
    elif args.skip_download:
        print(f"error: --skip-download given but {chemreps} does not exist")
        return None
    else:
        print(f"[1/3] downloading ChEMBL {args.chembl_version} chemreps -> {chemreps}")
        fetch.fetch_chemreps(args.chembl_version, chemreps.parent)
    return str(chemreps), args.out or str(
        paths.data_dir(create=True) / "mmpdb" / f"chembl{args.chembl_version}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group()
    src.add_argument("--chemreps", default=None,
                     help="existing chemreps file (skips the download)")
    src.add_argument("--smi", default=None,
                     help="plain SMILES<TAB>id file instead of ChEMBL (no download at all)")
    ap.add_argument("--out", default=None,
                    help="output prefix; default $ISOGROW_DATA/mmpdb/chembl<version>")
    ap.add_argument("--limit", type=int, default=50_000,
                    help="max compounds, 0 = all of ChEMBL (default 50000)")
    ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--num-cuts", type=int, default=3,
                    help="bonds cut at most; SwissBioisostere uses 3 (default), "
                         "1 is far cheaper and covers side-chain replacements only")
    ap.add_argument("--max-heavies", type=int, default=60)
    ap.add_argument("--chembl-version", default="37")
    ap.add_argument("--activity-table", default=None,
                    help="CSV with molecule_chembl_id,pchembl -> adds #Better/#Similar/#Worse")
    ap.add_argument("--chembl-db", default=None,
                    help="ChEMBL SQLite dump: index its stored alogp instead of recomputing")
    ap.add_argument("--skip-download", action="store_true",
                    help="use the chemreps file already in $ISOGROW_DATA/chembl")
    ap.add_argument("--verify", default=DEFAULT_QUERY,
                    help="fragment to test the new database with ('' to skip)")
    ap.add_argument("--mmpdb-exe", default=None)
    args = ap.parse_args()

    t0 = time.time()
    resolved = _resolve_source(args)
    if resolved is None:
        return 2
    source, out = resolved
    if out.endswith(".mmpdb"):
        out = out[: -len(".mmpdb")]
    exe = paths.mmpdb_exe(args.mmpdb_exe)

    scope = "all compounds" if args.limit == 0 else f"limit {args.limit:,}"
    print(f"[2/3] building {out}.mmpdb  ({scope}, {args.jobs} jobs, {args.num_cuts} cuts)")
    if args.activity_table:
        db = chembl.build_from_activity(args.activity_table, source, out,
                                        jobs=args.jobs, num_cuts=args.num_cuts, mmpdb_exe=exe)
    elif args.smi:
        db = chembl.build_from_smiles(source, out, limit=args.limit, jobs=args.jobs,
                                      num_cuts=args.num_cuts, max_heavies=args.max_heavies,
                                      mmpdb_exe=exe)
    else:
        db = chembl.build_from_chemreps(source, out, limit=args.limit, jobs=args.jobs,
                                        num_cuts=args.num_cuts, max_heavies=args.max_heavies,
                                        mmpdb_exe=exe, alogp_db=args.chembl_db)
    print(f"      -> {db}  ({_size(db)})")

    if args.verify:
        print(f"[3/3] verification query: {args.verify}")
        df = engine.search(args.verify, str(db), topn=5)
        if df.empty:
            print("      no replacements found — the database built, but this fragment is "
                  "not covered (normal for a tiny SMILES file)")
        else:
            with pd.option_context("display.width", 160):
                print(df.to_string(index=False))

    print(f"\ndone in {(time.time() - t0) / 60:.1f} min")
    print(f"search it with:\n  isogrow search --db {db} --query '{args.verify or DEFAULT_QUERY}'")
    print("   (add --activity if the database was built with --activity-table)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
