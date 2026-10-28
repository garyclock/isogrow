#!/usr/bin/env python
"""Download the full ChEMBL SQLite dump (~5.4 GB compressed).

    python scripts/data/download_chembl_db.py

Only needed for the features that require the relational database: ChEMBL's own
logP (ACD/CX) and the assay/target/PubMed metadata.  The structure layer and the
activity sample (REST API) work without it.

The download resumes from a partial file; delete the ``.aria2``/partial file to
start over.  Expect hours on a throttled EBI link — see ``docs/data.md``.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from isogrow import fetch, paths  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--version", default="37")
    ap.add_argument("--out-dir", default=None, help="default: $ISOGROW_DATA/chembl")
    args = ap.parse_args()

    out = Path(args.out_dir) if args.out_dir else paths.data_dir(create=True) / "chembl"
    path = fetch.fetch_sqlite(args.version, out)
    print(f"\nsqlite dump -> {path}")
    print("unpack with: tar -xzf <file>   # yields chembl_<version>_sqlite/chembl_<version>.db")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
