#!/usr/bin/env python
"""Download the ChEMBL ``chemreps`` dump (all structures, ~280 MB).

    python scripts/data/download_chemreps.py
    python scripts/data/download_chemreps.py --version 37 --out-dir /data/chembl

The file lands in ``$ISOGROW_DATA/chembl/`` (default ``<repo>/data/chembl/``)
and is the input for ``scripts/data/build_mmpdb.py``.

Note: EBI throttles single connections to ~10 KB/s in some networks. If the
transfer crawls, see ``docs/data.md`` for the parallel ``aria2c`` recipe.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from isogrow import fetch, paths  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--version", default="37", help="ChEMBL release (default 37)")
    ap.add_argument("--out-dir", default=None, help="default: $ISOGROW_DATA/chembl")
    args = ap.parse_args()

    out = Path(args.out_dir) if args.out_dir else paths.data_dir(create=True) / "chembl"
    path = fetch.fetch_chemreps(args.version, out)
    print(f"\nchemreps -> {path}")
    print(f"next:\n  python scripts/data/build_mmpdb.py --chemreps {path} --limit 50000")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
