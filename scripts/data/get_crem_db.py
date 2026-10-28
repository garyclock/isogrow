#!/usr/bin/env python
"""Obtain a CReM fragment database (needed by pipeline stage B).

Two supported routes:

1. **Precompiled ChEMBL database** (what the original project used, ~2.5 GB)::

       python scripts/data/get_crem_db.py --url <direct-link-to.db>

   Precompiled databases are published at <http://www.qsar4u.com/pages/crem.php>.
   The maintainers do not expose a stable direct URL, so pass the download link
   (or the local path of an already-downloaded file) explicitly.

2. **Build your own** from a SMILES file — smaller and faster, at the cost of a
   less diverse fragment set::

       cremdb_create -i data/chembl_sample/chembl.smi -o data/crem/fragments.db -s chembl

Point the pipeline at the result with ``"crem_db"`` in the config (or
``$ISOGROW_CREM_DB``).
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from isogrow import fetch, paths  # noqa: E402

REGISTRY = "http://www.qsar4u.com/pages/crem.php"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--url", default=None, help="direct link to a .db file")
    ap.add_argument("--local", default=None, help="copy an existing .db instead of downloading")
    ap.add_argument("--out", default=None, help="default: $ISOGROW_DATA/crem/chembl.db")
    args = ap.parse_args()

    if not args.url and not args.local:
        print("No --url / --local given.\n")
        print(f"Precompiled CReM databases: {REGISTRY}")
        print("Download one, then either re-run with --local <path> or export ISOGROW_CREM_DB.")
        print("\nOr build a small one yourself:")
        print("  cremdb_create -i <molecules.smi> -o data/crem/fragments.db -s chembl")
        return 2

    dest = Path(args.out) if args.out else paths.data_dir(create=True) / "crem" / "chembl.db"
    if args.local:
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(args.local, dest)
        print(f"[ok] copied -> {dest}")
    else:
        fetch.download(args.url, dest)
        print(f"[ok] downloaded -> {dest}")

    size_gb = dest.stat().st_size / 1e9
    print(f"size: {size_gb:.2f} GB")
    print('set "crem_db" in your pipeline config to this path')
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
