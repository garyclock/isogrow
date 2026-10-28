#!/usr/bin/env python
"""Download the AutoDock Vina static binary from its official release.

    python scripts/data/get_vina.py            # picks the asset for this platform
    python scripts/data/get_vina.py --version 1.2.5

The binary goes to ``$ISOGROW_DATA/bin/``.  Point the pipeline at it either by
adding that directory to ``PATH`` or by exporting ``ISOGROW_VINA``.
Vina is Apache-2.0; we download rather than vendor it so the repository stays
platform independent and small.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from isogrow import fetch, paths  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--version", default="1.2.5")
    ap.add_argument("--out-dir", default=None, help="default: $ISOGROW_DATA/bin")
    args = ap.parse_args()

    out = Path(args.out_dir) if args.out_dir else paths.data_dir(create=True) / "bin"
    path = fetch.fetch_vina(args.version, out)
    print(f"\nvina -> {path}")
    print("verify with:  isogrow doctor")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
