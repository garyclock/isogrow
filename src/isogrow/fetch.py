"""Dataset and tool acquisition.

Nothing large is stored in the repository: the ChEMBL releases, the CReM
fragment database and the docking binaries are fetched by these helpers into
``$ISOGROW_DATA`` (default ``<repo>/data``), with resume support and optional
SHA-256 verification.

Why not ``aria2``/``wget``?  Pure Python keeps the package platform independent
and avoids bundling a GPL binary.  EBI throttles single connections hard
(~10 KB/s in our measurements), so :func:`download` reports throughput and
tells you when to switch to ``aria2c -x16`` for the multi-GB SQLite dump.
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
import time
from pathlib import Path

import requests

from . import paths

CHEMBL_FTP = "https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb"
CHEMBL_API = "https://www.ebi.ac.uk/chembl/api/data"
VINA_RELEASE = "https://github.com/ccsb-scripps/AutoDock-Vina/releases/download"


def _human(n: float) -> str:
    for unit in ("B", "KB", "MB", "GB"):
        if abs(n) < 1024:
            return f"{n:.1f} {unit}"
        n /= 1024
    return f"{n:.1f} TB"


def sha256sum(path: os.PathLike | str, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def download(
    url: str,
    dest: os.PathLike | str,
    expected_sha256: str | None = None,
    resume: bool = True,
    retries: int = 5,
    timeout: int = 90,
) -> Path:
    """Download ``url`` to ``dest``, resuming a partial file when possible.

    Returns the destination path.  Raises :class:`RuntimeError` if all retries
    fail or the checksum does not match.
    """
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)

    if dest.is_file() and expected_sha256 and sha256sum(dest) == expected_sha256:
        print(f"[skip] {dest.name} already present (sha256 ok)")
        return dest

    for attempt in range(1, retries + 1):
        pos = dest.stat().st_size if (resume and dest.is_file()) else 0
        headers = {"Range": f"bytes={pos}-"} if pos else {}
        mode = "ab" if pos else "wb"
        try:
            with requests.get(url, headers=headers, stream=True, timeout=timeout) as r:
                if r.status_code == 416:  # already complete
                    print(f"[ok] {dest.name} already complete")
                    break
                r.raise_for_status()
                if pos and r.status_code == 200:
                    # server ignored the range request -> restart from scratch
                    pos, mode = 0, "wb"
                total = int(r.headers.get("Content-Length", 0)) + pos
                t0, done, last = time.time(), pos, 0.0
                with open(dest, mode) as f:
                    for chunk in r.iter_content(1 << 20):
                        f.write(chunk)
                        done += len(chunk)
                        if time.time() - last > 5:
                            last = time.time()
                            rate = (done - pos) / max(1e-6, time.time() - t0)
                            pct = f"{done / total * 100:5.1f}%" if total else "  ?  "
                            print(f"\r  {dest.name} {pct} {_human(done)}/{_human(total)} "
                                  f"@ {_human(rate)}/s   ", end="", flush=True)
                print()
            if expected_sha256:
                got = sha256sum(dest)
                if got != expected_sha256:
                    raise RuntimeError(
                        f"sha256 mismatch for {dest.name}: {got} != {expected_sha256}"
                    )
            return dest
        except Exception as exc:  # noqa: BLE001 - network errors are broad by nature
            print(f"\n[warn] attempt {attempt}/{retries} failed: {exc}")
            if attempt == retries:
                raise RuntimeError(f"download failed after {retries} attempts: {url}") from exc
            time.sleep(min(30, 2 ** attempt))
    return dest


def chemreps_url(version: str = "37") -> str:
    return f"{CHEMBL_FTP}/latest/chembl_{version}_chemreps.txt.gz"


def sqlite_url(version: str = "37") -> str:
    return f"{CHEMBL_FTP}/latest/chembl_{version}_sqlite.tar.gz"


def fetch_chemreps(version: str = "37", out_dir: os.PathLike | str | None = None) -> Path:
    """ChEMBL ``chemreps`` (all structures, ~280 MB) — the input for building an MMP DB."""
    out = Path(out_dir) if out_dir else paths.data_dir(create=True) / "chembl"
    return download(chemreps_url(version), out / f"chembl_{version}_chemreps.txt.gz")


def fetch_sqlite(version: str = "37", out_dir: os.PathLike | str | None = None) -> Path:
    """Full ChEMBL SQLite dump (~5.4 GB) — needed for ACD logP and activity metadata."""
    out = Path(out_dir) if out_dir else paths.data_dir(create=True) / "chembl"
    return download(sqlite_url(version), out / f"chembl_{version}_sqlite.tar.gz")


# --------------------------------------------------------------------------- #
# ChEMBL REST API sample (small, instant, enough for the activity features)
# --------------------------------------------------------------------------- #
def fetch_activities(target: str = "CHEMBL203", std_type: str = "IC50",
                     max_rows: int = 5000) -> list[dict]:
    """Page through the activity endpoint until ``max_rows`` records are collected."""
    rows: list[dict] = []
    offset = 0
    while len(rows) < max_rows:
        r = requests.get(
            f"{CHEMBL_API}/activity.json",
            params={
                "target_chembl_id": target,
                "standard_type": std_type,
                "pchembl_value__isnull": "false",
                "limit": 1000,
                "offset": offset,
            },
            timeout=90,
        )
        r.raise_for_status()
        payload = r.json()
        acts = payload.get("activities", [])
        if not acts:
            break
        rows.extend(acts)
        offset += len(acts)
        total = payload.get("page_meta", {}).get("total_count", 0)
        print(f"  fetched {len(rows)}/{min(total, max_rows)}")
        if offset >= total:
            break
    return rows[:max_rows]


def save_activity_sample(target: str = "CHEMBL203", std_type: str = "IC50",
                         max_rows: int = 5000,
                         out_dir: os.PathLike | str | None = None) -> dict:
    """Fetch and persist a compact activity sample; returns its metadata."""
    import csv
    import statistics
    from collections import Counter, defaultdict

    out = Path(out_dir) if out_dir else paths.data_dir(create=True) / "chembl_sample"
    out.mkdir(parents=True, exist_ok=True)

    print(f"fetching {std_type} activities for {target} ...")
    acts = fetch_activities(target, std_type, max_rows)
    print(f"total activities: {len(acts)}")

    try:
        tinfo = requests.get(f"{CHEMBL_API}/target/{target}.json", timeout=60).json()
    except Exception:  # noqa: BLE001 - metadata is best-effort
        tinfo = {}

    per_mol: dict[str, list[float]] = defaultdict(list)
    assays: Counter = Counter()
    for a in acts:
        mol, pch = a.get("molecule_chembl_id"), a.get("pchembl_value")
        if mol and pch is not None:
            per_mol[mol].append(float(pch))
            assays[a.get("assay_chembl_id")] += 1

    act_cols = ["molecule_chembl_id", "assay_chembl_id", "standard_type",
                "standard_value", "standard_units", "pchembl_value"]
    with open(out / f"activities_{target}.csv", "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=act_cols, extrasaction="ignore")
        w.writeheader()
        for a in acts:
            w.writerow(a)
    with open(out / f"molecules_{target}.csv", "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["molecule_chembl_id", "pchembl", "n_activities"])
        w.writeheader()
        for mol, vals in per_mol.items():
            w.writerow({"molecule_chembl_id": mol, "pchembl": statistics.median(vals),
                        "n_activities": len(vals)})

    meta = {
        "target": target,
        "target_name": tinfo.get("pref_name"),
        "target_type": tinfo.get("target_type"),
        "std_type": std_type,
        "n_activities": len(acts),
        "n_molecules": len(per_mol),
        "n_assays": len(assays),
        "top_assays": assays.most_common(5),
        "max_rows": max_rows,
        "source": f"{CHEMBL_API}/activity.json",
        "license": "ChEMBL data is CC BY-SA 4.0",
    }
    (out / f"meta_{target}.json").write_text(
        json.dumps(meta, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(meta, ensure_ascii=False, indent=2))
    return meta


# --------------------------------------------------------------------------- #
# External binaries
# --------------------------------------------------------------------------- #
def vina_asset(version: str = "1.2.5") -> str:
    """Return the release asset name for the current platform."""
    if sys.platform.startswith("win"):
        return f"vina_{version}_win.exe"
    if sys.platform == "darwin":
        return f"vina_{version}_mac_x86_64"
    return f"vina_{version}_linux_x86_64"


def fetch_vina(version: str = "1.2.5", out_dir: os.PathLike | str | None = None) -> Path:
    """Download the AutoDock Vina static binary from the official release."""
    out = Path(out_dir) if out_dir else paths.data_dir(create=True) / "bin"
    dest = out / ("vina.exe" if sys.platform.startswith("win") else "vina")
    url = f"{VINA_RELEASE}/v{version}/{vina_asset(version)}"
    download(url, dest)
    if not sys.platform.startswith("win"):
        dest.chmod(0o755)
    print(f"[ok] vina -> {dest}\n     export ISOGROW_VINA={dest}")
    return dest
