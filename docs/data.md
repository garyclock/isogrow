# Data acquisition

Nothing in this table is committed to the repository. Everything lands in
`$ISOGROW_DATA` (default `<repo>/data`), which is git-ignored.

| Artifact | Size | Command | Needed for |
|---|---|---|---|
| ChEMBL chemreps | ~279 MB | `python scripts/data/download_chemreps.py` | building the structure database |
| ChEMBL SQLite dump | ~5.4 GB (compressed) | `python scripts/data/download_chembl_db.py` | ChEMBL logP, assay/target/PubMed metadata |
| Activity sample (REST) | ~300 KB | `python scripts/data/fetch_activities.py --target CHEMBL203` | #Better/#Similar/#Worse, ΔpChEMBL |
| CReM fragment DB | ~2.5 GB | `python scripts/data/get_crem_db.py --url …` or `cremdb_create` | stage B variants |
| AutoDock Vina | ~1 MB | `python scripts/data/get_vina.py` | stage F docking |
| fkcombu (KCOMBU) | ~1.7 MB | `conda install -c michellab fkcombu` | stage D flexible alignment |

## KCOMBU / fkcombu

**Required by default**: the pipeline ships with `align_backend: kcombu`. If you cannot install
it, run with `--align-backend mcs` (pure-Python rigid fit, no external dependency).

There is **no Windows build** — it is a C program for Linux/macOS:

```bash
conda install -c michellab fkcombu          # or build from pdbj.org/kcombu/download_src.html
export ISOGROW_KCOMBU="$(which fkcombu)"
isogrow doctor                              # should print the resolved path
```

On Windows, install it inside WSL and point the variable at a **command string**; the pipeline
translates `D:\...` paths to `/mnt/d/...` for you:

```bat
set ISOGROW_KCOMBU=wsl.exe -d Ubuntu -- /home/<user>/miniforge3/bin/fkcombu
```

Cost: one process per molecule (roughly 0.2–1 s), so 10k molecules is a 30–90 min step; run it
on the molecules you actually need rather than on a whole screening deck.

## Build sizes and resources

Measured with `--num-cuts 3` (SwissBioisostere's setting), symmetric rules, 6 jobs:

| Compounds | Database | `.fragdb` | Index peak RSS | Wall time |
|---|---|---|---|---|
| 50,000 | 6.5 GB | 0.4 GB | ~9 GB | ~25 min |
| 150,000 | 20 GB+ (est.) | 1.2 GB | > 16 GB (aborted) | — |
| 2.4 M (full) | not practical on a workstation | — | — | — |

Practical guidance: 50k for iteration, 200k+ on a machine with ≥32 GB RAM, and expect the
indexer — not the fragmenter — to be the memory bottleneck. `--num-cuts 1` reduces both size
and time by roughly an order of magnitude but only indexes side-chain replacements.

## Which logP goes into the database?

`isogrow build` indexes RDKit's Crippen logP, tPSA and MW. `--chembl-db <chembl_37.db>` swaps
the logP for ChEMBL's stored `compound_properties.alogp` instead (batched lookups while
streaming chemreps, Crippen fallback for ids the dump does not cover; adds about a minute).

**This does not reproduce the SwissBioisostere ΔlogP column.** ChEMBL stopped using ACD/Labs
for `ALOGP` at release 24 — since then it is computed with RDKit — while the ACD value only
lives in fields that are not part of the public dump. We verified it empirically: the Crippen
build and the alogp build return the identical candidate set and identical deltas. The
`--chembl-db` option is therefore a no-op for ChEMBL ≥ 24 and only matters for older dumps.
See [benchmark.md](benchmark.md) for the measured offset.

## When the download crawls

EBI throttles single connections; we measured ~10 KB/s on one stream versus ~350 KB/s with 16
parallel range requests. The Python downloader resumes from a partial file but uses one
connection. For multi-GB files, either run it repeatedly (it resumes) or use a parallel
downloader:

```bash
aria2c -x16 -s16 -k1M -c \
  -o data/chembl/chembl_37_sqlite.tar.gz \
  "https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_37_sqlite.tar.gz"
```

Keep the `.aria2` control file: it is what makes the transfer resumable.

## Versions

`ChEMBL 37` is the default (`--chembl-version`). The URL pattern is
`…/ChEMBLdb/latest/chembl_<version>_*`, so older releases are reachable by changing the
version string. Databases are **not** portable across ChEMBL versions: rebuild rather than
mixing a 36 database with 37 structures.

## Licences

* ChEMBL — CC BY-SA 4.0. Attribute and do not redistribute the dumps or databases derived from
  them from this repository. Rebuild them locally.
* CReM databases — distributed by the CReM authors / qsar4u; check their terms.
* AutoDock Vina — Apache-2.0 (download from the official release, do not vendor).

Details: [THIRD_PARTY.md](../THIRD_PARTY.md).

## Offline / air-gapped use

Copy `data/` to the target machine and export the paths explicitly:

```bash
export ISOGROW_DATA=/mnt/datasets/isogrow
export ISOGROW_MMPDB=/opt/mmpdb/bin/mmpdb
export ISOGROW_VINA=/opt/vina/bin/vina
isogrow doctor
```

`isogrow doctor` prints exactly what was resolved, which is the fastest way to diagnose a
path problem.
