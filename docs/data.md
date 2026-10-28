# Data acquisition

Nothing in this table is committed to the repository. Everything lands in
`$ISOGROW_DATA` (default `<repo>/data`), which is git-ignored.

## Download links

| Artifact | Direct link | Size |
|---|---|---|
| ChEMBL 37 chemreps (all structures) | https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_37_chemreps.txt.gz | 279 MB |
| ChEMBL 37 SQLite dump | https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_37_sqlite.tar.gz | 5.4 GB |
| ChEMBL 37 SDF (not used here) | https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_37.sdf.gz | 888 MB |
| Official checksums (SHA-256) | https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/checksums.txt | — |
| Older/newer releases | https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/ | — |
| CReM precompiled fragment DB | http://www.qsar4u.com/pages/crem.php | ~2.5 GB |
| AutoDock Vina binaries | https://github.com/ccsb-scripps/AutoDock-Vina/releases | 1 MB |
| fkcombu (KCOMBU) conda channel | https://conda.anaconda.org/michellab/linux-64/ + https://pdbj.org/kcombu/download_src.html | 1.7 MB |
| SwissBioisostere (reference service) | https://www.swissbioisostere.ch | — |

The `latest/` path always points at the current release; replace `37` with the version you want
(e.g. `chembl_36_chemreps.txt.gz`) or browse the release index. **The scripts already know these
URLs** — you normally do not type them.

| Artifact | Command | Needed for |
|---|---|---|
| ChEMBL chemreps | `python scripts/data/download_chemreps.py` | building the structure database |
| ChEMBL SQLite dump | `python scripts/data/download_chembl_db.py` | stored `alogp`, assay/target/PubMed metadata |
| Activity sample (REST) | `python scripts/data/fetch_activities.py --target CHEMBL203` | #Better/#Similar/#Worse, ΔpChEMBL |
| CReM fragment DB | `python scripts/data/get_crem_db.py --url …` or `cremdb_create` | stage B variants |
| AutoDock Vina | `python scripts/data/get_vina.py` | stage F docking |
| fkcombu (KCOMBU) | `conda install -c michellab fkcombu` | stage D flexible alignment |

## From zero to a searchable database

One command does download → fragment → index → verify:

```bash
pip install -e ".[mmpdb,pipeline]"          # provides the `mmpdb` CLI

# 50k compounds, ~6.5 GB database, ~25 min on 6 cores
python scripts/data/build_database.py --limit 50000 --jobs 6 --out data/mmpdb/chembl50k
```

It prints what it did and how to query it:

```text
[1/3] downloading ChEMBL 37 chemreps -> data/chembl/chembl_37_chemreps.txt.gz
[2/3] building data/mmpdb/chembl50k.mmpdb  (limit 50,000, 6 jobs, 3 cuts)
      -> data/mmpdb/chembl50k.mmpdb  (6.5 GB)
[3/3] verification query: [*:1]c1cc2ccccc2o1
  Candidate Fragments  Frequency   Δ logP  Δ tPSA   Δ MW
      [*:1]c1cc2ccccc2s1          2     0.47  -13.14  16.07
      [*:1]c1ccc(OC)cc1           2    -0.74   -3.91 -10.00
      ...

search it with:
  isogrow search --db data/mmpdb/chembl50k.mmpdb --query '[*:1]c1cc2ccccc2o1'
```

Under the hood that is three steps you can also run separately:

```bash
python scripts/data/download_chemreps.py                       # 1. 279 MB, resumable
isogrow build --chemreps data/chembl/chembl_37_chemreps.txt.gz \
              --limit 50000 --jobs 6 --out data/mmpdb/chembl50k  # 2. fragment + index
isogrow search --db data/mmpdb/chembl50k.mmpdb \
               --query '[*:1]c1cc2ccccc2o1' --topn 25              # 3. query
```

Then use it for real work:

```bash
isogrow benchmark --official export.csv --db data/mmpdb/chembl50k.mmpdb --out runs/bench
isogrow run --config config/example_local.json        # full A→G pipeline
```

### Which size should I build?

| Goal | `--limit` | `--num-cuts` | Database | Wall time | RAM peak |
|---|---|---|---|---|---|
| Smoke test / CI | `--smi tests/data/tiny.smi` | 3 | 2 MB | seconds | <1 GB |
| Iterate on queries | 50,000 | 3 | 6.5 GB | ~25 min | ~9 GB |
| Production-ish | 200,000 | 3 | ~20 GB (est.) | hours | ≥32 GB |
| All of ChEMBL | 0 | 3 | not practical on a workstation | — | — |
| All of ChEMBL | 0 | 1 | much smaller, side-chain replacements only | hours | depends |

Cost scales steeply with `--num-cuts` because the number of matched pairs grows super-linearly:
50k compounds at 3 cuts already produce a 6.5 GB database, and the indexer — not the fragmenter —
is the memory bottleneck. Build a subset you can afford, and note that **coverage of a query's
replacements is proportional to the number of compounds indexed** (a 50k subset reproduced 6.1 %
of the official benzofuran candidate list; see [benchmark.md](benchmark.md)).

### If you only have a SMILES file

```bash
python scripts/data/build_database.py --smi my_series.smi --out data/mmpdb/my_series
```

No ChEMBL download at all — useful for a company series or a quick experiment.


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

## Measured build costs

`--num-cuts 3` (SwissBioisostere's setting), symmetric rules, 6 jobs, on a 24-core/64 GB box:

| Compounds | Database | `.fragdb` | Index peak RSS | Wall time |
|---|---|---|---|---|
| 50,000 | 6.5 GB | 0.4 GB | ~9 GB | ~25 min |
| 150,000 | 20 GB+ (est.) | 1.2 GB | > 16 GB (aborted) | — |
| 2.4 M (full) | not practical on a workstation | — | — | — |

`--num-cuts 1` cuts both size and time by roughly an order of magnitude but only indexes
side-chain (single attachment point) replacements.

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
