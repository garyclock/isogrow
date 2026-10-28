# IsoGrow

**A local, reproducible replacement for SwissBioisostere — plus the pipeline that turns
fragment replacements into grown, aligned, docked and ranked molecules.**

[![CI](https://github.com/garyclock/isogrow/actions/workflows/ci.yml/badge.svg)](https://github.com/garyclock/isogrow/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](pyproject.toml)

SwissBioisostere is a web service: you submit a fragment, you get replacements. IsoGrow
rebuilds that search locally from matched molecular pairs (MMP) mined from ChEMBL, keeps the
same output layout, and then continues where the web service stops — generating variants,
building 3D conformers, aligning them to a reference pose, docking and ranking.

```text
fragment / scaffold
        │  A. bioisostere search      (local mmpdb, or an official SwissBioisostere CSV)
        ▼
   N replacements
        │  B. CReM mutation           (context-aware, purchasable-space variants)
        ▼
   M variants
        │  C. ETKDG + MMFF            D. MCS alignment      E. Meeko → PDBQT
        ▼
   aligned 3D poses
        │  F. docking (Vina by default, Uni-Dock / external pluggable)
        ▼
   G. ranked table + HTML report
```

## What is in the box

| Layer | Module | What it does |
|---|---|---|
| **Search** | `isogrow.bioisostere` | mmpdb fragmentation/indexing of ChEMBL, fragment → replacements with `Frequency`, `ΔlogP`, `ΔtPSA`, `ΔMW`, `#Better/#Similar/#Worse` and a success score |
| **Growth** | `isogrow.pipeline` | stages A→G above, pure Python (RDKit + Meeko + Vina), no py2, no `sbatch`, no compiled alignment tool |
| **Data** | `scripts/data/` | download and preparation of ChEMBL chemreps, the SQLite dump, activity samples, CReM databases and the Vina binary |
| **CLI** | `isogrow` | `doctor`, `fetch`, `build`, `search`, `run`, `benchmark` |

## Install

```bash
git clone https://github.com/garyclock/isogrow
cd isogrow
python -m venv .venv && source .venv/bin/activate      # Windows: .venv\Scripts\activate
pip install -e ".[all]"
isogrow doctor
```

> **Gotcha worth knowing:** the PyPI package named `mmpdb` installs but exposes **no module**.
> The working distribution is in the [rdkit/mmpdb](https://github.com/rdkit/mmpdb) repository —
> module `mmpdblib`, CLI `mmpdb`. That is why the dependency is written as
> `mmpdb @ git+https://github.com/rdkit/mmpdb.git` and not as a version pin.
> `meeko` additionally needs `gemmi` to import.

## Quickstart — 30 seconds, no downloads

The bundled 15-molecule MMP series is enough to build a real (2 MB) database:

```bash
isogrow build --smi tests/data/tiny.smi --out data/mmpdb/tiny.mmpdb
isogrow search --db data/mmpdb/tiny.mmpdb --query '[*:1]c1cc2ccccc2o1'
```

```text
  Candidate Fragments  Frequency  Δ logP  Δ tPSA   Δ MW
      [*:1]c1ccccc1          3   -0.75   -13.14 -40.02
  [*:1]c1ccc(OC)cc1          2   -0.74    -3.91 -10.00
 [*:1]c1cc2ccccc2s1          2    0.47   -13.14  16.07   ← benzothiophene
        [*:1]c1cccnc1          1   -1.35    -0.25 -39.03
         [*:1]c1ccco1          1   -1.15     0.00 -50.06
```

## Build a searchable database from ChEMBL

One command downloads the structures, fragments them, indexes them and verifies the result:

```bash
pip install -e ".[mmpdb,pipeline]"
python scripts/data/build_database.py --limit 50000 --jobs 6 --out data/mmpdb/chembl50k
# -> [3/3] verification query: [*:1]c1cc2ccccc2o1  ...  search it with: isogrow search ...
```

The equivalent explicit steps:

```bash
python scripts/data/download_chemreps.py                 # 279 MB, all ChEMBL structures
isogrow build --chemreps data/chembl/chembl_37_chemreps.txt.gz \
              --out data/mmpdb/chembl50k --limit 50000 --jobs 6
isogrow search --db data/mmpdb/chembl50k.mmpdb --query '[*:1]c1cc2ccccc2o1' --topn 25
```

Source URLs (the scripts use these; listed so you can fetch manually or mirror them):

| Artifact | URL |
|---|---|
| ChEMBL 37 structures | `https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_37_chemreps.txt.gz` |
| ChEMBL 37 SQLite | `https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/chembl_37_sqlite.tar.gz` |
| Checksums / other releases | `https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/checksums.txt` · `https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/` |
| CReM fragment DB / Vina | `http://www.qsar4u.com/pages/crem.php` · `https://github.com/ccsb-scripps/AutoDock-Vina/releases` |

**No ChEMBL at all?** Build from any SMILES file: `python scripts/data/build_database.py --smi my_series.smi`.
Full details, size/time tables and the "which `--limit`/`--num-cuts`" decision:
[docs/data.md](docs/data.md).

`--num-cuts 3` (the default) matches SwissBioisostere: it indexes side-chain (1 cut),
linker (2 cuts) and scaffold (3 cuts) replacements. 2/3-cut rules are only reachable with
multi-attachment fragments such as `[*:1]C(=O)Oc1ccc([*:2])cc1`.

**About ΔlogP.** We index RDKit's Crippen value. ChEMBL's own `alogp` column has been
computed with RDKit since ChEMBL 24, so indexing it instead changes nothing — we built both
databases and the candidate sets and deltas are identical. SwissBioisostere reports ACD/Labs
logP, which is **not** in the public ChEMBL dump, so our ΔlogP sits **0.50 log units below**
the export on average (sd 0.29, r = 0.92). A +0.5 calibration halves the error (0.51 → 0.22)
but is not applied by default: [docs/benchmark.md](docs/benchmark.md) shows the evidence.

```bash
# optional: index ChEMBL's stored alogp instead of recomputing Crippen
# (identical for ChEMBL >= 24; useful for older dumps where alogp was ACD-computed)
isogrow build --chemreps data/chembl/chembl_37_chemreps.txt.gz \
              --chembl-db data/chembl/chembl_37/chembl_37_sqlite/chembl_37.db \
              --out data/mmpdb/chembl50k --limit 50000 --jobs 6
```

**Size warning.** The 3-cut database grows steeply: 50k compounds ≈ 6.5 GB, 150k ≈ 20 GB+.
Plan disk and RAM accordingly; the indexer peaks around 9 GB on 50k compounds.

### Activity layer (#Better / #Similar / #Worse)

```bash
python scripts/data/fetch_activities.py --target CHEMBL203 --max-rows 5000
isogrow build --chemreps data/chembl/chembl_37_chemreps.txt.gz \
              --activity-table data/chembl_sample/molecules_CHEMBL203.csv \
              --out data/mmpdb/egfr
isogrow search --db data/mmpdb/egfr.mmpdb --query '[*:1]F' --activity
```

This adds `Δ pchembl`, the 0.5-log three-way classification, a success score and a Tanimoto
similarity column. The full ChEMBL SQLite dump (`scripts/data/download_chembl_db.py`) is only
needed for ChEMBL's own logP model and for assay/target/PubMed metadata.

## Full pipeline

```bash
cp config/example_local.json config/my_run.json     # edit paths / query
isogrow run --config config/my_run.json
```

`config/example_official_csv.json` shows the other input mode: feed a SwissBioisostere CSV
export straight into the same downstream stages. Outputs are written per stage
(`1_bioisosteres.csv`, `2_3d/`, `3_aligned/`, `4_pdbqt/`, `5_dock/`, `6_ranked.csv|html`).

Variant budget is allocated **per candidate** (`variants_per_candidate`) rather than as a
global slice, so the ranking does not collapse onto the single most frequent replacement.

## Alignment backends

Stage D uses KCOMBU's flexible fit **by default** (`align_backend` in the config, or
`--align-backend` on `isogrow run`). Two engines are available:

| Backend | What it does | Needs |
|---|---|---|
| `kcombu` (default) | KCOMBU's `fkcombu`: 2D-MCS correspondence, then **flexible** fit — torsion angles stamped from the reference and refined by steepest descent over `nini` initial conformers; can also avoid the receptor (`-P`) | `fkcombu` binary (Linux/macOS, or WSL) |
| `mcs` | rigid-body least-squares on the RDKit MCS — rotates/translates, never changes internal geometry | nothing (pure Python) |

```bash
conda install -c michellab fkcombu            # Linux / macOS
export ISOGROW_KCOMBU=/path/to/fkcombu        # or a WSL command string, see docs/data.md
isogrow run --config config/my_run.json                    # kcombu
isogrow run --config config/my_run.json --align-backend mcs   # no external binary
```

If `fkcombu` is missing the pipeline **fails immediately with install instructions** rather than
silently aligning nothing. CI installs `fkcombu` and Vina and runs the docking/alignment tests,
so the default path is covered by the badge above.

On Windows, drive a WSL install with a command string:

```bash
set ISOGROW_KCOMBU=wsl.exe -d Ubuntu -- /home/<user>/miniforge3/bin/fkcombu
```

**Measured difference (same 6 molecules, YEATS2 case):** docking scores move by at most
**0.06 kcal/mol** between backends (−4.98…−3.98 vs −5.03…−3.92 kcal/mol), because Vina re-optimises
the pose anyway. KCOMBU reports `rmsd_match = 0.16 Å` over *its* correspondence; a Rigid/MCS
RMSD is not comparable across backends since the atom correspondences differ. Use KCOMBU when
you need pocket-compatible conformers (rigid docking, pharmacophores) rather than as a drop-in
accuracy upgrade.

## Benchmark against the official service

Query `[*:1]c1cc2ccccc2o1` (benzofuran), 50k-compound database:

| Metric | Value | Reading |
|---|---|---|
| Precision vs official | **72.3 %** | local candidates are real SwissBioisostere hits |
| Recall@50 | **64 %** | the official top-50 ranking is largely reproduced |
| Coverage of all 1663 candidates | 6.1 % | bounded by database size (50k ≈ 1/48 of ChEMBL) |
| Δ error, tPSA — median / mean | **0.00 / 0.13** | exact; three inconsistent official rows inflate the mean |
| Δ error, MW — median / mean | **0.00 / 1.81** | same three rows (e.g. an official ΔMW of −78.7 for a −4.6 change) |
| Δ error, logP — median / mean | **0.50 / 0.51** | systematic: the service reports ACD/Labs logP, absent from public ChEMBL |

Stage B was verified end-to-end against a previous run of the original notebook: the CReM
variant sets matched **100 %** (50/50 and 53/53 molecules), and the RDKit alignment RMSD was
0.06 Å. Details in [docs/benchmark.md](docs/benchmark.md).

## Limitations (please read before citing numbers)

* **Coverage scales with the database.** A 50k subset reproduces ~6 % of a 2.4M-compound
  resource. Build a bigger database for production use; the method is not the bottleneck.
* **ΔlogP sits 0.50 log units below the export** (sd 0.29, r = 0.92). SwissBioisostere reports
  ACD/Labs logP; ChEMBL's own `alogp` has been RDKit-computed since release 24, so this gap
  cannot be closed with public data (we built both variants to prove it). ΔtPSA and ΔMW agree
  to rounding on 99 of 102 candidates (median error 0.00) — always look at the median, since a
  few inconsistent official rows wreck the mean. `report.json` reports both.
* **The success score is a proxy** — net-improvement rate × log(pairs). SwissBioisostere also
  weights target/target-class/Murcko-scaffold diversity, which needs the full activity tables.
* **Alignment defaults to KCOMBU (flexible)**; `--align-backend mcs` switches to the
  dependency-free rigid fit. py2 AutoDockTools is replaced by Meeko. Molecule identity is
  preserved in every path (heavy-atom counts verified).
* **Docking scores are engine specific.** Vina is the default; Uni-Dock and external engines
  plug in through `isogrow.pipeline.dock.BACKENDS`.

## Repository layout

```text
src/isogrow/
  bioisostere/   engine, mmpdb backend, Δ properties, activity layer, IO, plots, benchmark
  pipeline/      A→G stages + orchestrator
  paths.py       $ISOGROW_DATA / executable resolution
  config.py      JSON config with ${VAR} expansion and validation
  fetch.py       dataset and binary acquisition (resume + checksum)
  cli.py         isogrow doctor|fetch|build|search|run|benchmark
scripts/data/    download + preparation entry points
scripts/dev/     benchmark runner, report renderer, reproduction check
tests/           fixtures + 45 tests (a real mmpdb database is built on the fly)
config/          example pipeline configs
```

## Data policy

No ChEMBL data, no built database and no official SwissBioisostere export is committed.
`tests/data/` holds only ~113 KB of fixtures (a 15-molecule SMILES series, a receptor and a
reference ligand). See [THIRD_PARTY.md](THIRD_PARTY.md) for licenses and attribution.

## Citation & license

Code is MIT licensed. If you use the method, please also cite the SwissBioisostere paper and
ChEMBL — see [CITATION.cff](CITATION.cff) and [THIRD_PARTY.md](THIRD_PARTY.md).

中文说明：[README.zh-CN.md](README.zh-CN.md)
