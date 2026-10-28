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

## Search real data

```bash
python scripts/data/download_chemreps.py                 # 279 MB, all ChEMBL structures
isogrow build --chemreps data/chembl/chembl_37_chemreps.txt.gz \
              --out data/mmpdb/chembl50k --limit 50000 --jobs 6
isogrow search --db data/mmpdb/chembl50k.mmpdb --query '[*:1]c1cc2ccccc2o1' --topn 25
```

`--num-cuts 3` (the default) matches SwissBioisostere: it indexes side-chain (1 cut),
linker (2 cuts) and scaffold (3 cuts) replacements. 2/3-cut rules are only reachable with
multi-attachment fragments such as `[*:1]C(=O)Oc1ccc([*:2])cc1`.

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

## Benchmark against the official service

Query `[*:1]c1cc2ccccc2o1` (benzofuran), 50k-compound database:

| Metric | Value | Reading |
|---|---|---|
| Precision vs official | **72.3 %** | local candidates are real SwissBioisostere hits |
| Recall@50 | **64 %** | the official top-50 ranking is largely reproduced |
| Coverage of all 1663 candidates | 6.1 % | bounded by database size (50k ≈ 1/48 of ChEMBL) |
| Δ error, tPSA — median / mean | **0.00 / 0.13** | exact; three inconsistent official rows inflate the mean |
| Δ error, MW — median / mean | **0.00 / 1.81** | same three rows (e.g. an official ΔMW of −78.7 for a −4.6 change) |
| Δ error, logP — median / mean | **0.50 / 0.51** | a systematic offset, not scatter: Crippen vs ChEMBL ACD |

Stage B was verified end-to-end against a previous run of the original notebook: the CReM
variant sets matched **100 %** (50/50 and 53/53 molecules), and the RDKit alignment RMSD was
0.06 Å. Details in [docs/benchmark.md](docs/benchmark.md).

## Limitations (please read before citing numbers)

* **Coverage scales with the database.** A 50k subset reproduces ~6 % of a 2.4M-compound
  resource. Build a bigger database for production use; the method is not the bottleneck.
* **ΔlogP carries a systematic +0.50 offset**, because we use RDKit Crippen while
  SwissBioisostere reports ChEMBL's ACD logP. ΔtPSA and ΔMW agree to rounding on 99 of 102
  candidates (median error 0.00) — always look at the median, since a few inconsistent
  official rows wreck the mean. `report.json` reports both.
* **The success score is a proxy** — net-improvement rate × log(pairs). SwissBioisostere also
  weights target/target-class/Murcko-scaffold diversity, which needs the full activity tables.
* **KCOMBU is replaced by RDKit MCS alignment**, and py2 AutoDockTools by Meeko. Molecule
  identity is preserved (heavy-atom counts verified); pose geometry can differ slightly from
  the original KCOMBU output.
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
