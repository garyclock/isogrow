# Benchmark and validation

All numbers below were produced on this codebase; the raw artefacts are not committed
(they depend on ChEMBL data). Reproduce with the commands at the end of each section.

## 1 · Search layer vs the official service

Query `[*:1]c1cc2ccccc2o1` (benzofuran), official export with **1663** candidates,
local database built from a **50,000-compound** ChEMBL 37 subset (`--num-cuts 3`,
`--symmetric`).

| Metric | Value | Definition |
|---|---|---|
| Official candidates | 1663 | rows in the export (footer text excluded) |
| Local candidates | 141 | rows returned locally |
| Overlap | 102 | same canonical fragment on both sides |
| **Precision vs official** | **72.3 %** | overlap / local — "are our candidates real?" |
| **Recall@50** | **64.0 %** | share of the official top-50 we reproduce |
| Coverage of all candidates | 6.1 % | overlap / official — "how much of the space do we cover?" |
| Missing | 1561 | official-only |
| Extra | 39 | local-only |
| ΔtPSA — median / p90 / mean | **0.00 / 0.00 / 0.13** | deterministic; agreement is rounding |
| ΔMW — median / p90 / mean | **0.00 / 0.01 / 1.81** | see the outlier note below |
| ΔlogP — median / p90 / mean | **0.50 / 0.86 / 0.51** | systematic model offset |
| Δ compared | 102 | candidates present on both sides |

### How to read it

* **Precision and recall@50 measure the method; coverage measures the database.** 72 % of
  our candidates are genuine SwissBioisostere hits and 64 % of the official top-50 ranking is
  reproduced from a database ~1/48 the size of ChEMBL. Coverage will rise with `--limit`;
  precision/recall will not change much.
* **The top of the ranking is right.** Official top replacements — benzothiophene (ΔMW 16.07,
  ΔtPSA −13.14), indole, phenyl, naphthalene, thiophene, anisole — all appear locally with
  property deltas matching to two decimals on tPSA and MW.
* **ΔtPSA and ΔMW are exact for 99 of 102 candidates** (median 0.00, p90 ≤ 0.01). The mean is
  distorted by three official rows that are internally inconsistent, e.g.
  `Clc1ncccc1[*:1]` reports ΔMW −78.71 although swapping a benzofuran for a chloropyridine
  changes MW by ≈ −4.6. **Report the median, not the mean** — this is why `report.json` now
  carries `delta_median` and `delta_p90` next to `delta_mae`.
* **ΔlogP is offset by a constant +0.50** (median exactly 0.500, p90 0.86 — a shift, not
  scatter). We compute RDKit Crippen logP; SwissBioisostere reports ChEMBL's ACD logP. A check
  over 3000 ChEMBL molecules finds `alogp − Crippen`: median 0.00, mean +0.06, |offset| 0.08
  — i.e. the two models agree for most chemotypes but diverge by ~0.5 on aromatic
  heteroatom swaps, which is exactly the class being benchmarked here. Rebuilding the
  database with the `alogp` column from the ChEMBL SQLite dump should remove the offset;
  that rebuild has not been run yet, so treat the ΔlogP column as *model dependent*.

## 2 · Stage B (CReM variants) vs a previous run of the original notebook

Input: the official export filtered with the notebook's Δ ranges (1663 → 1353 candidates).
Reference: 501 retained SDFs from the earlier run (`1_crem_sdf/batch_000`), grouped by ring
system.

| Candidate | Ours | Reference | Overlap | Recall |
|---|---|---|---|---|
| `[*:1]c1cc2ccccc2s1` (benzothiophene) | 50 | 50 | 50 | **100 %** |
| `[*:1]c1cc2ccccc2[nH]1` (indole) | 53 | 53 | 53 | **100 %** |

No missing, no extra molecules. Parameters match the notebook exactly: `radius=2`,
`min_size=0`, `max_size=0`, `min_inc=±2`, `replace_cycles=False`, `protect_H=True`.

## 3 · Stages C–E: does the format chain preserve the molecule?

20 molecules sampled; heavy-atom counts compared across SDF → aligned PDB → PDBQT:

| Check | Result |
|---|---|
| Heavy atoms identical across all three formats | **20/20** |
| RDKit MCS alignment RMSD to the reference pose | **0.06 Å** (0.05–0.07) |

So replacing KCOMBU (C++) with RDKit's MCS alignment, and py2 AutoDockTools with Meeko, does
not alter molecular identity. Pose geometry may still differ slightly from KCOMBU output; that
is a geometry difference, not a chemistry difference.

## 4 · Activity layer spot checks (EGFR sample, 5000 IC50 records)

| Query | Replacements | Example |
|---|---|---|
| `[*:1]F` | 145 | `[*:1]OCc1ccccc1`: 8 pairs, 6 better / 2 similar, ΔpChEMBL +0.79, success 0.716 |
| `[*:1]c1ccccc1` | 233 | 3-chlorophenyl: 21 pairs, 16 better, ΔpChEMBL +0.88, success 1.02 |
| `[*:1]C(=O)Oc1ccc([*:2])cc1` | 274 | multi-attachment (2 dummies) query works |

ΔpChEMBL is `pChEMBL(to) − pChEMBL(from)`; the classification threshold is 0.5 log units,
matching SwissBioisostere.

## Reproducing

```bash
# 1 · data
python scripts/data/download_chemreps.py
isogrow build --chemreps data/chembl/chembl_37_chemreps.txt.gz \
              --out data/mmpdb/chembl50k --limit 50000 --jobs 6

# 2 · compare against your own SwissBioisostere export
python scripts/dev/benchmark.py --official export.csv \
       --db data/mmpdb/chembl50k.mmpdb --out runs/benchmark
python scripts/dev/make_report.py --run runs/benchmark        # HTML summary

# 3 · stage B reproduction against a reference SDF set
python scripts/dev/verify_crem_reproduction.py --official export.csv \
       --crem-db data/crem/chembl.db --reference-dir /path/to/1_crem_sdf
```

Benchmark runs also write `per_candidate.csv`, which lists every official candidate with its
rank, frequency and both sides' Δ values — the fastest way to see where the two disagree.
