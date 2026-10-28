# Architecture

## Two independent layers

```text
isogrow.bioisostere                     isogrow.pipeline
┌───────────────────────────┐           ┌──────────────────────────────────────┐
│ mmpdb_backend             │           │ inputs      A  candidates            │
│   build / query / props   │           │ crem_variants B  variants            │
│ props    ΔlogP ΔtPSA ΔMW  │           │ embed3d     C  ETKDG + MMFF          │
│ activity #B/#S/#W + score │           │ align3d     D  MCS / KCOMBU align    │
│ engine   search()         │           │ pdb2pdbqt   E  Meeko                 │
│ io       SB-compatible CSV│           │ dock        F  Vina / Uni-Dock / …   │
│ chembl   dataset builders │           │ rank_final  G  CSV + HTML            │
│ benchmark                 │           │ orchestrator  run(config)            │
└───────────────────────────┘           └──────────────────────────────────────┘
            ▲                                          ▲
            └──────────── isogrow.cli ─────────────────┘
                    paths · config · fetch
```

The pipeline depends on the search layer only through `inputs.get_bioisosteres()`, which
returns `(query, DataFrame)` in the SwissBioisostere column layout for both input modes.
That is what makes "official CSV" and "local database" interchangeable.

## Stage contracts

| Stage | Input | Output | Failure behaviour |
|---|---|---|---|
| A | query string, db path **or** CSV | `1_bioisosteres.csv` | raises — nothing downstream is meaningful |
| B | candidate fragments, CReM db | list of variant SMILES | per-candidate, skips molecules CReM rejects |
| C | SMILES | `2_3d/batch_*/mol_*.sdf` | skips molecules that fail embedding |
| D | SDF + reference ligand | `3_aligned/…_aligned.sdf` | unaligned molecules keep their conformer || E | aligned SDF | `4_pdbqt/*.pdbqt` | skips molecules Meeko cannot type |
| F | PDBQT + receptor + box | affinity per ligand | records `error` instead of raising |
| G | docking records | `6_ranked.csv` / `.html` | always produced |

## Deliberate choices

**Fragment canonicalisation is rooted at the dummy atom.** mmpdb stores replacement
fragments as `[*:1]c1cc2ccccc2o1`, while RDKit's default canonical form is
`c1ccc2oc([*:1])cc2c1`. Querying with the default form silently returns nothing, so
`canonical_fragment()` uses `MolToSmiles(rootedAtAtom=dummy)` and `fragment_forms()`
additionally enumerates one rooted form per attachment point for 2–3 cut (linker/scaffold)
fragments.

**Δ properties come from mmpdb statistics, not from a second pass.** Indexing with
`mmpdb index --properties <file> --symmetric` fills `rule_environment_statistics.avg`, which is
exactly SwissBioisostere's "mean difference of the property over all pairs". The query layer
only aggregates; it does not recompute chemistry per candidate.

**No hard-coded paths.** `isogrow.paths` resolves the data directory
(`$ISOGROW_DATA` → `<repo>/data`) and executables (`ISOGROW_VINA`, `ISOGROW_MMPDB` → `PATH`).
`isogrow.config` expands `${VAR}` in JSON configs and rejects undefined variables instead of
silently producing empty strings.

**Docking is an interface, not a call.** `dock.BACKENDS` maps a name to a function with the
signature `(receptor, ligand_pdbqt, center, size, out_pdbqt=…) -> {affinity, …}`. Vina and
Uni-Dock ship in-tree; an external engine is a three-line wrapper.

**Alignment is an interface too.** `align3d.BACKENDS` offers `mcs` (rigid-body fit on the RDKit
MCS, pure Python) and `kcombu` (KCOMBU's flexible fit: 2D-MCS correspondence, torsion stamping,
steepest-descent refinement over `nini` initial conformers, optional receptor avoidance). The
kcombu backend shells out, so its `exe` accepts a **command string** — that is how a Linux-only
binary inside WSL is driven from Windows (`ISOGROW_KCOMBU="wsl.exe -d Ubuntu -- /path/fkcombu"`),
with `D:\…` paths translated to `/mnt/d/…` automatically. Two consequences we had to handle:
`fkcombu`'s SDF omits some hydrogens (they are re-added with coordinates, otherwise Meeko
refuses the molecule) and leaves the 3D flag unset. Its reported `rmsd_match` is defined over
*its own* correspondence, so it is not comparable with the MCS backend's RMSD — compare
backends by docking score.

**Variant sampling is seeded.** CReM draws replacements at random; `crem.seed` is passed through
so a config reproduces its variant set (without it, two runs of the same config diverge).

**Variant budget is per candidate.** A global `variants[:N]` slice starves every candidate but
the first (they are emitted grouped by candidate). The orchestrator instead allocates
`variants_per_candidate` to each, mirroring SwissBioisostere's "rank candidates, then take
top-k" behaviour.

## Test strategy

`tests/conftest.py` builds a real mmpdb database from `tests/data/tiny.smi` once per session
(~2 MB, a few seconds), so the suite exercises the actual indexer, the actual SQL and the
actual RDKit calls without shipping any dataset. A second, activity-aware database is built
with synthetic pChEMBL values to test the three-way classification. Docking is skipped when no
Vina binary is available, which keeps CI green on a clean runner.
