# Test fixtures

Small, licence-clean inputs. No ChEMBL dump and no SwissBioisostere export is
stored here — those are regenerated at test time or downloaded by
`scripts/data/*`.

| File | Size | Purpose | Provenance |
|---|---|---|---|
| `tiny.smi` | 461 B | 15-molecule MMP series around benzofuran / benzothiophene / indole / phenyl / naphthalene / thiophene / furan / pyridine / anisole | written for this project |
| `receptor.pdbqt` | 108 KB | docking receptor (AutoDock4 atom types) | prepared from a public PDB structure for the YEATS2 project |
| `reference_ligand.sdf` | 1.5 KB | reference pose: defines the docking box and the alignment target | same project (co-crystallised ligand) |

`tiny.smi` is enough to build a ~2 MB mmpdb database in a few seconds
(`tests/conftest.py`), which then exercises canonicalisation, querying,
Δ-property statistics and — with synthetic pChEMBL values — the
#Better/#Similar/#Worse classification.

The SwissBioisostere-format CSV used by the I/O and pipeline tests is **generated**
at test time from that database (`sb_export_csv` fixture), so the repository
contains no third-party export.
