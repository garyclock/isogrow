# Third-party software, data and attribution

IsoGrow is MIT licensed, but it depends on third-party packages and can consume
third-party datasets whose licenses differ. This file records them so users can
comply. **No third-party dataset is redistributed in this repository.**

## Software dependencies

| Package | License | Used for |
|---|---|---|
| [RDKit](https://github.com/rdkit/rdkit) | BSD-3-Clause | parsing, canonicalisation, properties, conformers, alignment |
| [mmpdb](https://github.com/rdkit/mmpdb) | MIT | matched molecular pair fragmentation, indexing, statistics |
| [CReM](https://github.com/DrrDom/crem) | BSD-3-Clause | context-aware fragment mutation (stage B) |
| [Meeko](https://github.com/forlilab/Meeko) | BSD-3-Clause | PDB → PDBQT conversion (stage E) |
| [gemmi](https://github.com/project-gemmi/gemmi) | MPL-2.0 | structure I/O dependency of Meeko |
| [AutoDock Vina](https://github.com/ccsb-scripps/AutoDock-Vina) | Apache-2.0 | docking (downloaded, not vendored) |
| [pandas](https://github.com/pandas-dev/pandas) / [NumPy](https://numpy.org) | BSD-3-Clause | tabular data |
| [requests](https://github.com/psf/requests) | Apache-2.0 | dataset downloads |

## Data

| Dataset | License | How it is obtained |
|---|---|---|
| [ChEMBL](https://www.ebi.ac.uk/chembl/) structures (`chemreps`) | CC BY-SA 4.0 | `scripts/data/download_chemreps.py` |
| ChEMBL SQLite dump | CC BY-SA 4.0 | `scripts/data/download_chembl_db.py` |
| ChEMBL activity records (REST API) | CC BY-SA 4.0 | `scripts/data/fetch_activities.py` |
| CReM precompiled fragment databases | see CReM / qsar4u terms | `scripts/data/get_crem_db.py` |

ChEMBL attribution: Gaulton et al., *Nucleic Acids Research* (2017, 2024 update),
"ChEMBL: a large-scale bioactivity database". If you publish results derived from
ChEMBL, cite ChEMBL and keep the CC BY-SA 4.0 attribution. Because CC BY-SA is a
share-alike license, **do not commit the ChEMBL dumps or databases derived from
them into this repository**; rebuild them locally with the scripts above.

## Methods

The local fragment-replacement layer reproduces the methodology of
[SwissBioisostere](https://www.swissbioisostere.ch):

> Wirth M., Zoete V., Michielin O., Sauer W.H.B. *SwissBioisostere: a database of
> molecular replacements for ligand design.* Nucleic Acids Research 41 (D1),
> D1137–D1143, 2013. doi:10.1093/nar/gkt951

IsoGrow is an independent reimplementation based on the same public idea
(matched molecular pairs mined from ChEMBL) and is not affiliated with or
endorsed by the SwissBioisostere authors.
