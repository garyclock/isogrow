"""IsoGrow — local SwissBioisostere-style fragment replacement + growth pipeline.

Two layers, importable independently:

* :mod:`isogrow.bioisostere` — matched-molecular-pair search (mmpdb + ChEMBL),
  reproducing the SwissBioisostere output layout plus the activity columns.
* :mod:`isogrow.pipeline` — stages A→G: bioisostere input → CReM variants →
  3D → alignment → PDBQT → docking → ranking.
"""

__version__ = "0.1.0"

__all__ = ["__version__"]
