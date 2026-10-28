"""Alignment backends: rigid MCS (always) and KCOMBU flexible fit (when installed)."""
from pathlib import Path

import pytest
from rdkit import Chem

from isogrow.pipeline import align3d, embed3d

LEAD = "Cc1cc2ccccc2o1"
EXPECTED_HEAVY = Chem.MolFromSmiles(LEAD).GetNumHeavyAtoms()


def test_unknown_backend_rejected(tmp_path: Path, reference_ligand: Path):
    sdf = embed3d.batch_embed([LEAD], str(tmp_path / "3d"))[0]
    with pytest.raises(ValueError, match="unknown alignment backend"):
        align3d.align_file(sdf, str(reference_ligand), str(tmp_path / "x.sdf"), backend="nope")


def test_mcs_backend_preserves_molecule(tmp_path: Path, reference_ligand: Path):
    sdf = embed3d.batch_embed([LEAD], str(tmp_path / "3d"))[0]
    out = tmp_path / "mcs.sdf"
    mol, rmsd = align3d.align_file(sdf, str(reference_ligand), str(out), backend="mcs")

    assert out.is_file()
    assert mol.GetNumHeavyAtoms() == EXPECTED_HEAVY
    assert rmsd is not None and 0.0 <= rmsd < 2.0
    # the written file is readable and keeps 3D coordinates
    reloaded = Chem.MolFromMolFile(str(out), removeHs=False)
    assert reloaded is not None
    assert reloaded.GetConformer().Is3D()


def test_mapped_rmsd_is_zero_for_identical_mols():
    mol = Chem.MolFromMolFile(str(Path(__file__).parent / "data" / "reference_ligand.sdf"),
                              removeHs=False)
    mapping = (tuple(range(mol.GetNumAtoms())), tuple(range(mol.GetNumAtoms())))
    assert align3d.mapped_rmsd(mol, mol, mapping) == pytest.approx(0.0, abs=1e-9)


@pytest.mark.slow
def test_kcombu_backend(tmp_path: Path, reference_ligand: Path, kcombu: str):
    sdf = embed3d.batch_embed([LEAD], str(tmp_path / "3d"))[0]
    out = tmp_path / "kcombu.sdf"

    mol, rmsd = align3d.align_file(
        sdf, str(reference_ligand), str(out),
        backend="kcombu", exe=kcombu, rigid=False, nini=1,
    )

    assert out.is_file()
    assert mol.GetNumHeavyAtoms() == EXPECTED_HEAVY
    assert mol.GetConformer().Is3D()
    assert rmsd is not None and rmsd < 2.0


@pytest.mark.slow
def test_kcombu_rigid_matches_mcs_closely(tmp_path: Path, reference_ligand: Path, kcombu: str):
    sdf = embed3d.batch_embed([LEAD], str(tmp_path / "3d"))[0]
    _m, rms_mcs = align3d.align_file(sdf, str(reference_ligand), str(tmp_path / "a.sdf"),
                                     backend="mcs")
    _k, rms_kc = align3d.align_file(sdf, str(reference_ligand), str(tmp_path / "b.sdf"),
                                    backend="kcombu", exe=kcombu, rigid=True, nini=1)
    assert rms_mcs is not None and rms_kc is not None
    # both are least-squares fits on the same MCS, so they must be comparable
    assert abs(rms_mcs - rms_kc) < 0.5
