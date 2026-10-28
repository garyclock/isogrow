"""Stages A→F on the shipped fixtures (no external datasets, no CReM database)."""
from pathlib import Path

import pytest
from rdkit import Chem

from isogrow.pipeline import align3d, dock, embed3d, inputs, pdb2pdbqt

LEADS = ["Cc1cc2ccccc2o1", "Cc1cc2ccccc2s1", "Cc1cc2ccccc2[nH]1"]


def test_input_layer_from_generated_export(sb_export_csv: Path):
    query, df = inputs.get_bioisosteres(
        source="official_csv", csv_path=str(sb_export_csv),
        logp_range=(-10, 10), tpsa_range=(-300, 300), mw_range=(-200, 80),
    )
    assert query == "[*:1]c1cc2ccccc2o1"
    assert not df.empty


def test_input_layer_requires_paths():
    with pytest.raises(ValueError):
        inputs.get_bioisosteres(source="official_csv", csv_path=None)
    with pytest.raises(ValueError):
        inputs.get_bioisosteres(source="local", query=None, mmpdb_path=None)


def test_embed_align_pdbqt(tmp_path: Path, reference_ligand: Path):
    sdfs = embed3d.batch_embed(LEADS, str(tmp_path / "3d"))
    assert len(sdfs) == len(LEADS)

    aligned = align3d.batch_align(sdfs, str(reference_ligand), str(tmp_path / "aligned"),
                              backend="mcs")   # dependency-free backend for the CI path
    assert len(aligned) >= 1
    # alignment must not change the molecule: heavy-atom count is preserved
    for original, path in zip(sdfs, aligned, strict=False):
        n_before = Chem.MolFromMolFile(original, removeHs=False).GetNumHeavyAtoms()
        n_after = Chem.MolFromMolFile(path, removeHs=False).GetNumHeavyAtoms()
        assert n_before == n_after

    pdbqts = pdb2pdbqt.batch_pdbqt(aligned, str(tmp_path / "pdbqt"))
    assert len(pdbqts) == len(aligned)
    for path in pdbqts:
        text = Path(path).read_text(encoding="utf-8", errors="replace")
        assert "ATOM" in text or "HETATM" in text


def test_box_from_reference_ligand(reference_ligand: Path):
    center, size = align3d.box_from_ligand(str(reference_ligand), padding=8.0)
    assert len(center) == 3 and len(size) == 3
    assert all(s >= 8.0 for s in size)


@pytest.mark.slow
def test_docking_scores_one_ligand(tmp_path: Path, receptor: Path, reference_ligand: Path, vina: str):
    sdf = embed3d.batch_embed(LEADS[:1], str(tmp_path / "3d"))[0]
    aligned = align3d.batch_align([sdf], str(reference_ligand), str(tmp_path / "aligned"),
                              backend="mcs")[0]
    pdbqt = pdb2pdbqt.batch_pdbqt([aligned], str(tmp_path / "pdbqt"))[0]

    center, size = align3d.box_from_ligand(str(reference_ligand), padding=8.0)
    res = dock.vina_dock(str(receptor), pdbqt, center, size,
                         out_pdbqt=str(tmp_path / "dock.pdbqt"),
                         exhaustiveness=1, cpu=1, vina_exe=vina)
    assert res["affinity"] is not None, res.get("error")
    assert res["affinity"] < 0  # kcal/mol, favourable
