"""Δ properties of a replacement, checked against known reference values.

The reference numbers come from the official SwissBioisostere export used during
development (benzofuran → benzothiophene: ΔtPSA = -13.14, ΔMW = 16.07).  tPSA and
MW are deterministic RDKit calculations, so they must match to rounding; logP is
model dependent (Crippen here vs. ChEMBL ACD upstream) and is only checked for
direction.
"""
import pytest

from isogrow.bioisostere.props import delta_props

BENZOFURAN = "[*:1]c1cc2ccccc2o1"
BENZOTHIOPHENE = "[*:1]c1cc2ccccc2s1"


def test_benzothiophene_delta_matches_reference():
    d = delta_props(BENZOFURAN, BENZOTHIOPHENE)
    assert d is not None
    assert d["dtpsa"] == pytest.approx(-13.14, abs=0.02)
    assert d["dmw"] == pytest.approx(16.07, abs=0.02)


def test_logp_direction():
    d = delta_props(BENZOFURAN, BENZOTHIOPHENE)
    assert d["dlogp"] > 0  # sulfur is more lipophilic than oxygen


def test_self_replacement_is_zero():
    d = delta_props(BENZOFURAN, BENZOFURAN)
    assert d["dlogp"] == pytest.approx(0.0, abs=1e-6)
    assert d["dmw"] == pytest.approx(0.0, abs=1e-6)


def test_extra_descriptors_present():
    d = delta_props(BENZOFURAN, "[*:1]c1ccccc1")
    for key in ("dhac", "drings", "dhbd", "dhba", "drotb"):
        assert key in d
    # benzofuran (2 rings) -> phenyl (1 ring)
    assert d["drings"] == -1
    # same ring system, so the ring count does not change
    assert delta_props(BENZOFURAN, BENZOTHIOPHENE)["drings"] == 0


def test_invalid_input():
    assert delta_props(BENZOFURAN, "not-a-smiles") is None
