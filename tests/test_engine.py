"""Querying a real (if tiny) mmpdb database."""
from pathlib import Path

import pytest

from isogrow.bioisostere import engine

BENZOFURAN = "[*:1]c1cc2ccccc2o1"


def test_benzofuran_replacements(tiny_mmpdb: Path):
    df = engine.search(BENZOFURAN, str(tiny_mmpdb))
    found = set(df["Candidate Fragments"])
    # the classic ring replacements that SwissBioisostere also reports first
    assert "[*:1]c1cc2ccccc2s1" in found      # benzothiophene
    assert "[*:1]c1ccccc1" in found           # phenyl
    assert len(found) >= 6


def test_self_replacement_excluded(tiny_mmpdb: Path):
    df = engine.search(BENZOFURAN, str(tiny_mmpdb))
    assert BENZOFURAN not in set(df["Candidate Fragments"])


def test_frequency_is_sorted_descending(tiny_mmpdb: Path):
    df = engine.search(BENZOFURAN, str(tiny_mmpdb))
    freqs = list(df["Frequency"])
    assert freqs == sorted(freqs, reverse=True)


def test_topn_limits_rows(tiny_mmpdb: Path):
    df = engine.search(BENZOFURAN, str(tiny_mmpdb), topn=3)
    assert len(df) == 3


def test_delta_columns_are_numeric(tiny_mmpdb: Path):
    df = engine.search(BENZOFURAN, str(tiny_mmpdb))
    row = df[df["Candidate Fragments"] == "[*:1]c1cc2ccccc2s1"].iloc[0]
    assert float(row["Δ MW"]) == pytest.approx(16.07, abs=0.05)
    assert float(row["Δ tPSA"]) == pytest.approx(-13.14, abs=0.05)


def test_mw_filter(tiny_mmpdb: Path):
    full = engine.search(BENZOFURAN, str(tiny_mmpdb))
    narrow = engine.search(BENZOFURAN, str(tiny_mmpdb), mw_range=(-5, 5))
    assert len(narrow) < len(full)
    assert all(abs(float(v)) <= 5 for v in narrow["Δ MW"])


def test_unknown_fragment_returns_empty_frame(tiny_mmpdb: Path):
    df = engine.search("[*:1]C#N", str(tiny_mmpdb))
    assert df.empty
    assert list(df.columns)[:3] == ["Candidate Fragments", "Frequency", "# Better"]


def test_l2_descriptors(tiny_mmpdb: Path):
    df = engine.search(BENZOFURAN, str(tiny_mmpdb), add_l2=True)
    for col in ("ΔHAC", "ΔRings", "ΔHBD", "ΔHBA", "ΔRotB"):
        assert col in df.columns
