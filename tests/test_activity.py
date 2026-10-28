"""Activity layer: classification, ΔpChEMBL, success score, specific pair, similarity."""
from pathlib import Path

import pytest

from isogrow.bioisostere.activity import add_similarity, query_activity, query_pair

BENZOFURAN = "[*:1]c1cc2ccccc2o1"
BENZOTHIOPHENE = "[*:1]c1cc2ccccc2s1"


def test_classification_and_delta(tiny_activity_db: Path):
    df = query_activity(str(tiny_activity_db), BENZOFURAN)
    assert not df.empty
    row = df[df["Candidate Fragments"] == BENZOTHIOPHENE].iloc[0]
    # the fixture ranks this replacement better in at least one matched pair
    assert row["# Better"] >= 1
    assert row["# Worse"] == 0
    # ΔpChEMBL is averaged over every matched pair of the replacement, so it sits
    # between the improved pair (+1) and the neutral one (0)
    assert 0 < float(row["Δ pchembl"]) <= 1.0
    assert 0 < row["success_score"] < 1.0


def test_results_sorted_by_success(tiny_activity_db: Path):
    df = query_activity(str(tiny_activity_db), BENZOFURAN)
    scores = list(df["success_score"])
    assert scores == sorted(scores, reverse=True)


def test_threshold_controls_classification(tiny_activity_db: Path):
    # with a 1.5 log threshold the +1.0 replacement is only "similar"
    df = query_activity(str(tiny_activity_db), BENZOFURAN, threshold=1.5)
    row = df[df["Candidate Fragments"] == BENZOTHIOPHENE].iloc[0]
    assert row["# Better"] == 0
    assert row["# Similar"] >= 1


def test_min_pairs_filter(tiny_activity_db: Path):
    df = query_activity(str(tiny_activity_db), BENZOFURAN, min_pairs=99)
    assert df.empty


def test_query_pair(tiny_activity_db: Path):
    res = query_pair(str(tiny_activity_db), BENZOFURAN, BENZOTHIOPHENE)
    assert res is not None
    assert res["Frequency"] >= 1
    assert res["# Better"] >= 1
    assert float(res["Δ pchembl"]) > 0


def test_query_pair_missing(tiny_activity_db: Path):
    assert query_pair(str(tiny_activity_db), BENZOFURAN, "[*:1]C#N") is None


def test_requires_pchembl_property(tiny_mmpdb: Path):
    with pytest.raises(ValueError, match="pchembl"):
        query_activity(str(tiny_mmpdb), BENZOFURAN)


def test_similarity(tiny_activity_db: Path):
    df = query_activity(str(tiny_activity_db), BENZOFURAN)
    out = add_similarity(df, BENZOFURAN)
    assert "similarity" in out.columns
    assert all(0.0 <= v <= 1.0 for v in out["similarity"])
