"""SwissBioisostere CSV layout: read, write and round-trip."""
from pathlib import Path

from isogrow.bioisostere import io as sbio


def test_query_line_and_columns(sb_export_csv: Path):
    query, df = sbio.read_sb_csv(str(sb_export_csv))
    assert query == "[*:1]c1cc2ccccc2o1"
    assert not df.empty
    for col in sbio.CANON_COLS:
        assert col in df.columns


def test_header_uses_non_breaking_space(sb_export_csv: Path):
    lines = sb_export_csv.read_text(encoding="utf-8").splitlines()
    assert lines[0].startswith("Query fragment:")
    assert "\xa0" in lines[1]          # official exports use NBSP in headers


def test_round_trip(sb_export_csv: Path, tmp_path: Path):
    query, df = sbio.read_sb_csv(str(sb_export_csv))
    out = tmp_path / "roundtrip.csv"
    sbio.write_sb_csv(df, str(out), query)

    query2, df2 = sbio.read_sb_csv(str(out))
    assert query2 == query
    assert len(df2) == len(df)
    assert list(df2.columns)[: len(sbio.CANON_COLS)] == sbio.CANON_COLS
    # the candidate order and the frequency values survive unchanged
    assert list(df2["Candidate Fragments"]) == list(df["Candidate Fragments"])


def test_write_requires_all_columns(tmp_path: Path):
    import pandas as pd
    import pytest

    with pytest.raises(ValueError):
        sbio.write_sb_csv(pd.DataFrame({"x": [1]}), str(tmp_path / "bad.csv"), "[*:1]C")
