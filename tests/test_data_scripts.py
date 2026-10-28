"""The one-command database builder, exercised offline on the bundled fixtures."""
import os
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts" / "data" / "build_database.py"
TINY = Path(__file__).parent / "data" / "tiny.smi"


@pytest.mark.slow
def test_build_database_from_smiles(tmp_path: Path, mmpdb_exe: str):
    """--smi path: no download, builds a real database and verifies it."""
    out = tmp_path / "tiny"
    env = {**os.environ, "ISOGROW_MMPDB": mmpdb_exe}
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), "--smi", str(TINY), "--out", str(out), "--jobs", "2"],
        capture_output=True, text=True, env=env, timeout=600, cwd=str(REPO),
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert (tmp_path / "tiny.mmpdb").is_file()
    # it must report both the build and a verification query
    assert "verification query" in proc.stdout
    assert "isogrow search --db" in proc.stdout
    # ... and the query must actually return rows for this fixture
    assert "Candidate Fragments" in proc.stdout


@pytest.mark.slow
def test_build_database_skip_download_without_file_fails(tmp_path: Path, mmpdb_exe: str):
    env = {**os.environ, "ISOGROW_MMPDB": mmpdb_exe, "ISOGROW_DATA": str(tmp_path / "data")}
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), "--skip-download", "--out", str(tmp_path / "x")],
        capture_output=True, text=True, env=env, timeout=120, cwd=str(REPO),
    )
    assert proc.returncode == 2
    assert "--skip-download" in proc.stdout
