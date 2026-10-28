"""Config loading: variable expansion, path resolution, validation."""
from pathlib import Path

import pytest

from isogrow import config as cfgmod
from isogrow import paths

EXAMPLE = paths.repo_root() / "config" / "example_local.json"


def test_expand_builtin_vars():
    out = cfgmod.expand_vars({"db": "${ISOGROW_DATA}/mmpdb/x.mmpdb"})
    assert Path(out["db"]).parts[-3:] == ("data", "mmpdb", "x.mmpdb")


def test_expand_env_var(monkeypatch):
    monkeypatch.setenv("MY_TEST_DIR", "/tmp/abc")
    assert cfgmod.expand_vars("${MY_TEST_DIR}/f.csv") == "/tmp/abc/f.csv"


def test_undefined_var_raises():
    with pytest.raises(cfgmod.ConfigError, match="undefined variable"):
        cfgmod.expand_vars("${NOT_DEFINED_ANYWHERE}/x")


def test_expand_is_recursive():
    out = cfgmod.expand_vars({"a": ["${ISOGROW_HOME}/x"], "b": {"c": "${ISOGROW_HOME}"}})
    assert out["a"][0] == f"{paths.repo_root()}/x"
    assert out["b"]["c"] == str(paths.repo_root())


def test_resolve_paths_makes_relative_absolute():
    out = cfgmod.resolve_paths({"out": "runs/demo", "receptor": "a.pdbqt"}, base="/tmp/base")
    assert out["out"] == str(Path("/tmp/base/runs/demo").resolve())
    assert Path(out["receptor"]).is_absolute()


def test_validate_requires_keys():
    with pytest.raises(cfgmod.ConfigError, match="missing required"):
        cfgmod.validate({"out": "x"})


def test_validate_source_specific():
    base = {"out": "x", "crem_db": "c", "ref_ligand": "r", "receptor": "p"}
    with pytest.raises(cfgmod.ConfigError, match="csv_path"):
        cfgmod.validate({**base, "source": "official_csv"})
    with pytest.raises(cfgmod.ConfigError, match="query"):
        cfgmod.validate({**base, "source": "local"})
    with pytest.raises(cfgmod.ConfigError, match="unknown source"):
        cfgmod.validate({**base, "source": "nope"})


def test_load_example_config():
    cfg = cfgmod.load_config(EXAMPLE)
    assert cfg["source"] == "local"
    assert Path(cfg["mmpdb_path"]).is_absolute()
    assert isinstance(cfg["logp_range"], list)


def test_overrides_win(tmp_path: Path):
    cfg = cfgmod.load_config(EXAMPLE, overrides={"out": str(tmp_path / "run")})
    assert cfg["out"] == str((tmp_path / "run").resolve())


def test_missing_file():
    with pytest.raises(cfgmod.ConfigError, match="not found"):
        cfgmod.load_config("does/not/exist.json")
