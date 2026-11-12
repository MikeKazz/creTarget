"""Tests for config handling.

Two of the bugs fixed in 0.2.0 lived here, so both are pinned with
regression tests:
  - configured paths never had ~ expanded, while the docs told users to
    write "~/..." (a literal tilde then reached bedtools)
  - the config writer was a line-level string replace that could not add a
    missing key and mangled values containing quotes or backslashes
"""

import importlib

import pytest


@pytest.fixture
def cfg(tmp_path, monkeypatch):
    """A cretarget.config bound to a throwaway config directory."""
    monkeypatch.setenv("CRETARGET_CONFIG_DIR", str(tmp_path / "config"))
    monkeypatch.delenv("CRETARGET_CACHE_DIR", raising=False)
    import cretarget.config as config
    importlib.reload(config)
    return config


class TestExpand:
    def test_expands_tilde(self, cfg, monkeypatch):
        monkeypatch.setenv("HOME", "/home/someone")
        assert cfg.expand("~/ref/hg38.fa") == "/home/someone/ref/hg38.fa"

    def test_expands_environment_variables(self, cfg, monkeypatch):
        monkeypatch.setenv("SCRATCH", "/scratch/proj")
        assert cfg.expand("$SCRATCH/hg38.fa") == "/scratch/proj/hg38.fa"

    @pytest.mark.parametrize("value", ["", None])
    def test_passes_empty_through(self, cfg, value):
        assert cfg.expand(value) == value

    def test_leaves_plain_paths_alone(self, cfg):
        assert cfg.expand("/data/hg38.fa") == "/data/hg38.fa"


class TestGetAndRequirePath:
    def test_get_path_expands(self, cfg, monkeypatch):
        monkeypatch.setenv("HOME", "/home/someone")
        config = {"reference": {"fasta": "~/hg38.fa"}}
        assert cfg.get_path(config, "reference", "fasta") == "/home/someone/hg38.fa"

    def test_get_path_returns_default_when_unset(self, cfg):
        assert cfg.get_path({}, "reference", "fasta", "/fallback") == "/fallback"

    def test_require_path_raises_with_guidance(self, cfg):
        with pytest.raises(SystemExit) as excinfo:
            cfg.require_path({"reference": {"fasta": ""}}, "reference", "fasta")
        assert "cretarget configure --set reference.fasta" in str(excinfo.value)


class TestSetValues:
    def test_round_trips_a_simple_value(self, cfg):
        cfg.set_values(["reference.fasta=/data/hg38.fa"])
        assert cfg.load_config()["reference"]["fasta"] == "/data/hg38.fa"

    def test_round_trips_quotes_and_backslashes(self, cfg):
        tricky = '/opt/my "weird" \\path'
        cfg.set_values([f"seed_finder.dir={tricky}"])
        assert cfg.load_config()["seed_finder"]["dir"] == tricky

    def test_preserves_literal_tilde_on_disk(self, cfg):
        # Expansion happens at read time, not write time: the config should
        # stay portable between machines with different home directories.
        cfg.set_values(["reference.fasta=~/hg38.fa"])
        assert cfg.load_config()["reference"]["fasta"] == "~/hg38.fa"

    def test_can_add_a_key_that_was_absent(self, cfg):
        cfg.set_values(["brand_new.key=value"])
        assert cfg.load_config()["brand_new"]["key"] == "value"

    def test_applies_multiple_pairs(self, cfg):
        cfg.set_values(["reference.fasta=/a.fa", "reference.flashfry_jar=/b.jar"])
        config = cfg.load_config()
        assert config["reference"]["fasta"] == "/a.fa"
        assert config["reference"]["flashfry_jar"] == "/b.jar"

    def test_does_not_clobber_other_values(self, cfg):
        cfg.set_values(["reference.fasta=/a.fa"])
        cfg.set_values(["seed_finder.dir=/sf"])
        config = cfg.load_config()
        assert config["reference"]["fasta"] == "/a.fa"
        assert config["seed_finder"]["dir"] == "/sf"

    @pytest.mark.parametrize("bad", ["noequals", "missingsection=value"])
    def test_rejects_malformed_pairs(self, cfg, bad):
        with pytest.raises(SystemExit):
            cfg.set_values([bad])

    def test_rejects_atomically(self, cfg):
        """A bad pair must not leave a half-applied config behind."""
        cfg.set_values(["reference.fasta=/original.fa"])
        with pytest.raises(SystemExit):
            cfg.set_values(["reference.fasta=/new.fa", "malformed"])
        assert cfg.load_config()["reference"]["fasta"] == "/original.fa"


class TestCacheDir:
    """Resolution order: $CRETARGET_CACHE_DIR -> [cache].dir -> ~/.cretarget"""

    def test_environment_overrides_the_config_file(self, cfg, monkeypatch):
        # Regression test. The order used to be the other way round, so
        # `CRETARGET_CACHE_DIR=/tmp/x cretarget setup-reference` silently
        # ignored the override and downloaded 3.2 GB to the configured
        # location instead.
        monkeypatch.setenv("CRETARGET_CACHE_DIR", "/from/env")
        resolved = cfg.cache_dir({"cache": {"dir": "/from/config"}})
        assert str(resolved) == "/from/env"

    def test_falls_back_to_the_config_file(self, cfg, monkeypatch):
        monkeypatch.delenv("CRETARGET_CACHE_DIR", raising=False)
        resolved = cfg.cache_dir({"cache": {"dir": "/from/config"}})
        assert str(resolved) == "/from/config"

    def test_uses_environment_when_no_config_value(self, cfg, monkeypatch):
        monkeypatch.setenv("CRETARGET_CACHE_DIR", "/from/env")
        assert str(cfg.cache_dir({})) == "/from/env"

    def test_environment_value_is_expanded(self, cfg, monkeypatch):
        monkeypatch.setenv("HOME", "/home/someone")
        monkeypatch.setenv("CRETARGET_CACHE_DIR", "~/scratch/cache")
        assert str(cfg.cache_dir({})) == "/home/someone/scratch/cache"

    def test_defaults_to_home(self, cfg, monkeypatch):
        monkeypatch.delenv("CRETARGET_CACHE_DIR", raising=False)
        monkeypatch.setenv("HOME", "/home/someone")
        assert str(cfg.cache_dir({})) == "/home/someone/.cretarget"
