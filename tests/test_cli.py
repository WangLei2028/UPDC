from pathlib import Path

from updc.cli import init_config, validate


def test_init_and_validate(tmp_path: Path):
    cfg = tmp_path / "updc.toml"
    init_config(cfg)
    text = cfg.read_text()
    assert "schema_version = 1" in text
    # The template uses relative example directories. Create them so filesystem
    # role validation is exercised as part of the public validate workflow.
    (tmp_path / "data" / "level2").mkdir(parents=True)
    (tmp_path / "data" / "psf").mkdir(parents=True)
    result = validate(cfg, quiet=True)
    assert result["instrument"]["adapter"] == "jwst"
