from pathlib import Path

import pytest

from updc.adapters import get_adapter
from updc.config import engine_argv, load_config


def _minimal(tmp_path: Path, adapter="jwst") -> Path:
    cfg = tmp_path / "updc.toml"
    cfg.write_text(f'''\n[updc]\nschema_version = 1\n\n[instrument]\nadapter = "{adapter}"\n\n[input]\ninput_dir = "data"\n\n[psf]\npsf_dir = "psf"\n\n[output]\noutput_dir = "out"\n''')
    return cfg


def test_schema_v1_resolves_relative_paths(tmp_path: Path):
    cfg = load_config(_minimal(tmp_path))
    assert cfg["input"]["input_dir"] == (tmp_path / "data").resolve()
    assert cfg["psf"]["psf_dir"] == (tmp_path / "psf").resolve()
    assert cfg["output"]["output_dir"] == (tmp_path / "out").resolve()
    assert cfg["updc"]["schema_version"] == 1


def test_cli_override_translation(tmp_path: Path):
    cfg = load_config(_minimal(tmp_path), ["solver.iterations=31", "performance.fft_workers=4"])
    cfg = get_adapter("jwst").apply_defaults(cfg)
    argv = engine_argv(cfg)
    assert argv[argv.index("--iterations") + 1] == "31"
    assert argv[argv.index("--fft-workers") + 1] == "4"


def test_cli_override_survives_adapter_defaults(tmp_path: Path):
    cfg = load_config(_minimal(tmp_path, "fits"), [
        'input.input_glob="*_sci.fits"',
        'input.pairing="index"',
        'psf.psf_orientation="stored-transposed"',
        "psf.psf_crop_size=9",
    ])
    cfg = get_adapter("fits").apply_defaults(cfg)
    assert cfg["input"]["input_glob"] == "*_sci.fits"
    assert cfg["input"]["pairing"] == "index"
    assert cfg["psf"]["psf_orientation"] == "stored-transposed"
    assert cfg["psf"]["psf_crop_size"] == 9
    jwst = load_config(_minimal(tmp_path), ['input.wcs_mode="auto"'])
    assert get_adapter("jwst").apply_defaults(jwst)["input"]["wcs_mode"] == "auto"


def test_fits_adapter_defaults_are_not_jwst_specific(tmp_path: Path):
    cfg = load_config(_minimal(tmp_path, "fits"))
    cfg = get_adapter("fits").apply_defaults(cfg)
    assert cfg["input"]["wcs_mode"] == "fits"
    assert cfg["input"]["pairing"] == "sorted"
    assert cfg["input"]["input_glob"] == "*.fits"
    assert cfg["psf"]["psf_orientation"] == "as-is"
    assert cfg["psf"]["psf_crop_size"] == 0


def test_unknown_key_and_schema_version_rejected(tmp_path: Path):
    bad = tmp_path / "bad.toml"
    bad.write_text('''[updc]\nschema_version=1\n[instrument]\nadapter="jwst"\n[input]\ninput_dir="data"\nmade_up=1\n[psf]\npsf_dir="psf"\n''')
    with pytest.raises(ValueError, match="made_up"):
        load_config(bad)
    bad.write_text('''[updc]\nschema_version=2\n[instrument]\nadapter="jwst"\n[input]\ninput_dir="data"\n[psf]\npsf_dir="psf"\n''')
    with pytest.raises(ValueError, match="Unsupported"):
        load_config(bad)


def test_validate_paths_rejects_directory_used_as_manifest(tmp_path: Path):
    from updc.config import validate_paths
    cfg_path = _minimal(tmp_path)
    (tmp_path / "data").mkdir()
    (tmp_path / "psf").mkdir()
    cfg = load_config(cfg_path)
    validate_paths(cfg)
    cfg["input"]["manifest"] = tmp_path / "data"
    with pytest.raises(ValueError, match="manifest must be a file"):
        validate_paths(cfg)


def test_cache_lifecycle_translation(tmp_path: Path):
    cfg_path = _minimal(tmp_path)
    text = cfg_path.read_text() + '\n[cache]\nenabled=true\nreset="current"\nprune_after_days=30\n'
    cfg_path.write_text(text)
    cfg = load_config(cfg_path)
    argv = engine_argv(cfg)
    assert argv[argv.index("--cache-reset") + 1] == "current"
    assert argv[argv.index("--cache-prune-after-days") + 1] == "30"


@pytest.mark.parametrize("override, error", [
    ("updc.schema_version=true", "schema_version"),
    ("solver.iterations=true", "solver.iterations"),
    ("solver.save_every=false", "solver.save_every"),
    ("psf.psf_crop_size=false", "psf.psf_crop_size"),
    ("psf.psf_oversample=0", "psf.psf_oversample"),
    ("input.bad_dq_mask=-1", "input.bad_dq_mask"),
    ("input.poisson_pedestal_sigma=nan", "input.poisson_pedestal_sigma"),
    ("performance.fft_autotune=1", "performance.fft_autotune"),
    ("cache.prune_after_days=nan", "cache.prune_after_days"),
])
def test_invalid_scalar_rejected_by_validate(tmp_path: Path, override: str, error: str):
    with pytest.raises(ValueError, match=error):
        load_config(_minimal(tmp_path), [override])
