from pathlib import Path
import json

import pytest

pytest.importorskip("astropy")
from updc import _engine


def test_option_value_is_not_legacy_config():
    assert _engine._legacy_config_candidate(["--input-dir", "/tmp/level2"]) is None


def test_first_positional_remains_legacy_config():
    assert _engine._legacy_config_candidate(["legacy.toml", "--iterations", "3"]) == Path("legacy.toml")


def test_nonfinite_diagnostics_serialize_as_json_null():
    payload = _engine._json_finite({"psf_center": float("nan"), "minimum_mu": float("inf")})
    encoded = json.dumps(payload, allow_nan=False)
    assert json.loads(encoded) == {"psf_center": None, "minimum_mu": None}
