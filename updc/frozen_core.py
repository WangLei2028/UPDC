"""Integrity metadata for the frozen v17 C numerical core."""
from __future__ import annotations

import hashlib
from pathlib import Path

FROZEN_CORE_SHA256 = {
    "updc_core.c": "df2cfffce5cf9b81749ae41dff8574c1fd8f2d4d1b409890a1e23062ad5807bb",
    "overlap.c": "c32b447dff787f599207ea35b676be6c1cfab9f76d8a4e95966bb9db38b68d71",
    "overlap_geometry.h": "f297ccf2c7eff07e5e8f73fc9a901c41d33a59a6b147ed264995ed47ef7d43f8",
    "psf_projector.c": "a6f6baba9cb2bba5c9761b1590918a4f45aa82f0791863c88a7cfa4ff3b5bfd1",
}


def core_directory() -> Path:
    return Path(__file__).resolve().parent / "c_core"


def verify_frozen_core() -> dict[str, bool]:
    root = core_directory()
    result: dict[str, bool] = {}
    for name, expected in FROZEN_CORE_SHA256.items():
        actual = hashlib.sha256((root / name).read_bytes()).hexdigest()
        result[name] = actual == expected
    return result


def assert_frozen_core() -> None:
    failures = [name for name, ok in verify_frozen_core().items() if not ok]
    if failures:
        raise RuntimeError("Frozen UPDC v17 C core integrity check failed: " + ", ".join(failures))
