from __future__ import annotations

from copy import deepcopy
from typing import Any, Mapping

from .base import InstrumentAdapter, ConfigDict


class JWSTAdapter(InstrumentAdapter):
    name = "jwst"
    description = "JWST Level-2 adapter using authoritative meta.wcs/GWCS"

    def apply_defaults(self, config: ConfigDict) -> ConfigDict:
        cfg = deepcopy(config)
        inp = cfg.setdefault("input", {})
        explicit = cfg.get("_explicit", set())
        if ("input", "input_glob") not in explicit: inp["input_glob"] = "*_cal*.fits"
        if ("input", "pairing") not in explicit: inp["pairing"] = "index"
        if ("input", "wcs_mode") not in explicit: inp["wcs_mode"] = "gwcs"
        return cfg

    def validate(self, config: Mapping[str, Any]) -> list[str]:
        errors: list[str] = []
        inp = config.get("input", {})
        if inp.get("wcs_mode") not in ("gwcs", "auto"):
            errors.append("JWST adapter expects input.wcs_mode='gwcs' or 'auto'")
        return errors


class GenericFITSAdapter(InstrumentAdapter):
    name = "fits"
    description = "Generic normalized FITS adapter using linear celestial FITS WCS"

    def apply_defaults(self, config: ConfigDict) -> ConfigDict:
        cfg = deepcopy(config)
        inp = cfg.setdefault("input", {})
        psf = cfg.setdefault("psf", {})
        explicit = cfg.get("_explicit", set())
        if ("input", "input_glob") not in explicit: inp["input_glob"] = "*.fits"
        if ("input", "pairing") not in explicit: inp["pairing"] = "sorted"
        if ("input", "wcs_mode") not in explicit: inp["wcs_mode"] = "fits"
        if ("psf", "psf_orientation") not in explicit: psf["psf_orientation"] = "as-is"
        if ("psf", "psf_crop_size") not in explicit: psf["psf_crop_size"] = 0
        return cfg

    def validate(self, config: Mapping[str, Any]) -> list[str]:
        errors: list[str] = []
        inp = config.get("input", {})
        if inp.get("wcs_mode") != "fits":
            errors.append("FITS adapter requires input.wcs_mode='fits'")
        return errors
