"""Stable TOML configuration schema v1 for UPDC."""
from __future__ import annotations

from copy import deepcopy
import math
from pathlib import Path
from typing import Any

from .version import CONFIG_SCHEMA_VERSION

try:
    import tomllib  # py311+
except ModuleNotFoundError:  # pragma: no cover
    import tomli as tomllib  # type: ignore


SECTION_KEYS: dict[str, set[str]] = {
    "updc": {"schema_version"},
    "instrument": {"adapter", "options"},
    "input": {"input_dir", "input_glob", "pairing", "manifest", "jwst_dir", "wcs_mode",
              "weight_mode", "bad_dq_mask", "poisson_pedestal_sigma"},
    "psf": {"psf_dir", "psf_glob", "psf_orientation", "psf_frame", "psf_oversample",
            "psf_crop_size", "psf_anchor_x", "psf_anchor_y", "write_transformed_psf"},
    "solver": {"iterations", "upsample", "save_every"},
    "performance": {"workers", "fft_workers", "exposure_workers", "preprocess_workers",
                    "fft_plan", "fft_autotune", "retune_fft"},
    "cache": {"enabled", "cache_dir", "fftw_wisdom", "reset", "prune_after_days"},
    "output": {"output_dir", "output_stem", "metrics_json", "keep_work", "work_dir"},
    "build": {"cc", "c_core_binary", "psf_projector_binary"},
}

PATH_KEYS = {
    ("input", "input_dir"), ("input", "manifest"), ("input", "jwst_dir"),
    ("psf", "psf_dir"), ("cache", "cache_dir"), ("cache", "fftw_wisdom"),
    ("output", "output_dir"), ("output", "metrics_json"), ("output", "work_dir"),
    ("build", "c_core_binary"), ("build", "psf_projector_binary"),
}

DEFAULTS: dict[str, dict[str, Any]] = {
    "updc": {"schema_version": CONFIG_SCHEMA_VERSION},
    "instrument": {"adapter": "jwst", "options": {}},
    "input": {
        "input_glob": "*_cal*.fits", "pairing": "index", "wcs_mode": "gwcs",
        "weight_mode": "rnoise-sky", "bad_dq_mask": 531, "poisson_pedestal_sigma": 0.0,
    },
    "psf": {
        "psf_glob": "*.fits", "psf_orientation": "stored-transposed", "psf_frame": "detector",
        "psf_oversample": 2, "psf_crop_size": 1431, "write_transformed_psf": False,
    },
    "solver": {"iterations": 19, "upsample": 2, "save_every": 5},
    "performance": {
        "workers": 1, "fft_workers": 2, "exposure_workers": 2, "preprocess_workers": 4,
        "fft_plan": "measure", "fft_autotune": True, "retune_fft": False,
    },
    "cache": {"enabled": True, "reset": "none", "prune_after_days": 0.0},
    "output": {"keep_work": False},
    "build": {"cc": "cc"},
}


def _merge(base: dict[str, Any], update: dict[str, Any]) -> dict[str, Any]:
    out = deepcopy(base)
    for key, value in update.items():
        if isinstance(value, dict) and isinstance(out.get(key), dict):
            out[key] = _merge(out[key], value)
        else:
            out[key] = value
    return out


def _parse_scalar(value: str) -> Any:
    # TOML itself is a convenient safe scalar parser.
    doc = tomllib.loads("x = " + value)
    return doc["x"]


def apply_overrides(config: dict[str, Any], overrides: list[str]) -> dict[str, Any]:
    out = deepcopy(config)
    for item in overrides:
        if "=" not in item or "." not in item.split("=", 1)[0]:
            raise ValueError(f"Invalid --set override {item!r}; expected section.key=value")
        lhs, rhs = item.split("=", 1)
        section, key = lhs.split(".", 1)
        if section not in SECTION_KEYS or key not in SECTION_KEYS[section]:
            raise ValueError(f"Unknown configuration key {section}.{key}")
        out.setdefault(section, {})[key] = _parse_scalar(rhs)
        # Adapter defaults must not replace a value explicitly supplied via --set.
        out.setdefault("_explicit", set()).add((section, key))
    return out


def load_config(path: Path, overrides: list[str] | None = None) -> dict[str, Any]:
    path = path.expanduser().resolve()
    with path.open("rb") as f:
        raw = tomllib.load(f)
    unknown_sections = set(raw) - set(SECTION_KEYS)
    if unknown_sections:
        raise ValueError("Unknown TOML section(s): " + ", ".join(sorted(unknown_sections)))
    for section, values in raw.items():
        if not isinstance(values, dict):
            raise ValueError(f"[{section}] must be a table")
        unknown = set(values) - SECTION_KEYS[section]
        if unknown:
            raise ValueError(f"Unknown key(s) in [{section}]: " + ", ".join(sorted(unknown)))
    cfg = _merge(DEFAULTS, raw)
    cfg["_explicit"] = {(section, key) for section, values in raw.items() for key in values}
    if overrides:
        cfg = apply_overrides(cfg, overrides)
    cfg["_config_path"] = path
    cfg["_config_dir"] = path.parent
    _resolve_paths(cfg, path.parent)
    validate_config(cfg)
    return cfg


def _resolve_paths(cfg: dict[str, Any], base: Path) -> None:
    for section, key in PATH_KEYS:
        value = cfg.get(section, {}).get(key)
        if value in (None, ""):
            continue
        p = Path(str(value)).expanduser()
        if not p.is_absolute():
            p = base / p
        cfg[section][key] = p.resolve()


def validate_config(cfg: dict[str, Any]) -> None:
    schema_version = cfg.get("updc", {}).get("schema_version")
    if type(schema_version) is not int or schema_version != CONFIG_SCHEMA_VERSION:
        raise ValueError(
            f"Unsupported updc.schema_version={schema_version!r}; supported version is {CONFIG_SCHEMA_VERSION}"
        )
    adapter = cfg.get("instrument", {}).get("adapter")
    if not isinstance(adapter, str) or not adapter.strip():
        raise ValueError("instrument.adapter must be a non-empty string")
    if not isinstance(cfg["instrument"].get("options", {}), dict):
        raise ValueError("instrument.options must be a table")
    inp = cfg["input"]; psf = cfg["psf"]; solver = cfg["solver"]; perf = cfg["performance"]
    if not any(inp.get(k) for k in ("input_dir", "manifest", "jwst_dir")):
        raise ValueError("Set one of input.input_dir, input.manifest, or input.jwst_dir")
    if not inp.get("manifest") and not (psf.get("psf_dir") or inp.get("jwst_dir")):
        raise ValueError("Set psf.psf_dir unless input.manifest or input.jwst_dir supplies PSFs")
    if inp.get("pairing") not in ("index", "sorted"):
        raise ValueError("input.pairing must be 'index' or 'sorted'")
    if inp.get("wcs_mode") not in ("auto", "gwcs", "fits"):
        raise ValueError("input.wcs_mode must be auto, gwcs, or fits")
    if inp.get("weight_mode") not in ("rnoise-sky", "totalvar", "mask"):
        raise ValueError("input.weight_mode must be rnoise-sky, totalvar, or mask")
    dq_mask = inp.get("bad_dq_mask")
    if type(dq_mask) is not int or not 0 <= dq_mask <= 0xFFFFFFFF:
        raise ValueError("input.bad_dq_mask must be a 32-bit non-negative integer")
    pedestal_sigma = inp.get("poisson_pedestal_sigma")
    if (type(pedestal_sigma) not in (int, float) or
            not math.isfinite(pedestal_sigma) or pedestal_sigma < 0):
        raise ValueError("input.poisson_pedestal_sigma must be a finite non-negative number")
    if psf.get("psf_orientation") not in ("stored-transposed", "as-is"):
        raise ValueError("psf.psf_orientation must be stored-transposed or as-is")
    if psf.get("psf_frame") not in ("detector", "target"):
        raise ValueError("psf.psf_frame must be detector or target")
    for section, key in (("solver", "iterations"), ("solver", "upsample"),
                         ("performance", "workers"), ("performance", "fft_workers"),
                         ("performance", "exposure_workers"), ("performance", "preprocess_workers")):
        value = cfg[section].get(key)
        if type(value) is not int or value < 1:
            raise ValueError(f"{section}.{key} must be a positive integer")
    if type(solver.get("save_every")) is not int or solver["save_every"] < 0:
        raise ValueError("solver.save_every must be >= 0")
    oversample = psf.get("psf_oversample")
    if type(oversample) is not int or oversample < 1:
        raise ValueError("psf.psf_oversample must be a positive integer")
    crop = psf.get("psf_crop_size")
    if type(crop) is not int or crop < 0 or (crop and crop % 2 == 0):
        raise ValueError("psf.psf_crop_size must be 0 or a positive odd integer")
    for key in ("psf_anchor_x", "psf_anchor_y"):
        value = psf.get(key)
        if value is not None and (type(value) not in (int, float) or not math.isfinite(value)):
            raise ValueError(f"psf.{key} must be a finite number")
    if (psf.get("psf_anchor_x") is None) != (psf.get("psf_anchor_y") is None):
        raise ValueError("psf.psf_anchor_x and psf.psf_anchor_y must be supplied together")
    for section, key in (("psf", "write_transformed_psf"),
                         ("performance", "fft_autotune"),
                         ("performance", "retune_fft"),
                         ("output", "keep_work")):
        if not isinstance(cfg[section].get(key), bool):
            raise ValueError(f"{section}.{key} must be true or false")
    if perf.get("fft_plan") not in ("estimate", "measure"):
        raise ValueError("performance.fft_plan must be estimate or measure")
    cache = cfg["cache"]
    if not isinstance(cache.get("enabled"), bool):
        raise ValueError("cache.enabled must be true or false")
    if cache.get("reset") not in ("none", "current", "all"):
        raise ValueError("cache.reset must be none, current, or all")
    prune_days = cache.get("prune_after_days")
    if (type(prune_days) not in (int, float) or
            not math.isfinite(prune_days) or prune_days < 0):
        raise ValueError("cache.prune_after_days must be a finite non-negative number")


def validate_paths(cfg: dict[str, Any]) -> None:
    """Validate filesystem roles without changing or creating data products."""
    inp = cfg["input"]
    psf = cfg["psf"]
    out = cfg["output"]
    cache = cfg["cache"]
    build = cfg["build"]

    def require_dir(section: str, key: str, value: Any | None) -> None:
        if value in (None, ""):
            return
        p = Path(value)
        if not p.exists():
            raise ValueError(f"{section}.{key} directory does not exist: {p}")
        if not p.is_dir():
            raise ValueError(f"{section}.{key} must be a directory: {p}")

    def require_file(section: str, key: str, value: Any | None) -> None:
        if value in (None, ""):
            return
        p = Path(value)
        if not p.exists():
            raise ValueError(f"{section}.{key} file does not exist: {p}")
        if not p.is_file():
            raise ValueError(f"{section}.{key} must be a file: {p}")

    require_dir("input", "input_dir", inp.get("input_dir"))
    require_dir("input", "jwst_dir", inp.get("jwst_dir"))
    require_file("input", "manifest", inp.get("manifest"))
    require_dir("psf", "psf_dir", psf.get("psf_dir"))

    # Output/cache/work directories may not exist yet, but if they do they must
    # be directories.  Their parents are created by the engine when needed.
    for section, key, value in (
        ("output", "output_dir", out.get("output_dir")),
        ("output", "work_dir", out.get("work_dir")),
        ("cache", "cache_dir", cache.get("cache_dir")),
    ):
        if value not in (None, ""):
            q = Path(value)
            if q.exists() and not q.is_dir():
                raise ValueError(f"{section}.{key} must be a directory: {q}")

    require_file("build", "c_core_binary", build.get("c_core_binary"))
    require_file("build", "psf_projector_binary", build.get("psf_projector_binary"))


def engine_argv(cfg: dict[str, Any]) -> list[str]:
    """Translate stable schema v1 into the frozen v17 engine CLI contract."""
    args: list[str] = []
    def add(flag: str, value: Any | None):
        if value is not None:
            args.extend([flag, str(value)])
    inp=cfg["input"]; psf=cfg["psf"]; sol=cfg["solver"]; perf=cfg["performance"]
    cache=cfg["cache"]; out=cfg["output"]; build=cfg["build"]
    add("--input-dir", inp.get("input_dir")); add("--input-glob", inp.get("input_glob"))
    add("--pairing", inp.get("pairing")); add("--manifest", inp.get("manifest")); add("--jwst-dir", inp.get("jwst_dir"))
    add("--wcs-mode", inp.get("wcs_mode")); add("--weight-mode", inp.get("weight_mode")); add("--bad-dq-mask", inp.get("bad_dq_mask"))
    add("--poisson-pedestal-sigma", inp.get("poisson_pedestal_sigma"))
    add("--psf-dir", psf.get("psf_dir")); add("--psf-glob", psf.get("psf_glob")); add("--psf-orientation", psf.get("psf_orientation"))
    add("--psf-frame", psf.get("psf_frame")); add("--psf-oversample", psf.get("psf_oversample")); add("--psf-crop-size", psf.get("psf_crop_size"))
    add("--psf-anchor-x", psf.get("psf_anchor_x")); add("--psf-anchor-y", psf.get("psf_anchor_y"))
    if psf.get("write_transformed_psf"): args.append("--write-transformed-psf")
    add("--iterations", sol.get("iterations")); add("--upsample", sol.get("upsample")); add("--save-every", sol.get("save_every"))
    add("--workers", perf.get("workers")); add("--fft-workers", perf.get("fft_workers")); add("--exposure-workers", perf.get("exposure_workers")); add("--preprocess-workers", perf.get("preprocess_workers")); add("--fft-plan", perf.get("fft_plan"))
    args.append("--fft-autotune" if perf.get("fft_autotune", True) else "--no-fft-autotune")
    args.append("--retune-fft" if perf.get("retune_fft", False) else "--no-retune-fft")
    add("--cache-dir", cache.get("cache_dir")); add("--fftw-wisdom", cache.get("fftw_wisdom"))
    add("--cache-reset", cache.get("reset", "none")); add("--cache-prune-after-days", cache.get("prune_after_days", 0.0))
    args.append("--preprocess-cache" if cache.get("enabled", True) else "--no-preprocess-cache")
    add("--output-dir", out.get("output_dir")); add("--output-stem", out.get("output_stem")); add("--metrics-json", out.get("metrics_json")); add("--work-dir", out.get("work_dir"))
    if out.get("keep_work"): args.append("--keep-work")
    add("--cc", build.get("cc")); add("--c-core-binary", build.get("c_core_binary")); add("--psf-projector-binary", build.get("psf_projector_binary"))
    return args
