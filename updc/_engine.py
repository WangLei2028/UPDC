#!/usr/bin/env python3
"""Run signed-x multiplicative UPDC on astronomy FITS exposures.

The bundled frontend supports JWST Level-2/GWCS workflows and a FITS-WCS fallback.

The SCI extension is *not* background-subtracted by this driver.  A robust scalar
background is estimated per exposure and enters the forward model explicitly:

    y_eff = SCI + pedestal
    raw_mu = H x + background + pedestal
    ratio = y_eff/raw_mu, permanently 1 after raw_mu <= 0 at that pixel
    x <- x * H^T[ W * ratio ] / H^T W

The non-negative pedestal adds only a constant mean to numerator and model;
it never adds a noise realization. It retains otherwise valid SCI<=0 pixels.
"""
from __future__ import annotations

import argparse
import atexit
import csv
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import json
import hashlib
import os
import resource
import re
import shutil
import platform
import subprocess
import sys
import tempfile
import time
import warnings
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from .cache_policy import prune_dataset_caches, reset_dataset_caches, touch_dataset_cache
from .build import (
    build_signature, compiler_identity, discover_fftw, native_arch_flag,
    probe_c11_atomics, sanitize_c_environment, signature_matches, write_signature,
)

PACKAGE = Path(__file__).resolve().parent


class Stage:
    __slots__ = ("name", "wall", "self_cpu", "child_cpu", "cpu")
    def __init__(self, name):
        self.name = name; self.wall = self.self_cpu = self.child_cpu = self.cpu = 0.0

STAGES: list[Stage] = []


class timer:
    def __init__(self, name: str, quiet: bool = False):
        self.stage = Stage(name); self.quiet = quiet
    def __enter__(self):
        self._t0 = time.perf_counter(); self._c0 = time.process_time()
        ru = resource.getrusage(resource.RUSAGE_CHILDREN)
        self._ch0 = ru.ru_utime + ru.ru_stime
        return self.stage
    def __exit__(self, exc_type, exc, tb):
        s = self.stage
        s.wall = time.perf_counter() - self._t0
        s.self_cpu = time.process_time() - self._c0
        ru = resource.getrusage(resource.RUSAGE_CHILDREN)
        s.child_cpu = (ru.ru_utime + ru.ru_stime) - self._ch0
        s.cpu = s.self_cpu + s.child_cpu
        STAGES.append(s)
        if not self.quiet:
            print(f"[TIME] {s.name:<28} CPU={s.cpu:9.3f}s  wall={s.wall:9.3f}s", flush=True)
        return False


def report_total(label: str = "TOTAL") -> None:
    cpu = sum(s.cpu for s in STAGES); wall = sum(s.wall for s in STAGES)
    print("-" * 72)
    for s in STAGES:
        print(f"{s.name:<28} CPU={s.cpu:9.3f}s   wall={s.wall:9.3f}s")
    print("-" * 72)
    print(f"{label:<28} CPU={cpu:9.3f}s   wall={wall:9.3f}s")


def _json_finite(value):
    """Keep optional non-finite diagnostics valid in standards-compliant JSON."""
    if isinstance(value, np.ndarray):
        return _json_finite(value.tolist())
    if isinstance(value, np.generic):
        return _json_finite(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    if isinstance(value, dict):
        return {key: _json_finite(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_finite(item) for item in value]
    return value





_CONFIG_PATH_KEYS = {
    "input_dir", "psf_dir", "manifest", "jwst_dir", "output_dir", "work_dir",
    "cache_dir", "metrics_json", "fftw_wisdom", "c_core_binary", "psf_projector_binary",
}
_CONFIG_SECTIONS = {
    "input": {"input_dir", "input_glob", "pairing", "manifest", "jwst_dir", "wcs_mode",
              "weight_mode", "bad_dq_mask", "poisson_pedestal_sigma"},
    "psf": {"psf_dir", "psf_glob", "psf_orientation", "psf_frame", "psf_oversample",
            "psf_crop_size", "psf_anchor_x", "psf_anchor_y", "write_transformed_psf"},
    "solver": {"iterations", "upsample", "save_every"},
    "performance": {"workers", "fft_workers", "exposure_workers", "preprocess_workers",
                    "fft_plan", "fft_autotune", "retune_fft"},
    "output": {"output_dir", "output_stem", "metrics_json", "keep_work", "work_dir"},
    "cache": {"enabled", "cache_dir", "fftw_wisdom", "reset", "prune_after_days"},
    "build": {"cc", "c_core_binary", "psf_projector_binary"},
}


def _load_toml_module():
    try:
        import tomllib  # type: ignore
        return tomllib
    except ModuleNotFoundError:  # Python 3.10
        import tomli  # type: ignore
        return tomli


def load_config_defaults(path: Path) -> dict:
    """Load a strict TOML config. Relative paths are relative to the config file."""
    toml = _load_toml_module()
    path = path.expanduser().resolve()
    with path.open("rb") as f:
        data = toml.load(f)
    unknown_sections = set(data) - set(_CONFIG_SECTIONS)
    if unknown_sections:
        raise ValueError(f"Unknown config section(s): {', '.join(sorted(unknown_sections))}")
    defaults: dict = {}
    for section, allowed in _CONFIG_SECTIONS.items():
        values = data.get(section, {})
        if not isinstance(values, dict):
            raise ValueError(f"Config section [{section}] must be a table")
        unknown = set(values) - allowed
        if unknown:
            raise ValueError(f"Unknown key(s) in [{section}]: {', '.join(sorted(unknown))}")
        for key, value in values.items():
            if section == "cache" and key == "enabled":
                if not isinstance(value, bool):
                    raise ValueError("cache.enabled must be true or false")
                defaults["no_preprocess_cache"] = not value
                continue
            if section == "cache" and key == "reset":
                defaults["cache_reset"] = value
                continue
            if section == "cache" and key == "prune_after_days":
                defaults["cache_prune_after_days"] = value
                continue
            if key in _CONFIG_PATH_KEYS and value is not None:
                q = Path(str(value)).expanduser()
                if not q.is_absolute():
                    q = path.parent / q
                value = q.resolve()
            defaults[key] = value
    defaults["config_file"] = path
    return defaults

def _thread_map_ordered(func, items, workers: int):
    """Run independent exposure tasks concurrently while preserving input order."""
    seq = list(items)
    if not seq:
        return []
    nworkers = max(1, min(int(workers), len(seq)))
    if nworkers == 1:
        return [func(item) for item in seq]
    with ThreadPoolExecutor(max_workers=nworkers, thread_name_prefix="updc-pre") as pool:
        return list(pool.map(func, seq))

def sci_header(path: Path) -> fits.Header:
    with fits.open(path, memmap=True) as hdul:
        if "SCI" not in hdul:
            raise ValueError(f"Missing SCI extension in {path}")
        return hdul["SCI"].header.copy()


def celestial_wcs(path: Path) -> WCS:
    """Linear FITS WCS, used for target-grid orientation/scale and as fallback."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        w = WCS(sci_header(path)).celestial
    if not w.has_celestial:
        raise ValueError(f"SCI extension has no usable celestial FITS WCS: {path}")
    return w


def pixel_to_world(path: Path, x: np.ndarray, y: np.ndarray, mode: str):
    """Use authoritative JWST GWCS when available; FITS WCS is the fallback."""
    if mode != "fits":
        try:
            try:
                from jwst import datamodels
            except ImportError:
                from stdatamodels.jwst import datamodels
            with datamodels.open(path) as model:
                if model.meta.wcs is None:
                    raise RuntimeError("meta.wcs is absent")
                world = model.meta.wcs(x, y)
                return np.asarray(world[0]), np.asarray(world[1]), "gwcs"
        except Exception as exc:
            if mode == "gwcs":
                raise RuntimeError(f"GWCS requested but unavailable for {path}: {exc}") from exc
            if not getattr(pixel_to_world, "_warned", False):
                print(f"WARNING: JWST GWCS unavailable ({exc}); falling back to SCI FITS WCS", flush=True)
                pixel_to_world._warned = True
    w = celestial_wcs(path)
    ra, dec = w.all_pix2world(x, y, 0)
    return np.asarray(ra), np.asarray(dec), "fits"


def make_target(images: list[Path], wcs_mode: str, upsample: int,
                input_shape: tuple[int, int], preprocess_workers: int = 1) -> tuple[WCS, fits.Header, tuple[int, int]]:
    first = celestial_wcs(images[0])
    h0 = sci_header(images[0])
    target = WCS(naxis=2)
    target.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    target.wcs.cunit = ["deg", "deg"]
    target.wcs.crval = first.wcs.crval.copy()
    target.wcs.cd = first.pixel_scale_matrix / upsample
    target.wcs.crpix = first.wcs.crpix * upsample - (upsample - 1) / 2

    input_height, input_width = input_shape
    edge_x = np.linspace(-0.5, input_width - 0.5, 65)
    edge_y = np.linspace(-0.5, input_height - 0.5, 65)
    boundary = np.concatenate([
        np.column_stack((edge_x, np.full_like(edge_x, -0.5))),
        np.column_stack((edge_x, np.full_like(edge_x, input_height - 0.5))),
        np.column_stack((np.full_like(edge_y, -0.5), edge_y)),
        np.column_stack((np.full_like(edge_y, input_width - 0.5), edge_y)),
    ])

    def map_boundary(item):
        path = item
        ra, dec, used = pixel_to_world(path, boundary[:, 0], boundary[:, 1], wcs_mode)
        return np.asarray(ra), np.asarray(dec), used

    mapped = _thread_map_ordered(map_boundary, images, preprocess_workers)
    xy = []
    used_modes = set()
    for ra, dec, used in mapped:
        used_modes.add(used)
        xy.append(target.wcs_world2pix(np.column_stack((ra, dec)), 0))
    xy = np.concatenate(xy)
    low = np.floor(np.nanmin(xy, axis=0) - 4).astype(int)
    high = np.ceil(np.nanmax(xy, axis=0) + 4).astype(int)
    target.wcs.crpix -= low
    width, height = (high - low + 1).tolist()
    header = target.to_header()
    header["BUNIT"] = h0.get("BUNIT", "MJy/sr")
    header["FILTER"] = h0.get("FILTER", fits.getheader(images[0], 0).get("FILTER", "UNKNOWN"))
    header["DETECTOR"] = h0.get("DETECTOR", fits.getheader(images[0], 0).get("DETECTOR", "UNKNOWN"))
    header["HISTORY"] = f"UPDC target grid from first Level-2 WCS; pixel scale divided by {upsample}."
    header["WCSMAP"] = "+".join(sorted(used_modes)).upper()
    return target, header, (width, height)

def prepare_maps(images: list[Path], target: WCS, directory: Path, wcs_mode: str,
                 input_shape: tuple[int, int], preprocess_workers: int = 1) -> None:
    input_height, input_width = input_shape
    yy, xx = np.indices((input_height + 1, input_width + 1), dtype=np.float64)
    xx -= 0.5; yy -= 0.5

    def prepare_one(item):
        index, path = item
        ra, dec, used = pixel_to_world(path, xx, yy, wcs_mode)
        local_target = WCS(target.to_header())
        px, py = local_target.wcs_world2pix(ra, dec, 0)
        np.stack((px, py), axis=-1).astype("<f8").tofile(directory / f"vertices_{index}.bin")
        return index, used

    results = _thread_map_ordered(prepare_one, list(enumerate(images, 1)), preprocess_workers)
    for index, used in results:
        print(f"WCS vertices {index}/{len(images)} prepared with {used.upper()}", flush=True)

def robust_background_mask(science: np.ndarray, base_valid: np.ndarray,
                           max_iter: int = 6) -> tuple[np.ndarray, float, float]:
    """Asymmetric sigma clipping to reject positive astronomical sources."""
    mask = base_valid.copy()
    for _ in range(max_iter):
        vals = science[mask]
        if vals.size < max(32, min(1000, base_valid.size // 100)):
            raise ValueError("Too few valid pixels for robust background estimation")
        med = float(np.median(vals))
        mad = float(np.median(np.abs(vals - med)))
        sig = 1.4826 * mad
        if not np.isfinite(sig) or sig <= 0:
            sig = float(np.std(vals))
        if not np.isfinite(sig) or sig <= 0:
            raise ValueError("Cannot estimate a positive background RMS")
        new = base_valid & (science > med - 4.5 * sig) & (science < med + 2.5 * sig)
        if np.array_equal(new, mask):
            mask = new; break
        mask = new
    vals = science[mask]
    med = float(np.median(vals))
    sig = float(1.4826 * np.median(np.abs(vals - med)))
    return mask, med, sig


def build_weight(mode: str, var_p: np.ndarray, var_r: np.ndarray, var_f: np.ndarray,
                 sky_mask: np.ndarray) -> tuple[np.ndarray, float]:
    sky_var_p = float(np.median(var_p[sky_mask & np.isfinite(var_p) & (var_p >= 0)]))
    if not np.isfinite(sky_var_p) or sky_var_p < 0:
        raise ValueError("Cannot estimate sky VAR_POISSON")
    if mode == "mask":
        w = np.ones_like(var_r, dtype=np.float64)
    elif mode == "rnoise-sky":
        # Recommended default: retain exposure/pixel instrumental-noise information
        # without using source-dependent VAR_POISSON at bright-object pixels.
        variance = var_r + sky_var_p
        with np.errstate(divide="ignore", invalid="ignore"):
            w = 1.0 / variance
    elif mode == "totalvar":
        variance = var_p + var_r + var_f
        with np.errstate(divide="ignore", invalid="ignore"):
            w = 1.0 / variance
    else:
        raise ValueError(mode)
    return w, sky_var_p


def prepare_weighted_inputs(images: list[Path], directory: Path, weight_mode: str,
                            bad_dq_mask: int, pedestal_sigma: float,
                            input_shape: tuple[int, int], preprocess_workers: int = 1):
    """Create raw float32 SCI/WHT working planes for the C numerical core.

    Exposure files are independent and are prepared concurrently.  Results are
    returned and printed in exposure order so runs remain easy to compare.
    """
    def prepare_one(item):
        index, path = item
        with fits.open(path, memmap=False) as h:
            for name in ("SCI", "DQ", "VAR_POISSON", "VAR_RNOISE", "VAR_FLAT"):
                if name not in h:
                    raise ValueError(f"Missing {name} in {path}")
            science = np.asarray(h["SCI"].data, dtype=np.float32)
            dq = np.asarray(h["DQ"].data, dtype=np.uint32)
            var_p = np.asarray(h["VAR_POISSON"].data, dtype=np.float64)
            var_r = np.asarray(h["VAR_RNOISE"].data, dtype=np.float64)
            var_f = np.asarray(h["VAR_FLAT"].data, dtype=np.float64)
            if science.shape != input_shape:
                raise ValueError(f"Unexpected SCI shape {science.shape} in {path}")
            if any(a.shape != science.shape for a in (dq, var_p, var_r, var_f)):
                raise ValueError(f"Variance/DQ shape mismatch in {path}")

            base_valid = np.isfinite(science) & ((dq & np.uint32(bad_dq_mask)) == 0)
            base_valid &= np.isfinite(var_p) & np.isfinite(var_r)
            if weight_mode == "totalvar":
                base_valid &= np.isfinite(var_f)
            sky_mask, background, sky_sigma = robust_background_mask(science, base_valid)

            raw_weight, sky_var_p = build_weight(weight_mode, var_p, var_r, var_f, sky_mask)
            usable = base_valid & np.isfinite(raw_weight) & (raw_weight > 0)
            if not np.any(usable):
                raise ValueError(f"No usable pixels remain in {path}")
            min_valid_science = float(np.min(science[usable].astype(np.float64)))
            positive_floor = max(1.0e-12, 1.0e-6 * float(sky_sigma))
            auto_pedestal = max(0.0, positive_floor - min_valid_science)
            requested_pedestal = max(0.0, pedestal_sigma * sky_sigma - background)
            pedestal = max(auto_pedestal, requested_pedestal)
            y_eff = science.astype(np.float64) + pedestal
            nonpositive_before = usable & (science <= 0)
            nonpositive_after = usable & (y_eff <= 0)
            if np.any(nonpositive_after):
                raise ValueError(f"Automatic pedestal failed to make y_eff positive in {path}")
            wht = np.where(usable, raw_weight, 0.0).astype(np.float32)
            if not np.any(wht > 0):
                raise ValueError(f"No usable pixels remain in {path}")

            science_raw = directory / f"level2_science_{index:03d}.f32"
            weight_raw = directory / f"level2_weight_{index:03d}.f32"
            np.asarray(science, dtype="<f4").tofile(science_raw)
            np.asarray(wht, dtype="<f4").tofile(weight_raw)

            metric = {
                "exposure": index, "source": str(path), "background": background,
                "background_sigma": sky_sigma, "poisson_pedestal": pedestal,
                "auto_positivity_pedestal": auto_pedestal,
                "requested_pedestal": requested_pedestal,
                "min_valid_science": min_valid_science,
                "sky_var_poisson": sky_var_p, "weight_mode": weight_mode,
                "base_invalid": int((~base_valid).sum()),
                "nonpositive_valid_before_pedestal": int(nonpositive_before.sum()),
                "nonpositive_after_pedestal": int(nonpositive_after.sum()),
                "wht_zero": int((wht == 0).sum()),
                "wht_valid_median": float(np.median(wht[wht > 0])),
            }
            log = (
                f"Input {index}/{len(images)}: B={background:.8g}, sigma={sky_sigma:.8g}, "
                f"pedestal={pedestal:.8g} (auto={auto_pedestal:.8g}), "
                f"SCI<=0 retained={metric['nonpositive_valid_before_pedestal']}, "
                f"zero-WHT={metric['wht_zero']}"
            )
            return (science_raw, weight_raw), metric, log

    results = _thread_map_ordered(prepare_one, list(enumerate(images, 1)), preprocess_workers)
    prepared = []
    metrics = []
    for pair, metric, log in results:
        prepared.append(pair)
        metrics.append(metric)
        print(log, flush=True)
    return prepared, metrics

def read_psf_image(path: Path) -> np.ndarray:
    """Read the first 2-D image from a PSF FITS file as float32."""
    with fits.open(path, memmap=False) as hdul:
        for hdu in hdul:
            data = getattr(hdu, "data", None)
            if data is not None and np.ndim(data) == 2:
                image = np.asarray(data, dtype=np.float32)
                if image.size == 0:
                    break
                return image
    raise ValueError(f"No readable 2-D PSF image in {path}")


def _canonical_psf(image: np.ndarray, orientation: str) -> np.ndarray:
    """Return PSF axes in logical detector x/y order before WCS transformation."""
    if orientation == "stored-transposed":
        return np.ascontiguousarray(image.T, dtype=np.float32)
    return np.ascontiguousarray(image, dtype=np.float32)


def compile_psf_projector(args, directory: Path) -> Path:
    """Build the standalone polygon PSF projector (no FFTW dependency)."""
    if args.psf_projector_binary:
        binary = args.psf_projector_binary.expanduser().resolve()
        if not binary.is_file():
            raise ValueError(f"PSF projector binary does not exist: {binary}")
        return binary
    c_core = PACKAGE / "c_core"
    binary = directory / ("updc_psf_projector.exe" if os.name == "nt" else "updc_psf_projector")
    sources = [c_core / "psf_projector.c", c_core / "overlap.c", c_core / "overlap_geometry.h"]
    source_digest = _source_digest(sources)
    build_env, removed = sanitize_c_environment()
    if removed:
        print("Build environment: ignored C++ include paths while compiling C sources", flush=True)
    cc = shutil.which(args.cc, path=build_env.get("PATH")) or args.cc
    compiler = compiler_identity(cc, build_env)
    native = native_arch_flag(cc, build_env)
    flags = ["-O3", "-std=c11", "-DNDEBUG", "-Wall", "-Wextra", "-pthread"]
    if native:
        flags.insert(0, native)
    signature = Path(str(binary) + ".build.json")
    expected = build_signature(source_digest, compiler, flags, libraries=["-lm"])
    if binary.is_file() and signature_matches(signature, expected):
        subprocess.run([str(binary), "--self-test"], check=True, env=build_env)
        return binary
    base = [cc, *flags, "-o", str(binary),
            str(c_core / "psf_projector.c"), str(c_core / "overlap.c"), "-lm"]
    subprocess.run(base, check=True, env=build_env)
    subprocess.run([str(binary), "--self-test"], check=True, env=build_env)
    write_signature(signature, expected)
    return binary


def _psf_detector_vertex_grid(shape: tuple[int, int], anchor: tuple[float, float],
                              oversample: int) -> tuple[np.ndarray, np.ndarray]:
    """Detector coordinates of all oversampled PSF pixel corners.

    An odd-sized HybPSF pixel at its array center is centered exactly on the
    detector anchor.  A PSF pixel is 1/oversample native detector pixel wide.
    """
    ph, pw = shape
    ax, ay = anchor
    x = ax + (np.arange(pw + 1, dtype=np.float64) - pw / 2.0) / float(oversample)
    y = ay + (np.arange(ph + 1, dtype=np.float64) - ph / 2.0) / float(oversample)
    xx, yy = np.meshgrid(x, y)
    return xx, yy


def _map_psf_vertices(path: Path, target: WCS, shape: tuple[int, int],
                      anchor: tuple[float, float], oversample: int,
                      input_shape: tuple[int, int], wcs_mode: str,
                      out_path: Path) -> dict:
    """Map a detector-center oversampled PSF grid through full GWCS to target.

    The target WCS is the distortion-free first-order TAN grid built by
    ``make_target``.  Dither translation is removed by subtracting the mapped
    detector anchor, while rotation/scale/shear and all higher-order detector
    distortion encoded in the JWST GWCS are retained across the PSF support.
    """
    ph, pw = shape
    xx, yy = _psf_detector_vertex_grid(shape, anchor, oversample)
    hdet, wdet = input_shape
    if (float(np.min(xx)) < -0.5 or float(np.max(xx)) > wdet - 0.5 or
            float(np.min(yy)) < -0.5 or float(np.max(yy)) > hdet - 0.5):
        raise ValueError(
            f"PSF support {pw}x{ph} at anchor {anchor} and oversample={oversample} "
            "extends outside the detector; choose a smaller PSF or a valid anchor"
        )
    ra, dec, used = pixel_to_world(path, xx, yy, wcs_mode)
    local_target = WCS(target.to_header())
    px, py = local_target.wcs_world2pix(ra, dec, 0)
    cra, cdec, used_center = pixel_to_world(
        path, np.asarray([anchor[0]], dtype=np.float64),
        np.asarray([anchor[1]], dtype=np.float64), wcs_mode
    )
    cx, cy = local_target.wcs_world2pix(cra, cdec, 0)
    center_x, center_y = float(cx[0]), float(cy[0])
    px = np.asarray(px, dtype=np.float64) - center_x
    py = np.asarray(py, dtype=np.float64) - center_y
    if not (np.isfinite(px).all() and np.isfinite(py).all()):
        raise ValueError(f"Non-finite full-GWCS PSF mapping for {path}")
    np.stack((px, py), axis=-1).astype("<f8").tofile(out_path)
    return {
        "anchor_x": float(anchor[0]), "anchor_y": float(anchor[1]),
        "center_target_x": center_x, "center_target_y": center_y,
        "extent_x": float(np.max(np.abs(px))),
        "extent_y": float(np.max(np.abs(py))),
        "wcs_used": used if used == used_center else f"{used}+{used_center}",
    }


def prepare_psf_raw(psfs: list[Path], images: list[Path], target: WCS,
                    anchors: list[tuple[float, float]], input_shape: tuple[int, int],
                    directory: Path, args, output_dir: Path):
    """Prepare UPDC-ready PSFs with exposure-parallel GWCS preparation.

    Detector-frame HybPSFs are treated as oversampled detector images and are
    pushed through each exposure's full detector->sky GWCS into the linear
    distortion-free target grid using exact polygon-overlap flux fractions.
    Independent FITS reads/GWCS vertex maps/crops are parallelized with
    --preprocess-workers; the C polygon projector uses the same worker count.
    """
    if len(psfs) != len(images) or len(anchors) != len(images):
        raise ValueError("PSF/exposure/anchor count mismatch")
    oversample = int(args.psf_oversample or args.upsample)

    if args.psf_frame == "target":
        def bypass_one(item):
            index, path, anchor = item
            image = read_psf_image(path)
            raw_h, raw_w = image.shape
            raw = directory / f"psf_{index:03d}.f32"
            np.asarray(image, dtype="<f4").tofile(raw)
            total = float(np.sum(image, dtype=np.float64))
            record = {
                "exposure": index, "anchor_x": float(anchor[0]), "anchor_y": float(anchor[1]),
                "oversample": oversample, "center_target_x": float("nan"),
                "center_target_y": float("nan"), "extent_x": float((raw_w - 1) / 2),
                "extent_y": float((raw_h - 1) / 2), "input_width": raw_w,
                "input_height": raw_h, "mapped_width": raw_w, "mapped_height": raw_h,
                "output_width": raw_w, "output_height": raw_h,
                "input_sum": total, "mapped_sum": total, "output_sum": total,
                "crop_retained_flux": 1.0, "crop_discarded_flux": 0.0, "flux_rel_error": 0.0,
                "overlaps": 0, "wcs_used": "none", "method": "NONE",
            }
            return (raw, int(raw_w), int(raw_h)), record, (
                f"PSF raw {index}/{len(psfs)}: target-frame bypass {raw_w}x{raw_h}"
            )

        items = [(i, p, a) for i, (p, a) in enumerate(zip(psfs, anchors), 1)]
        results = _thread_map_ordered(bypass_one, items, args.preprocess_workers)
        prepared, records = [], []
        for entry, record, log in results:
            prepared.append(entry); records.append(record); print(log, flush=True)
        return prepared, records

    projector = compile_psf_projector(args, directory)

    def geometry_one(item):
        index, path, exposure, anchor = item
        image = _canonical_psf(read_psf_image(path), args.psf_orientation)
        if image.shape[0] % 2 != 1 or image.shape[1] % 2 != 1:
            raise ValueError(f"Detector-frame HybPSF must have odd dimensions: {path} has {image.shape}")
        if not np.isfinite(image).all() or np.any(image < 0):
            raise ValueError(f"Detector-frame HybPSF must be finite and non-negative: {path}")
        ph, pw = image.shape
        source_raw = directory / f"psf_source_{index:03d}.f32"
        vertex_raw = directory / f"psf_vertices_{index:03d}.bin"
        output_raw = directory / f"psf_mapped_{index:03d}.f32"
        np.asarray(image, dtype="<f4").tofile(source_raw)
        record = {
            "exposure": index, "oversample": oversample,
            "input_width": int(pw), "input_height": int(ph), "method": "GWCS_POLYGON",
        }
        record.update(_map_psf_vertices(
            exposure, target, image.shape, anchor, oversample,
            input_shape, args.wcs_mode, vertex_raw
        ))
        row = f"{source_raw}\t{pw}\t{ph}\t{vertex_raw}\t{output_raw}\t{index}\n"
        log = (
            f"PSF geometry {index}/{len(psfs)}: detector anchor=({anchor[0]:.3f},{anchor[1]:.3f}), "
            f"oversample={oversample}, mapped half-extent=({record['extent_x']:.3f},"
            f"{record['extent_y']:.3f}) target pix, WCS={record['wcs_used'].upper()}"
        )
        return row, record, log

    geom_items = [
        (i, p, e, a) for i, (p, e, a) in enumerate(zip(psfs, images, anchors), 1)
    ]
    geom_results = _thread_map_ordered(geometry_one, geom_items, args.preprocess_workers)
    rows, records = [], []
    for row, record, log in geom_results:
        rows.append(row); records.append(record); print(log, flush=True)

    manifest = directory / "psf_project_manifest.txt"
    results = directory / "psf_project_results.csv"
    manifest.write_text("".join(rows))
    subprocess.run([str(projector), str(manifest), str(args.preprocess_workers), str(results)], check=True)
    with results.open() as f:
        projected = list(csv.DictReader(f))
    if len(projected) != len(psfs):
        raise ValueError("PSF polygon projector result count mismatch")

    def finalize_one(item):
        index, path, row, base_record = item
        record = dict(base_record)
        if int(row["exposure"]) != index:
            raise ValueError("PSF polygon projector result ordering mismatch")
        mapped_w, mapped_h = int(row["width"]), int(row["height"])
        mapped_raw = directory / f"psf_mapped_{index:03d}.f32"
        solver_raw = directory / f"psf_{index:03d}.f32"
        mapped_sum = float(row["output_sum"])
        crop = int(args.psf_crop_size)
        warning = None
        if crop:
            if crop > mapped_w or crop > mapped_h:
                raise ValueError(
                    f"--psf-crop-size={crop} exceeds mapped PSF {mapped_w}x{mapped_h} for exposure {index}"
                )
            mapped = np.memmap(mapped_raw, dtype="<f4", mode="r", shape=(mapped_h, mapped_w))
            x0 = (mapped_w - crop) // 2
            y0 = (mapped_h - crop) // 2
            cropped = np.asarray(mapped[y0:y0 + crop, x0:x0 + crop], dtype="<f4")
            cropped.tofile(solver_raw)
            crop_sum = float(np.sum(cropped, dtype=np.float64))
            retained = crop_sum / mapped_sum if mapped_sum > 0 else float("nan")
            solver_w = solver_h = crop
            if 1.0 - retained > 1.0e-3:
                warning = f"WARNING: PSF {index} crop discards {(1.0-retained)*100:.4f}% of mapped flux"
        else:
            solver_raw.write_bytes(mapped_raw.read_bytes())
            crop_sum = mapped_sum
            retained = 1.0
            solver_w, solver_h = mapped_w, mapped_h
        record.update({
            "mapped_width": mapped_w, "mapped_height": mapped_h,
            "output_width": solver_w, "output_height": solver_h,
            "input_sum": float(row["input_sum"]), "mapped_sum": mapped_sum,
            "output_sum": crop_sum, "flux_rel_error": float(row["flux_rel_error"]),
            "crop_retained_flux": retained, "crop_discarded_flux": 1.0 - retained,
            "overlaps": int(row["overlaps"]),
        })
        if args.write_transformed_psf:
            data = np.memmap(solver_raw, dtype="<f4", mode="r", shape=(solver_h, solver_w))
            h = fits.Header()
            h["PSFFRAME"] = "TARGET"
            h["PSFXFORM"] = "GWCS-POLY"
            h["PSFOVR"] = (oversample, "HybPSF samples per native detector pixel")
            h["PSFANCHX"] = float(record["anchor_x"])
            h["PSFANCHY"] = float(record["anchor_y"])
            h["PSFCENX"] = float(record["center_target_x"])
            h["PSFCENY"] = float(record["center_target_y"])
            h["FLUXERR"] = float(record["flux_rel_error"])
            h["PSFCROP"] = (int(args.psf_crop_size), "Odd square crop after GWCS polygon mapping; 0=off")
            h["CROPRETN"] = float(record["crop_retained_flux"])
            fits.PrimaryHDU(np.asarray(data, dtype=np.float32), header=h).writeto(
                output_dir / f"{args.output_stem}_psf{index:02d}_target.fits", overwrite=True
            )
        log = (
            f"PSF raw {index}/{len(psfs)}: input {record['input_width']}x{record['input_height']} "
            f"-> mapped {mapped_w}x{mapped_h} -> solver {solver_w}x{solver_h}; "
            f"polygon flux rel={record['flux_rel_error']:.3e}, crop retained={retained:.9f}"
        )
        return (solver_raw, solver_w, solver_h), record, warning, log

    fin_items = [
        (i, p, row, rec) for i, (p, row, rec) in enumerate(zip(psfs, projected, records), 1)
    ]
    fin_results = _thread_map_ordered(finalize_one, fin_items, args.preprocess_workers)
    prepared = []
    final_records = []
    for entry, record, warning, log in fin_results:
        prepared.append(entry); final_records.append(record)
        if warning: print(warning, flush=True)
        print(log, flush=True)
    return prepared, final_records

def coadd_metadata(images: list[Path], psfs: list[Path], args,
                   target: WCS, psf_records: list[dict]) -> fits.Header:
    """Record inputs, time span, weights, geometry, and solver configuration."""
    headers = [fits.getheader(path, ext=0) for path in images]
    first = headers[0]
    meta = fits.Header()
    for key in ("TELESCOP", "INSTRUME", "DETECTOR", "FILTER", "PUPIL",
                "TARGNAME", "PROGRAM", "OBSERVTN", "VISIT"):
        if key in first:
            values = {str(h.get(key, "")) for h in headers}
            meta[key] = first[key] if len(values) == 1 else "MULTI"
    begin = [h["DATE-BEG"] for h in headers if h.get("DATE-BEG")]
    end = [h["DATE-END"] for h in headers if h.get("DATE-END")]
    if len(begin) == len(end) == len(images):
        meta["DATE-BEG"] = min(begin)
        meta["DATE-END"] = max(end)
        meta["DATE-OBS"] = min(begin)[:10]
        meta["TIME-OBS"] = min(begin)[11:]
    starts = [float(h["EXPSTART"]) for h in headers if h.get("EXPSTART") is not None]
    ends = [float(h["EXPEND"]) for h in headers if h.get("EXPEND") is not None]
    if len(starts) == len(ends) == len(images):
        meta["EXPSTART"] = (min(starts), "First exposure start, MJD")
        meta["EXPEND"] = (max(ends), "Last exposure end, MJD")
        meta["DURATION"] = ((max(ends) - min(starts)) * 86400.0, "Observation span, seconds")
    exposure_times = [float(h.get("EFFEXPTM", h.get("EXPTIME", 0.0))) for h in headers]
    if all(t > 0 for t in exposure_times):
        meta["TEXPTIME"] = (sum(exposure_times), "Sum of effective exposures, seconds")
    meta["DATE"] = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S")
    meta["ORIGIN"] = "UPDC SignedX v17-generic-final"
    meta["NCOMBINE"] = len(images)
    meta["UPDCALG"] = ("BG-WRL-SX", "Signed-x UPDC; persistent mu freeze")
    meta["MUMARKER"] = (1.0e-12, "Diagnostic mu marker for frozen detector pixels")
    meta["ITERMAX"] = args.iterations
    meta["WHTMODE"] = args.weight_mode
    meta["WHTFORM"] = "SUM(W*A/J)"
    meta["BADDQMSK"] = args.bad_dq_mask
    meta["BKGMODE"] = "SCALAR"
    meta["POISPSIG"] = args.poisson_pedestal_sigma
    meta["WCSMODE"] = args.wcs_mode.upper()
    meta["TGTDIST"] = ("NONE", "No target high-order distortion")
    meta["GEOMALG"] = ("OVERLAP", "Polygon overlap geometry")
    meta["GEOAREA"] = ("FULLQUAD", "Full mapped pixel area J")
    meta["LINKFORM"] = ("U32_CSR", "u32 offsets + f32 A/J links")
    meta["LINKCACH"] = ("PERSIST", "Persistent mmap LinkMaps")
    meta["SCICACHE"] = ("RAM", "SCI and WHT loaded once per exposure")
    meta["PSFCACHE"] = ("PERSIST", "Persistent mmap PSF spectra")
    meta["PIXFRAC"] = (1.0, "Full detector pixel footprint")
    meta["KERNEL"] = "square"
    meta["PSR"] = (1.0 / args.upsample, "Target/input pixel scale ratio")
    meta["SAMPFACT"] = (args.upsample, "Target samples per input pixel axis")
    scales = target.proj_plane_pixel_scales()
    meta["PIXSCALE"] = (float(np.mean([s.value for s in scales]) * 3600), "Target arcsec/pixel")
    meta["PSFTRNS"] = (args.psf_orientation == "stored-transposed",
                       "Input PSF transposed into detector x/y order")
    meta["PSFNORM"] = ("UNIT_SUM", "Each PSF normalized to sum one")
    meta["PSFADJ"] = ("FFT_CONJ", "Adjoint PSF via Fourier conjugate")
    meta["PSFFRAME"] = (args.psf_frame.upper(), "Input PSF coordinate frame")
    meta["PSFXFORM"] = (("GWCS-POLY" if args.psf_frame == "detector" else "NONE"),
                        "Full-WCS polygon PSF mapping")
    meta["PSFRSMP"] = (("POLYGON" if args.psf_frame == "detector" else "NONE"),
                        "PSF resampling method")
    meta["PSFOVR"] = (int(args.psf_oversample or args.upsample),
                      "PSF samples/native pixel")
    meta["PSFANCH"] = ("PER-EXP", "Detector location assigned to empirical PSF")
    meta["FFTWORK"] = (args.fft_workers, "Threads inside each FFTW transform")
    meta["EXPWORK"] = (args.exposure_workers, "Concurrent iterative exposure workers")
    meta["PREWORK"] = (args.preprocess_workers, "Preprocess exposure workers")
    meta["CORELANG"] = ("C", "Numerical core implementation")
    meta["PRECISN"] = ("F64XY_F32", "f64 geometry coords; f32 area/solver")
    meta["FFTPLAN"] = (args.fft_plan.upper(), "FFTW final plan mode")
    meta["FFTAUTO"] = (bool(args.fft_autotune), "FFT padding autotune")
    meta["SAVEEVRY"] = (args.save_every, "Iteration save cadence; 0=final")
    meta["ITERPIPE"] = ("FUSED-FREQ", "Fused detector + freq reduction")
    meta["PREPCACH"] = (not args.no_preprocess_cache, "Persistent prep cache")
    meta["PSFCROP"] = (args.psf_crop_size, "Odd PSF crop; 0=off")
    for i, (image, psf) in enumerate(zip(images, psfs), 1):
        meta[f"INFILE{i:02d}"] = (image.name, f"Input exposure {i}")
        meta[f"PSFNAME{i}"] = (psf.name, f"PSF used for exposure {i}")
    meta.add_history("Input Level-2 SCI is not background-subtracted before solving.")
    meta.add_history("Signed x is updated multiplicatively without positivity clipping.")
    meta.add_history("SCI<=0 is retained via a constant pedestal; no noise is added.")
    meta.add_history("If raw mu=H*x+b+pedestal <= 0, the detector pixel freezes at ratio=1.")
    meta.add_history("Frozen pixels remain frozen; 1e-12 is a marker, never a ratio denominator.")
    meta.add_history("Sampling S=A/J uses the full mapped detector pixel area J.")
    meta.add_history("WHT is a coverage sum and omits post-deconvolution covariance.")
    if args.psf_frame == "detector":
        meta.add_history("Each oversampled detector-frame HybPSF is mapped through the full exposure GWCS.")
        meta.add_history("Target PSF resampling uses flux-conserving polygon overlap; no Lanczos/affine approximation.")
        meta.add_history("The mapped detector-center translation is removed; distortion/rotation/scale/shear are retained.")
        meta.add_history("Exposure 1 is resampled like every other exposure to remove detector field distortion.")
    return meta


def psf_geometry_hdu(psf_records: list[dict], psfs: list[Path]) -> fits.BinTableHDU:
    names = [p.name for p in psfs]
    name_width = max(1, max(map(len, names)))
    def arr(key, dtype=np.float64):
        return np.asarray([row[key] for row in psf_records], dtype=dtype)
    wcs_used = [str(row.get("wcs_used", "")) for row in psf_records]
    methods = [str(row.get("method", "")) for row in psf_records]
    wcs_width = max(1, max(map(len, wcs_used)))
    method_width = max(1, max(map(len, methods)))
    cols = [
        fits.Column(name="EXPOSURE", format="J", array=np.arange(1, len(psf_records) + 1, dtype=np.int32)),
        fits.Column(name="PSFNAME", format=f"{name_width}A", array=np.asarray(names)),
        fits.Column(name="ANCHOR_X", format="D", array=arr("anchor_x")),
        fits.Column(name="ANCHOR_Y", format="D", array=arr("anchor_y")),
        fits.Column(name="OVERSAMP", format="J", array=arr("oversample", np.int32)),
        fits.Column(name="CENTER_X", format="D", array=arr("center_target_x")),
        fits.Column(name="CENTER_Y", format="D", array=arr("center_target_y")),
        fits.Column(name="EXTENT_X", format="D", array=arr("extent_x")),
        fits.Column(name="EXTENT_Y", format="D", array=arr("extent_y")),
        fits.Column(name="IN_W", format="J", array=arr("input_width", np.int32)),
        fits.Column(name="IN_H", format="J", array=arr("input_height", np.int32)),
        fits.Column(name="MAP_W", format="J", array=arr("mapped_width", np.int32)),
        fits.Column(name="MAP_H", format="J", array=arr("mapped_height", np.int32)),
        fits.Column(name="OUT_W", format="J", array=arr("output_width", np.int32)),
        fits.Column(name="OUT_H", format="J", array=arr("output_height", np.int32)),
        fits.Column(name="INPUT_SUM", format="D", array=arr("input_sum")),
        fits.Column(name="OUTPUT_SUM", format="D", array=arr("output_sum")),
        fits.Column(name="CROP_RETN", format="D", array=arr("crop_retained_flux")),
        fits.Column(name="CROP_LOSS", format="D", array=arr("crop_discarded_flux")),
        fits.Column(name="FLUX_ERR", format="D", array=arr("flux_rel_error")),
        fits.Column(name="OVERLAPS", format="K", array=arr("overlaps", np.int64)),
        fits.Column(name="METHOD", format=f"{method_width}A", array=np.asarray(methods)),
        fits.Column(name="WCS_USED", format=f"{wcs_width}A", array=np.asarray(wcs_used)),
    ]
    return fits.BinTableHDU.from_columns(cols, name="UPDC_PSFGEOM")


def output_hdul(data: np.ndarray, wht_data: np.ndarray, target_header: fits.Header,
                metadata: fits.Header, input_metrics: list[dict], iteration: int,
                filename: str, weight_mode: str, freeze_rows: list[dict],
                psf_records: list[dict], psfs: list[Path]) -> fits.HDUList:
    primary = fits.PrimaryHDU()
    primary.header.extend(metadata, update=True)
    primary.header["FILENAME"] = filename
    primary.header["UPDCITER"] = iteration
    primary.header["MUFRTOT"] = (sum(r["n_total_frozen"] for r in freeze_rows), "Frozen detector pixels at this iteration")
    primary.header["MUFRNEW"] = (sum(r["n_new_frozen"] for r in freeze_rows), "Newly frozen detector pixels")
    sci = fits.ImageHDU(data=data, header=target_header.copy(), name="SCI")
    sci.header.extend(metadata, update=True)
    sci.header["FILENAME"] = filename
    sci.header["UPDCITER"] = iteration
    sci.header["MUFRTOT"] = primary.header["MUFRTOT"]
    sci.header["MUFRNEW"] = primary.header["MUFRNEW"]
    sci.header["SCIMODE"] = ("SIGNED-X", "Signed multiplicative reconstruction")
    wht_hdu = fits.ImageHDU(data=wht_data, header=target_header.copy(), name="WHT")
    wht_hdu.header.extend(metadata, update=True)
    wht_hdu.header["FILENAME"] = filename
    wht_hdu.header["UPDCITER"] = iteration
    science_unit = str(target_header.get("BUNIT", "MJy/sr"))
    wht_hdu.header["BUNIT"] = "1" if weight_mode == "mask" else f"({science_unit})^-2"
    wht_hdu.header["WHTTYPE"] = "SUM(W*A/J)"
    filenames = [Path(m["source"]).name for m in input_metrics]
    name_width = max(len(name) for name in filenames)
    bkg_tab = fits.BinTableHDU.from_columns([
        fits.Column(name="EXPOSURE", format="J", array=np.arange(1, len(input_metrics) + 1, dtype=np.int32)),
        fits.Column(name="INPUTFILE", format=f"{name_width}A", array=np.array(filenames)),
        fits.Column(name="BACKGROUND", format="D", array=np.array([m["background"] for m in input_metrics])),
        fits.Column(name="BKG_SIGMA", format="D", array=np.array([m["background_sigma"] for m in input_metrics])),
        fits.Column(name="POIS_PEDESTAL", format="D", array=np.array([m["poisson_pedestal"] for m in input_metrics])),
        fits.Column(name="AUTO_PEDESTAL", format="D", array=np.array([m["auto_positivity_pedestal"] for m in input_metrics])),
        fits.Column(name="N_NONPOS", format="K", array=np.array([m["nonpositive_valid_before_pedestal"] for m in input_metrics], dtype=np.int64)),
        fits.Column(name="SKY_VAR_POISSON", format="D", array=np.array([m["sky_var_poisson"] for m in input_metrics])),
    ], name="UPDC_BKG")
    psf_tab = psf_geometry_hdu(psf_records, psfs)
    return fits.HDUList([primary, sci, wht_hdu, bkg_tab, psf_tab])


def resolve_pairs(args, parser: argparse.ArgumentParser) -> tuple[list[Path], list[Path], Path]:
    """Pair exposures and PSFs explicitly or by a shared numeric exposure ID."""
    if args.jwst_dir and args.input_dir:
        parser.error("Specify either --jwst-dir or --input-dir")
    if args.manifest:
        with args.manifest.open(newline="") as f:
            rows = list(csv.DictReader(f))
        if not rows or not {"exposure", "psf"}.issubset(rows[0]):
            parser.error("--manifest needs CSV columns exposure,psf")
        def path_from_manifest(value: str) -> Path:
            path = Path(value)
            return (path if path.is_absolute() else args.manifest.parent / path).resolve()
        images = [path_from_manifest(row["exposure"]) for row in rows]
        psfs = [path_from_manifest(row["psf"]) for row in rows]
        input_dir = args.input_dir or (args.jwst_dir / "NIRCam" / "level2" if args.jwst_dir else images[0].parent)
    else:
        if not args.input_dir and not args.jwst_dir:
            parser.error("Provide --input-dir, --jwst-dir, or --manifest")
        input_dir = args.input_dir or args.jwst_dir / "NIRCam" / "level2"
        if not input_dir.is_dir():
            parser.error(f"Input directory does not exist: {input_dir}")
        images = sorted(input_dir.glob(args.input_glob))
        if args.psf_dir:
            psf_dir = args.psf_dir
        elif args.jwst_dir:
            psf_dir = args.jwst_dir / "PSF"
        else:
            candidates = [input_dir / "PSF", input_dir.parent / "PSF", input_dir]
            psf_dir = next((p for p in candidates if p.is_dir()), None)
        if psf_dir is None or not psf_dir.is_dir():
            parser.error("Provide --psf-dir or place PSFs in an adjacent PSF directory")
        image_paths = {p.resolve() for p in images}
        psfs = sorted(p for p in psf_dir.glob(args.psf_glob)
                      if p.resolve() not in image_paths)
        if len(images) != len(psfs):
            parser.error(f"Found {len(images)} exposures but {len(psfs)} PSFs")
        if args.pairing == "index":
            def exposure_index(path: Path):
                match = re.search(r"_(\d{5})_[A-Za-z0-9]+_(?:cal|crf)", path.name, re.I)
                return int(match.group(1)) if match else None
            def psf_index(path: Path):
                match = re.search(r"_(\d+)\.fits(?:\.gz)?$", path.name, re.I)
                return int(match.group(1)) if match else None
            iidx = [exposure_index(p) for p in images]
            pidx = [psf_index(p) for p in psfs]
            if None in iidx or None in pidx or len(set(iidx)) != len(iidx) or set(iidx) != set(pidx):
                parser.error("Cannot match exposure/PSF numeric IDs; use --manifest or --pairing sorted")
            psf_by_id = dict(zip(pidx, psfs))
            psfs = [psf_by_id[i] for i in iidx]
    if not images or any(not p.is_file() for p in [*images, *psfs]):
        parser.error("Input exposure or PSF list is empty or contains a missing file")
    if len({p.resolve() for p in images}) != len(images):
        parser.error("Exposure paths must be unique")
    return [p.resolve() for p in images], [p.resolve() for p in psfs], input_dir.resolve()



CACHE_SCHEMA = "updc-v16-preprocess-1"

def _file_signature(path: Path) -> dict:
    st = path.stat()
    return {"path": str(path.resolve()), "size": int(st.st_size), "mtime_ns": int(st.st_mtime_ns)}

def _source_digest(paths) -> str:
    h = hashlib.sha256()
    for path in paths:
        p = Path(path)
        h.update(p.name.encode()); h.update(b"\0"); h.update(p.read_bytes())
    return h.hexdigest()

def _preprocess_cache_key(args, images: list[Path], psfs: list[Path],
                          input_shape: tuple[int, int], anchors: list[tuple[float, float]]) -> str:
    payload = {
        "schema": CACHE_SCHEMA,
        "images": [_file_signature(p) for p in images],
        "psfs": [_file_signature(p) for p in psfs],
        "input_shape": list(input_shape),
        "anchors": [[float(x), float(y)] for x, y in anchors],
        "wcs_mode": args.wcs_mode,
        "upsample": int(args.upsample),
        "psf_oversample": int(args.psf_oversample or args.upsample),
        "psf_orientation": args.psf_orientation,
        "psf_frame": args.psf_frame,
        "psf_crop_size": int(args.psf_crop_size),
        "weight_mode": args.weight_mode,
        "bad_dq_mask": int(args.bad_dq_mask),
        "poisson_pedestal_sigma": float(args.poisson_pedestal_sigma),
    }
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(blob).hexdigest()[:24]

def _cache_marker_valid(directory: Path, marker: dict, nexp: int, input_shape: tuple[int, int]) -> bool:
    if marker.get("schema") != CACHE_SCHEMA or len(marker.get("prepared", [])) != nexp:
        return False
    ih, iw = input_shape
    plane_bytes = ih * iw * 4
    try:
        for sci, wht in marker["prepared"]:
            if (directory / sci).stat().st_size != plane_bytes or (directory / wht).stat().st_size != plane_bytes:
                return False
        for psf, pw, ph in marker["prepared_psfs"]:
            if (directory / psf).stat().st_size != int(pw) * int(ph) * 4:
                return False
        for i in range(1, nexp + 1):
            if (directory / f"vertices_{i}.bin").stat().st_size != 2 * (iw + 1) * (ih + 1) * 8:
                return False
    except (KeyError, OSError, TypeError, ValueError):
        return False
    return True

def _load_python_preprocess_cache(directory: Path, nexp: int, input_shape: tuple[int, int]):
    marker_path = directory / "preprocess_cache.json"
    if not marker_path.is_file():
        return None
    try:
        marker = json.loads(marker_path.read_text())
        if not _cache_marker_valid(directory, marker, nexp, input_shape):
            return None
        header = fits.Header.fromstring(marker["target_header"], sep="\n")
        target = WCS(header)
        prepared = [(directory / a, directory / b) for a, b in marker["prepared"]]
        prepared_psfs = [(directory / a, int(w), int(h)) for a, w, h in marker["prepared_psfs"]]
        return target, header, (int(marker["width"]), int(marker["height"])), prepared, marker["input_metrics"], prepared_psfs, marker["psf_records"]
    except Exception as exc:
        print(f"WARNING: ignoring invalid preprocessing cache: {exc}", flush=True)
        return None

def _save_python_preprocess_cache(directory: Path, target: WCS, header: fits.Header,
                                  shape: tuple[int, int], prepared, input_metrics,
                                  prepared_psfs, psf_records) -> None:
    width, height = shape
    marker = {
        "schema": CACHE_SCHEMA,
        "width": int(width), "height": int(height),
        "target_header": header.tostring(sep="\n", endcard=True, padding=False),
        "prepared": [[Path(a).name, Path(b).name] for a, b in prepared],
        "input_metrics": input_metrics,
        "prepared_psfs": [[Path(a).name, int(w), int(h)] for a, w, h in prepared_psfs],
        "psf_records": psf_records,
    }
    tmp = directory / f"preprocess_cache.json.tmp.{os.getpid()}"
    tmp.write_text(json.dumps(marker, indent=2, sort_keys=True))
    os.replace(tmp, directory / "preprocess_cache.json")

def resolve_psf_anchors(args, images: list[Path], input_shape: tuple[int, int]) -> list[tuple[float, float]]:
    """Resolve per-exposure detector coordinates where each PSF was measured."""
    h, w = input_shape
    default = ((w - 1) / 2.0, (h - 1) / 2.0)
    if (args.psf_anchor_x is None) ^ (args.psf_anchor_y is None):
        raise ValueError("--psf-anchor-x and --psf-anchor-y must be supplied together")
    if args.psf_anchor_x is not None:
        anchors = [(float(args.psf_anchor_x), float(args.psf_anchor_y)) for _ in images]
    elif args.manifest:
        with args.manifest.open(newline="") as f:
            rows = list(csv.DictReader(f))
        has_x = rows and "psf_x" in rows[0]
        has_y = rows and "psf_y" in rows[0]
        if has_x != has_y:
            raise ValueError("Manifest must contain both psf_x and psf_y, or neither")
        if has_x:
            anchors = []
            for i, row in enumerate(rows, 1):
                try:
                    anchors.append((float(row["psf_x"]), float(row["psf_y"])))
                except (TypeError, ValueError) as exc:
                    raise ValueError(f"Invalid psf_x/psf_y in manifest row {i}") from exc
        else:
            anchors = [default for _ in images]
    else:
        anchors = [default for _ in images]
    for i, (x, y) in enumerate(anchors, 1):
        if not (np.isfinite(x) and np.isfinite(y)):
            raise ValueError(f"Non-finite PSF anchor for exposure {i}")
        if not (-0.5 <= x <= w - 0.5 and -0.5 <= y <= h - 0.5):
            raise ValueError(f"PSF anchor ({x},{y}) for exposure {i} lies outside the detector")
    return anchors


def _legacy_config_candidate(argv) -> Path | None:
    """Return a legacy positional TOML only when it is argv[0].

    The packaged CLI calls this frozen v17 engine with option-only argv such as
    ``--input-dir /path``.  ``argparse.parse_known_args`` cannot safely discover
    an optional positional config in that situation because the value of an
    unknown option may be consumed as the positional.  Only the first token is
    therefore eligible for the legacy config role.
    """
    seq = list(sys.argv[1:] if argv is None else argv)
    if seq and not str(seq[0]).startswith("-"):
        return Path(seq[0])
    return None


def main(argv=None) -> None:
    argv = list(sys.argv[1:] if argv is None else argv)
    legacy_config = _legacy_config_candidate(argv)
    config_defaults = {}
    if legacy_config is not None:
        try:
            config_defaults = load_config_defaults(legacy_config)
        except (OSError, ValueError, ModuleNotFoundError) as exc:
            pre = argparse.ArgumentParser(add_help=False)
            pre.add_argument("config", nargs="?", type=Path, help="TOML configuration file")
            pre.error(f"cannot load config: {exc}")

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", nargs="?", type=Path, help="TOML configuration file; command-line options override it")
    parser.add_argument("--iterations", type=int, default=9)
    parser.add_argument("--workers", type=int, default=1,
                        help="Backward-compatible FFT worker default; overridden by --fft-workers")
    parser.add_argument("--fft-workers", type=int,
                        help="Threads used inside each FFTW transform; default equals --workers")
    parser.add_argument("--exposure-workers", type=int, default=1,
                        help="Concurrent exposure workers for iterative forward/backprojection; default 1")
    parser.add_argument("--preprocess-workers", type=int,
                        help="Concurrent exposure workers before iteration (Level-2 prep, GWCS/PSF mapping, LinkMaps, PSF spectra, drizzle, sensitivity); default min(4, CPU count, exposure count)")
    parser.add_argument("--fft-plan", choices=("estimate", "measure"), default="measure",
                        help="FFTW planning mode; measure is slower to plan but usually faster over many iterations")
    parser.add_argument("--fft-autotune", dest="fft_autotune", action="store_true",
                        help="Enable runtime benchmarking/caching of valid FFT padding shapes")
    parser.add_argument("--no-fft-autotune", dest="fft_autotune", action="store_false",
                        help="Disable runtime benchmarking of several valid FFT padding shapes")
    parser.add_argument("--retune-fft", dest="retune_fft", action="store_true",
                        help="Ignore a cached FFT autotune choice and benchmark candidate shapes again")
    parser.add_argument("--no-retune-fft", dest="retune_fft", action="store_false",
                        help="Allow reuse of a cached FFT autotune choice")
    parser.set_defaults(fft_autotune=True, retune_fft=False)
    parser.add_argument("--fftw-wisdom", type=Path,
                        help="Persistent FFTW wisdom file; default ~/.cache/updc-jwst/fftwf_wisdom")
    parser.add_argument("--save-every", type=int,
                        help="Save every Nth iteration; 0 saves only the final iteration (default)")
    parser.add_argument("--save-iterations", choices=("all", "final"),
                        help="Deprecated compatibility alias: all -> --save-every 1, final -> --save-every 0")
    parser.add_argument("--input-dir", type=Path, help="Directory containing Level-2 exposure FITS files")
    parser.add_argument("--psf-dir", type=Path, help="Directory containing matching per-exposure PSFs")
    parser.add_argument("--input-glob", default="*_cal*.fits", help="Exposure filename glob within input-dir")
    parser.add_argument("--psf-glob", default="*.fits", help="PSF filename glob within psf-dir")
    parser.add_argument("--pairing", choices=("index", "sorted"), default="index",
                        help="Match numeric exposure/PSF IDs, or explicitly pair sorted files")
    parser.add_argument("--manifest", type=Path,
                        help="CSV with exposure,psf columns; overrides automatic pairing")
    parser.add_argument("--jwst-dir", type=Path,
                        help="Convenience alias for JWST/NIRCam/level2 and JWST/PSF")
    parser.add_argument("--upsample", type=int, default=2,
                        help="Target samples per native detector pixel axis; default 2")
    parser.add_argument("--psf-orientation", choices=("stored-transposed", "as-is"),
                        default="stored-transposed",
                        help="Storage orientation of the input PSF relative to detector x/y axes")
    parser.add_argument("--psf-frame", choices=("detector", "target"), default="detector",
                        help="Input PSF frame; detector maps HybPSFs through full GWCS with polygon overlap")
    parser.add_argument("--psf-oversample", type=int,
                        help="HybPSF samples per native detector pixel; default equals --upsample")
    parser.add_argument("--psf-crop-size", type=int, default=1431,
                        help="Odd square PSF size after GWCS polygon resampling; 0 disables cropping; default 1431")
    parser.add_argument("--psf-anchor-x", type=float,
                        help="Global detector x coordinate where PSFs were measured; default detector center")
    parser.add_argument("--psf-anchor-y", type=float,
                        help="Global detector y coordinate where PSFs were measured; default detector center")
    parser.add_argument("--write-transformed-psf", action="store_true",
                        help="Write per-exposure full-GWCS polygon-mapped target PSFs as diagnostic FITS files")
    parser.add_argument("--weight-mode", choices=("rnoise-sky", "totalvar", "mask"),
                        default="rnoise-sky",
                        help="Default avoids source-dependent per-pixel VAR_POISSON weighting")
    parser.add_argument("--bad-dq-mask", type=int, default=531,
                        help="Bits rejected by default: DO_NOT_USE(1)+SATURATED(2)+OUTLIER(16)+NON_SCIENCE(512)")
    parser.add_argument("--poisson-pedestal-sigma", type=float, default=0.0,
                        help="Optional effective-background floor in sky sigmas; automatic SCI<=0 pedestal remains active")
    parser.add_argument("--wcs-mode", choices=("auto", "gwcs", "fits"), default="auto",
                        help="Use JWST meta.wcs when available; FITS is a fallback/diagnostic mode")
    parser.add_argument("--output-dir", type=Path,
                        help="Output directory; defaults beside input-dir")
    parser.add_argument("--output-stem", help="Filename stem; default derived from filter and detector")
    parser.add_argument("--keep-work", action="store_true",
                        help="Retain intermediate FITS, vertices, links, raw iterations and solver timing")
    parser.add_argument("--work-dir", type=Path,
                        help="Explicit persistent intermediate directory (implies --keep-work)")
    parser.add_argument("--cache-dir", type=Path,
                        help="Persistent content-addressed preprocessing cache; default ~/.cache/updc-jwst/data")
    parser.add_argument("--cache-reset", choices=("none", "current", "all"), default="none",
                        help="Reset persistent preprocessing cache: current dataset, all dataset caches, or none")
    parser.add_argument("--cache-prune-after-days", type=float, default=0.0,
                        help="On startup, delete dataset caches unused for this many days; 0 disables automatic pruning")
    parser.add_argument("--preprocess-cache", dest="no_preprocess_cache", action="store_false",
                        help="Enable persistent Python/C preprocessing cache")
    parser.add_argument("--no-preprocess-cache", dest="no_preprocess_cache", action="store_true",
                        help="Disable persistent Python/C preprocessing cache and use a temporary work directory")
    parser.add_argument("--metrics-json", type=Path,
                        help="Optional output path for detailed metrics JSON")
    parser.add_argument("--cc", default="cc", help="C compiler used for the numerical core")
    parser.add_argument("--c-core-binary", type=Path,
                        help="Use a prebuilt updc_jwst_c binary instead of compiling the bundled C core")
    parser.add_argument("--psf-projector-binary", type=Path,
                        help="Use a prebuilt updc_psf_projector binary for detector-frame HybPSF mapping")
    if config_defaults:
        parser.set_defaults(**{k: v for k, v in config_defaults.items() if k != "config_file"})
    args = parser.parse_args(argv)
    if legacy_config is not None:
        args.config = legacy_config.expanduser().resolve()
        print(f"Configuration: {args.config}", flush=True)
    if args.iterations < 1:
        parser.error("iterations must be positive")
    if args.workers < 1:
        parser.error("workers must be positive")
    if args.fft_workers is None:
        args.fft_workers = args.workers
    if args.fft_workers < 1:
        parser.error("fft-workers must be positive")
    if args.exposure_workers < 1:
        parser.error("exposure-workers must be positive")
    if args.preprocess_workers is not None and args.preprocess_workers < 1:
        parser.error("preprocess-workers must be positive")
    if args.cache_prune_after_days < 0:
        parser.error("cache-prune-after-days must be non-negative")
    if args.no_preprocess_cache and (args.cache_reset != "none" or args.cache_prune_after_days > 0):
        print("WARNING: cache reset/prune policy ignored because persistent preprocessing cache is disabled", flush=True)
    if args.save_every is None:
        args.save_every = 1 if args.save_iterations == "all" else 0
    elif args.save_iterations is not None:
        parser.error("use either --save-every or --save-iterations, not both")
    if args.save_every < 0:
        parser.error("save-every must be >= 0")
    if args.poisson_pedestal_sigma < 0:
        parser.error("poisson pedestal sigma cannot be negative")
    if args.upsample < 1:
        parser.error("upsample must be at least one")
    if args.psf_oversample is not None and args.psf_oversample < 1:
        parser.error("psf-oversample must be at least one")
    if args.psf_crop_size < 0 or (args.psf_crop_size and args.psf_crop_size % 2 != 1):
        parser.error("psf-crop-size must be 0 or a positive odd integer")
    if args.output_stem and Path(args.output_stem).name != args.output_stem:
        parser.error("output-stem must be one filename component")

    images, psfs, input_dir = resolve_pairs(args, parser)
    if args.preprocess_workers is None:
        args.preprocess_workers = min(4, len(images), os.cpu_count() or 1)
    first_header = sci_header(images[0])
    input_shape = (int(first_header["NAXIS2"]), int(first_header["NAXIS1"]))
    first_primary = fits.getheader(images[0], 0)
    for path in images[1:]:
        h = sci_header(path)
        if (int(h["NAXIS2"]), int(h["NAXIS1"])) != input_shape:
            parser.error("All input SCI images must have the same dimensions")
        if h.get("BUNIT") != first_header.get("BUNIT"):
            parser.error("All input SCI images must have the same BUNIT")
        primary = fits.getheader(path, 0)
        for key in ("FILTER", "PUPIL"):
            if primary.get(key) != first_primary.get(key):
                parser.error(f"Mixed {key} values require separate runs")
    if args.output_stem is None:
        detectors = {fits.getheader(path, 0).get("DETECTOR") for path in images}
        detector = first_primary.get("DETECTOR", "mosaic") if len(detectors) == 1 else "mosaic"
        tag = f"updc_{first_primary.get('FILTER', 'nircam')}_{detector}"
        args.output_stem = re.sub(r"[^A-Za-z0-9_]+", "_", tag).strip("_").lower()
    try:
        psf_anchors = resolve_psf_anchors(args, images, input_shape)
    except ValueError as exc:
        parser.error(str(exc))
    output = (args.output_dir.resolve() if args.output_dir else
              (args.jwst_dir.resolve() / "output" if args.jwst_dir else input_dir.parent / "updc_output"))
    if output == input_dir:
        parser.error("Output directory must differ from input directory")

    output.mkdir(parents=True, exist_ok=True)
    temp_work = None
    preprocess_cache_hit = False
    cache_key = None
    if args.work_dir:
        work = args.work_dir.resolve(); work.mkdir(parents=True, exist_ok=True)
    elif args.keep_work:
        work = (output / f"{args.output_stem}_work").resolve(); work.mkdir(parents=True, exist_ok=True)
    elif not args.no_preprocess_cache:
        cache_root = (args.cache_dir.expanduser() if args.cache_dir else Path.home() / ".cache" / "updc-jwst" / "data")
        try:
            cache_root.mkdir(parents=True, exist_ok=True)
            cache_key = _preprocess_cache_key(args, images, psfs, input_shape, psf_anchors)
            removed = reset_dataset_caches(cache_root, args.cache_reset, cache_key)
            if removed:
                print(f"Cache reset ({args.cache_reset}): removed {len(removed)} dataset cache(s)", flush=True)
            pruned = prune_dataset_caches(cache_root, args.cache_prune_after_days, keep_key=cache_key)
            if pruned:
                print(f"Cache auto-prune: removed {len(pruned)} dataset cache(s) unused for >= {args.cache_prune_after_days:g} days", flush=True)
            work = (cache_root / cache_key).resolve(); work.mkdir(parents=True, exist_ok=True)
            touch_dataset_cache(work)
            print(f"Preprocessing cache: {work}", flush=True)
        except (OSError, ValueError) as exc:
            print(f"WARNING: preprocessing cache unavailable ({exc}); using temporary work directory", flush=True)
            temp_work = tempfile.TemporaryDirectory(prefix="updc-jwst-"); atexit.register(temp_work.cleanup); work = Path(temp_work.name)
    else:
        temp_work = tempfile.TemporaryDirectory(prefix="updc-jwst-"); atexit.register(temp_work.cleanup); work = Path(temp_work.name)

    wisdom_path = (args.fftw_wisdom.expanduser() if args.fftw_wisdom else
                   Path.home() / ".cache" / "updc-jwst" / "fftwf_wisdom")
    try:
        wisdom_path.parent.mkdir(parents=True, exist_ok=True)
        wisdom_path = wisdom_path.resolve()
    except OSError as exc:
        wisdom_path = work / "fftwf_wisdom"
        print(f"WARNING: persistent FFTW wisdom unavailable ({exc}); using {wisdom_path}", flush=True)

    print(f"Preprocessing exposure workers: {args.preprocess_workers}", flush=True)
    cached = None if args.no_preprocess_cache or args.work_dir or args.keep_work else _load_python_preprocess_cache(work, len(images), input_shape)
    if cached is not None:
        with timer("load preprocessing cache"):
            target, header, (width, height), prepared, input_metrics, prepared_psfs, psf_records = cached
        preprocess_cache_hit = True
        print(f"Persistent preprocessing cache hit: target {width} x {height}", flush=True)
        if args.write_transformed_psf:
            for index, ((raw, pw, ph), record) in enumerate(zip(prepared_psfs, psf_records), 1):
                data = np.memmap(raw, dtype="<f4", mode="r", shape=(ph, pw))
                h = fits.Header(); h["PSFFRAME"]="TARGET"; h["PSFXFORM"]="GWCS-POLY"; h["PSFCROP"]=int(args.psf_crop_size)
                h["CROPRETN"]=float(record.get("crop_retained_flux",1.0))
                fits.PrimaryHDU(np.asarray(data,dtype=np.float32),header=h).writeto(
                    output / f"{args.output_stem}_psf{index:02d}_target.fits", overwrite=True)
    else:
        with timer("prepare Level-2 SCI/WHT/BKG"):
            prepared, input_metrics = prepare_weighted_inputs(
                images, work, args.weight_mode, args.bad_dq_mask, args.poisson_pedestal_sigma,
                input_shape, args.preprocess_workers
            )
        with timer("make target WCS"):
            target, header, (width, height) = make_target(images, args.wcs_mode,
                                                         args.upsample, input_shape,
                                                         args.preprocess_workers)
        print(f"Target grid: {width} x {height}", flush=True)
        with timer("GWCS polygon-map HybPSFs"):
            prepared_psfs, psf_records = prepare_psf_raw(
                psfs, images, target, psf_anchors, input_shape, work, args, output
            )
        with timer("prepare polygon vertices"):
            prepare_maps(images, target, work, args.wcs_mode, input_shape, args.preprocess_workers)
        if not args.no_preprocess_cache and not args.work_dir and not args.keep_work:
            with timer("save preprocessing cache"):
                _save_python_preprocess_cache(work, target, header, (width, height), prepared,
                                              input_metrics, prepared_psfs, psf_records)
    metadata = coadd_metadata(images, psfs, args, target, psf_records)

    manifest = work / "manifest.txt"
    with timer("write manifest"):
        rows = []
        for i, ((science_raw, weight_raw), (psf_raw, psf_width, psf_height), metric) in enumerate(
                zip(prepared, prepared_psfs, input_metrics), 1):
            rows.append(
                f"{science_raw}\t{weight_raw}\t{psf_raw}\t{psf_width}\t{psf_height}\t"
                f"{work / f'vertices_{i}.bin'}\t{metric['background']:.17g}\t"
                f"{metric['poisson_pedestal']:.17g}\n"
            )
        manifest.write_text("".join(rows))

    c_core = PACKAGE / "c_core"
    if args.c_core_binary:
        binary = args.c_core_binary.expanduser().resolve()
        if not binary.is_file():
            parser.error(f"C core binary does not exist: {binary}")
    else:
        binary = work / ("updc_jwst_c.exe" if os.name == "nt" else "updc_jwst_c")
        core_sources = [c_core / "updc_core.c", c_core / "overlap.c", c_core / "overlap_geometry.h"]
        core_digest = _source_digest(core_sources)
        build_env, removed = sanitize_c_environment()
        if removed:
            details = "; ".join(f"{k}: {', '.join(v)}" for k, v in sorted(removed.items()))
            print(f"Build environment: ignored C++ include paths for C compilation ({details})", flush=True)
        cc = shutil.which(args.cc, path=build_env.get("PATH")) or args.cc
        compiler = compiler_identity(cc, build_env)

        atomic_ok, atomic_log = probe_c11_atomics(cc, build_env)
        if not atomic_ok:
            raise RuntimeError(
                "The selected C compiler cannot build a C11 <stdatomic.h> probe. "
                "On HPC systems, check for C++ header shadowing in C_INCLUDE_PATH/CPATH. "
                "UPDC already filters include/c++ paths, so a remaining failure usually means the compiler "
                "installation itself is incomplete.\n" + atomic_log[-4000:]
            )

        fftw = discover_fftw(cc, build_env)
        native = native_arch_flag(cc, build_env)
        flags = ["-O3", "-std=c11", "-DNDEBUG", "-Wall", "-Wextra",
                 "-DUPDC_FFTW_THREADS", "-DUPDC_PARALLEL", "-pthread"]
        if native:
            flags.insert(0, native)
        else:
            print("WARNING: native CPU flag probe failed; building portable C core without it", flush=True)
        compile_flags = [*flags, *fftw["cflags"]]
        libs = list(fftw["libs"])
        compiler_sig = Path(str(binary) + ".build.json")
        expected = build_signature(core_digest, compiler, compile_flags, libraries=libs, fftw=fftw)
        compile_needed = not (binary.is_file() and signature_matches(compiler_sig, expected))
        if not compile_needed:
            print(f"C core binary persistent-cache hit: {binary}", flush=True)
        if compile_needed:
          with timer("compile C core"):
            cmd = [cc, *compile_flags, "-o", str(binary),
                   str(c_core / "updc_core.c"), str(c_core / "overlap.c"), *libs]
            subprocess.run(cmd, check=True, env=build_env)
            write_signature(compiler_sig, expected)
            # Remove legacy source-only signatures so they cannot be mistaken for the new cache contract.
            legacy_sig = Path(str(binary) + ".sha256")
            try:
                legacy_sig.unlink()
            except FileNotFoundError:
                pass


    subprocess.run([str(binary), "--self-test"], check=True,
                   env=(build_env if 'build_env' in locals() else None))

    persistent_cache_safe = bool(cache_key is not None and not args.no_preprocess_cache and not args.work_dir and not args.keep_work)
    run_temp = None
    if persistent_cache_safe:
        run_temp = tempfile.TemporaryDirectory(prefix="updc-run-")
        run_work = Path(run_temp.name)
    else:
        run_work = work
    raw_prefix = run_work / "iter"
    raw_wht = run_work / "target_wht.bin"
    solver_env = os.environ.copy()
    solver_env["UPDC_PERSIST_CACHE"] = "1" if persistent_cache_safe else "0"
    with timer("run signed-x multiplicative UPDC"):
        subprocess.run([
            str(binary), str(manifest), str(width), str(height),
            str(input_shape[1]), str(input_shape[0]), str(len(images)),
            "1" if (args.psf_frame == "target" and args.psf_orientation == "stored-transposed") else "0",
            str(args.iterations), str(raw_prefix), str(raw_wht),
            str(args.fft_workers), str(args.exposure_workers), str(args.preprocess_workers),
            "1" if args.fft_plan == "measure" else "0",
            "1" if args.fft_autotune else "0", "1" if args.retune_fft else "0", str(wisdom_path),
            str(args.save_every)
        ], check=True, env=solver_env)

    freeze_path = run_work / "iter_mufreeze.csv"
    with freeze_path.open() as f:
        freeze_rows = [
            {
                "iteration": int(row["iteration"]),
                "exposure": int(row["exposure"]),
                "n_new_frozen": int(row["n_new_frozen"]),
                "n_total_frozen": int(row["n_total_frozen"]),
                "min_raw_mu": float(row["min_raw_mu"]),
                "mu_marker": float(row["mu_marker"]),
                "frozen_ratio": float(row["frozen_ratio"]),
            }
            for row in csv.DictReader(f)
        ]
    if len(freeze_rows) != args.iterations * len(images):
        raise ValueError("Solver mu-freeze CSV row mismatch")

    iteration_files = []
    if args.save_every > 0:
        output_iterations = list(range(args.save_every, args.iterations + 1, args.save_every))
        if not output_iterations or output_iterations[-1] != args.iterations:
            output_iterations.append(args.iterations)
    else:
        output_iterations = [args.iterations]
    with timer("write iteration FITS"):
        wht_data = np.memmap(raw_wht, dtype="<f4", mode="r", shape=(height, width))
        for iteration in output_iterations:
            raw = run_work / f"iter_{iteration:02d}.bin"
            data = np.memmap(raw, dtype="<f4", mode="r", shape=(height, width))
            iteration_file = output / f"{args.output_stem}_iter{iteration:02d}_i2d.fits"
            rows = [row for row in freeze_rows if row["iteration"] == iteration]
            with output_hdul(data, wht_data, header, metadata, input_metrics,
                             iteration, iteration_file.name, args.weight_mode, rows,
                             psf_records, psfs) as hdul:
                hdul.writeto(iteration_file, overwrite=True, checksum=True)
            iteration_files.append(iteration_file)
            print(f"Iteration {iteration}/{args.iterations} FITS: {iteration_file}", flush=True)

    with timer("validate output FITS"):
        if len(iteration_files) != len(output_iterations):
            raise ValueError("Missing requested iteration FITS files")
        for path in iteration_files:
            with fits.open(path, memmap=True, checksum=True) as check:
                if not np.isfinite(check["SCI"].data).all():
                    raise ValueError(f"SCI contains nonfinite values: {path}")
                if not np.isfinite(check["WHT"].data).all() or np.any(check["WHT"].data < 0):
                    raise ValueError(f"WHT contains invalid values: {path}")
                if len(check["UPDC_BKG"].data) != len(images):
                    raise ValueError(f"UPDC_BKG row mismatch: {path}")
                if len(check["UPDC_PSFGEOM"].data) != len(images):
                    raise ValueError(f"UPDC_PSFGEOM row mismatch: {path}")
                for name in ("SCI", "WHT", "UPDC_BKG", "UPDC_PSFGEOM"):
                    if check[name].verify_checksum() != 1 or check[name].verify_datasum() != 1:
                        raise ValueError(f"{name} FITS checksum failure: {path}")

    solver_timing_path = run_work / "iter_timing.csv"
    with solver_timing_path.open() as f:
        solver_stages = [{"stage": row["stage"], "wall_seconds": float(row["wall_seconds"])}
                         for row in csv.DictReader(f)]
    if args.metrics_json:
        metrics_path = args.metrics_json.resolve()
        metrics_path.parent.mkdir(parents=True, exist_ok=True)
        metrics_path.write_text(json.dumps(_json_finite({
            "weight_mode": args.weight_mode,
            "bad_dq_mask": args.bad_dq_mask,
            "poisson_pedestal_sigma": args.poisson_pedestal_sigma,
            "wcs_mode_requested": args.wcs_mode,
            "input_shape": list(input_shape),
            "target_shape": [height, width],
            "upsample": args.upsample,
            "psf_orientation": args.psf_orientation,
            "psf_frame": args.psf_frame,
            "psf_oversample": int(args.psf_oversample or args.upsample),
            "psf_mapping": "full_gwcs_polygon_overlap" if args.psf_frame == "detector" else "target_frame_bypass",
            "psf_geometry": [
                {k: (v.tolist() if isinstance(v, np.ndarray) else v)
                 for k, v in row.items()}
                for row in psf_records
            ],
            "performance": {
                "workers": args.workers,
                "fft_workers": args.fft_workers,
                "exposure_workers": args.exposure_workers,
                "preprocess_workers": args.preprocess_workers,
                "save_every": args.save_every,
                "psf_crop_size": args.psf_crop_size,
                "fft_plan": args.fft_plan,
                "fft_autotune": args.fft_autotune,
                "retune_fft": args.retune_fft,
                "fftw_wisdom": str(wisdom_path),
                "save_iterations": args.save_iterations,
                "core_language": "c",
                "c_core_binary": str(binary),
                "compiler": args.cc if not args.c_core_binary else "prebuilt",
                "target_cpu": "native" if not args.c_core_binary else "prebuilt",
                "preprocess_cache_enabled": not args.no_preprocess_cache,
                "preprocess_cache_hit": preprocess_cache_hit,
                "preprocess_cache_key": cache_key,
                "preprocess_cache_dir": str(work) if not args.no_preprocess_cache else None,
                "runtime_work_dir": str(run_work),
                "config_file": str(args.config) if args.config else None,
                "iteration_pipeline": "fused_forward_freezemu_scatter_frequency_reduction",
            },
            "inputs": input_metrics,
            "iterations": args.iterations,
            "geometry_backend": "C Sutherland-Hodgman adaptive DELTA16/ABS32 mmap CSR with precomputed A/J",
            "output_semantics": {
                "SCI": "signed multiplicative reconstruction x; no positivity clipping",
                "mu_freeze": "raw_mu<=0 permanently freezes detector pixel at ratio=1",
                "mu_marker": 1.0e-12,
                "source_residual_split": False,
            },
            "mu_freeze_diagnostics": freeze_rows,
            "iteration_files": [str(p) for p in iteration_files],
            "solver_stages": solver_stages,
            "stage_timing": [{"name": s.name, "cpu_seconds": s.cpu, "wall_seconds": s.wall} for s in STAGES],
        }), indent=2, allow_nan=False))
        print(f"Metrics {metrics_path}")
    report_total()
    print(f"Wrote {len(iteration_files)} iteration FITS files to {output}")
    # Release memmaps before removing the per-run staging directory.
    try:
        del data, wht_data
    except NameError:
        pass
    if run_temp is not None:
        run_temp.cleanup()
    if temp_work is not None:
        temp_work.cleanup()


if __name__ == "__main__":
    main()
