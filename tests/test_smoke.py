"""End-to-end test of variable exposure/PSF sizes and temporary cleanup."""
from __future__ import annotations

import shutil
import subprocess
import sys
import json
import csv
from pathlib import Path

import numpy as np
import pytest
fits = pytest.importorskip("astropy.io.fits")


def test_two_exposures_produce_only_two_iterations(tmp_path: Path) -> None:
    if not shutil.which("pkg-config"):
        pytest.skip("pkg-config is not installed")
    if subprocess.run(["pkg-config", "--exists", "fftw3f"]).returncode:
        pytest.skip("FFTW3f development libraries are not installed")

    images, psfs, output = (tmp_path / name for name in ("level2", "PSF", "results/output"))
    images.mkdir()
    psfs.mkdir()
    yy, xx = np.indices((64, 64), dtype=float)
    for index in (1, 2):
        primary = fits.PrimaryHDU()
        primary.header["FILTER"] = "F200W"
        primary.header["PUPIL"] = "CLEAR"
        primary.header["DETECTOR"] = "NRCA1"
        primary.header["EFFEXPTM"] = 10.0
        primary.header["DATE-BEG"] = f"2022-01-01T00:0{index}:00.000"
        primary.header["DATE-END"] = f"2022-01-01T00:0{index}:10.000"
        sci = (-0.1 + 0.02 * np.sin(xx / 2) + 2 * np.exp(
            -((xx - 31.5 - 0.2 * index) ** 2
              + (yy - 30.5) ** 2) / (2 * 2.2 ** 2))).astype("f4")
        h = fits.Header()
        h["CTYPE1"], h["CTYPE2"] = "RA---TAN", "DEC--TAN"
        h["CRVAL1"], h["CRVAL2"] = 100.0, 20.0
        h["CRPIX1"], h["CRPIX2"] = 32.5 + 0.2 * index, 32.5
        h["CD1_1"], h["CD1_2"] = -1e-5, 0.0
        h["CD2_1"], h["CD2_2"] = 0.0, 1e-5
        h["BUNIT"] = "Jy/pixel"
        hdus = [
            primary,
            fits.ImageHDU(sci, header=h, name="SCI"),
            fits.ImageHDU(np.zeros((64, 64), dtype="u4"), name="DQ"),
            fits.ImageHDU(np.full((64, 64), 0.001, dtype="f4"), name="VAR_POISSON"),
            fits.ImageHDU(np.full((64, 64), 0.0005, dtype="f4"), name="VAR_RNOISE"),
            fits.ImageHDU(np.full((64, 64), 0.0001, dtype="f4"), name="VAR_FLAT"),
        ]
        fits.HDUList(hdus).writeto(
            images / f"jw00001001001_02101_{index:05d}_nrca1_cal.fits"
        )
        sy, sx = (7, 9) if index == 1 else (9, 7)
        py, px = np.indices((sy, sx))
        psf = np.exp(-((px - (sx - 1) / 2) ** 2
                       + (py - (sy - 1) / 2) ** 2) / (2 * 1.2 ** 2))
        fits.PrimaryHDU((psf / psf.sum()).astype("f4")).writeto(psfs / f"PSF_{index}.fits")

    script = Path(__file__).resolve().parents[1] / "run_updc.py"
    subprocess.run([
        sys.executable, str(script), "--input-dir", str(images),
        "--psf-dir", str(psfs), "--output-dir", str(output),
        "--wcs-mode", "fits", "--psf-crop-size", "0",
        "--iterations", "2", "--save-every", "1",
        "--metrics-json", str(tmp_path / "metrics.json"),
    ], check=True, capture_output=True, text=True, timeout=120)

    files = sorted(output.iterdir())
    assert [p.name for p in files] == [
        "updc_f200w_nrca1_iter01_i2d.fits",
        "updc_f200w_nrca1_iter02_i2d.fits",
    ]
    for iteration, path in enumerate(files, 1):
        with fits.open(path, checksum=True) as hdul:
            height, width = hdul["SCI"].data.shape
            assert hdul[0].header["NCOMBINE"] == 2
            assert hdul[0].header["UPDCITER"] == iteration
            assert len(hdul["UPDC_BKG"].data) == 2
            assert len(hdul["UPDC_PSFGEOM"].data) == 2
            assert np.isfinite(hdul["SCI"].data).all()
            assert np.any(hdul["SCI"].data < 0)
            assert hdul["UPDC_BKG"].data["N_NONPOS"].sum() > 0
            assert hdul[0].header["MUFRTOT"] >= 0
            assert hdul["WHT"].header["BUNIT"] == "(Jy/pixel)^-2"
            assert hdul["SCI"].verify_checksum() == 1
    metrics = json.loads((tmp_path / "metrics.json").read_text())
    assert len(metrics["mu_freeze_diagnostics"]) == 4
    assert metrics["performance"]["workers"] == 1
    assert metrics["performance"]["fft_autotune"] is True

    # Final-only mode should retain all solver diagnostics while emitting only
    # the converged FITS product.
    final_output = tmp_path / "final_only_output"
    subprocess.run([
        sys.executable, str(script), "--input-dir", str(images),
        "--psf-dir", str(psfs), "--output-dir", str(final_output),
        "--wcs-mode", "fits", "--psf-crop-size", "0", "--iterations", "2",
        "--fft-plan", "estimate", "--no-fft-autotune",
        "--save-iterations", "final",
    ], check=True, capture_output=True, text=True, timeout=120)
    assert [p.name for p in sorted(final_output.iterdir())] == [
        "updc_f200w_nrca1_iter02_i2d.fits"
    ]

    # Force a non-positive forward mean at the detector footprint edge and
    # verify the C solver records permanent pixel freezing.
    work = tmp_path / "work"
    subprocess.run([
        sys.executable, str(script), "--input-dir", str(images),
        "--psf-dir", str(psfs), "--output-dir", str(tmp_path / "second_output"),
        "--work-dir", str(work), "--wcs-mode", "fits",
        "--psf-crop-size", "0", "--iterations", "1",
    ], check=True, capture_output=True, text=True, timeout=120)
    manifest = work / "manifest.txt"
    lines = []
    for line in manifest.read_text().splitlines():
        fields = line.split("\t")
        fields[6] = "-100"
        lines.append("\t".join(fields))
    manifest.write_text("\n".join(lines) + "\n")
    binary = work / "updc_jwst_c"
    subprocess.run([
        str(binary), str(manifest), str(width), str(height), "64", "64",
        "2", "0", "2", str(work / "forced_iter"), str(work / "forced_wht.bin"),
        "1", "1", "1", "0", "0", "0", "", "1",
    ], check=True, capture_output=True, text=True, timeout=120)
    with (work / "forced_iter_mufreeze.csv").open() as f:
        rows = list(csv.DictReader(f))
    assert any(int(row["n_new_frozen"]) > 0 for row in rows)
    for exposure in ("1", "2"):
        counts = [int(row["n_total_frozen"]) for row in rows if row["exposure"] == exposure]
        new = [int(row["n_new_frozen"]) for row in rows if row["exposure"] == exposure]
        assert counts[0] == new[0]
        assert counts[1] == counts[0] + new[1]
