#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
from astropy.io import fits


def stats(a: np.ndarray, b: np.ndarray) -> dict[str, float | int]:
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    if a.shape != b.shape:
        raise ValueError(f"shape mismatch: {a.shape} vs {b.shape}")
    d = b - a
    scale = max(float(np.nanmax(np.abs(a))), 1.0e-30)
    mask = np.abs(a) > 1.0e-8 * scale
    rel = np.abs(d[mask] / a[mask]) if np.any(mask) else np.array([0.0])
    return {
        "max_abs": float(np.nanmax(np.abs(d))),
        "rms": float(np.sqrt(np.nanmean(d * d))),
        "max_rel_masked": float(np.nanmax(rel)),
        "ref_sum": float(np.nansum(a)),
        "new_sum": float(np.nansum(b)),
        "ref_negative": int(np.sum(a < 0)),
        "new_negative": int(np.sum(b < 0)),
        "nonfinite_ref": int(np.sum(~np.isfinite(a))),
        "nonfinite_new": int(np.sum(~np.isfinite(b))),
    }


def main() -> None:
    p = argparse.ArgumentParser(description="Compare UPDC reference and candidate FITS SCI/WHT arrays")
    p.add_argument("reference", type=Path)
    p.add_argument("new", type=Path)
    args = p.parse_args()
    with fits.open(args.reference, memmap=True) as ref, fits.open(args.new, memmap=True) as new:
        for ext in ("SCI", "WHT"):
            print(f"[{ext}]")
            for key, value in stats(ref[ext].data, new[ext].data).items():
                print(f"  {key:18s} {value}")
        for key in ("MUFRTOT", "MUFRNEW", "UPDCITER"):
            print(f"[HEADER] {key}: ref={ref[0].header.get(key)} new={new[0].header.get(key)}")


if __name__ == "__main__":
    main()
