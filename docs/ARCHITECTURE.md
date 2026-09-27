# Architecture

UPDC 1.0.0 separates the validated reconstruction operator from telescope-specific ingestion.

```text
native telescope products
        │
        ▼
instrument adapter
(validate / normalize / stage)
        │
        ▼
normalized FITS exposure contract + per-exposure PSFs
        │
        ▼
Python orchestration
GWCS/FITS geometry, PSF transport, cache management
        │
        ▼
frozen v17 C numerical core
LinkMaps + FFTW + signed-x UPDC iterations
        │
        ▼
FITS reconstruction products + metrics
```

## Frozen layer

`updc/c_core/` is the frozen v17 numerical implementation. The files `updc_core.c`, `overlap.c`, `overlap_geometry.h`, and `psf_projector.c` are protected by SHA-256 values in `updc/frozen_core.py` and by CI tests.

Numerical changes require an explicit new core freeze version and updated validation; ordinary packaging or adapter work must not modify these files.

## Orchestration layer

`updc/_engine.py` is the v17 Python orchestration implementation with only an invocation-level change allowing an explicit argument vector. It still performs:

- robust background and WHT preparation;
- target WCS construction;
- detector-frame PSF full-WCS polygon transport;
- detector-pixel vertex generation;
- persistent preprocessing/cache orchestration;
- C-core build and invocation;
- FITS product creation and validation.

## Stable facade

`updc/cli.py` and `updc/config.py` form the public stable interface. They validate schema-v1 TOML, resolve an instrument adapter, translate configuration into the frozen engine's CLI contract, and execute it.

## Adapter boundary

Adapters are intentionally outside the numerical core. A future HST, Roman, ground-based, or simulated-data adapter may convert native products into the normalized exposure contract without changing the solver.

## Compatibility layer

`updc_jwst` and `run_updc.py` preserve the v17 flag-based interface. New workflows should use the `updc` command and schema-v1 TOML.
