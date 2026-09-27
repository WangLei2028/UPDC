# Configuration schema v1

UPDC production runs use strict TOML configuration. Unknown top-level sections and unknown keys are rejected.

Every file begins with:

```toml
[updc]
schema_version = 1

[instrument]
adapter = "jwst"
```

The machine-readable JSON Schema is packaged as `updc/schemas/updc.schema.json`.

## Sections

### `[instrument]`

- `adapter`: built-in `jwst` or `fits`, or a third-party entry-point name.
- `options`: optional adapter-specific table. The frozen engine does not interpret it; adapters do.

### `[input]`

Input discovery, WCS mode, weighting, DQ mask, and pedestal settings. One of `input_dir`, `manifest`, or `jwst_dir` must be supplied.

### `[psf]`

PSF discovery, coordinate frame, sampling, storage orientation, crop size, and optional detector anchor.

### `[solver]`

- `iterations`
- `upsample`
- `save_every`: `0` means final iteration only; otherwise save every Nth iteration plus the final iteration.

### `[performance]`

Independent thread controls for FFTW, iterative exposure parallelism, and preprocessing. FFT padding autotuning and its cache are also controlled here.

### `[cache]`

Persistent cache enable/disable and optional locations. The default path is deliberately compatible with the validated v17 cache layout.

### `[output]`

Output directory/stem, optional metrics JSON, and temporary-work retention.

### `[build]`

C compiler or prebuilt C-core/projector paths.

## Relative paths

All relative filesystem paths are resolved relative to the TOML file, not the current shell directory.

## Command-line override

Use a typed TOML scalar expression:

```bash
updc updc.toml --set solver.iterations=39
updc updc.toml --set performance.fft_workers=4
updc updc.toml --set cache.enabled=false
```

The override key must already be part of schema v1.

## Schema evolution policy

A backward-incompatible configuration change increments `updc.schema_version`. Schema v1 files will never be silently reinterpreted under a different meaning. A future release may provide an explicit migration command, but it must not mutate production configs without user action.

## Cache lifecycle policy (1.0.2+)

The stable schema v1 accepts two optional cache lifecycle controls:

```toml
[cache]
enabled = true
reset = "none"          # none | current | all
prune_after_days = 0    # 0 disables automatic pruning
```

`reset = "current"` deletes only the dataset-keyed preprocessing cache for the current dataset before rebuilding it. `reset = "all"` deletes all UPDC dataset-cache directories under `cache_dir`; unrelated files/directories and the global FFTW wisdom file are left untouched. Because reset is evaluated every run, return it to `"none"` after the intended rebuild. The key hashes input/PSF paths, sizes, modification times, and processing settings; it does not hash every input byte.

`prune_after_days > 0` removes dataset caches whose last-use marker is older than that age. Pruning occurs when UPDC next starts; UPDC does not run a background daemon. Existing pre-1.1.0 caches without a marker fall back to their preprocessing-marker/directory modification time.

Process RAM, mmap mappings, and FFTW workspaces are released when the UPDC process exits. Linux may retain file pages in the kernel page cache, but those pages are reclaimable automatically under memory pressure. Persistent files under `~/.cache/updc-jwst/data` remain on disk until reset/pruning/manual deletion.
