# UPDC 1.1.0

**English-only documentation edition.** All user-facing manuals, release notes, GitHub templates, and documentation in this package are provided in English.

UPDC is a high-performance, signed-x, background-aware astronomical image reconstruction package. The **v17 numerical C core is frozen byte-for-byte**; this release completes the production packaging around it: stable TOML schema, instrument adapters, versioned cache ABI, cache lifecycle controls, portable Linux/HPC build probes, tests, and a complete manual.

## Citation for publications

UPDC is released under the MIT License. If UPDC code, its algorithmic implementation, or scientific results produced with UPDC are used in a journal paper, conference proceeding, thesis/dissertation, technical report, data-release paper, or other publication, the UPDC project requires citation of **both** papers:

1. Wang, L., Li, G., & Kang, X. 2022, *Towards super-resolution via iterative multi-exposure coaddition*, **MNRAS, 517, 787**, doi: `10.1093/mnras/stac2664`.
2. Wang, L., Shan, H., Nie, L., et al. 2025, *UPdec-Webb: A Data Set for Coaddition of JWST NIRCam Images*, **ApJS, 276, 36**, doi: `10.3847/1538-4365/ad9566`.

See [`CITATION.cff`](CITATION.cff), [`CITATION.md`](CITATION.md), and [`NOTICE`](NOTICE). The publication citation policy is a scientific-attribution policy separate from the MIT software license.

## Scientific operator

For exposure `k`, UPDC uses

```text
H_k = S_k P_k
```

where `P_k` is target-grid PSF convolution and `S_k` is exact mapped detector-pixel area integration. The Level-2 science plane is not pre-background-subtracted. The frozen solver uses

```text
y_eff  = SCI + pedestal
raw_mu = H_k x + background + pedestal
ratio  = y_eff/raw_mu   while raw_mu > 0
ratio  = 1              permanently after raw_mu <= 0 (FreezeMu)

x <- x * sum H_k^T(W_k * ratio) / sum H_k^T W_k
```

The reconstruction `x` is signed; negative target pixels are not clipped away.

## Geometry and PSF handling

- Detector pixels are integrated as mapped quadrilaterals, not point samples.
- JWST detector-frame PSFs are transported through each exposure's full detector->sky GWCS, then into a fresh linear TAN target grid.
- Dither translation is removed from the convolution kernel while distortion/rotation/scale/shear are retained.
- Exposure 1 is transformed too; there is no identity exception.
- PSF transport uses exact polygon-overlap flux fractions, not Lanczos/bilinear interpolation.
- The validated JWST configuration center-crops the mapped PSF to 1431x1431 after transport; the solver then normalizes the retained kernel.

See the [complete English manual](docs/UPDC_COMPLETE_MANUAL_EN.md), the focused [algorithm and numerical implementation manual](docs/UPDC_ALGORITHM_MANUAL_EN.md), and the shorter [technical manual](docs/UPDC_MANUAL_EN.md). Detailed PSF/GWCS and precision notes are also available in [PSF_GWCS_POLYGON.md](docs/PSF_GWCS_POLYGON.md) and [PRECISION.md](docs/PRECISION.md).

## Major performance design

The frozen v17 implementation includes:

- fused `forward -> FreezeMu -> weighted scatter` detector loop;
- frequency-domain adjoint reduction, giving `2N+2` FFT transforms per iteration for `N` exposures;
- the same reduction for the sensitivity denominator;
- persistent exposure worker pool;
- exposure-level, FFT-level, and preprocessing parallelism;
- adaptive LinkMap v2 (`DELTA16` when safe, `ABS32` fallback otherwise);
- FFT-shape autotuning + persistent result cache + FFTW wisdom;
- dataset-keyed preprocessing cache (paths, sizes, mtimes, and processing settings);
- per-run temporary iteration staging, separate from persistent cache.

## Installation

Python >= 3.10 is required. UPDC runtime-compiles its frozen C core unless prebuilt binaries are configured.

### Ubuntu / Debian

```bash
sudo apt-get update
sudo apt-get install -y libfftw3-dev pkg-config cmake build-essential
python -m pip install updc_reconstruction-1.1.0-py3-none-any.whl
```

### macOS

```bash
brew install fftw pkg-config cmake
python -m pip install updc_reconstruction-1.1.0-py3-none-any.whl
```

### JWST support

```bash
python -m pip install 'updc-reconstruction[jwst]'
```

Then verify the environment:

```bash
updc doctor
```

On HPC systems `updc doctor --cc gcc` also diagnoses C11 atomics, native CPU flags, FFTW discovery, and C++-header shadowing. See [Linux/HPC guide](docs/HPC_LINUX.md).

## Run from one TOML file

```bash
updc init updc.toml
# edit input/PSF/output paths and settings
updc validate updc.toml
updc updc.toml
```

Equivalent:

```bash
python -m updc updc.toml
```

Temporary override without editing the file:

```bash
updc updc.toml --set solver.iterations=39 --set performance.fft_workers=4
```

## Cache reset and automatic pruning

```toml
[cache]
enabled = true
reset = "none"          # none | current | all
prune_after_days = 30   # 0 disables automatic pruning
```

- `current`: delete/rebuild only the current dataset cache. Set it back to `none` after the intended run.
- `all`: delete all UPDC dataset-keyed caches under the cache root; unrelated paths and global FFTW wisdom are preserved.
- `prune_after_days`: stale dataset caches are removed when UPDC next starts; UPDC does not run a background daemon.

Process RAM, mmap mappings, and FFTW workspaces are automatically released at process exit. Linux filesystem page cache is reclaimable by the kernel under memory pressure.

## Portable build behavior

UPDC 1.1.0 hardens runtime compilation without changing the frozen core:

- filters `.../include/c++/...` entries from the C build environment to prevent `stdatomic.h` shadowing on HPC modules;
- independently probes C11 atomics;
- independently probes `-march=native` / `-mcpu=native`, so unrelated compile errors are not misreported as CPU-flag failures;
- finds FFTW via pkg-config first, then explicit/Conda/common filesystem locations;
- compile/link probes FFTW thread support;
- includes compiler path/version, compile flags, FFTW metadata, OS/architecture, and frozen-source hash in the compiled-binary cache signature.

## Instrument adapters

Built-ins:

- `jwst`: JWST Level-2 products using authoritative `meta.wcs`/GWCS when requested;
- `fits`: normalized astronomical FITS products with celestial FITS WCS.

Third-party adapters register through:

```text
updc.instrument_adapters
```

Telescope-specific ingestion belongs in adapters; changing the numerical operator requires a new core freeze.

## Reference acceptance benchmark

A validated Linux/HPC run on AMD EPYC 7763 (Zen 3), GCC 13.3 (`-march=native -> znver3`), FFTW3f 3.3.5:

- 9 x 2048^2 exposures
- target: 4565 x 4530
- PSF: 1431^2
- 19 iterations
- FFT: 6000 x 6000
- `fft_workers=9`, `exposure_workers=9`, `preprocess_workers=9`
- warm-cache solver wall: **21.207 s**
- warm-cache total wall: **23.196 s**
- peak RSS: **9.189 GiB**
- iteration FFT count: **380** (`19 x (2*9+2)`)
- adaptive DELTA16 LinkMaps: **2.190 GiB** mapped

This is a regression/reference measurement, not a runtime guarantee for other hosts or datasets.

## Documentation

- [Complete English manual](docs/UPDC_COMPLETE_MANUAL_EN.md) and [PDF](docs/UPDC_COMPLETE_MANUAL_EN.pdf)
- [Algorithm and numerical implementation manual](docs/UPDC_ALGORITHM_MANUAL_EN.md) and [PDF](docs/UPDC_ALGORITHM_MANUAL_EN.pdf)
- [English technical manual](docs/UPDC_MANUAL_EN.md) and [PDF](docs/UPDC_MANUAL_EN.pdf)
- [Algorithm and optimization summary](docs/ALGORITHM_AND_OPTIMIZATION.md)
- [Installation](docs/INSTALLATION.md)
- [Configuration schema](docs/CONFIGURATION.md)
- [Architecture](docs/ARCHITECTURE.md)
- [Instrument adapters](docs/ADAPTERS.md)
- [Cache ABI](docs/CACHE_ABI.md)
- [Frozen core policy](docs/CORE_FREEZE.md)
- [PSF/GWCS polygon mapping](docs/PSF_GWCS_POLYGON.md)
- [Precision](docs/PRECISION.md)
- [Linux/HPC guide](docs/HPC_LINUX.md)
- [Release checklist](docs/GITHUB_RELEASE.md)
- [License and citation policy](docs/LICENSE_AND_CITATION.md)

## Testing

```bash
python -m pip install '.[test]'
pytest -q
```

C core:

```bash
cmake -S updc/c_core -B build -DUPDC_NATIVE_ARCH=ON
cmake --build build --config Release
ctest --test-dir build --output-on-failure
```

## Frozen numerical core

The following source SHA-256 values are enforced by CI/tests:

```text
updc_core.c          df2cfffce5cf9b81749ae41dff8574c1fd8f2d4d1b409890a1e23062ad5807bb
overlap.c            c32b447dff787f599207ea35b676be6c1cfab9f76d8a4e95966bb9db38b68d71
overlap_geometry.h   f297ccf2c7eff07e5e8f73fc9a901c41d33a59a6b147ed264995ed47ef7d43f8
psf_projector.c      a6f6baba9cb2bba5c9761b1590918a4f45aa82f0791863c88a7cfa4ff3b5bfd1
```

Packaging, CLI, adapter, cache-management, and build-detection changes must not silently modify these files.

## License

UPDC is distributed under the standard [MIT License](LICENSE). Preserve the copyright and permission notice as required by that license. For scholarly or other publication use, also follow the project citation policy in [CITATION.md](CITATION.md) and [NOTICE](NOTICE), which requires citation of both UPDC papers listed above. The citation policy is separate from and does not modify the MIT License.
