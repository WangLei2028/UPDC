# UPDC 1.1.0 release notes

UPDC 1.1.0 is the GitHub-ready release of the frozen v17 numerical core.

## Scientific core

- Signed-x, background-aware multi-exposure reconstruction.
- Exact mapped detector-pixel polygon integration on a common target grid.
- Per-exposure detector-frame PSF transport through full WCS/GWCS geometry.
- Frozen v17 C numerical core protected by SHA-256 regression checks.

## Performance / engineering

- Fused forward -> FreezeMu -> weighted-scatter detector loop.
- Frequency-domain adjoint reduction: `2N+2` FFT transforms per iteration.
- Persistent exposure workers and independent FFT/exposure/preprocess worker controls.
- Adaptive LinkMap v2 (`DELTA16` with safe `ABS32` fallback).
- FFT-shape autotuning, FFTW wisdom, preprocessing cache, cache reset/pruning.
- Linux/HPC build probes for C11 atomics, FFTW, native CPU flags, and C++ header shadowing.

## Packaging

- Stable TOML schema v1 and extensible instrument-adapter entry point.
- Built-in `jwst` and normalized `fits` adapters.
- MIT License.
- `CITATION.cff`, `CITATION.md`, and `NOTICE` define the scientific publication citation policy.
- English-only documentation set: complete Markdown/PDF manual, focused algorithm Markdown/PDF manual, and technical manual.
- GitHub Actions, issue templates, PR template, security policy, tests, and release manifests.

## Publication citations

Publications using UPDC code, its algorithmic implementation, or UPDC-derived scientific results are required by the UPDC project citation policy to cite both:

1. Wang, L., Li, G., & Kang, X. 2022, *Towards super-resolution via iterative multi-exposure coaddition*, MNRAS, 517, 787, doi:10.1093/mnras/stac2664.
2. Wang, L., Shan, H., Nie, L., et al. 2025, *UPdec-Webb: A Data Set for Coaddition of JWST NIRCam Images*, ApJS, 276, 36, doi:10.3847/1538-4365/ad9566.

The citation policy is separate from the MIT software license.
