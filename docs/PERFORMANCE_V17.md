# Performance notes - frozen v17 core

The frozen v17 core keeps the v16 `2N+2` FFT-per-iteration algorithm and adds adaptive LinkMap compression and FFT-autotune caching.

## Per-iteration transform count

For `N` exposures:

- 1 reconstruction forward FFT;
- N forward-model inverse FFTs;
- N scatter forward FFTs;
- 1 combined adjoint inverse FFT.

Total: `2N+2` transforms/iteration. For 9 exposures this is 20; for 19 iterations it is 380.

## Fused detector path

Prediction, FreezeMu, weighted ratio, and scatter reuse the same LinkMap stencil in one pass, reducing sparse memory traffic.

## Frequency-domain adjoint reduction

Instead of one adjoint inverse FFT per exposure, all `conj(PSF_k) * FFT(scatter_k)` spectra are accumulated first and transformed back once. The sensitivity denominator uses the same linearity.

## Persistent exposure workers

Iteration pthreads are created once and reused. `exposure_workers` and `fft_workers` form a two-level parallel topology; benchmark combinations on each host.

## FFT autotune result cache

Candidate padding shapes are benchmarked and the selected shape is cached by minimum FFT shape, FFT worker count, host/OS/architecture, and FFTW runtime signature. `retune_fft=true` forces a new measurement for the current topology.

## Adaptive LinkMap v2

Each exposure selects one of:

- `DELTA16`: uint32 CSR offsets, uint32 per-pixel base target, uint16 link deltas, float32 fractions;
- `ABS32`: uint32 CSR offsets, uint32 absolute link targets, float32 fractions.

DELTA16 is used only if all actual detector-pixel target-index spans fit `UINT16_MAX`; otherwise that exposure automatically falls back to ABS32. No telescope or regular-stencil assumption enters the numerical path.

## Runtime staging

Iteration raw binaries and diagnostic CSVs are per-run temporary files rather than preprocessing cache entries.

## EPYC 7763 acceptance

On the validated 64-core AMD EPYC 7763 host with GCC 13.3 and FFTW3f 3.3.5, the 9-exposure JWST acceptance dataset with a 9x9x9 worker configuration reached 21.207 s solver / 23.196 s total warm-cache wall time, 9.189 GiB peak RSS, 2.190 GiB mapped LinkMaps, and the expected 380 iteration transforms.
