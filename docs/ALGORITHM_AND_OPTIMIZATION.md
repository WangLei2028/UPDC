# Algorithm and optimization summary

## Forward model

For exposure `k`, the frozen v17 core uses `H_k = S_k P_k`, where `P_k` is target-grid PSF convolution and `S_k` is exact detector-pixel area integration derived from mapped detector-pixel quadrilaterals.

The driver does not pre-subtract the scalar background from Level-2 SCI. Instead:

```text
y_eff  = SCI + pedestal
raw_mu = H_k x + background + pedestal
ratio  = y_eff/raw_mu while raw_mu > 0
ratio  = 1 permanently after a pixel first reaches raw_mu <= 0
x <- x * sum H_k^T(W_k ratio) / sum H_k^T W_k
```

The reconstruction `x` remains signed; negative target pixels are not clipped.

## Cases explicitly handled

- positive and negative reconstruction pixels;
- valid detector pixels with non-positive raw SCI through a deterministic positivity pedestal;
- persistent invalid forward means through FreezeMu;
- per-exposure scalar background and weight planes;
- detector distortion through exact mapped-pixel polygons;
- detector-frame oversampled PSFs transformed through full per-exposure GWCS;
- exposure-to-exposure rotation/scale/shear/distortion differences;
- arbitrary sparse detector-to-target stencils via generic LinkMaps;
- LinkMap spans that do not fit 16 bits through automatic ABS32 fallback;
- multiple telescope formats through adapters, while preserving one frozen numerical operator.

The validated JWST workflow assumes one spatially invariant PSF kernel per exposure. Intra-exposure spatially varying PSFs are not part of the frozen v17 operator.

## Main computational optimizations

1. FFTW3f convolution on an autotuned padded grid.
2. `2N+2` FFT transforms/iteration through frequency-domain adjoint reduction.
3. Same one-inverse-FFT reduction for the precomputed sensitivity denominator.
4. Fused prediction -> FreezeMu -> weighted scatter detector loop.
5. Persistent exposure worker pool.
6. Independent FFT/exposure/preprocessing worker counts.
7. Adaptive LinkMap v2 (`DELTA16` or `ABS32`).
8. Persistent PSF-spectrum, LinkMap, initial-estimate, denominator, FFT-autotune, and geometry caches.
9. FFTW wisdom persistence.
10. Per-run iteration staging outside persistent cache.
11. Mixed precision: float64 geometry topology, float32 area and numerical solver.
12. Toolchain-aware C binary cache so compiler/flags/FFTW changes trigger rebuilds.

See `UPDC_ALGORITHM_MANUAL_EN.md` for the focused algorithm description and `UPDC_COMPLETE_MANUAL_EN.md` for the full scientific and engineering manual.
