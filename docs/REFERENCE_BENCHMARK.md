# Reference acceptance benchmarks

These measurements are regression references, not performance guarantees.

## Linux / AMD EPYC 7763 acceptance

Validated host/toolchain:

- AMD EPYC 7763 (Zen 3)
- GCC 13.3.0; `-march=native` resolves to `znver3`
- FFTW3f 3.3.5
- 64-core server, 256 GiB RAM

Dataset/configuration:

- 9 x 2048 x 2048 exposures
- target 4565 x 4530
- mapped/cropped PSF 1431 x 1431
- 19 iterations
- FFT autotune selected 6000 x 6000 for `fft_workers=9`
- `fft_workers=9`
- `exposure_workers=9`
- `preprocess_workers=9`
- adaptive LinkMap encoding: DELTA16 for all 9 exposures

Warm-cache acceptance:

- solver wall: **21.207 s**
- total wall: **23.196 s**
- peak process RSS: **9.189 GiB**
- LinkMaps mapped: **2.190 GiB**
- iteration transforms: **190 forward + 190 inverse = 380**

A preceding rebuild/autotune run selected the same 6000 x 6000 FFT shape and demonstrated successful GCC/C11-atomics compilation after removing HPC C++ header shadowing.

## Earlier macOS engineering reference

The v17 development sequence also used a 9-exposure JWST/NIRCam dataset on macOS. A representative warm-cache 2x2 topology measured about 33.4 s total / 32.1 s solver with 5.3 GiB peak RSS. This remains useful for cross-host regression context but should not be compared directly to the EPYC server as a hardware-normalized benchmark.

## Regression interpretation

A packaging or adapter change should not cause major changes in:

- FFT transform count for the same numerical configuration;
- LinkMap encoding/size when geometry is unchanged;
- reconstruction/WHT/FreezeMu diagnostics beyond expected float32 reduction-order effects;
- warm-cache wall time beyond normal system/scheduler variability.
