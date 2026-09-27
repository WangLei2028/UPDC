# Validation — v17

## C core build/self-test

The v17 C core was compiled with `-std=c11 -Wall -Wextra -Werror` and linked against a real FFTW3f runtime. `updc_jwst_c --self-test` passed.

## Adaptive LinkMap branch tests

Two synthetic geometries were executed through the complete solver:

1. A small-span footprint selected `encoding=DELTA16` and completed normally.
2. A deliberately extreme target with width 70000 made one detector pixel span target indices separated by 70000. The same core automatically selected `encoding=ABS32` and completed normally.

This validates the telescope/data-independent fallback rather than assuming JWST-like geometry.

## regular-3x3 statistic

A synthetic detector pixel mapped exactly onto a 3x3 target block. The diagnostic reported `9-link=1 regular-3x3=1 (100%)`. This statistic does not alter LinkMap encoding or solver arithmetic.

## FFT autotune cache

A first persistent-cache run benchmarked candidate FFT shapes and saved the selected shape. The second identical run reported `FFT autotune persistent-cache hit`. A third run with retuning enabled benchmarked candidates again.

## v16 numerical inheritance

v17 preserves the v16 fused detector pipeline and frequency-domain adjoint reduction. v16 serial/parallel fixtures were bitwise identical for the tested worker topology, and the v16 operator differed from the pre-frequency-reduction implementation only at expected float32 summation-order level.

## Production acceptance recommendation

For a new telescope/data adapter, validate:

- target-grid geometry and detector vertex maps;
- PSF coordinate convention and normalization;
- `DELTA16` versus `ABS32` encoding logs;
- total flux, max/RMS SCI differences, WHT, and FreezeMu diagnostics against a trusted run;
- multiple `fft_workers × exposure_workers` topologies.
