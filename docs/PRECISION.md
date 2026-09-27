# Precision policy

The release uses deliberate mixed precision.

## f64

Only geometric coordinates and operations needed to determine polygon topology stay in double precision: WCS-produced vertices, clipping tests, and segment/boundary intersections.

## f32

Polygon area and the numerical solver use single precision. Area is not evaluated from absolute ~thousands-pixel coordinates directly. Each clipped polygon is translated to its first vertex in f64, those small local coordinate differences are cast to f32, and the shoelace cross products are accumulated with `fmaf`.

The rest of the scientific path is f32: overlap area/fraction, SCI/WHT, PSF normalization, FFTW3f, forward prediction, raw_mu, FreezeMu ratio, scatter, sensitivity, multiplicative correction, and estimate.

PSF normalization uses compensated single-precision summation to reduce long-sum error without reintroducing f64 into that path.

## Expected differences

Changing FFTW plan, FFT padding shape, worker count, or operation ordering can produce normal float32-level differences. Validation should therefore compare max/RMS differences, total flux, negative-pixel count, WHT, and FreezeMu diagnostics rather than require bitwise identity across all configurations.
