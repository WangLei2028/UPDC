# Detector-frame HybPSF -> UPDC target-frame PSF

## Assumption

The input HybPSF is an oversampled, spatially invariant PSF expressed in detector coordinates. In the intended NIRCam use case it is sampled at 2x native detector sampling and centered on the detector center.

## Coordinate construction

For an input PSF of width `PW`, height `PH`, oversampling `s`, and detector anchor `(u0,v0)`, PSF pixel-boundary coordinates are

```text
u_j = u0 + (j - PW/2) / s,  j=0..PW
v_i = v0 + (i - PH/2) / s,  i=0..PH
```

Thus an odd-sized PSF has its central pixel centered exactly on `(u0,v0)`, and each PSF pixel has native-detector width `1/s`.

Every corner is evaluated with the exposure's authoritative detector-to-sky GWCS and then with the inverse of the common target WCS. The target WCS is a fresh linear TAN grid containing only first-order scale/orientation terms.

The target coordinate of `(u0,v0)` is subtracted from every mapped corner. This removes pointing/dither translation from the convolution kernel while retaining the PSF shape transformation.

## Flux push-forward

Each mapped input PSF pixel is a quadrilateral in target coordinates. Its input flux is distributed to overlapping target pixels by

```text
fraction = overlap_area / full_mapped_quadrilateral_area
output[target] += input_flux * fraction
```

Coordinates/intersections use float64; polygon areas and flux arithmetic use float32. The output support is automatically chosen as an odd, symmetric grid around the translated PSF center.

This procedure is non-negative for non-negative input PSFs and does not use Lanczos, bilinear, spline, or any other interpolation kernel.

## Exposure 1

Exposure 1 is not an identity special case. Its detector-frame PSF is also mapped through its full GWCS to remove detector field distortion before FFT convolution on the distortion-free target grid.

## Source and target distortion roles

- Source detector -> sky: full JWST GWCS, including higher-order field distortion.
- Target sky -> target pixels: linear TAN WCS only (CRVAL/CRPIX plus first-order CD/PC-equivalent matrix).

Therefore field distortion is corrected during the detector-to-target push-forward but is not encoded again in the target WCS.

## v14 centered support crop

After flux-conserving polygon mapping, detector-frame HybPSFs are center-cropped to a common odd square before the solver FFT. The default is 1431x1431 (`--psf-crop-size 1431`); `0` disables the crop. Cropping occurs only after the full GWCS mapping, so the distortion correction itself is unchanged. The solver subsequently normalizes the retained kernel to unit sum. `UPDC_PSFGEOM` records the mapped dimensions, solver dimensions, and retained/discarded flux fraction.
