# UPDC 1.1.0 Algorithm and Numerical Implementation Manual

**Numerical core:** frozen v17  
**Scientific model:** signed-x, background-aware, multi-exposure multiplicative reconstruction  
**Geometry:** exact overlap of mapped detector-pixel quadrilaterals  
**PSF handling:** per-exposure detector-frame PSF transported through the full exposure WCS/GWCS  
**Precision:** float64 geometry decisions + float32 solver/FFTW3f

This document focuses on the scientific operator and the numerical/computational choices that implement it. Installation, configuration, adapters, and release procedures are covered in `UPDC_COMPLETE_MANUAL_EN.md`.

---

## 1. Reconstruction objective and scope

UPDC reconstructs one signed target image `x` from multiple astronomical exposures that may differ in pointing, geometric distortion, detector sampling, scalar sky background, weights, and PSF.

For exposure `k`, the forward model is

$$
H_k = S_k P_k.
$$

`P_k` performs convolution on the common target grid. `S_k` integrates the convolved target image over the mapped detector pixels of exposure `k`.

The frozen core explicitly handles:

- dithered and rotated exposures;
- detector distortion represented by each exposure WCS/GWCS;
- exact detector-pixel quadrilateral integration on the target plane;
- one detector-frame PSF per exposure;
- signed reconstruction values;
- per-exposure scalar background;
- masked pixels and nonuniform weights;
- persistent protection against non-positive modeled detector means;
- high-performance multi-exposure FFT and sparse-geometry execution.

The validated end-to-end instrument workflow is JWST/NIRCam. Other instruments can use the normalized FITS/adaptor architecture but require their own scientific acceptance.

---

## 2. Detector-pixel forward operator

After PSF convolution, a detector pixel is not evaluated at its center. Its four detector-space corners are mapped to the common target plane and form a quadrilateral.

For detector pixel `p` and target pixel `i`, define

$$
S_{k,pi} = \frac{A_{k,pi}}{J_{k,p}},
$$

where

- $A_{k,pi}$ is the polygon-overlap area between mapped detector pixel `p` and target pixel `i`;
- $J_{k,p}$ is the full mapped area of detector pixel `p`.

Thus

$$
(H_kx)_p = \sum_i S_{k,pi}(P_kx)_i.
$$

The full mapped detector-pixel area remains the denominator even when the quadrilateral intersects the edge of the target grid. The model therefore does not renormalize partial edge coverage.

The geometry is exact for the straight-sided quadrilateral defined by the four mapped detector corners. It is not an analytic integration of continuously curved WCS behavior inside one detector pixel.

---

## 3. Signed-x background-aware multiplicative update

Let

- $y_k=SCI_k$;
- $b_k$ be the scalar background of exposure `k`;
- $p_k\ge0$ be the pedestal;
- $W_k$ be the detector weight map.

UPDC uses

$$
y_{\mathrm{eff},k}=y_k+p_k,
$$

$$
\mu_k=H_kx+b_k+p_k.
$$

For an unfrozen detector pixel with $\mu_k>0$,

$$
r_k=\frac{y_{\mathrm{eff},k}}{\mu_k}.
$$

The update is

$$
D=\sum_k H_k^TW_k,
$$

$$
x^{(n+1)}=x^{(n)}\odot
\frac{\sum_k H_k^T(W_k\odot r_k^{(n)})}{D}.
$$

`D` is computed once before the iterative loop.

The reconstruction is signed: negative target pixels are preserved. This distinguishes the implementation from a positivity-constrained Richardson-Lucy interpretation.

---

## 4. Pedestal and FreezeMu

### 4.1 Pedestal

The pedestal is a deterministic constant added to both observation and model:

$$
SCI+p,
\qquad
Hx+b+p.
$$

It does not alter the detector weights and is not a simulated Poisson draw.

The automatic pedestal ensures that all usable `SCI + pedestal` values are positive. A configurable RMS-based floor can increase the pedestal if requested.

### 4.2 FreezeMu

If a detector pixel reaches

$$
\mu_k\le0,
$$

its ratio is permanently fixed to

$$
r_k=1.
$$

The pixel remains in the weighted scatter; it is not discarded. FreezeMu is persistent across all subsequent iterations.

The implementation records newly frozen counts, total frozen counts, and the minimum raw modeled mean for numerical diagnostics.

---

## 5. Background and weights

Each exposure has an independently estimated scalar sky background. Robust asymmetric clipping suppresses positive source contamination.

The main validated weight choice is

$$
W=\frac{1}{VAR_{RNOISE}+\mathrm{median}_{sky}(VAR_{POISSON})}.
$$

Other supported modes are full pixelwise variance and unit-mask weighting.

DQ masking is applied before the solve according to the instrument-appropriate `bad_dq_mask`.

---

## 6. Target WCS and distortion treatment

The common output grid is a new linear TAN WCS. It does not repeat detector SIP/high-order distortion in the output WCS.

Instead, detector distortion is represented by mapping each exposure's detector geometry through its own WCS/GWCS into the common target grid. This separates the source detector geometry from the target coordinate system.

This design is used for both science detector pixels and detector-frame PSF transport.

---

## 7. Detector-frame PSF transport

For a detector-frame PSF with width `PW`, height `PH`, oversampling `s`, and detector anchor $(u_0,v_0)$,

$$
u_j=u_0+\frac{j-PW/2}{s},
$$

$$
v_i=v_0+\frac{i-PH/2}{s}.
$$

Every PSF-sample corner is mapped through the exposure's full detector-to-sky transformation and then into the target WCS. The mapped detector-anchor position is subtracted so the convolution kernel does not contain the exposure's dither translation.

This preserves exposure-dependent:

- rotation;
- scale;
- shear;
- detector distortion over the PSF support.

Exposure 1 follows the same procedure; there is no identity shortcut.

Flux is pushed forward by polygon overlap:

$$
f_{out,j} \mathrel{+}=f_{in}\frac{A_{overlap,j}}{A_{quad}}.
$$

A centered odd-square support crop can be applied after the full transport. The retained kernel is then normalized before FFT spectrum generation.

---

## 8. Adjoint operator and sparse LinkMaps

The transpose of detector integration is implemented by scattering detector values through the same stored overlap fractions. Each exposure's sparse overlap structure is precomputed as a LinkMap.

LinkMap v2 selects between two target-index formats:

- `DELTA16`: `base_target[p] + uint16 delta[z]` with float32 overlap fraction;
- `ABS32`: full `uint32 target[z]` with float32 fraction.

DELTA16 is used only when safe for that exposure. Otherwise the core automatically uses ABS32. The choice affects memory traffic, not the mathematical coefficients.

---

## 9. FFT implementation

PSF convolution is performed with single-precision FFTW3f on a padded grid.

One reconstruction FFT is reused for all exposure PSFs. For each exposure, multiplication by its PSF spectrum followed by inverse FFT gives the convolved target image.

The adjoint uses the conjugate PSF spectrum.

---

## 10. Frequency-domain adjoint reduction

A direct implementation could compute one adjoint inverse FFT per exposure. UPDC uses linearity:

$$
\sum_k P_k^Ts_k
=F^{-1}\left[\sum_k P_k^*F(s_k)\right].
$$

Each exposure produces its detector-scatter spectrum; the PSF-conjugated spectra are summed, and one final inverse transform produces the total correction.

For `N` exposures, one iteration therefore requires:

- 1 reconstruction forward FFT;
- `N` forward-model inverse FFTs;
- `N` detector-scatter forward FFTs;
- 1 combined adjoint inverse FFT.

Hence

$$
N_{FFT}=2N+2.
$$

For 9 exposures this is 20 FFTs per iteration, or 380 transforms over 19 iterations.

The sensitivity denominator uses the same combined-adjoint idea during precomputation.

---

## 11. Fused detector hot path

The main detector loop is fused so one local LinkMap traversal performs:

```text
forward prediction
-> add background/pedestal
-> FreezeMu decision
-> ratio
-> multiply by detector weight
-> immediate sparse scatter
```

This avoids writing a large intermediate detector-state array and avoids rereading the LinkMap in separate forward, FreezeMu, and backprojection phases.

---

## 12. Parallel execution

UPDC separates parallelism into three controls:

- `exposure_workers`
- `fft_workers`
- `preprocess_workers`

Exposure workers persist through the solve rather than being created for every phase.

Nested FFT parallelism can exceed the number of physical cores if both `fft_workers` and `exposure_workers` are large. Performance must therefore be measured on the target machine. Large FFTs frequently become memory-bandwidth or NUMA limited before all cores are computationally saturated.

The AMD EPYC 7763 acceptance run used 9 exposure workers, 9 FFT threads, and 9 preprocessing workers, but that setting is a machine-specific reference, not a universal recommendation.

---

## 13. Mixed precision and reproducibility

Geometry decisions use float64. The numerical solve uses float32, including FFTW3f.

Because floating-point addition is not associative, changing FFT padding, FFTW plans, thread scheduling, or reduction order can change low-order bits. Scientific regression therefore uses tolerances and diagnostics rather than universal cross-platform bitwise identity.

Relevant comparison quantities include:

- maximum and RMS SCI difference;
- reconstruction sum/flux diagnostics;
- negative-pixel count;
- WHT;
- FreezeMu counts;
- minimum raw mu;
- WCS/header consistency.

---

## 14. Precomputation and persistent cache

The following expensive quantities can be reused across runs when the scientific inputs and relevant settings are unchanged:

- normalized input working planes;
- target WCS;
- mapped detector-frame PSFs;
- detector vertices;
- LinkMaps;
- PSF spectra;
- signed initial estimate;
- sensitivity denominator;
- compiled C binary keyed by toolchain signature.

Cache reuse changes performance only. It must not change the operator or the final scientific result beyond expected floating-point ordering differences.

---

## 15. Scientific and implementation limitations

The frozen v17 operator should not be interpreted beyond the following boundaries.

- It consumes externally generated PSFs; it is not a PSF modeling package.
- It supports one spatially invariant PSF per exposure, not a continuously varying PSF field within one exposure.
- Its target WCS is linear TAN; detector high-order distortion is handled through source-to-target mapping.
- Detector integration is exact polygon overlap for mapped corner quadrilaterals, not continuous integration of WCS curvature inside a pixel.
- WHT is a coverage/weight product, not a covariance-aware inverse variance after deconvolution.
- JWST/NIRCam has real-data end-to-end acceptance; other instruments require their own adapter validation and science regression.

---

## 16. Frozen-core policy

The following C files are frozen for the 1.1.x release line:

- `updc_core.c`
- `overlap.c`
- `overlap_geometry.h`
- `psf_projector.c`

Their SHA-256 values are enforced by tests and `updc doctor`.

Changes to CLI, adapters, cache lifecycle, build-system probes, packaging, or documentation must not silently modify these files. A numerical-core modification requires a new core identifier and an explicit scientific validation campaign.

---

## 17. Publication citation

UPDC is distributed under the MIT License. For publications using UPDC code, its algorithmic implementation, or scientific results produced with UPDC, the UPDC project requires citation of both:

1. Wang, L., Li, G., & Kang, X. 2022, *Towards super-resolution via iterative multi-exposure coaddition*, **MNRAS, 517, 787**, doi: `10.1093/mnras/stac2664`.
2. Wang, L., Shan, H., Nie, L., et al. 2025, *UPdec-Webb: A Data Set for Coaddition of JWST NIRCam Images*, **ApJS, 276, 36**, doi: `10.3847/1538-4365/ad9566`.

See `CITATION.cff`, `CITATION.md`, and `NOTICE` for machine-readable metadata and the project publication-citation policy.
