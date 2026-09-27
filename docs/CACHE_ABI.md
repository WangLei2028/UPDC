# Cache ABI policy

Persistent cache is a performance feature, not part of the scientific data product. It nevertheless has explicit ABI versions so old bytes are never silently reinterpreted.

## Current ABIs

| Cache | ABI | Identifier |
|---|---:|---|
| Python preprocessing | 1 | `updc-v16-preprocess-1` |
| C LinkMap | 2 | header magic `ULM2` / version 2 |
| FFT autotune | 1 | `fft_autotune_v1_*` |

The Python marker string is intentionally retained from v16/v17 so validated caches remain reusable after packaging.

## LinkMap v2

Each exposure stores a versioned header and independently chooses one of:

- `ABS32`: general 32-bit absolute target indices;
- `DELTA16`: per-detector-pixel 32-bit base target plus 16-bit deltas when every actual span fits.

This is capability-detected per exposure. No telescope, 3×3 stencil, upsample factor, or detector geometry is hard-coded.

## Compatibility rules

1. A reader must validate magic/version/shape/count/encoding before mapping a cache file.
2. A format change increments the corresponding ABI.
3. Unknown/newer ABI is a cache miss, never a best-effort read.
4. Scientific output must not depend on whether a cache was hit; only floating-point reduction-order differences already qualified by the frozen core are permitted.
5. Cache directories may be safely deleted; UPDC rebuilds them.

## Lifecycle metadata

UPDC 1.0.2 adds a `.last_used` marker inside each dataset-keyed cache. This marker is lifecycle metadata only; it is not part of the numerical cache ABI and does not alter cache keys or numerical results. Automatic pruning only considers directories whose names match the 24-hex-character UPDC dataset key format. The key uses path, size, and modification time for each input/PSF rather than hashing entire file contents; force a rebuild if an external operation changes bytes without changing this metadata.
