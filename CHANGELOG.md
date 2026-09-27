# Changelog

## 1.1.0 - frozen-core GitHub release

Numerical core remains byte-for-byte frozen v17.

### Post-review corrections

- Add English technical documentation based on a source-level review, including complete and algorithm-focused manuals.
- Preserve explicit `--set` values when built-in adapters apply defaults.
- Disable the JWST-specific 1431-pixel PSF crop by default for generic FITS runs; update the generic example.
- Create missing parent directories for output paths and derive WHT `BUNIT` from input SCI units.
- Reject invalid Boolean/numeric TOML values during stable configuration validation.
- Serialize non-finite optional diagnostics as JSON `null` rather than nonstandard `NaN`/`Infinity` tokens.
- Restore the small-PSF end-to-end smoke test by explicitly disabling crop in the legacy CLI workflow.

### Open-source publication metadata

- Replace the pre-publication restrictive license placeholder with the standard MIT License.
- Add `CITATION.cff`, `CITATION.md`, and `NOTICE` with the UPDC scientific publication citation policy.
- Publications using UPDC code, its algorithmic implementation, or UPDC-derived scientific results are required by the project citation policy to cite both the 2022 MNRAS method paper and the 2025 ApJS UPdec-Webb paper.
- Add GitHub issue/PR templates and a security-reporting policy.
- Remove macOS/Pytest development artefacts from the release tree.

### Packaging / operations

- Complete GitHub-ready English documentation, including complete and algorithm-focused Markdown/PDF manuals.
- Stable TOML schema v1 and instrument adapter API retained.
- Cache lifecycle controls retained and documented:
  - `cache.reset = "none" | "current" | "all"`
  - `cache.prune_after_days = N`
- Runtime iteration staging remains outside persistent preprocessing cache.

### Linux / HPC portability

- Filter C++ include directories from the C compilation environment to prevent HPC `C_INCLUDE_PATH` from shadowing GCC's C11 `<stdatomic.h>`.
- Add independent C11 atomics compile probe.
- Probe native CPU flags independently; unrelated source/link failures are no longer mislabeled as `-march=native` failures.
- FFTW discovery now prefers pkg-config but can fall back to explicit/Conda/common filesystem locations and verifies thread support with a real compile/link probe.
- Compiled C binary cache signature now includes compiler path/version, compile flags, FFTW metadata/link flags, OS/architecture, and frozen source hash.
- `updc doctor --cc <compiler>` reports build-environment diagnostics.

### Reference acceptance

- Linux AMD EPYC 7763 / GCC 13.3 / FFTW3f 3.3.5 acceptance run documented.
- 9 exposures, 19 iterations, 9x9x9 workers, 6000x6000 FFT: warm-cache total 23.196 s, solver 21.207 s, peak RSS 9.189 GiB.

## 1.0.2

- Added safe persistent-cache reset and stale-cache pruning policy.

## 1.0.1

- Fixed stable CLI -> frozen-engine config handoff so option values such as `--input-dir` are not misinterpreted as legacy TOML paths.

## 1.0.0

- First production packaging of frozen v17 numerical core with stable CLI/TOML schema, instrument adapters, cache ABI, tests, and documentation.
