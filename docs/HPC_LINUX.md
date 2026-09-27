# Linux / HPC build and performance guide

## Recommended AMD EPYC toolchain

For AMD EPYC Zen 3 systems, GCC 13+ with `-O3 -march=native` and FFTW3f is a suitable production baseline. On the validated EPYC 7763 host, GCC 13.3 resolves `-march=native` to `znver3`.

## C11 atomic header shadowing

Some HPC compiler modules export C++ header directories through `C_INCLUDE_PATH`, e.g. `.../include/c++/13`. When a C compiler processes `<stdatomic.h>`, this can shadow GCC's real C11 header and produce errors such as:

```text
unknown type name 'atomic_uint_fast64_t'
memory_order_relaxed undeclared
```

UPDC 1.1.0 sanitizes the C build environment by removing C++ include directories from `CPATH`/`C_INCLUDE_PATH` and dropping `CPLUS_INCLUDE_PATH` for C compilation only. Other C include directories are preserved.

Diagnose manually:

```bash
echo '#include <stdatomic.h>' | gcc -std=c11 -H -E -x c - >/dev/null
gcc -print-file-name=include/stdatomic.h
```

A minimal probe should compile and return zero:

```c
#include <stdint.h>
#include <stdatomic.h>
static atomic_uint_fast64_t x;
int main(void) {
    atomic_init(&x, 0);
    atomic_fetch_add_explicit(&x, 1, memory_order_relaxed);
    return atomic_load_explicit(&x, memory_order_relaxed) == 1 ? 0 : 1;
}
```

## `updc doctor`

Run after loading modules/activating the environment:

```bash
updc doctor --cc gcc
```

It reports compiler identity, detected C++ header shadowing, C11 atomics, native CPU flag support, FFTW3f discovery, frozen-core integrity, cache ABIs, and installed adapters.

## FFTW discovery

Discovery order:

1. `pkg-config fftw3f` + optional `fftw3f_threads` metadata;
2. `UPDC_FFTW_INCLUDE` / `UPDC_FFTW_LIB`;
3. Conda/micromamba prefix;
4. common `/usr` and `/usr/local` locations.

A real compile/link probe verifies `fftwf_init_threads()` before the core is built.

## Build cache safety

Compiled C binaries are reused only when the signature matches source hash, compiler path/version, flags, FFTW metadata/link flags, OS, and architecture.

## Threading

`fft_workers * exposure_workers` is a nominal nested maximum. Benchmark rather than assuming all physical cores should be used. On a 64-core EPYC 7763 acceptance host with 9 exposures, a 9x9 topology gave the best measured warm-cache wall time in the tested configurations, but a scheduler-limited job may prefer 7x9 to keep the nominal maximum at 63.
