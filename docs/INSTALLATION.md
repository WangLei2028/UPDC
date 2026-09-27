# Installation

UPDC requires Python >=3.10 and the single-precision FFTW3 development interface (`fftw3f` plus thread wrapper).

## macOS (Homebrew)

```bash
brew install fftw pkg-config cmake
python -m pip install updc_reconstruction-1.1.0-py3-none-any.whl
```

## Ubuntu / Debian

```bash
sudo apt-get update
sudo apt-get install -y libfftw3-dev pkg-config cmake build-essential
python -m pip install updc_reconstruction-1.1.0-py3-none-any.whl
```

## Conda / micromamba / HPC

```bash
conda install -c conda-forge fftw pkg-config
python -m pip install updc_reconstruction-1.1.0-py3-none-any.whl
```

For JWST GWCS support:

```bash
python -m pip install 'updc-reconstruction[jwst]'
```

## Verify the toolchain

```bash
updc doctor
# or, on a module-based HPC system:
updc doctor --cc gcc
```

The doctor reports frozen-core hashes, cache ABIs, compiler identity, C11 atomics, native CPU flag support, FFTW3f discovery, and adapters.

## HPC header-shadowing protection

Some compiler modules inject C++ directories such as `.../include/c++/13` into `C_INCLUDE_PATH`. That can cause a C compilation to read libstdc++'s compatibility `stdatomic.h` instead of GCC's C11 header. UPDC 1.1.0 filters these C++ paths for its C build automatically.

Manual diagnosis:

```bash
echo '#include <stdatomic.h>' | gcc -std=c11 -H -E -x c - >/dev/null
gcc -print-file-name=include/stdatomic.h
```

If needed, a shell-level workaround is:

```bash
unset C_INCLUDE_PATH
unset CPLUS_INCLUDE_PATH
```

See `docs/HPC_LINUX.md`.

## FFTW discovery fallback

UPDC tries pkg-config first. If metadata is absent it can use:

```bash
export UPDC_FFTW_INCLUDE=/path/to/include
export UPDC_FFTW_LIB=/path/to/lib
```

It also checks Conda/micromamba and common system locations, then compile/links a small FFTW-thread probe.

## Prebuilt C binaries

Sites that forbid runtime compilation can set `build.c_core_binary` and `build.psf_projector_binary` in TOML. They should be built from the same frozen core version and scientifically qualified for the deployment environment.
