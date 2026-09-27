# Development guide

## Local checks

```bash
python -m pip install -e '.[test]'
pytest -q
updc doctor
```

Build the C core directly:

```bash
cmake -S updc/c_core -B build -DUPDC_NATIVE_ARCH=OFF
cmake --build build --config Release
ctest --test-dir build --output-on-failure
```

## Change boundaries

- `updc/c_core/`: frozen numerical code; protected by hashes.
- `updc/_engine.py`: validated v17 orchestration; avoid refactors unless required for correctness.
- `updc/config.py`, `updc/cli.py`: stable user interface.
- `updc/adapters/`: telescope-specific extension boundary.

## Adding an adapter

Prefer a separate package using the `updc.instrument_adapters` entry-point group. Keep telescope SDK dependencies out of base UPDC when possible.

## Regression policy

Every release should run:

- frozen-core integrity test;
- config/schema/CLI tests;
- adapter discovery tests;
- C core self-test on Linux and macOS;
- end-to-end normalized FITS smoke test with FFTW;
- wheel build/install test.

For numerical-core changes, follow `CORE_FREEZE.md` rather than the normal packaging release process.
