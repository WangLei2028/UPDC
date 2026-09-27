# Numerical core freeze

UPDC 1.0.0 freezes the validated v17 C numerical core.

## Protected files

| File | SHA-256 |
|---|---|
| `updc_core.c` | `df2cfffce5cf9b81749ae41dff8574c1fd8f2d4d1b409890a1e23062ad5807bb` |
| `overlap.c` | `c32b447dff787f599207ea35b676be6c1cfab9f76d8a4e95966bb9db38b68d71` |
| `overlap_geometry.h` | `f297ccf2c7eff07e5e8f73fc9a901c41d33a59a6b147ed264995ed47ef7d43f8` |
| `psf_projector.c` | `a6f6baba9cb2bba5c9761b1590918a4f45aa82f0791863c88a7cfa4ff3b5bfd1` |

These values are enforced by `tests/test_core_freeze.py` and visible through `updc doctor`.

## What “frozen” means

Packaging, documentation, adapters, CLI, and configuration can evolve without changing these files. Any intentional numerical-core modification requires:

1. a new core freeze version;
2. explicit source-hash changes;
3. numerical comparison against the v17 reference implementation;
4. synthetic forward/adjoint tests;
5. real-data acceptance tests;
6. performance regression measurements;
7. cache ABI review.

A maintainer must not simply update the stored hashes to make CI pass.
