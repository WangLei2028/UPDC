# Contributing

The v17 numerical core is frozen. Normal contributions should target adapters, configuration, packaging, tests, or documentation. Do not modify files under `updc/c_core/` unless the change is explicitly being developed and validated as a new numerical-core freeze.

Before submitting changes:

```bash
pytest -q
updc doctor
```

See `docs/DEVELOPMENT.md` and `docs/CORE_FREEZE.md`.
