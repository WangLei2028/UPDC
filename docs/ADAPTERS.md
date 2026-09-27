# Instrument adapters

Instrument adapters isolate telescope-specific ingestion from the frozen reconstruction core.

## Normalized exposure contract

The frozen v17 engine currently consumes same-sized FITS exposures containing:

- `SCI`
- `DQ`
- `VAR_POISSON`
- `VAR_RNOISE`
- `VAR_FLAT`

For the `fits` adapter, the `SCI` header supplies celestial FITS WCS. For `jwst`, `meta.wcs` is used when requested. Each exposure is paired with one 2-D PSF FITS image.

An adapter for another telescope may either point UPDC at already-normalized files or materialize normalized files in its temporary staging directory.

## API

Subclass `updc.adapters.InstrumentAdapter`:

```python
from updc.adapters import InstrumentAdapter, PreparedRun

class MyAdapter(InstrumentAdapter):
    name = "my-scope"
    description = "My telescope adapter"

    def apply_defaults(self, config):
        return config

    def validate(self, config):
        return []

    def prepare(self, config, context):
        # Optional: write normalized FITS files under context.staging_dir
        # and return a modified config that points to them.
        return PreparedRun(config=config)
```

`prepare` is called after stable-schema validation and before the frozen engine starts. The staging directory remains alive for the duration of the run and is deleted afterwards.

## Plugin registration

In a third-party `pyproject.toml`:

```toml
[project.entry-points."updc.instrument_adapters"]
my-scope = "mypackage.adapter:MyAdapter"
```

Then:

```bash
updc adapters
```

will discover it, and configuration may use:

```toml
[instrument]
adapter = "my-scope"
```

## Generalization rule

Adapters may change ingestion and normalization but must not reinterpret the frozen core's numerical semantics. If a telescope needs a genuinely different numerical operator, that is a new core capability and requires a new core freeze/validation, not an adapter shortcut.
