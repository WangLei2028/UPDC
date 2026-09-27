"""Minimal third-party adapter example; package this in a separate distribution."""
from copy import deepcopy
from updc.adapters import InstrumentAdapter, PreparedRun


class ExampleAdapter(InstrumentAdapter):
    name = "example"
    description = "Example telescope adapter"

    def apply_defaults(self, config):
        cfg = deepcopy(config)
        cfg["input"]["wcs_mode"] = "fits"
        cfg["input"]["pairing"] = "sorted"
        return cfg

    def validate(self, config):
        return []

    def prepare(self, config, context):
        # Convert native telescope files into the normalized FITS contract under
        # context.staging_dir here, then update config['input']['input_dir'] and
        # config['psf']['psf_dir'] before returning.
        return PreparedRun(config=config, metadata={"example": True})
