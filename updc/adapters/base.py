"""Public instrument-adapter API.

Adapters are deliberately outside the frozen numerical core. They may validate,
normalize, or stage telescope-specific products into the FITS exposure contract
consumed by the engine, then return configuration overrides.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping

ConfigDict = dict[str, Any]


@dataclass(frozen=True)
class AdapterContext:
    config_path: Path
    staging_dir: Path


@dataclass
class PreparedRun:
    config: ConfigDict
    metadata: dict[str, Any] = field(default_factory=dict)


class InstrumentAdapter:
    """Base class for built-in and third-party instrument adapters.

    Third-party distributions register instances/classes in the entry-point group
    ``updc.instrument_adapters``. ``prepare`` may materialize normalized FITS files
    under ``context.staging_dir`` and return a modified configuration.
    """

    name = "base"
    description = "abstract UPDC instrument adapter"

    def apply_defaults(self, config: ConfigDict) -> ConfigDict:
        return config

    def validate(self, config: Mapping[str, Any]) -> list[str]:
        return []

    def prepare(self, config: ConfigDict, context: AdapterContext) -> PreparedRun:
        return PreparedRun(config=config)
