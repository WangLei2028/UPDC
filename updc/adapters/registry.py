from __future__ import annotations

from importlib import metadata
from typing import Iterable

from .base import InstrumentAdapter
from .builtin import GenericFITSAdapter, JWSTAdapter

_BUILTINS = {"jwst": JWSTAdapter, "fits": GenericFITSAdapter}
ENTRYPOINT_GROUP = "updc.instrument_adapters"


def _coerce_adapter(obj) -> InstrumentAdapter:
    if isinstance(obj, InstrumentAdapter):
        return obj
    if isinstance(obj, type) and issubclass(obj, InstrumentAdapter):
        return obj()
    if callable(obj):
        candidate = obj()
        if isinstance(candidate, InstrumentAdapter):
            return candidate
    raise TypeError("instrument adapter entry point must provide InstrumentAdapter")


def get_adapter(name: str) -> InstrumentAdapter:
    key = name.strip().lower()
    if key in _BUILTINS:
        return _BUILTINS[key]()
    try:
        eps = metadata.entry_points(group=ENTRYPOINT_GROUP)
    except TypeError:  # compatibility with older importlib.metadata APIs
        eps = metadata.entry_points().get(ENTRYPOINT_GROUP, ())
    for ep in eps:
        if ep.name.lower() == key:
            adapter = _coerce_adapter(ep.load())
            if adapter.name.lower() != key:
                adapter.name = key
            return adapter
    raise KeyError(f"Unknown instrument adapter {name!r}")


def adapter_descriptions() -> list[tuple[str, str, str]]:
    rows = [(name, cls().description, "builtin") for name, cls in sorted(_BUILTINS.items())]
    try:
        eps: Iterable = metadata.entry_points(group=ENTRYPOINT_GROUP)
    except TypeError:
        eps = metadata.entry_points().get(ENTRYPOINT_GROUP, ())
    for ep in eps:
        if ep.name.lower() not in _BUILTINS:
            rows.append((ep.name, f"plugin: {ep.value}", "plugin"))
    return rows
