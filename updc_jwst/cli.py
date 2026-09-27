"""Backward-compatible module facade for the validated v17 interface."""
from updc import _engine as _legacy

main = _legacy.main


def __getattr__(name):
    return getattr(_legacy, name)


def __dir__():
    return sorted(set(globals()) | set(dir(_legacy)))

__all__ = ["main"]
