"""Safe lifecycle management for persistent UPDC preprocessing caches.

Only content-addressed dataset directories (24 lowercase hex characters) are
eligible for automatic deletion.  FFTW wisdom outside the dataset cache root is
never deleted by these helpers.
"""
from __future__ import annotations

import re
import shutil
import time
from pathlib import Path

_DATASET_KEY_RE = re.compile(r"^[0-9a-f]{24}$")
_LAST_USED = ".last_used"


def _dataset_dirs(cache_root: Path):
    if not cache_root.is_dir():
        return []
    out = []
    for child in cache_root.iterdir():
        if child.is_dir() and _DATASET_KEY_RE.fullmatch(child.name):
            out.append(child)
    return out


def touch_dataset_cache(directory: Path) -> None:
    """Record use of a dataset cache without changing its cache key."""
    directory.mkdir(parents=True, exist_ok=True)
    marker = directory / _LAST_USED
    marker.touch(exist_ok=True)


def _last_used_time(directory: Path) -> float:
    marker = directory / _LAST_USED
    try:
        return marker.stat().st_mtime
    except OSError:
        pass
    # Backward compatibility for pre-1.1.0 caches.
    for candidate in (directory / "preprocess_cache.json", directory):
        try:
            return candidate.stat().st_mtime
        except OSError:
            continue
    return time.time()


def reset_dataset_caches(cache_root: Path, mode: str, current_key: str | None = None) -> list[Path]:
    """Delete persistent preprocessing caches according to ``mode``.

    ``current`` deletes only ``current_key``. ``all`` deletes all directories
    whose names match UPDC's content-addressed dataset-key format.  Other files
    and directories under ``cache_root`` are intentionally left untouched.
    """
    if mode == "none":
        return []
    if mode not in {"current", "all"}:
        raise ValueError(f"unsupported cache reset mode: {mode}")
    removed: list[Path] = []
    if mode == "current":
        if not current_key or not _DATASET_KEY_RE.fullmatch(current_key):
            raise ValueError("current cache reset requires a valid dataset cache key")
        targets = [cache_root / current_key]
    else:
        targets = _dataset_dirs(cache_root)
    for path in targets:
        if path.is_dir():
            shutil.rmtree(path)
            removed.append(path)
    return removed


def prune_dataset_caches(cache_root: Path, after_days: float, keep_key: str | None = None,
                         now: float | None = None) -> list[Path]:
    """Delete dataset caches unused for at least ``after_days``.

    Pruning happens when UPDC starts; it is not a background daemon.  A value
    <= 0 disables pruning. ``keep_key`` is never deleted during the current run.
    """
    if after_days <= 0:
        return []
    now = time.time() if now is None else float(now)
    cutoff = now - float(after_days) * 86400.0
    removed: list[Path] = []
    for path in _dataset_dirs(cache_root):
        if keep_key and path.name == keep_key:
            continue
        if _last_used_time(path) < cutoff:
            shutil.rmtree(path)
            removed.append(path)
    return removed
