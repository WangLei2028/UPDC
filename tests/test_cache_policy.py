from pathlib import Path

from updc.cache_policy import prune_dataset_caches, reset_dataset_caches, touch_dataset_cache


def test_reset_current_and_all_are_scoped(tmp_path: Path):
    a = tmp_path / ("a" * 24)
    b = tmp_path / ("b" * 24)
    unrelated = tmp_path / "do-not-delete"
    for p in (a, b, unrelated):
        p.mkdir()
        (p / "x").write_text("x")
    removed = reset_dataset_caches(tmp_path, "current", a.name)
    assert removed == [a]
    assert not a.exists() and b.exists() and unrelated.exists()
    removed = reset_dataset_caches(tmp_path, "all")
    assert b in removed
    assert not b.exists() and unrelated.exists()


def test_prune_uses_last_used_and_keeps_current(tmp_path: Path):
    old = tmp_path / ("1" * 24)
    current = tmp_path / ("2" * 24)
    old.mkdir(); current.mkdir()
    touch_dataset_cache(old); touch_dataset_cache(current)
    old_marker = old / ".last_used"
    old_time = 1_000_000.0
    import os
    os.utime(old_marker, (old_time, old_time))
    removed = prune_dataset_caches(tmp_path, after_days=30, keep_key=current.name,
                                   now=old_time + 31 * 86400)
    assert removed == [old]
    assert not old.exists() and current.exists()
