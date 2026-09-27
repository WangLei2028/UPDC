from updc.cache_abi import FFT_AUTOTUNE_CACHE_ABI, LINKMAP_ABI_VERSION, PYTHON_PREPROCESS_CACHE_ABI
from updc.frozen_core import verify_frozen_core


def test_frozen_v17_core_hashes_are_unchanged():
    result = verify_frozen_core()
    assert result
    assert all(result.values()), result


def test_cache_abi_versions_are_explicit():
    assert PYTHON_PREPROCESS_CACHE_ABI == 1
    assert LINKMAP_ABI_VERSION == 2
    assert FFT_AUTOTUNE_CACHE_ABI == 1
