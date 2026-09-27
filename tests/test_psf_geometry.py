import numpy as np

pytest_plugins = []


def test_psf_vertex_geometry_if_runtime_dependencies_available():
    try:
        from updc._engine import _psf_detector_vertex_grid
    except ModuleNotFoundError as exc:
        import pytest
        pytest.skip(f"optional runtime dependency unavailable: {exc}")
    xx, yy = _psf_detector_vertex_grid((9, 9), (1023.5, 1023.5), 2)
    assert xx.shape == (10, 10)
    assert yy.shape == (10, 10)
    assert np.allclose(np.diff(xx[0]), 0.5)
    assert np.allclose(np.diff(yy[:, 0]), 0.5)
    assert np.isclose(xx[4, 4], 1023.25)
    assert np.isclose(xx[4, 5], 1023.75)
