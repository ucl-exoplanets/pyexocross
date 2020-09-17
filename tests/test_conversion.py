import pytest
import numpy as np

def test_spectral_conv():
    from pyexocross.util import convert_to_wavenumber

    assert convert_to_wavenumber(0.1, 'um') == 10000/0.1
    assert convert_to_wavenumber(100, 'k') == 100.0
    np.testing.assert_almost_equal(convert_to_wavenumber((0.1, 10), 'um'), 
                                   np.array([100000.0, 1000]))