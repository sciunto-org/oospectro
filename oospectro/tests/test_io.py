import pytest
from numpy.testing import assert_equal

from oospectro import io

path = 'spectra/sample1/000266.xy'

def test_load_file():
    data = io.load_spectrum(path)
    assert(len(data) == 2)
    assert(data[0].shape == data[1].shape)
