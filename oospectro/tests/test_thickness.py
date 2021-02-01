import pytest

import numpy as np
from numpy.testing import assert_almost_equal
import os.path

from oospectro import io
from oospectro.thickness import thickness_from_minmax, thickness_from_fft
from spectra import known_values



def test_minmax_sample1():
    min_peak_prominence = 0.02
    min_peak_distance = 10
    index = known_values.sample1['refractive_index']
    skipped = ('sample1/011137.xy',
               'sample1/012426.xy',
               'sample1/012795.xy',
               'sample1/012979.xy',
               'sample1/011321.xy', #Insufficient number of data points
               )
    for name, value in known_values.sample1['known_thicknesses'].items():
        if name in skipped:
            continue
        path = os.path.join('spectra', name)
        lambdas, intensities = io.load_spectrum(path, lambda_min=450)
        result = thickness_from_minmax(lambdas, intensities, refractive_index=index,
                                       min_peak_prominence=min_peak_prominence,
                                       min_peak_distance=min_peak_distance,
                                       method='ransac', debug=False)
        error = 100 * np.abs(result.thickness - value) / value
        assert(error < 20)


def test_fft_sample1():
    for name, value in known_values.sample1['known_thicknesses'].items():
        index = known_values.sample1['refractive_index']
        # this method works only for large thicknesses
        if value < 3000:
            continue
        path = os.path.join('spectra', name)
        lambdas, intensities = io.load_spectrum(path, lambda_min=450)
        result = thickness_from_fft(lambdas, intensities, refractive_index=index,)
        error = 100 * np.abs(result.thickness - value) / value
        assert(error < 20)
