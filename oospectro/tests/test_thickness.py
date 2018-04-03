import pytest

import numpy as np
from numpy.testing import assert_almost_equal
import os.path

from oospectro import io
from oospectro.thickness import thickness_from_minmax
from spectra import known_values



def test_sample1():
    min_peak_prominence = 0.02
    index = known_values.sample1['refractive_index']
    for name, value in known_values.sample1['known_thicknesses'].items():
        skipped = ('sample1/011137.xy',
                   'sample1/012426.xy',
                   'sample1/012795.xy',
                   'sample1/012979.xy',
                   )
        if name in skipped:
            continue
        path = os.path.join('spectra', name)
        lambdas, intensities = io.load_spectrum(path, lambda_min=450)
        result = thickness_from_minmax(lambdas, intensities, refractive_index=index,
                                       min_peak_prominence=min_peak_prominence,
                                       method='ransac', debug=False)
        error = 100 * np.abs(result - value) / value
        assert(error < 20)
