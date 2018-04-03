Example
=======

.. code-block:: python

    from oospectro import load_data, thickness_from_minmax

    spectrum = 'spectra/sample1/007084.xy'
    lambdas, intensities = load_data(spectrum, lambda_min=450)
    result = thickness_from_minmax(lambdas, intensities,
                                   refractive_index=1.33,
                                   min_peak_prominence=0.02,
                                   method='ransac', debug=True)
    print(result)
