Example
=======

.. code-block:: python

    from oospectro import load_spectrum, thickness_from_minmax

    spectrum = 'spectra/sample1/007084.xy'
    lambdas, intensities = load_spectrum(spectrum, lambda_min=450)
    result = thickness_from_minmax(lambdas, intensities,
                                   refractive_index=1.33,
                                   min_peak_prominence=0.02,
                                   method='ransac', debug=True)
    # result is a class with various attributes
    # related to the algorithm outputs.
    # Print the thickness
    print(result.thickness)



Large thicknesses can be determined by FFT (requires numerous oscillations)

.. code-block:: python

    from oospectro import load_spectrum, thickness_from_fft

    spectrum = 'spectra/sample1/003582.xy'
    lambdas, intensities = load_spectrum(spectrum, lambda_min=450)
    result = thickness_from_fft(lambdas, intensities,
                                refractive_index=1.33)
    print(result.thickness)
