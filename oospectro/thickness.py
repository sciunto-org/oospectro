import numpy as np
#from sklearn import linear_model, datasets
from skimage.measure import ransac, LineModelND
from scipy import stats
from scipy.signal import find_peaks
from scipy.interpolate import interp1d
from scipy.fftpack import fft, fftfreq

import matplotlib.pyplot as plt


class OptimizeResult(dict):
    """ Represents the optimization result.

    Notes
    -----
    This class has been copied from scipy.optimize

    """
    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError:
            raise AttributeError(name)

    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

    def __repr__(self):
        if self.keys():
            m = max(map(len, list(self.keys()))) + 1
            return '\n'.join([k.rjust(m) + ': ' + repr(v)
                              for k, v in sorted(self.items())])
        else:
            return self.__class__.__name__ + "()"

    def __dir__(self):
        return list(self.keys())


def thickness_from_minmax(lambdas, intensities, refractive_index=1.,
                          min_peak_prominence=0.01, min_peak_distance=10,
                          method='linreg', debug=False):
    """
    Return the thickness from a min-max detection.

    Parameters
    ----------
    lambdas : array
        Wavelength values in nm.
    intensities : array
        Intensity values.
    refractive_index : scalar, optional
        Value of the refractive index of the medium.
    min_peak_prominence : scalar, optional
        Required prominence of peaks.
    min_peak_distance : scalar, optional
        Minimum distance between peaks.
    method : string, optional
        Either 'linreg' for linear regression or 'ransac'
        for Randon Sampling Consensus.
    debug : boolean, optional
        Show plots of peak detection and lin regression.

    Returns
    -------
    results : Instance of `OptimizeResult` class.
        The attribute `thickness` gives the thickness value in nm.

    Notes
    -----
    For more details about `min_peak_prominence` and `min_peak_distance`,
    see the documentation of `scipy.signal.find_peaks`. This function
    is used to find extrema.
    """
    peaks_max, _ = find_peaks(intensities, prominence=min_peak_prominence, distance=min_peak_distance)
    peaks_min, _ = find_peaks(-intensities, prominence=min_peak_prominence, distance=min_peak_distance)
    peaks = np.concatenate((peaks_min, peaks_max))
    peaks.sort()

    k_values = np.arange(len(1/lambdas[peaks]))

    if k_values.size < 2:
        # Can't fit if less than two points.
        return np.nan

    if method.lower() == 'ransac':
        residual_threshold = 4e-5
        min_samples = 2
        # Scikit-image
        data = np.column_stack([k_values, 1/lambdas[peaks][::-1]])
        model_robust, inliers = ransac(data, LineModelND, min_samples=min_samples,
                                       residual_threshold=residual_threshold,
                                       max_trials=100)
        outliers = (inliers == False)
        slope = model_robust.params[1][1]
        thickness_minmax = 1 / slope / refractive_index / 4

        # Scikit-learn
        #X, y = k_values.reshape(-1, 1), 1/lambdas[peaks][::-1]

        ## Fit line using all data
        #lr = linear_model.LinearRegression()
        #lr.fit(X, y)

        #slransac = linear_model.RANSACRegressor(min_samples=min_samples,
        #                                        residual_threshold=residual_threshold)
        #slransac.fit(X, y)
        #inlier_mask = slransac.inlier_mask_
        #outlier_mask = np.logical_not(inlier_mask)

        ## Predict data of estimated models
        #line_X = np.arange(X.min(), X.max())[:, np.newaxis]
        #line_y = lr.predict(line_X)
        #line_y_ransac = slransac.predict(line_X)

        #slope = slransac.estimator_.coef_[0]

        if debug:
            fig, ax = plt.subplots(ncols=2, figsize=(15, 6))

            ax[0].set_xlabel('lambda')
            ax[0].set_ylabel('Intensity')
            ax[0].plot(lambdas, intensities)
            ax[0].plot(lambdas[peaks_min], intensities[peaks_min], "s")
            ax[0].plot(lambdas[peaks_max], intensities[peaks_max], "o")

            ax[1].set_xlabel('1 / lambda')
            ax[1].set_ylabel('min & max')
            ax[1].plot(data[inliers, 0], data[inliers, 1], 'xb', alpha=0.6, label='Inlier data')
            ax[1].plot(data[outliers, 0], data[outliers, 1], '+r', alpha=0.6, label='Outlier data')
            ax[1].plot(k_values, model_robust.predict_y(k_values), '-g', label='Robust line model')

            plt.legend(loc='lower right')
            plt.show()

        return OptimizeResult(thickness=thickness_minmax,
                              num_inliers=inliers.sum(),
                              num_outliers=outliers.sum(),
                              peaks_max=peaks_max,
                              peaks_min=peaks_min,
                              )

    elif method.lower() == 'linreg':
        slope, intercept, r_value, p_value, std_err = stats.linregress(k_values,
                                                                       1/lambdas[peaks][::-1])

        thickness_minmax = 1 / slope / refractive_index / 4

        if debug:
            fig, axes = plt.subplots(ncols=2, figsize=(15, 6))
            ax = axes.ravel()

            ax[0].set_xlabel('lambda')
            ax[0].set_ylabel('Intensity')
            ax[0].plot(lambdas, intensities)
            ax[0].plot(lambdas[peaks_min], intensities[peaks_min], "s")
            ax[0].plot(lambdas[peaks_max], intensities[peaks_max], "o")

            ax[1].set_xlabel('1 / lambda')
            ax[1].set_ylabel('min & max')
            ax[1].plot(1 / lambdas[peaks][::-1], 's')
            ax[1].plot(k_values, intercept + k_values * slope)
            plt.show()

        return OptimizeResult(thickness=thickness_minmax,
                              peaks_max=peaks_max,
                              peaks_min=peaks_min,
                              stderr=std_err)

    else:
        raise ValueError('Wrong method')


def thickness_from_fft(lambdas, intensities,
                       refractive_index=1.,
                       num_half_space=None,
                       debug=False):
    """
    Determine the tickness by Fast Fourier Transform.

    Parameters
    ----------
    lambdas : array
        Wavelength values in nm.
    intensities : array
        Intensity values.
    refractive_index : scalar, optional
        Value of the refractive index of the medium.
    num_half_space : scalar, optional
        Number of points to compute FFT's half space.
        If `None`, default corresponds to `10*len(lambdas)`.
    debug : boolean, optional
        Show plot of the transformed signal and the peak detection.

    Returns
    -------
    results : Instance of `OptimizeResult` class.
        The attribute `thickness` gives the thickness value in nm.
    """
    if num_half_space is None:
        num_half_space = 10 * len(lambdas)

    # FFT requires evenly spaced data.
    # So, we interpolate the signal.
    # Interpolate to get a linear increase of 1 / lambdas.
    inverse_lambdas = 1 / lambdas
    f = interp1d(inverse_lambdas, intensities)

    inverse_lambdas_linspace = np.linspace(inverse_lambdas[0],
                                           inverse_lambdas[-1],
                                           2*num_half_space)
    intensities_linspace = f(inverse_lambdas_linspace)

    # Perform FFT
    density = (inverse_lambdas[-1]-inverse_lambdas[0]) / (2*num_half_space)
    inverse_lambdas_fft = fftfreq(2*num_half_space, density)
    intensities_fft = fft(intensities_linspace)

    # The FFT is symetrical over [0:N] and [N:2N].
    # Keep over [N:2N], ie for positive freq.
    intensities_fft = intensities_fft[num_half_space:2*num_half_space]
    inverse_lambdas_fft = inverse_lambdas_fft[num_half_space:2*num_half_space]

    idx_max_fft = np.argmax(abs(intensities_fft))
    freq_max = inverse_lambdas_fft[idx_max_fft]
    thickness_fft = freq_max / 2. / refractive_index

    if debug:
        plt.loglog(inverse_lambdas_fft, np.abs(intensities_fft))
        plt.loglog(freq_max, np.abs(intensities_fft[idx_max_fft]), 'o')
        plt.xlabel('Frequency')
        plt.ylabel('FFT(Intensity)')
        plt.title(f'Thickness={thickness_fft:.2f}')

    return OptimizeResult(thickness=thickness_fft,)
