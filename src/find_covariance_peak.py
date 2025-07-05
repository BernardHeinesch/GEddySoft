"""Module for determining time lags between correlated atmospheric signals.

This module provides functionality to estimate time lags between two signals by
analyzing their cross-covariance function.

The implementation uses smoothed cross-covariance analysis with configurable
window sizes and filtering to handle noisy signals robustly.

Author
------
B. Heinesch
University of Liege, Gembloux Agro-Bio Tech
"""

import numpy as np
from scipy.signal import filtfilt
import matplotlib.pyplot as plt


def find_covariance_peak(covariance, outer_window_size, inner_window_size, window_center, filter_length, positive_peak_only, plot=False, filename=''):
    """Estimate time lag between two signals using cross-covariance analysis.

    For two time-shifted signals that are either positively or negatively
    correlated, this function estimates their relative time lag by finding
    the maximum/minimum of their smoothed cross-covariance function within
    a specified search window.

    Parameters
    ----------
    covariance : array_like
        Cross-covariance function between the two signals
    outer_window_size : int
        Size (±samples) of the outer window used for visualization.
        The total window size will be 2*outer_window_size + 1.
    inner_window_size : int
        Size (±samples) of the inner window used for lag search.
        The peak will be searched within window_center ± inner_window_size.
    window_center : int
        Expected lag in samples around which to center the search windows
    filter_length : int
        Length of the Hamming window used for smoothing the covariance
        function before peak detection. Set to 1 for no filtering.
    positive_peak_only : bool
        If True, only search for positive peaks (maximum).
        If False, determine peak sign based on mean covariance.
    plot : bool, optional
        If True, generate diagnostic plots showing the covariance
        function and detected peak
    filename : str, optional
        Name of the file being processed, used in plot titles

    Returns
    -------
    int
        Estimated lag in samples between the two signals

    Notes
    -----
    The function implements the following steps:
    1. Optional smoothing using a Hamming window and filtfilt
    2. Peak detection within the inner search window
    3. Sign handling:
       - If positive_peak_only=True: always find maximum
       - If positive_peak_only=False: find maximum for positive
         correlation, minimum for negative correlation

    The smoothing helps reduce false peak detection in noisy data.
    The filter_length parameter should be chosen carefully:
    - Too small: susceptible to noise
    - Too large: may cause errors in filtfilt and lose true peaks

    This implementation differs from some others (e.g., InnFlux)
    by not using derivative-based methods, which can be unreliable
    with low signal-to-noise ratio data.

    Examples
    --------
    >>> # Generate sample data with known lag
    >>> import numpy as np
    >>> x = np.random.randn(1000)
    >>> y = np.roll(x, 5)  # Create copy of x shifted by 5 samples
    >>> xcov = np.correlate(x, y, mode='full')
    >>> lag = find_covariance_peak(
    ...     xcov, 50, 20, 0, 5, True, plot=True
    ... )
    >>> print(f'Detected lag: {lag} samples')

    See Also
    --------
    scipy.signal.filtfilt : Zero-phase filtering used for smoothing
    numpy.correlate : Correlation calculation
    xcov : Cross-covariance calculation
    """

    if positive_peak_only:
        b = np.hamming(filter_length)/sum(np.hamming(filter_length))
        a = 1
        filtered_cov = filtfilt(b, a, covariance.astype(float), axis=0, padtype='odd', padlen=3*(max(len(b), len(a))-1))
        idx = np.argmax(filtered_cov[(outer_window_size - inner_window_size):(outer_window_size + inner_window_size)])
        lag_samples = window_center - inner_window_size + idx
    else:
        if filter_length != 1:
            filtered_cov = filtfilt(np.hamming(filter_length)/sum(np.hamming(filter_length)), 1, covariance.astype(float))
        else:
            # no filtering
            filtered_cov = covariance

        if np.mean(filtered_cov[(outer_window_size - inner_window_size):(outer_window_size + inner_window_size)]) > 0:
            # covariance is positive
            idx = np.argmax(filtered_cov[(outer_window_size - inner_window_size):(outer_window_size + inner_window_size)])
        else:
            # covariance is negative
            idx = np.argmin(filtered_cov[(outer_window_size - inner_window_size):(outer_window_size + inner_window_size)])

        lag_samples = window_center - inner_window_size + idx

    # plot results if requested
    if plot:

        fig, ax = plt.subplots(2, 1, figsize=(8, 6))

        # plot whole cov function and found maximum
        xrange = np.linspace(-outer_window_size, outer_window_size, num=outer_window_size*2+1) + window_center
        ax[0].plot(xrange, covariance)
        ax[0].plot(xrange, filtered_cov)
        ax[0].set_ylabel('Cov (uncal units)')
        ax[0].set_xlabel('Samples')

        ax[0].set_title('Lag determination for file :' + filename)

        ax[0].axvline(x=window_center, color='black', label='axvline - full height')
        ax[0].axvline(x=-inner_window_size+window_center, color='black', linestyle='dashed')
        ax[0].axvline(x=inner_window_size + window_center, color='black', linestyle='dashed')
        ax[0].plot(lag_samples, filtered_cov[outer_window_size + lag_samples - window_center], marker="o", markersize=8, markeredgecolor="red", markerfacecolor="red")

        # zoom and search range
        xrange = np.linspace(-inner_window_size, inner_window_size, num=(2*inner_window_size)) + window_center
        ax[1].plot(xrange, covariance[(outer_window_size - inner_window_size):(outer_window_size + inner_window_size)])
        ax[1].plot(xrange, filtered_cov[(outer_window_size - inner_window_size):(outer_window_size + inner_window_size)])
        ax[1].set_ylabel('Cov (uncal units)')
        ax[1].set_xlabel('Samples')

        ax[1].axvline(x=window_center, color='black')
        ax[1].axvline(x=-inner_window_size + window_center, color='black', linestyle='dashed')
        ax[1].axvline(x=inner_window_size + window_center, color='black', linestyle='dashed')
        ax[1].plot(lag_samples, filtered_cov[outer_window_size + lag_samples - window_center], marker="o", markersize=8, markeredgecolor="red", markerfacecolor="red")

        fig.subplots_adjust(hspace=0.3)

        plt.show(block=False)
        # plt.pause(1)  # pauses one second before disappearing
        # plt.close('all')  # close all previous plots

    return int(lag_samples)
