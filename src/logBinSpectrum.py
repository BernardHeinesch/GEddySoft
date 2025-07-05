"""Module for logarithmic binning of spectral data.

This module provides functionality for logarithmic binning of
spectral data, particularly useful for analyzing eddy covariance
spectra and cospectra. It implements:

- Logarithmic frequency binning
- Special handling of zero frequency
- Automatic bin averaging
- Flexible bin count control

Based on the InnFlux code for spectral analysis.

Author
------
B. Heinesch
University of Liege, Gembloux Agro-Bio Tech
November 2, 2022
"""

import math
import numpy as np
import warnings


def logBinSpectrum(f, y, N, f_min, f_max):
    """Perform logarithmic binning of spectral data.

    This function bins spectral data using logarithmically spaced
    frequency intervals. The zero frequency component is handled
    separately in the first bin.

    Parameters
    ----------
    f : numpy.ndarray
        Frequency array [Hz]
        Must be monotonically increasing
        First element should be 0 Hz
    y : numpy.ndarray
        Spectral values array
        Must be same length as f
    N : int
        Number of logarithmic bins
        Must be > 1
    f_min : float
        Minimum frequency for binning [Hz]
        Must be > 0
    f_max : float
        Maximum frequency for binning [Hz]
        Must be > f_min

    Returns
    -------
    list
        [f_out, y_out] where:
        f_out : numpy.ndarray
            Logarithmically binned frequencies [Hz]
            Length = N
            f_out[0] = 0 for zero frequency
        y_out : numpy.ndarray
            Binned spectral values
            Length = N
            y_out[0] = y[0] for zero frequency

    Notes
    -----
    1. Bin edges are calculated as:
       exp(linspace(log(f_min), log(f_max), N))
    2. Bin centers are geometric means of edges
    3. Values within each bin are arithmetically averaged
    4. Warnings are suppressed during computation
    5. NaN is used for empty bins

    Examples
    --------
    >>> f = np.linspace(0, 10, 1000)
    >>> y = np.sin(2*np.pi*f)
    >>> f_binned, y_binned = logBinSpectrum(f, y, 20, 0.1, 10)
    """

    warnings.filterwarnings("ignore")

    bounds = np.concatenate((np.array([0]), np.exp(np.linspace( math.log(f_min), math.log(f_max), N ))), axis=0)

    f_out = np.full([N], np.nan)
    y_out = np.full([N], np.nan)

    f_out[0] = 0  # first bin contains only the f=0 value
    y_out[0] = y[0]

    for i, name in enumerate(f_out[0:N-1], start=1):
        f_out[i] = np.exp( (np.log(bounds[i]) + np.log(bounds[i+1]) )/2 )

        idx = np.logical_and(f > bounds[i], f <= bounds[i+1])
        y_out[i] = np.mean(y[idx])

    return [f_out, y_out]
