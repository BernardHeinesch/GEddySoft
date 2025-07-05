"""Module for computing cross-covariance between two time series.

This module provides a function to calculate the cross-covariance between two signals
with specified time lags, similar to MATLAB's xcov function but with important
differences in mean handling. The implementation avoids circular effects from numpy's
roll operation and handles NaN values appropriately.

Author
------
B. Heinesch
University of Liege, Gembloux Agro-Bio Tech
"""

import numpy as np


def xcov(x, y, lag):
    """Calculate cross-covariance between two time series.

    This function computes the cross-covariance between two signals for specified
    time lags. It handles NaN values and avoids circular effects from numpy's roll
    operation. The implementation differs from MATLAB's xcov in how means are
    computed for non-zero lags.

    Parameters
    ----------
    x : array_like
        First time series signal
    y : array_like
        Second time series signal
    lag : array_like
        Two-element list/array specifying [start, end] lag interval in samples.
        Both positive and negative lags are supported.

    Returns
    -------
    array_like or float
        Cross-covariance values for each lag. Returns a float if only one lag
        is specified, otherwise returns an array.

    Notes
    -----
    Key differences from MATLAB's xcov:
    1. For non-zero lags, this implementation uses means of the truncated arrays,
       while MATLAB uses means of the full arrays.
    2. Results will match MATLAB for zero lag but differ for non-zero lags.
    3. NaN values are properly handled using numpy's nansum and count_nonzero.

    The cross-covariance is computed as:
    cov(x,y) = Σ((x - mean(x))(y - mean(y))) / (n-1)
    where the summation and means are computed over the overlapping portions
    of the signals at each lag.

    Examples
    --------
    >>> x = np.array([1, 2, 3, 4, 5])
    >>> y = np.array([2, 3, 4, 5, 6])
    >>> result = xcov(x, y, [0, 2])  # Compute for lags 0, 1, and 2
    >>> print(result)  # Shows covariance at each lag

    See Also
    --------
    numpy.correlate : NumPy's correlation function
    numpy.cov : NumPy's covariance function
    """

    lag = np.linspace(lag[0], lag[1], num=lag[1]-lag[0]+1, dtype='int16')

    # create output array
    crosscov = np.full((len(lag)), np.nan)

    for i in range(0, len(lag), 1):

        yshifted = np.roll(y, lag[i])
        if lag[i] > 0:

            crosscov[i] = (
                (np.nansum(x[lag[i]:] * yshifted[lag[i]:]) -
                 np.nansum(x[lag[i]:]) * np.nansum(yshifted[lag[i]:]) / np.count_nonzero(~np.isnan(x[lag[i]:])))
                / (np.count_nonzero(~np.isnan(x[lag[i]:])) - 1)
            )
        else:
            crosscov[i] = (
                (np.nansum(x[0:len(x)-1-lag[i]] * yshifted[0:len(x)-1-lag[i]]) -
                 np.nansum(x[0:len(x)-1-lag[i]]) * np.nansum(yshifted[0:len(x)-1-lag[i]]) /
                 np.count_nonzero(~np.isnan(x[0:len(x)-1-lag[i]])))
                / (np.count_nonzero(~np.isnan(x[0:len(x)-1-lag[i]])) - 1)
            )

    # reduce to float if unique lag sent
    if len(crosscov) == 1:
        crosscov = crosscov[0]

    return crosscov
