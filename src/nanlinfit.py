"""Module for linear regression with NaN handling.

This module provides functionality for fitting linear trends to
data that may contain NaN values. It implements:

- Automatic NaN removal
- Linear regression using numpy.polyfit
- Index-based time axis generation

Primarily used for detrending time series in eddy covariance
processing.

Author
------
B. Heinesch
University of Liege, Gembloux Agro-Bio Tech
"""

import numpy as np


def nanlinfit(x):
    """Fit linear trend to data with NaN handling.

    This function fits a linear trend to data by first removing
    NaN values and then using polynomial fitting. The time axis
    is generated as sequential indices.

    Parameters
    ----------
    x : numpy.ndarray
        1-dimensional array of values
        May contain NaN values
        Must have length > 1

    Returns
    -------
    list
        [slope, offset] where:
        slope : float
            Rate of change per index unit
        offset : float
            Y-intercept of the fitted line

    Notes
    -----
    1. NaN values are removed before fitting
    2. Time axis is 0-based sequential indices
    3. Uses numpy.polyfit with degree=1
    4. Returns parameters in descending order (slope, offset)

    Examples
    --------
    >>> x = np.array([1, 2, np.nan, 4, 5])
    >>> slope, offset = nanlinfit(x)
    >>> print(f'Trend: {slope:.2f}x + {offset:.2f}')
    """


    x = np.delete(x, np.where(np.isnan(x)))
    t = np.transpose(np.arange(0, len(x)))
    coeff = np.polyfit(t, x, 1)
    slope = coeff[0]
    offset = coeff[1]

    return [slope, offset]
