"""Module for detrending data arrays containing NaN values.

This module provides functionality to remove linear trends from
data arrays while properly handling NaN values. It supports:

- Automatic detection of NaN values
- Linear trend removal using robust fitting
- Fallback to scipy.signal.detrend for NaN-free data

Author
------
B. Heinesch
University of Liege, Gembloux Agro-Bio Tech
"""

import numpy as np
from nanlinfit import nanlinfit
from scipy import signal


def nandetrend(x):
    """Remove linear trend from array while handling NaN values.

    This function removes linear trends from data arrays that may
    contain NaN values. It uses a robust fitting method when NaNs
    are present, and falls back to scipy.signal.detrend for
    NaN-free data.

    Parameters
    ----------
    x : ndarray
        1-dimensional numpy array to detrend. May contain NaN values.

    Returns
    -------
    ndarray
        Array of same shape as input with linear trend removed.
        NaN values in input remain NaN in output.

    Notes
    -----
    The detrending process:
    1. Check for presence of NaN values
    2. If NaNs present:
       - Fit linear trend using nanlinfit
       - Subtract trend from data
    3. If no NaNs:
       - Use scipy.signal.detrend
    """

    nanidx = np.where(np.isnan(x))
    if len(nanidx) > 0:
        coeff = nanlinfit(x)
        y = x - np.transpose(np.polyval(coeff, np.arange(0, len(x))))
    else:
        y = signal.detrend(x)

    return y
