import numpy as np


def nanlinfit(x):
    """
    Calculate linear fit for a vector X. X may contain non-finite (NaN/Inf)
    entries: the fit is performed on the original sample indices of the
    finite values.

    parameters
    ----------
    x : 1 dim np array

    returns
    -------
    slope and offset of the linear fit

    comments
    --------
    Written by B. Heinesch.
    University of Liege, Gembloux Agro-Bio Tech.
    """

    t = np.transpose(np.arange(0, len(x)))
    ii = np.where(np.isfinite(x))[0]
    x = x[ii]
    t = t[ii]
    coeff = np.polyfit(t, x, 1)
    slope = coeff[0]
    offset = coeff[1]

    return [slope, offset]
