import numpy as np
from nanlinfit import nanlinfit
from scipy import signal


def nandetrend(x):
    """
    remove linear trend from a vector that may contain NaN

    Parameters
    ----------
    x : 1 dim np array

    Returns
    -------
    the x array with linear trend removed

    comments
    --------
    Written by B. Heinesch.
    University of Liege, Gembloux Agro-Bio Tech.
    """

    nanidx = np.where(np.isnan(x))
    if len(nanidx) > 0:
        coeff = nanlinfit(x)
        y = x - np.transpose(np.polyval(coeff, np.arange(0, len(x))))
    else:
        y = signal.detrend(x)

    return y
