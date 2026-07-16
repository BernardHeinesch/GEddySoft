import numpy as np


def nanlinfit(x):
    """
    Calculate linear fit for a vector X. X may contain NaN entries.

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

    x = np.delete(x, np.where(np.isnan(x)))
    t = np.transpose(np.arange(0, len(x)))
    coeff = np.polyfit(t, x, 1)
    slope = coeff[0]
    offset = coeff[1]

    return [slope, offset]
