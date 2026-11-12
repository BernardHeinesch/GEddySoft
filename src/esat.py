import numpy as np


def esat(T):
    """
    Saturation water vapour pressure over liquid water, using the Magnus
    formula with the improved coefficients of Alduchov & Eskridge (1996).

    parameters
    ----------
    T : float or np array, temperature (K)

    returns
    -------
    saturation vapour pressure (hPa)
    """

    Tc = T - 273.15
    return 6.112 * np.exp(17.62 * Tc / (243.12 + Tc))
