import numpy as np
import math
from xcov import xcov


def test_steady_state_FW96(w_prime, c_prime):
    """
    Performs the steady state test for w'c' according to Foken & Wichura (1996)
    see Chapter 9 in X. Lee et al. (eds.), Handbook of Micrometeorology (2004)

    parameters
    ----------
    w_prime :  raw high-frequency time series of vertical wind component (i.e. W)
    c_prime :  raw high-frequency time series of scalar atmospheric variable (e.g. CO2)

    returns
    -------
    S_FW96: Returns the non-stationary test statistic.

    Comments
    --------

    Written by B. Heinesch, 29 October, 2022.
    University of Liege, Gembloux Agro-Bio Tech.

    """

    subset_size = math.floor(len(w_prime)/5)
    wc_subset = np.zeros(5)
    for j in range(1, 6):
        wc_subset[j-1] = xcov(w_prime[((j-1)*subset_size):(j*subset_size)],
                            c_prime[((j-1)*subset_size):(j*subset_size)],
                            [0, 0])
    wc_total = xcov(w_prime, c_prime, [0, 0])
    S_FW96 = abs((np.nanmean(wc_subset) - wc_total)/wc_total)

    return S_FW96
