import math
import numpy as np
import warnings


def logBinSpectrum(f, y, N, f_min, f_max):
    """
    log-bin average y(f), from InnFlux code

    parameters
    ----------
    f: numpy arr of float, frequency axis [s-1]
    y: numpy arr of float, cospectrum
    N: number of bins
    f_min: lower frequency to be considered
    f_max: highest frequency to be considered

    returns
    -------
    f_out: numpy arr of float, log-binned frequency
    y_out: numpy arr of float, log-binned (co)spectrum

    comments
    --------
    Written by B. Heinesch, 2 November, 2022.
    University of Liege, Gembloux Agro-Bio Tech.

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
