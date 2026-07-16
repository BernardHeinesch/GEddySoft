import numpy as np


def xcov(x, y, lag):
    """
    Perform Cross-Correlation on x and y

    parameters
    ----------
    x    : 1st signal
    y    : 2nd signal
    lag  : lag interval [start,end] in samples

    returns
    -------
    corr : lagged covariance

    comments:
        - not influenced by circular behaviour of the roll operation 
        - when translating from MATLAB (InnFlux), no xcov function was found in python packages
        - MATLAB computes the mean of the initial array while python (and excel) uses the mean
          of the truncated arrays. So gives similar results for a zero lag but different results
          for a non-zero lag !

    Written by B. Heinesch.
    University of Liege, Gembloux Agro-Bio Tech.

    """

    lag = np.linspace(lag[0], lag[1], num=lag[1]-lag[0]+1, dtype='int16')

    # create output array
    crosscov = np.full((len(lag)), np.nan)

    x = np.asarray(x)
    y = np.asarray(y)

    for i in range(0, len(lag), 1):

        l = int(lag[i])
        yshifted = np.roll(y, l)

        if l > 0:
            xs = x[l:]
            ys = yshifted[l:]
        else:
            xs = x[0:len(x)-1-l]
            ys = yshifted[0:len(x)-1-l]

        mask = np.isfinite(xs) & np.isfinite(ys)
        n = int(mask.sum())
        if n < 2:
            crosscov[i] = np.nan
            continue

        xs = xs[mask]
        ys = ys[mask]
        crosscov[i] = (np.sum(xs * ys) - (np.sum(xs) * np.sum(ys) / n)) / (n - 1)

    # reduce to float if unique lag sent
    if len(crosscov) == 1:
        crosscov = crosscov[0]

    return crosscov
