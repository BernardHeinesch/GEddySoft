import numpy as np
import math
from xcov import xcov


def test_steady_state_M98(x, y):
    """
    Performs the non-stationary ratio test described by Mahrt (1998)
    after Vitale et al. 2020, Biogeosciences, "Robust data cleaning procedure
    for eddy covariance flux measurements"" and associated inst_prob_test
    routine found on Gitlab

    parameters
    ----------
    x :  raw high-frequency time series of vertical wind component (i.e. W)
    y :  raw high-frequency time series of scalar atmospheric variable (e.g. CO2)

    returns
    -------
    S_M98: Returns the non-stationary ratio test statistic.

    comments
    --------
    - computation of covariances might be more efficient using rolling and cov functions
      on df  (e.g. https://www.geeksforgeeks.org/how-to-calculate-rolling-correlation-in-python/),
      or the np.roll function

    Written by B. Heinesch, 29 October, 2022.
    University of Liege, Gembloux Agro-Bio Tech.
    """

    # compute COV, covariance on the whole dataset
    COV = xcov(x, y, [0, 0])

    # compute COVs, covariances over 6 intervals
    subset_size = math.floor(len(x)/6)
    COVs = np.zeros(6)
    for j in range(1, 7):
        COVs[j-1] = xcov(x[((j-1)*subset_size):(j*subset_size-1)],
                         y[((j-1)*subset_size):(j*subset_size-1)],
                         [0, 0])

    # compute COVw, covariances over 6*6 intervals
    subset_size = math.floor(len(x)/36)
    COVw = np.zeros(36)
    for j in range(1, 37):
        COVw[j-1] = xcov(x[((j-1)*subset_size):(j*subset_size-1)],
                         y[((j-1)*subset_size):(j*subset_size-1)],
                         [0, 0])

    sigmaB = np.sqrt(np.nansum((COVs - COV)**2)/(len(COVs)-1))

    sigmaW1 = np.sqrt(1/5*np.nansum((COVw[0:5] - COVs[0])**2)) if len(COVw) >= 6 else np.na
    sigmaW2 = np.sqrt(1/5*np.nansum((COVw[6:11] - COVs[1])**2)) if len(COVw) >= 12 else  np.na
    sigmaW3 = np.sqrt(1/5*np.nansum((COVw[12:17] - COVs[2])**2)) if len(COVw) >= 18 else np.na
    sigmaW4 = np.sqrt(1/5*np.nansum((COVw[18:23] - COVs[3])**2)) if len(COVw) >= 24 else np.na
    sigmaW5 = np.sqrt(1/5*np.nansum((COVw[24:29] - COVs[4])**2)) if len(COVw) >= 30 else np.na
    sigmaW6 = np.sqrt(1/len(COVw)*np.nansum((COVw[30:len(COVw)]-COVs[5])**2)) if (len(COVw) > 30 and len(COVw) <= 36) else np.nan
    sigmaWi = [sigmaW1, sigmaW2, sigmaW3, sigmaW4, sigmaW5, sigmaW6]
    sigmaW = np.nanmean(sigmaWi)

    # compute the non-stationarity ratio
    S_M98 = sigmaB/(sigmaW/np.sqrt(np.nansum(~np.isnan(sigmaWi))))

    return S_M98
