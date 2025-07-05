import numpy as np
import matplotlib.pyplot as plt
from scipy import stats


def test_steady_state_D99(w_prime, c_prime, plot=False):
    """
    Performs the non-stationary test described by Dutaur (1999) and Nemitz (2002)

    parameters
    ----------
    x :  raw high-frequency time series of vertical wind component (i.e. W)
    y :  raw high-frequency time series of scalar atmospheric variable (e.g. CO2)

    returns
    -------
    RCSC : Returns the relative covariance stationarity criterion .

    comments
    --------

    Written by D. Tzvetkov, April, 2023.
    University of Liege, Gembloux Agro-Bio Tech.

    """

    # cumulated sum of the product
    F_T = np.cumsum(w_prime*c_prime)/len(w_prime)

    t = np.linspace(0, len(w_prime), num=len(w_prime), dtype=int)

    # linear fit
    slope, intercept, r_value, p_value, std_mean = stats.linregress(t, F_T)

    # compute the residuals
    fitted_values = slope * t + intercept
    residuals = F_T - slope * t + intercept
    
    # compute the standard deviation of the residuals
    stddev = residuals.std()

    # compute the relative covariance stat criterion
    RCSC = abs(2*stddev/np.nanmean(w_prime*c_prime))

    # plot results if requested
    if plot:

        plt.figure(dpi=400)
        plt.title('Dutaur 1999 stationarity test')
        plt.xlabel('time')
        plt.ylabel('cum covar')
        plt.plot(t, F_T, label='data')
        plt.plot(t, t*slope + intercept, label='lin reg')
        plt.legend(title=("sigma residuals = " + str('{:.2e}'.format(np.round(stddev, 7))) + '\n' +
                   'RCSC = '+ str('{:.2e}'.format(np.round(RCSC, 6)))))

    return RCSC