import numpy as np
from xcov import xcov


def flux_uncertainties(ini, w_prime, c_prime):
    """
    Computes different forms of turbulent flux uncertainties

    parameters
    ----------
    ini: dict, initialisation information
    w_prime :  np 1D array, raw high-frequency time series of vertical wind component (i.e. W)
    c_prime :  np 1D array, raw high-frequency time series of scalar atmospheric variable (e.g. CO2)

    returns
    -------
    flux_noise_mean: float, flux detection limit using mean of covariance
                     between +/- 160-180s (Spirig et al. 2005)
    flux_noise_std: float, flux detection limit using STD of covariance
                     between +/- 160-180s (Spirig et al. 2005)
    flux_noise_rmse: float, flux detection limit: flux noise criterium
                     using RMSE noise of covariance between +/- 160-180s
    random_error_FS: float, random error as described by Finkelstein and Sims 2001
    random_flux: float, flux detection limit: random shuffle criterium as described by Billesbach 2011
    random_error_noise: float, random error from white noise using
                        autocovariance (Lenschow 2000, Mauder 2013, Langford 2015)

    Comments
    --------
    - for the white noise computation, modification compared to MATLAB Innflux :
      sqrt of the difference instead of difference of the sqrt (cfr. Langford 2015 p 4200 and 4201)
      warning : in the MATLAB Innflux version, interpolated variance at zero lag can be negative and the MATLAB sqrt gives an imaginary number

    Written by B. Heinesch.
    University of Liege, Gembloux Agro-Bio Tech.
    """

    # flux detection limit: flux noise criterium using STD noise of covariance between +/- 160-180s (Spirig et al. 2005)
    cov_wc_left = xcov(w_prime, c_prime, [- 180*ini['param']['SAMPLING_RATE_FINAL'], -160*ini['param']['SAMPLING_RATE_FINAL']])
    cov_wc_right = xcov(w_prime, c_prime, [160*ini['param']['SAMPLING_RATE_FINAL'], 180*ini['param']['SAMPLING_RATE_FINAL']])
    flux_noise_mean = np.nanmean(np.concatenate((cov_wc_left, cov_wc_right)))
    flux_noise_std = np.nanstd(np.concatenate((cov_wc_left, cov_wc_right)), ddof=1)
    flux_noise_rmse = np.sqrt(0.5*((np.nanstd(cov_wc_left, ddof=1))**2 +
                                      (np.nanmean(cov_wc_left))**2 +
                                      (np.nanstd(cov_wc_right, ddof=1))**2 +
                                      (np.nanmean(cov_wc_right))**2
                                    ))

    # estimate random error as described by Finkelstein and Sims 2001
    random_error_FS = np.sqrt(1/(ini['param']['WINDOW_LENGTH']*ini['param']['SAMPLING_RATE_FINAL'])
                              * (sum(xcov(w_prime, w_prime, [0, 99]) * xcov(c_prime, c_prime, [0, 99]))
                                 + sum(xcov(w_prime, c_prime, [0, 99])*xcov(c_prime, w_prime, [0, 99]))
                               )
                              )

    # flux detection limit: random shuffle criterium as described by Billesbach 2011
    xcov_rand = np.empty(10)
    for irand in range(0, 10):
        w_rand = w_prime[np.random.permutation(len(w_prime))]
        xcov_rand[irand] = xcov(w_rand, c_prime, [0, 0])
    random_flux = np.nanstd(xcov_rand)

    # estimate white noise using autocovariance (Lenschow 2000, Mauder 2013, Langford 2015)
    autocov_c = xcov(c_prime, c_prime, [0, 4])
    x_fit = np.arange(0, 4, dtype=float)
    y_fit = np.asarray(autocov_c[0:4], dtype=float)
    mask = np.isfinite(x_fit) & np.isfinite(y_fit)
    if np.isfinite(autocov_c[4]) and int(mask.sum()) >= 2:
        try:
            p = np.polyfit(x_fit[mask], y_fit[mask], 1)
            y4 = np.polyval(p, 4)
            white_noise_c = np.sqrt(autocov_c[4] - y4)
        except np.linalg.LinAlgError:
            white_noise_c = np.nan
    else:
        white_noise_c = np.nan    
    random_error_noise = white_noise_c * np.nanstd(w_prime, ddof=1) / np.sqrt((ini['param']['WINDOW_LENGTH']*ini['param']['SAMPLING_RATE_FINAL']))

    return flux_noise_mean, flux_noise_std, flux_noise_rmse, random_error_FS, random_flux, random_error_noise
