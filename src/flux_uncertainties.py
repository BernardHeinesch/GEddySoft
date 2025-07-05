import numpy as np
from xcov import xcov


def flux_uncertainties(ini, w_prime, c_prime):
    """
    Computes different forms of turbulent flux uncertainties.

    This function implements various methods to estimate flux uncertainties and detection limits
    in eddy covariance measurements, including methods from Spirig et al. (2005),
    Finkelstein and Sims (2001), and others.

    Parameters
    ----------
    ini : dict
        Initialization dictionary containing:
        
        - param.SAMPLING_RATE_FINAL : float
            Final sampling rate in Hz
        - param.WINDOW_LENGTH : float
            Length of the averaging window in seconds
    w_prime : numpy.ndarray
        1D array of high-frequency fluctuations of vertical wind component (W)
        Units: m s⁻¹
    c_prime : numpy.ndarray
        1D array of high-frequency fluctuations of scalar atmospheric variable
        (e.g., CO₂, CH₄, or VOCs)
        Units: Dependent on the scalar (e.g., ppb for VOCs)

    Returns
    -------
    flux_noise_mean : float
        Flux detection limit using mean of covariance between ±160-180s.
        Method from Spirig et al. (2005)
    flux_noise_std : float
        Flux detection limit using standard deviation of covariance
        between ±160-180s. Method from Spirig et al. (2005)
    flux_noise_rmse : float
        Flux detection limit using RMSE of covariance between ±160-180s
    random_error_FS : float
        Random error calculated using the Finkelstein and Sims (2001) method
    random_flux : float
        Flux detection limit using random shuffle method from Billesbach (2011)
    random_error_noise : float
        Random error from white noise using autocovariance methods from
        Lenschow (2000), Mauder (2013), and Langford (2015)

    Notes
    -----
    The white noise computation uses sqrt(difference) instead of difference(sqrt)
    as described in Langford (2015) pp. 4200-4201. This differs from the MATLAB
    InnFLUX implementation where interpolated variance at zero lag could be negative,
    leading to imaginary numbers in the MATLAB sqrt function.

    References
    ----------
    .. [1] Spirig, C., et al. (2005). Eddy covariance flux measurements of
           biogenic VOCs during ECHO 2003 using proton transfer reaction mass
           spectrometry. Atmospheric Chemistry and Physics, 5(2), 465-481.
    .. [2] Finkelstein, P. L., & Sims, P. F. (2001). Sampling error in eddy
           correlation flux measurements. Journal of Geophysical Research:
           Atmospheres, 106(D4), 3503-3509.
    .. [3] Billesbach, D. P. (2011). Estimating uncertainties in individual eddy
           covariance flux measurements: A comparison of methods and a proposed new
           method. Agricultural and Forest Meteorology, 151(3), 394-405.
    .. [4] Langford, B., et al. (2015). Eddy-covariance data with low signal-to-noise
           ratio: time-lag determination, uncertainties and limit of detection.
           Atmospheric Measurement Techniques, 8(10), 4197-4213.

    See Also
    --------
    xcov : Cross-covariance calculation
    cospectrum : Co-spectrum analysis

    Examples
    --------
    >>> import numpy as np
    >>> # Setup example data (30 minutes at 10 Hz)
    >>> ini = {'param': {'SAMPLING_RATE_FINAL': 10, 'WINDOW_LENGTH': 1800}}
    >>> w = np.random.randn(18000)  # Vertical wind fluctuations
    >>> c = np.random.randn(18000)  # Scalar concentration fluctuations
    >>> results = flux_uncertainties(ini, w, c)
    >>> print(f'Random error (FS method): {results[3]:.2e}')

    Author
    ------
    Written by Bernard Heinesch
    University of Liege, Gembloux Agro-Bio Tech
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
    white_noise_c = np.sqrt(autocov_c[4] - np.polyval(np.polyfit(list(range(0, 4)), np.transpose(autocov_c[0: 4]), 1), 4))
    random_error_noise = white_noise_c * np.nanstd(w_prime, ddof=1) / np.sqrt((ini['param']['WINDOW_LENGTH']*ini['param']['SAMPLING_RATE_FINAL']))

    return flux_noise_mean, flux_noise_std, flux_noise_rmse, random_error_FS, random_flux, random_error_noise
