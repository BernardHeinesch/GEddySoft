import numpy as np
from logBinSpectrum import logBinSpectrum

from theor_cosp_Kaimal import Kaimal_cosp
import matplotlib.pyplot as plt


def cospectrum(ini, x, y, f, mean_wind_speed, plot=False, filename='', spectrum_type=''):
    """
    Compute power spectrum or co-spectrum of atmospheric time series.

    This function calculates either the power spectrum (when x=y) or the co-spectrum
    (when x≠y) of atmospheric time series data. It implements standard eddy covariance
    spectral analysis techniques, including normalization and scaling according to
    Kaimal et al. (1972) conventions.

    Parameters
    ----------
    ini : dict
        Initialization dictionary containing:

        - param.WINDOW_LENGTH : float
            Length of the averaging window in seconds
        - param.SAMPLING_RATE_FINAL : float
            Final sampling rate in Hz
        - param.NUM_FREQ_BINS : int
            Number of frequency bins for log-averaging
        - param.SENSOR_HEIGHT : float
            Height of the sensor above ground [m]
    x : numpy.ndarray
        Fluctuations of the first time series. For power spectrum,
        this is the only relevant time series.
        Units: Variable-specific (e.g., m s⁻¹ for wind, ppb for concentrations)
    y : numpy.ndarray
        Fluctuations of the second time series. For power spectrum,
        should be identical to x.
        Units: Variable-specific
    f : numpy.ndarray
        Frequency axis of spectra [Hz]
    mean_wind_speed : float
        Mean wind speed used for normalization [m s⁻¹]
    plot : bool, optional
        If True, creates diagnostic plots of the (co)spectra
    filename : str, optional
        Name of the file being processed, used in plot titles
    spectrum_type : {'spec', 'cospec'}, optional
        Type of analysis being performed:
        - 'spec': power spectrum
        - 'cospec': co-spectrum

    Returns
    -------
    cospec_xy : numpy.ndarray
        Log-bin averaged (co)spectrum Co(x',y',f)
        Units: Product of input variable units per Hz
    cospec_xy_scaled : numpy.ndarray
        Normalized and frequency-scaled (co)spectrum f·Co(x',y',f)/xy
        Units: Dimensionless
    cospec_xy_integral : float
        Integral of the (co)spectrum, equals covariance for co-spectrum
        or variance for power spectrum
        Units: Product of input variable units

    Notes
    -----
    The function follows standard procedures for spectral analysis in
    micrometeorology:
    
    1. Computes two-sided spectrum using FFT
    2. Extracts and scales the one-sided spectrum
    3. Log-bins the results for clearer visualization
    4. Normalizes by integral and scales by frequency
    
    For co-spectra, the results can be compared to the Kaimal formulations
    for unstable conditions when properly normalized.

    See Also
    --------
    logBinSpectrum : Log-binning of spectral data
    theor_cosp_Kaimal : Theoretical Kaimal co-spectrum

    References
    ----------
    .. [1] Kaimal, J. C., et al. (1972). Spectral characteristics of surface‐layer
           turbulence. Q.J.R. Meteorol. Soc., 98: 563-589.
    .. [2] Stull, R. B. (1988). An Introduction to Boundary Layer Meteorology.
           Kluwer Academic Publishers.

    Examples
    --------
    >>> # Calculate power spectrum of vertical wind speed
    >>> w = np.random.randn(18000)  # 30 minutes at 10 Hz
    >>> f = np.fft.fftfreq(18000, d=0.1)[:9000]  # Positive frequencies only
    >>> ini = {'param': {
    ...     'WINDOW_LENGTH': 1800,
    ...     'SAMPLING_RATE_FINAL': 10,
    ...     'NUM_FREQ_BINS': 50,
    ...     'SENSOR_HEIGHT': 23.5
    ... }}
    >>> spec, spec_scaled, variance = cospectrum(
    ...     ini, w, w, f, 2.5, plot=True,
    ...     spectrum_type='spec'
    ... )

    Author
    ------
    Written by Bernard Heinesch (December 2022)
    University of Liege, Gembloux Agro-Bio Tech
    Updated July 2024: Added power spectrum functionality
    """

    # two-sided cospec
    cospec2 = np.real(np.conj(np.fft.fft(x))/(ini['param']['WINDOW_LENGTH']*ini['param']['SAMPLING_RATE_FINAL']) *
                      np.fft.fft(y)/(ini['param']['WINDOW_LENGTH']*ini['param']['SAMPLING_RATE_FINAL']))
    # extract the lower ~1/2 of the spec
    cospec = np.copy(cospec2[0:len(f)])
    # multiply all but the first and the last spectral entry by 2
    cospec[1:len(cospec)-2] = 2*cospec[1:len(cospec)-2]
    # scale by delta_f
    cospec = cospec / (f[1]-f[0])
    # compute integral (for normalisation purposes)
    cospec_xy_integral = np.trapz(cospec, f)

    # cospectrum log-bin-averaged in regular f domain
    f_bin, cospec_xy = logBinSpectrum(f, cospec, ini['param']['NUM_FREQ_BINS'], f[1], f[-1])
    # normalized, f-scaled log-bin-averaged cospectrum
    cospec_xy_scaled = f_bin*cospec_xy/cospec_xy_integral

    # additional comments:
    # - first value of cospec2 is the mean. In this case =zero because fluctuations are provided as input
    # - multiplication by 2 is done for retreiving the energy content (Parsival theorem)
    # - with a discrete approach, the covariance can also be calculated as the
    #   sum of cospec2 and the sum of cospec before scaling by delta f)

    # plot results if requested
    if plot:
        fig, ax = plt.subplots(2, 1, figsize=(8, 10))

        # plot f*cospec in log-log scales
        ax[0].scatter(f, f*cospec, s=10, color='grey', label='all')
        ax[0].plot(f_bin, f_bin*cospec_xy, color='blue', label='binned')
        ax[0].set_xscale('log')
        ax[0].set_yscale('log')
        ax[0].set_ylabel(f'f*{spectrum_type} (uncal units)')
        ax[0].set_xlabel('natural f [Hz]')

        ax[0].set_title(f'{spectrum_type} for file :' + filename)

        # creation of the normalized frequency scale, for ploting purposes only
        fnorm = f*ini['param']['SENSOR_HEIGHT']/mean_wind_speed
        fnorm_bin = f_bin*ini['param']['SENSOR_HEIGHT']/mean_wind_speed

        # plot f*(co)spec/cov in normalized frequency scale and log-log scales
        ax[1].scatter(fnorm, f*cospec/cospec_xy_integral, s=10, color='grey', label='all')
        ax[1].plot(fnorm_bin, cospec_xy_scaled, color='red', label='binned')
        ax[1].set_xscale('log')
        ax[1].set_yscale('log')
        if spectrum_type == 'cospec':
            ax[1].set_ylabel(f'f*{spectrum_type}/cov (uncal units)')
        else:
            ax[1].set_ylabel(f'f*{spectrum_type}/var (uncal units)')
        ax[1].set_xlabel('normalized f [-]')

        # add (unstable) Kaimal in case of cospectrum
        if spectrum_type == 'cospec':
            # compute and plot Kaimal cospectrum
            kaimal = np.zeros(len(f_bin))
            for x in range(1, len(f_bin)):
                kaimal[x] = Kaimal_cosp(f_bin[x], f_bin[x]*ini['param']['SENSOR_HEIGHT']/mean_wind_speed, -0.001)
            ax[1].plot(fnorm_bin, f_bin*kaimal, color='black', label='Kaimal)')
            ax[1].legend(loc="lower left")

        plt.show(block=False)

    return [cospec_xy, cospec_xy_scaled, cospec_xy_integral]
