"""Module for detecting instrumental problems in eddy covariance data.

This module implements the instrumental problem detection methods from
Vitale et al. 2020 (Biogeosciences), providing comprehensive quality
checks for eddy covariance measurements. It includes:

- Statistical tests for data quality
- Visualization of test results
- Configurable thresholds for different variables
- Comprehensive error detection metrics

References
----------
Vitale et al. 2020, Biogeosciences
    "Robust data cleaning procedure for eddy covariance flux measurements"

Author
------
B. Heinesch
University of Liege, Gembloux Agro-Bio Tech
October 17, 2022
"""

import numpy as np
from scipy.stats import kurtosis
from scipy.stats import skew as skewness
from statsmodels.tsa.stattools import acf
from statsmodels.robust.scale import qn_scale
# from qn_scale_test import qn_scale
import matplotlib.pyplot as plt
from diptest import diptest

from nandetrend import nandetrend


def inst_prob_test(x, detrend=False, hz=20, plot=False, var_name='w'):
    """Detect instrumental problems in eddy covariance time series.

    This function implements multiple statistical tests to identify
    potential instrumental problems in high-frequency eddy covariance
    data, following the methodology of Vitale et al. 2020.

    Parameters
    ----------
    x : array_like
        Raw high-frequency eddy covariance time series
    detrend : bool, optional
        If True, removes linear trend from data. Default is False.
    hz : int, optional
        Data acquisition frequency in Hz (e.g., 10 or 20).
        Default is 20 Hz.
    plot : bool, optional
        If True, generates diagnostic plots. Default is False.
    var_name : {'U', 'V', 'W', 'T_SONIC', 'CO2', 'H2O'}, optional
        Variable name for OOR test limits and plot labels.
        Default is 'w'.

    Returns
    -------
    list
        Quality metrics in the following order:
        - S_VM97 : float
            Skewness (Vickers and Mahrt 1997)
        - K_VM97 : float
            Kurtosis (Vickers and Mahrt 1997)
        - KID0 : float
            Kurtosis index on differenced data (with zeros)
        - KID : float
            Kurtosis index on differenced data (no zeros)
        - HF5 : float
            Fluctuation homogeneity (5σ limit) [%]
        - HF10 : float
            Fluctuation homogeneity (10σ limit) [%]
        - HD5 : float
            Differenced data homogeneity (5σ limit) [%]
        - HD10 : float
            Differenced data homogeneity (10σ limit) [%]
        - AL1 : float
            Autocorrelation at lag 1
        - DDI : float
            Data distribution integrity
        - DIP : float
            Hartigan's dip test p-value
        - OOR : float
            Out of range percentage [%]

    Notes
    -----
    Implementation details:
    1. Uses Qn scale estimator for robust standard deviation
       Qn was found to be extremely slow. And far worse if using the robustbase package.
    2. Removes repeated values for differenced data (for HD5, HD10 and KID) using 1E-3 threshold,
       which might not be suitable for all tracers and might affect .
    3. Assumes 30-minute averaging periods
    4. OOR limits (in m/s and °C):
       - U, V: [-30, 30]
       - W: [-10, 10]
       - T_SONIC: [-50, 50]
       - CO2: [300, 1000]
       - H2O: [0, 30]

    See Vitale 2020, Table 1 for recommended thresholds. 
    See also the associated Rflux manual and Rflux vignette.

    References
    ----------
    .. [1] Vitale et al. 2020, Biogeosciences
           "Robust data cleaning procedure for eddy covariance"
    .. [2] Vickers and Mahrt 1997
           "Quality Control and Flux Sampling Problems for Tower"
    """


    nl = 36000 if len(x) > 20000 else 18000

    # compute statistics

    if not detrend:
        flucts = x - np.nanmean(x)
    else:
        nans, xx = np.isnan(x), lambda z: z.nonzero()[0]
        x[nans] = np.interp(xx(nans), xx(~nans), x[~nans])
        flucts = nandetrend(x)

    sigma_f = max(0.01, qn_scale(flucts[~(np.isnan(flucts))]))

    HF5 = np.count_nonzero(abs(flucts) > 5 * sigma_f) / len(flucts) * 100.
    HF10 = np.count_nonzero(abs(flucts) > 10 * sigma_f) / len(flucts) * 100.

    d0 = np.diff(x)
    x_r = np.copy(x)
    x_r[np.where(abs(d0) < 1E-3)[0] + 1] = np.nan
    d1 = np.diff(x_r)
    if np.sum(~np.isnan(d1)) > 1000:
        sigma_d = max(0.01, qn_scale(d1[~(np.isnan(d1))]))
        HD5 = np.count_nonzero(abs(d0) > 5*sigma_d) / len(d0) * 100.
        HD10 = np.count_nonzero(abs(d0) > 10*sigma_d) / len(d0) * 100.
    else:
        # too much repeated values
        HD5 = np.nan
        HD10 = np.nan

    detrended_x = nandetrend(x[~np.isnan(x)])
    S_VM97 = skewness(detrended_x)
    K_VM97 = 3+kurtosis(detrended_x)
    KID0 = 3+kurtosis(d0, nan_policy='omit')
    KID = 3+kurtosis(np.diff(x_r[~np.isnan(x_r)]), nan_policy='omit')

    # Check if `x` has variation
    x_min, x_max = np.nanmin(x), np.nanmax(x)
    has_variation = x_min != x_max

    AL1 = acf(x, missing='conservative')[1] if has_variation else np.nan
    breaks = np.histogram_bin_edges(x, bins='fd', range=(x_min, x_max))
    DDI = max(np.histogram(x, bins=breaks)[0]) if has_variation else nl

    dip, pval = diptest(flucts)
    DIP = pval if has_variation else np.nan

    # Out of range test (with limits hard-coded, in ms-1 and °C)
    if var_name in {'U', 'V', 'W', 'T_SONIC'}:
        min_val = {'U': [-30], 'V': [-30], 'W': [-10], 'T_SONIC': [-50], 'CO2': [300], 'H2O': [0]}
        max_val = {'U': [30], 'V': [30], 'W': [10], 'T_SONIC': [50], 'CO2': [1000], 'H2O': [30]}
        OOR = (np.sum((x < min_val[var_name]) | (x > max_val[var_name]))) / len(x) * 100
    else:
        OOR = np.nan

    # plot results if requested
    if plot:

        # treshold definitions (based on the RFlux vignette)
        # be careful taht for DIP test, the severe threshold is sometimes set to 0.05
        # and sometimes grounded (like in the L2 levels of the Carbon Portal!)
        mod_thr = [np.nan, np.nan, np.nan, 30, 2, 0.5, 2, 0.5, 0.75, hz*150, 0.1]
        sev_thr = [np.nan, np.nan, np.nan, 50, 4, 1, 4, 1, 0.5, hz*300, 0.1]

        if (var_name == 'U' or var_name == 'V' or var_name == 'W'): lab = var_name + '\n(ms-1)'
        if (var_name == 'T_SONIC'): lab = var_name + '\n(K)'
        if (var_name == 'c or h'): lab = var_name + '\n(\u03BCmol/mol or (mmol/mol))'
        if (var_name == 'tracer'): lab = var_name + '\n(mmol/mol)'

        fig, ax = plt.subplots(4, 2, figsize=(8, 6), gridspec_kw={'width_ratios': [4, 1]})

        # plot color meaning
        ax[0, 1].text(0.6, 1, 'SevEr', c='red')
        ax[0, 1].text(0.6, 0.8, 'ModEr', c='orange')
        ax[0, 1].text(0.6, 0.6, 'NoEr', c='gray')

        # plot fluctuations data test
        ax[0, 0].plot(x)
        ax[0, 0].set_ylabel(lab)
        ax[0, 0].set_xlabel('Samples')

        ax[0, 1].set_axis_off()
        color = 'red' if HF5 > sev_thr[4] else ('orange' if (HF5 > mod_thr[4] and HF5 < sev_thr[4]) else 'gray')
        ax[0, 1].text(0, 0.4, 'HF5 : ' + str(HF5), c=color)
        color = 'red' if HF10 > sev_thr[5] else ('orange' if (HF10 > mod_thr[5] and HF10 < sev_thr[5]) else 'gray')
        ax[0, 1].text(0, 0.2, 'HF10: ' + str(HF10), c=color)

        # plot differenced data test
        ax[1, 0].plot(d0)
        ax[1, 0].set_ylabel('delta' + lab)
        ax[1, 0].set_xlabel('Samples')

        ax[1, 1].set_axis_off()
        color = 'red' if KID > sev_thr[3] else ('orange' if (KID > mod_thr[3] and KID < sev_thr[3]) else 'gray')
        ax[1, 1].text(0.03, 0.7, 'KID : ' + f'{KID:.{3}f}', c=color)
        color = 'red' if HD5 > sev_thr[6] else ('orange' if (HD5 > mod_thr[6] and HD5 < sev_thr[6]) else 'gray')
        ax[1, 1].text(0, 0.5, 'HD5: ' + str(HD5), c=color)
        color = 'red' if HD10 > sev_thr[7] else ('orange' if (HD10 > mod_thr[7] and HD10 < sev_thr[7]) else 'gray')
        ax[1, 1].text(0, 0.3, 'HD10: ' + str(HD10), c=color)

        # plot autocorrelation at lag1 test
        ax[2, 0].plot(acf(x, missing='conservative')[0:100])
        ax[2, 0].set_ylabel('ACF')
        ax[2, 0].set_xlabel('Lags')

        ax[2, 1].set_axis_off()
        color = 'red' if AL1 < sev_thr[8] else ('orange' if (AL1 < mod_thr[8]) else 'gray')
        ax[2, 1].text(0, 0.5, 'AL1 : ' + f'{AL1:.{3}f}', c=color)

        # plot data distribution integrity and Hartigan's Dip tests
        ax[3, 0].hist(x, breaks)
        ax[3, 0].set_ylabel('Frequency')
        ax[3, 0].set_xlabel(lab)  

        ax[3, 1].set_axis_off()
        color = 'red' if DDI > sev_thr[9] else ('orange' if (DDI > mod_thr[9] and DDI < sev_thr[9]) else 'gray')
        ax[3, 1].text(0, 0.6, 'DDI : ' + str(DDI), c=color)
        color = 'red' if DIP < sev_thr[10] else ('orange' if (DIP < mod_thr[10] and DIP > sev_thr[10]) else 'gray')
        ax[3, 1].text(0, 0.4, 'DIP :' + f'{DIP:.{3}f}', c=color)

        fig.subplots_adjust(hspace=0.9)

        plt.show(block=False)
        # plt.pause(1)  # pauses one second before disappearing
#        plt.close('all')  # close all previous plots

    return [S_VM97, K_VM97, KID0, KID, HF5, HF10, HD5, HD10, AL1, DDI, DIP, OOR]
