"""Module for computing low-pass filtering correction factors in eddy covariance.

This module provides functions to calculate correction factors for high-frequency
losses in eddy covariance measurements. It implements two methods:

1. Theoretical approach using Massman (2004) cospectra and a transfer function
   based on a half-power cutoff frequency.
2. Empirical approach using pre-calculated correction factors as a function of
   wind speed and atmospheric stability.

References
----------
.. [1] Massman, W. J. (2004). Concerning the measurement of atmospheric trace
       gas fluxes with open- and closed-path eddy covariance systems: The WPL
       terms and spectral attenuation. In Handbook of Micrometeorology (pp. 133-160).
.. [2] Peltola, O., et al. (2021). Insights into updating algorithms for spectral
       corrections in eddy covariance flux measurements. Atmospheric Measurement
       Techniques, 14(8), 5071-5088.

Author
------
B. Heinesch
University of Liege, Gembloux Agro-Bio Tech

Created
-------
2025-02-20
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import trapz

# %% Massman cospectrum


def theor_cospectra_massman(zoL, nf, kf, a0_st, kf0_st, mu_st,
                            a0_un, kf0_un, mu_un):
    """Calculate reference cospectra using the Massman model.

    This function implements the Massman (2004) model for scalar flux cospectra.
    The model uses different parameters for stable and unstable conditions to
    account for the effects of atmospheric stability on turbulent transport.

    Parameters
    ----------
    zoL : float
        Stability parameter (z-d)/L, where z is measurement height,
        d is displacement height, and L is Obukhov length [-]
    nf : ndarray
        Natural frequency array [Hz]
    kf : ndarray
        Normalized frequency array, f*(z-d)/U, where U is wind speed [-]
    a0_st : float
        Amplitude parameter for stable conditions [-]
    kf0_st : float
        Peak frequency parameter for stable conditions [-]
    mu_st : float
        Shape parameter for stable conditions [-]
    a0_un : float
        Amplitude parameter for unstable conditions [-]
    kf0_un : float
        Peak frequency parameter for unstable conditions [-]
    mu_un : float
        Shape parameter for unstable conditions [-]

    Returns
    -------
    ndarray
        Normalized cospectrum values for each input frequency

    Notes
    -----
    The model follows Eq. 4.2 from Massman (2004) with the form:
    Co(f) = a0 * (kf/kf0) / (1 + (kf/kf0)^(2μ))^(7/6μ) / f

    Different parameters are used for stable (zoL > 0) and unstable
    (zoL ≤ 0) conditions to better match observed cospectra.
    """
    # Initialize cospectra array
    cosp = np.zeros_like(nf)

    # Calculate cospectra using Massman formulation
    if zoL > 0:  # Stable conditions
        cosp = a0_st * (kf/kf0_st) / ((1.0 + (kf/kf0_st)**(2*mu_st))**(1.1667/mu_st)) / nf
    else:  # Unstable conditions
        cosp = a0_un * (kf/kf0_un) / ((1.0 + (kf/kf0_un)**(2*mu_un))**(1.1667/mu_un)) / nf

    return cosp


# %% correction factor for the low-pas filtering


def correction_factor_lpf(u, zoL, ini, df_lpfc=False, ctrplot=0):
    """Calculate correction factors for high-frequency flux losses.

    This function implements two methods to correct for high-frequency losses
    in eddy covariance measurements:

    1. Theoretical approach (ini['param']['LPFC']=1):
       Uses Massman model cospectra and a transfer function based on
       half-power cutoff frequency.

    2. Empirical approach (ini['param']['LPFC']=2):
       Uses pre-calculated correction factors based on wind speed
       and atmospheric stability class.

    Parameters
    ----------
    u : array_like
        Wind speed [m s⁻¹]
    zoL : array_like
        Stability parameter (z-d)/L [-]
    ini : dict
        Configuration dictionary containing:
        - param.LPFC : int
            Method selection (1 or 2)
        - param.SAMPLING_RATE_TRACER : float
            Sampling frequency [Hz]
        - param.SENSOR_HEIGHT : float
            Measurement height [m]
    df_lpfc : pandas.DataFrame, optional
        Lookup table for correction factors containing columns:
        - stability_class : str
            'stable', 'unstable', or 'all'
        - name : str
            Parameter name ('cof', 'A0', 'kf0', 'mu')
        - value : float
            Parameter value
        - ws_max : float
            Maximum wind speed for empirical method
        - CF_L : float
            Correction factor for empirical method
    ctrplot : bool, optional
        If True, generate diagnostic plots for method 1

    Returns
    -------
    float or array_like
        Correction factor(s) for low-pass filtering losses

    Notes
    -----
    For method 1, the correction factor is calculated as:
    CF = ∫(Co(f)df) / ∫(TF(f)·Co(f)df)
    where Co(f) is the Massman model cospectrum and TF(f) is
    a Lorentzian transfer function.

    For method 2, the correction factor is interpolated from
    a lookup table based on wind speed and stability class.

    See Also
    --------
    theor_cospectra_massman : Massman model implementation

    Author
    ------
    Written by B. Heinesch, 20 February, 2025.
    University of Liege, Gembloux Agro-Bio Tech.
    """

    if ini['param']['LPFC'] == 1:

        cof_all = df_lpfc[(df_lpfc['stability_class'] == 'all') & (df_lpfc['name'] == 'cof')]['value'].values[0]
        A0_unstable = df_lpfc[(df_lpfc['stability_class'] == 'unstable') & (df_lpfc['name'] == 'A0')]['value'].values[0]
        kf0_unstable = df_lpfc[(df_lpfc['stability_class'] == 'unstable') & (df_lpfc['name'] == 'kf0')]['value'].values[0]
        mu_unstable = df_lpfc[(df_lpfc['stability_class'] == 'unstable') & (df_lpfc['name'] == 'mu')]['value'].values[0]
        A0_stable = df_lpfc[(df_lpfc['stability_class'] == 'stable') & (df_lpfc['name'] == 'A0')]['value'].values[0]
        kf0_stable = df_lpfc[(df_lpfc['stability_class'] == 'stable') & (df_lpfc['name'] == 'kf0')]['value'].values[0]
        mu_stable = df_lpfc[(df_lpfc['stability_class'] == 'stable') & (df_lpfc['name'] == 'mu')]['value'].values[0]

        num_spec = 6000

        # Initialize natural frequency array
        nf = np.zeros(num_spec)
        for i in range(num_spec):
            nf[i] = float(i) * ini['param']['SAMPLING_RATE_TRACER'] / (2.0 * float(num_spec))
        nf[0] = 0.0001

        # Normalized frequencies
        kf = nf * abs((ini['param']['SENSOR_HEIGHT']) / u)

        # Get reference normalised cospectrum based on Massman model
        ncosp = theor_cospectra_massman(zoL, nf, kf, A0_stable, kf0_stable, mu_stable,
                                       A0_unstable, kf0_unstable, mu_unstable)

        # Compute low-pass transfer function using a Lorentzian
        # if cof calculated with sqrt(H), then it has to be applied as sqrt(H) as well (Peltola 2021)
        TF_LPF = (1/(1+(nf/cof_all)**2))**0.5

        # Compute correction factor
        cf_lpf = trapz(ncosp, nf) / trapz(TF_LPF * ncosp, nf)

        if ctrplot:
            # Create figure and subplots (now with only 2 plots)
            fig, (ax1) = plt.subplots(1, 1, figsize=(8, 8))

            # --- Plot 1: Natural frequency (Log-Log Scale) ---
            ax1.loglog(nf, nf * ncosp, 'g--', linewidth=2, label='Massman nat freq')  # cosp curve
            ax1.loglog(nf, nf * ncosp * TF_LPF, 'g-.', linewidth=2, label='cosp * TF_LPF')  # Product curve
            ax1.set_xlabel('Natural frequency [Hz]', fontsize=12, labelpad=10)
            ax1.set_ylabel('n*Cowc/cov', fontsize=12, labelpad=10)
            ax1.set_title(f'ref. cospectrum and its attenuation for (z-d)/L={zoL:.5f}', fontsize=14)
            ax1.grid(True)
            ax1.set_xlim(0.001, 10)
            ax1.set_ylim(0.001, 1)

            # Fill area between cosp and cosp * TF_LPF
            ax1.fill_between(nf, nf * ncosp, nf * ncosp * TF_LPF, color='green', alpha=0.3, label='Shaded Area')

            # Create a twin y-axis for TF_LPF (linear scale)
            ax1_right = ax1.twinx()
            l1 = ax1.plot(nf, nf * ncosp, 'g--', linewidth=2, label='Massman nat freq')  # Main plot
            l2 = ax1.plot(nf, nf * ncosp * TF_LPF, 'g-.', linewidth=2, label='cosp * TF_LPF')  # Product plot
            l3 = ax1_right.plot(nf, TF_LPF, 'b-', linewidth=2, label='TF_LPF')  # Secondary plot

            ax1_right.set_ylabel('TF_LPF', fontsize=12, labelpad=10, color='b')
            ax1_right.tick_params(axis='y', labelcolor='b')

            # Merge legends from both axes
            lines = l1 + l2 + l3
            labels = [line.get_label() for line in lines]
            ax1.legend(lines, labels, loc='upper right')

            # Adjust layout to prevent labels from being cut off
            plt.tight_layout(pad=3.0)
            plt.subplots_adjust(bottom=0.1, hspace=0.3)

            plt.show()

        return cf_lpf

    elif ini['param']['LPFC'] == 2:

        # Choose stability case
        stability_class = 'unstable' if zoL < 0.01 else 'stable'
        subset = df_lpfc[df_lpfc['stability_class'] == stability_class]

        # Find the first ws_max >= wind_speed, or the last if wind_speed is too high
        subset_sorted = subset.sort_values('ws_max')
        row = subset_sorted[subset_sorted['ws_max'] >= u].head(1)
        if row.empty:
            row = subset_sorted.tail(1)
        return float(row['CF_L'].values[0])
