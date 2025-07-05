"""Module for converting eddy covariance flux units.

This module provides functionality to convert kinematic flux units (w'c')
to final physical units (e.g., μg m⁻² s⁻¹ or μmol m⁻² s⁻¹) in eddy
covariance measurements. The conversions account for air density and,
when needed, molecular mass of the measured scalar.

The module assumes standard units for input variables:
- Vertical wind (w): m s⁻¹
- IRGA concentrations: ppm (μmol mol⁻¹)
- TRACER concentrations: ppb (nmol mol⁻¹)

Author
------
B. Heinesch
University of Liege, Gembloux Agro-Bio Tech
"""

import sys
import numpy as np


def flux_unit_conversion(flux_individual, flux_noise_mean, flux_noise_std, flux_noise_rmse, random_error_FS, random_flux, random_error_noise,
                         conv_type, air_mol_conc, MM=np.nan):
    """Convert eddy covariance flux and uncertainty estimates to physical units.

    This function converts kinematic fluxes (w'c') and their associated
    uncertainty metrics from raw covariance units to physical flux units.
    The conversion requires air molar concentration and, for mass-based
    units, the molecular mass of the measured scalar.

    Parameters
    ----------
    flux_individual : float or array_like
        Individual flux values [raw kinematic units]
    flux_noise_mean : float or array_like
        Mean of covariance noise [raw kinematic units]
    flux_noise_std : float or array_like
        Standard deviation of covariance noise [raw kinematic units]
    flux_noise_rmse : float or array_like
        RMSE of covariance noise [raw kinematic units]
    random_error_FS : float or array_like
        Random error from Finkelstein & Sims method [raw kinematic units]
    random_flux : float or array_like
        Random error from flux randomization [raw kinematic units]
    random_error_noise : float or array_like
        Random error from white noise [raw kinematic units]
    conv_type : {'to ug m-2 s-1', 'to umol m-2 s-1'}
        Type of unit conversion to perform:
        - 'to ug m-2 s-1': Convert to micrograms per square meter per second
        - 'to umol m-2 s-1': Convert to micromoles per square meter per second
    air_mol_conc : float
        Air molar concentration [mol m⁻³]
    MM : float, optional
        Molar mass of measured scalar [g mol⁻¹]
        Required only for mass-based conversions (μg m⁻² s⁻¹)

    Returns
    -------
    tuple
        7-element tuple containing the converted values in the same order
        as the input parameters, now in physical units [μg m⁻² s⁻¹ or
        μmol m⁻² s⁻¹]

    Notes
    -----
    The conversion factors are:
    1. For molar fluxes (μmol m⁻² s⁻¹):
       flux * air_mol_conc
    2. For mass fluxes (μg m⁻² s⁻¹):
       flux * air_mol_conc * MM / 1000

    The division by 1000 in the mass flux conversion accounts for the
    conversion from ng to μg.

    Examples
    --------
    >>> # Convert CO₂ flux to μmol m⁻² s⁻¹
    >>> flux_conv = flux_unit_conversion(
    ...     0.1, 0.01, 0.02, 0.015, 0.03, 0.025, 0.02,
    ...     'to umol m-2 s-1', 41.6
    ... )
    >>> print(f'Converted flux: {flux_conv[0]:.2f} μmol m⁻² s⁻¹')

    >>> # Convert CH₄ flux to μg m⁻² s⁻¹
    >>> flux_conv = flux_unit_conversion(
    ...     0.1, 0.01, 0.02, 0.015, 0.03, 0.025, 0.02,
    ...     'to ug m-2 s-1', 41.6, MM=16.04
    ... )
    >>> print(f'Converted flux: {flux_conv[0]:.2f} μg m⁻² s⁻¹')
    """

    if conv_type == 'to ug m-2 s-1':
        factor = air_mol_conc * MM / 1000.
    elif conv_type == 'to umol m-2 s-1':
        factor = air_mol_conc
    else:
        sys.exit('this flux unit conversion is not allowed')

    return (
        flux_individual * factor,
        flux_noise_mean * factor,
        flux_noise_std * factor,
        flux_noise_rmse * factor,
        random_error_FS * factor,
        random_flux * factor,
        random_error_noise * factor
    )