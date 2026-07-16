import sys
import numpy as np


def flux_unit_conversion(flux_individual, flux_noise_mean, flux_noise_std, flux_noise_rmse, random_error_FS, random_flux, random_error_noise,
                         conv_type, air_mol_conc, MM=np.nan):
    """
    Convert flux units from kinematics to final ones
    ! postulates that w is in m s-1, IRGA concentration is in ppm and TRACER concentration in ppb

    parameters
    ----------
    air_mol_conc: float, air molar concentration (mol m-3)
    MM: float, molar mass (-)
    conv_type: string, convertion type ('to ug m-2 s-1', ...)

    returns
    -------
    all arguments but the 3 ones described above, converted to the desired units

    Comments
    --------

    Written by B. Heinesch.
    University of Liege, Gembloux Agro-Bio Tech.
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