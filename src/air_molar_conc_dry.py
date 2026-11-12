import numpy as np
from esat import esat


def air_molar_conc_dry(T_meteo, T_mean, rh_meteo, p_mean, chi_H2O=np.nan):
    """
    Dry air molar concentration (mol m-3):
        rho_dry = (p - ea) / (R * T)
    with ea the water vapour pressure (hPa).

    Fallback chain for ea:
        1) ea = rh_meteo/100 * esat(T)   (meteo file relative humidity, in %)
        2) ea = chi_H2O * p_mean         (measured H2O mole fraction)
        3) ea = 0                        (dry air assumption)

    parameters
    ----------
    T_meteo : float, air temperature from meteo file (K); NaN -> T_mean used
    T_mean  : float, mean sonic temperature (K)
    rh_meteo : float, relative humidity from meteo file (%); NaN allowed
    p_mean  : float, mean air pressure (hPa)
    chi_H2O : float, measured H2O mole fraction (-); NaN if not available

    returns
    -------
    air_mol_conc : float, dry air molar concentration (mol m-3)
    """

    Rgas = 8.314462618  # J mol-1 K-1 (CODATA 2018)

    T_use = T_meteo if np.isfinite(T_meteo) else T_mean

    ea = np.nan  # hPa
    if np.isfinite(rh_meteo):
        ea = (rh_meteo/100.) * esat(T_use)
    elif np.isfinite(chi_H2O):
        ea = chi_H2O * p_mean
    if not np.isfinite(ea) or ea < 0:
        ea = 0.
    ea = min(ea, p_mean)

    return (p_mean - ea)*1e2 / (Rgas * T_use)  # p_mean in hPa -> Pa
