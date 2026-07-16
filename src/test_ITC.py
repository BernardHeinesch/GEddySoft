import numpy as np


def test_ITC(w_prime, u_prime, T_prime, wT, zoL, u_star, lat):
    """
    test on developed turbulent conditions
    Foken and Wichura 1996
    calculates integral turbulence characteristics based on flux-variance similarity

    Parameters
    ----------
    w_prime, u_prime, T_prime: high frequency time series of fluctuations for w, u and T
    wT : temperature flux <w'T'>
    zoL: stability parameter
    u_star: friction velocity

    returns
    -------
    ITC_w, ITC_u, ITC_T: relative model deviation of integral turbulence characteristics test for w, u, and T

    Comments
    ----------
    Created on Sun May 29 12:00 2022
    @author: Bernard HEINESCH

    """

    # implementation from EddyPro v7.0.4

    # Coriolis factor
    omega = 2 * np.pi / (24 * 60 * 60)  # Earth's angular velocity in rad/s
    Fcor = abs(2 * omega * np.sin(lat * np.pi / 180))

    # z+, Thomas & Foken (2002)
    zplus = 1.0

    # Modeled characteristics
    if zoL < -0.2:
        # Unstable conditions
        norm_sigma_w_model = 1.3 * (1 - 2 * zoL)**(1/3)
        norm_sigma_u_model = 4.15 * abs(zoL)**(1/8)
    elif zoL >= -0.2:
        # Neutral/stable conditions
        norm_sigma_w_model = 0.21 * np.log(Fcor * zplus / u_star) + 3.1
        norm_sigma_u_model = 0.44 * np.log(Fcor * zplus / u_star) + 6.3

    if zoL < -1:
        norm_sigma_T_model = abs(zoL)**(-1/3)
    elif -1 <= zoL < -0.0625:
        norm_sigma_T_model = abs(zoL)**(-1/4)
    elif -0.0625 <= zoL < 0.02:
        norm_sigma_T_model = 0.5 * (abs(zoL)**(-0.5))
    elif zoL > 0.02:
        norm_sigma_T_model = 1.4 * abs(zoL)**(-1/4)

    T_star = -wT/u_star
    norm_sigma_w_meas = np.nanstd(w_prime, ddof=1)/u_star
    norm_sigma_u_meas = np.nanstd(u_prime, ddof=1)/u_star
    norm_sigma_T_meas = np.nanstd(T_prime, ddof=1)/T_star
    ITC_w = abs((norm_sigma_w_model - norm_sigma_w_meas)/norm_sigma_w_model)
    ITC_u = abs((norm_sigma_u_model - norm_sigma_u_meas)/norm_sigma_u_model)
    ITC_T = abs((norm_sigma_T_model - norm_sigma_T_meas)/norm_sigma_T_model)

    return ITC_w, ITC_u, ITC_T

    # code from InnFlux, does not cover the whole zoL range !?
    # if (zoL < -1):
    #     norm_sigma_w_model = 2*(-zoL)**(1/6)
    #     norm_sigma_u_model = 2.83*(-zoL)**(1/6)
    #     norm_sigma_T_model = (-zoL)**(-1/3)
    # elif (zoL < -0.0625):
    #     norm_sigma_w_model = 2*(-zoL)**(1/8)
    #     norm_sigma_u_model = 2.83*(-zoL)**(1/8)
    #     norm_sigma_T_model = (-zoL)**(-1/4)
    # elif (zoL < 0):
    #     norm_sigma_w_model = 1.41
    #     norm_sigma_u_model = 1.99
    #     norm_sigma_T_model = 0.5*(-zoL)**(-1/2)
    # else:
    #     norm_sigma_w_model = np.NaN
    #     norm_sigma_u_model = np.NaN
    #     norm_sigma_T_model = np.NaN
