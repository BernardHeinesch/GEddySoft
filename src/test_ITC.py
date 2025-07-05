import numpy as np


def test_ITC(w_prime, u_prime, T_prime, wT, zoL, u_star, lat):
    """
    Test for developed turbulent conditions using Integral Turbulence Characteristics (ITC).

    This function implements the ITC test based on flux-variance similarity as described
    in Foken and Wichura (1996) and Thomas & Foken (2002). It evaluates the degree of
    developed turbulence by comparing measured and modeled standard deviations of wind
    components and temperature, normalized by appropriate scaling parameters.

    Parameters
    ----------
    w_prime : numpy.ndarray
        High-frequency fluctuations of vertical wind velocity [m s⁻¹]
    u_prime : numpy.ndarray
        High-frequency fluctuations of horizontal wind velocity [m s⁻¹]
    T_prime : numpy.ndarray
        High-frequency fluctuations of temperature [K]
    wT : float
        Kinematic temperature flux <w'T'> [K m s⁻¹]
    zoL : float
        Monin-Obukhov stability parameter (z/L)
        z: measurement height [m]
        L: Obukhov length [m]
    u_star : float
        Friction velocity [m s⁻¹]
    lat : float
        Site latitude [degrees]

    Returns
    -------
    ITC_w : float
        Relative model deviation for vertical wind component
        ITC_w = |(σ_w/u_*)_model - (σ_w/u_*)_meas| / (σ_w/u_*)_model
    ITC_u : float
        Relative model deviation for horizontal wind component
        ITC_u = |(σ_u/u_*)_model - (σ_u/u_*)_meas| / (σ_u/u_*)_model
    ITC_T : float
        Relative model deviation for temperature
        ITC_T = |(σ_T/T_*)_model - (σ_T/T_*)_meas| / (σ_T/T_*)_model

    Notes
    -----

    References
    ----------
    .. [1] Foken, T. and Wichura, B. (1996). Tools for quality assessment of
           surface-based flux measurements. Agricultural and Forest Meteorology,
           78(1-2), 83-105.
    .. [2] Thomas, C. and Foken, T. (2002). Re-evaluation of integral turbulence
           characteristics and their parameterizations. 15th Conference on Boundary
           Layer and Turbulence. 15-19 July 2002, Wageningen, The Netherlands.

    See Also
    --------
    test_steady_state_FW96 : Test for steady state conditions
    test_steady_state_M98 : Alternative steady state test

    Examples
    --------
    >>> # Example with unstable conditions
    >>> w = np.random.normal(0, 0.3, 18000)  # 30 min at 10 Hz
    >>> u = np.random.normal(0, 0.5, 18000)
    >>> T = np.random.normal(0, 0.1, 18000)
    >>> wT = -0.1  # Upward heat flux
    >>> zoL = -0.5  # Unstable
    >>> u_star = 0.4
    >>> lat = 50.0
    >>> ITC_w, ITC_u, ITC_T = test_ITC(w, u, T, wT, zoL, u_star, lat)

    Author
    ------
    Written by Bernard Heinesch (May 2022)
    University of Liege, Gembloux Agro-Bio Tech
    Based on EddyPro v7.0.4 implementation
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
