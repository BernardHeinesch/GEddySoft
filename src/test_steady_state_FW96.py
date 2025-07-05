import numpy as np
import math
from xcov import xcov


def test_steady_state_FW96(w_prime, c_prime):
    """
    Test for stationarity of turbulent fluxes using the Foken & Wichura (1996) method.

    This function implements the steady state test for vertical fluxes by comparing
    the covariance calculated over the entire averaging period to the mean of
    covariances calculated over shorter subintervals. The test is particularly
    important for eddy covariance measurements where stationarity is a key assumption.

    Parameters
    ----------
    w_prime : numpy.ndarray
        High-frequency fluctuations of vertical wind velocity [m s⁻¹]
    c_prime : numpy.ndarray
        High-frequency fluctuations of scalar quantity
        (e.g., CO₂, H₂O, temperature)
        Units depend on scalar: [ppm], [ppb], [g m⁻³], or [K]

    Returns
    -------
    S_FW96 : float
        Non-stationarity parameter, calculated as:
        S_FW96 = |(<w'c'>_5 - <w'c'>_1) / <w'c'>_1|
        where:
        - <w'c'>_5 is mean covariance from five 6-min subintervals
        - <w'c'>_1 is covariance over full 30-min interval

    Notes
    -----
    The stationarity test follows these steps:
    1. Divide the 30-min period into 5 equal subintervals
    2. Calculate covariance for each subinterval
    3. Calculate covariance for the full period
    4. Compare mean of subinterval covariances to full period

    Quality classes based on S_FW96:
    - 0-30%: High quality data
    - 30-50%: Moderate quality, usable for general analysis
    - 50-100%: Low quality, not recommended for fundamental research
    - >100%: Data should be rejected

    References
    ----------
    .. [1] Foken, T. and Wichura, B. (1996). Tools for quality assessment of
           surface-based flux measurements. Agricultural and Forest Meteorology,
           78(1-2), 83-105.
    .. [2] Lee, X., Massman, W., and Law, B. (2004). Handbook of
           Micrometeorology: A Guide for Surface Flux Measurement and Analysis.
           Chapter 9. Kluwer Academic Publishers.

    See Also
    --------
    test_steady_state_M98 : Alternative steady state test by Mahrt (1998)
    test_ITC : Test for developed turbulent conditions
    xcov : Cross-covariance calculation used internally

    Examples
    --------
    >>> # Test stationarity of vertical CO2 flux
    >>> import numpy as np
    >>> # Generate sample data (30 min at 10 Hz)
    >>> w = np.random.normal(0, 0.3, 18000)
    >>> co2 = np.random.normal(400, 1, 18000)
    >>> # Add a trend to make it non-stationary
    >>> co2 += np.linspace(0, 2, 18000)
    >>> s = test_steady_state_FW96(w, co2)
    >>> print(f'Non-stationarity: {s*100:.1f}%')

    Author
    ------
    Written by Bernard Heinesch (October 2022)
    University of Liege, Gembloux Agro-Bio Tech
    """

    subset_size = math.floor(len(w_prime)/5)
    wc_subset = np.zeros(5)
    for j in range(1, 6):
        wc_subset[j-1] = xcov(w_prime[((j-1)*subset_size):(j*subset_size)],
                            c_prime[((j-1)*subset_size):(j*subset_size)],
                            [0, 0])
    wc_total = xcov(w_prime, c_prime, [0, 0])
    S_FW96 = abs((np.nanmean(wc_subset) - wc_total)/wc_total)

    return S_FW96
