import numpy as np

from xcov import xcov
from find_covariance_peak import find_covariance_peak


def compute_time_lag(ini, w_prime, c_prime, lag_center, prescribed_lag_samples, cov_wc, plot_label):
    """
    Compute the optimal time lag between vertical wind fluctuations and a scalar signal.

    parameters
    ----------
    ini: dict
        Main configuration dictionary.
        Required keys:
            - param.LAG_DETECT_METHOD
            - param.LAG_OUTER_WINDOW_SIZE
            - param.LAG_INNER_WINDOW_SIZE
            - param.LAG_COVPEAK_FILTER_LENGTH
            - run_param.PLOT_FIND_COVARIANCE_PEAK

    w_prime: np.ndarray
        Vertical wind fluctuation time-series for the averaging interval.

    c_prime: np.ndarray
        Scalar fluctuation time-series for the averaging interval.

    lag_center: int
        Center of the lag search window in samples (typically physical lag + clock drift).

    prescribed_lag_samples: int or float
        Prescribed lag in samples (used when LAG_DETECT_METHOD == 'PRESCRIBED').

    cov_wc: np.ndarray
        Preallocated cross-covariance window, shape (2*LAG_OUTER_WINDOW_SIZE + 1,).
        For CONST/PRESCRIBED methods, it is left unchanged (typically filled with NaNs).

    plot_label: str
        Label used when plotting the covariance peak search.

    returns
    -------
    lag_samples: int
        Detected lag in samples.

    cov_wc: np.ndarray
        Cross-covariance window used for peak detection. For CONST/PRESCRIBED methods,
        this remains the input array (typically NaNs).

    comments
    --------
    The function implements the lag logic used in GEddySoft:
        - CONST: use lag_center directly
        - PRESCRIBED: use the provided prescribed_lag_samples
        - MAX: compute a cross-covariance window and pick the peak
        - MAX_WITH_DEFAULT: same as MAX, but falls back to lag_center if the detected lag
          lies on the edge of the inner search window.

    The cross-covariance window is only computed for MAX/MAX_WITH_DEFAULT to avoid
    unnecessary computations for CONST/PRESCRIBED while keeping cov_wc defined.
    """

    if ini['param']['LAG_DETECT_METHOD'] == 'CONST':
        lag_samples = lag_center
    elif ini['param']['LAG_DETECT_METHOD'] == 'PRESCRIBED':
        lag_samples = int(prescribed_lag_samples)
    elif ini['param']['LAG_DETECT_METHOD'] in ('MAX', 'MAX_WITH_DEFAULT'):
        cov_wc = xcov(w_prime, c_prime, [
            lag_center - ini['param']['LAG_OUTER_WINDOW_SIZE'],
            lag_center + ini['param']['LAG_OUTER_WINDOW_SIZE']
        ])

        lag_samples = find_covariance_peak(
            cov_wc,
            ini['param']['LAG_OUTER_WINDOW_SIZE'],
            ini['param']['LAG_INNER_WINDOW_SIZE'],
            lag_center,
            ini['param']['LAG_COVPEAK_FILTER_LENGTH'],
            0,
            ini['run_param']['PLOT_FIND_COVARIANCE_PEAK'],
            plot_label,
        )  # lag in samples

        if ini['param']['LAG_DETECT_METHOD'] == 'MAX_WITH_DEFAULT':
            if (lag_samples == lag_center - ini['param']['LAG_INNER_WINDOW_SIZE'] + 1 or
                lag_samples == lag_center + ini['param']['LAG_INNER_WINDOW_SIZE']):
                lag_samples = lag_center
    else:
        raise ValueError(f"Unknown LAG_DETECT_METHOD: {ini['param']['LAG_DETECT_METHOD']}")

    return lag_samples, cov_wc
