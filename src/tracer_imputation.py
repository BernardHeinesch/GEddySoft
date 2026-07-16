import numpy as np
import scipy.signal
from scipy.interpolate import interp1d


def build_interp_c(tracer_time, tracer_conc, interp_time, sampling_rate_final, sampling_rate_tracer, imputation_method):
    """
    Build TRACER concentration time-series on the final time grid.

    parameters
    ----------
    tracer_time : TRACER timestamp array (seconds since epoch)
    tracer_conc : TRACER concentration array (same length as tracer_time)
    interp_time : target time grid (seconds since epoch), typically at SAMPLING_RATE_FINAL
    sampling_rate_final : sampling rate of interp_time (Hz)
    sampling_rate_tracer : sampling rate of the raw TRACER signal (Hz)
    imputation_method : integer
        0 : no imputation, only map measured timestamps to the final grid (implicit tolerance via rounding)
        1 : FFT resampling to sampling_rate_final then nearest interpolation
        2 : plateau (zero-order hold) upsampling to sampling_rate_final then nearest interpolation

    returns
    -------
    interp_c : 1 dim np array
        Concentration time-series on interp_time.

    comments
    --------
    This function is designed for the TRACER branch. It is used to build the concentration
    series on the SAMPLING_RATE_FINAL grid without modifying the raw TRACER data arrays.

    Written by B. Heinesch.
    University of Liege, Gembloux Agro-Bio Tech.

    """
    if imputation_method == 0:
        interp_c = np.full(len(interp_time), np.nan)
        _t_tr = np.asarray(tracer_time).reshape(-1)
        _c_tr = np.asarray(tracer_conc).reshape(-1)
        _idx = np.rint((_t_tr - interp_time[0]) * sampling_rate_final).astype(int)
        _mask = (_idx >= 0) & (_idx < len(interp_time)) & np.isfinite(_c_tr) & np.isfinite(_t_tr)
        interp_c[_idx[_mask]] = _c_tr[_mask]
        return interp_c

    if imputation_method not in (1, 2):
        raise ValueError(f"Unsupported IMPUTATION_METHOD: {imputation_method}")

    _t_raw = np.asarray(tracer_time).reshape(-1)
    _c_raw = np.asarray(tracer_conc).reshape(-1)

    if sampling_rate_tracer < sampling_rate_final:
        upsampling_factor = sampling_rate_final / sampling_rate_tracer
        upsampling_factor_int = int(round(upsampling_factor))
        if abs(upsampling_factor - upsampling_factor_int) > 1e-12:
            raise ValueError('SAMPLING_RATE_FINAL/SAMPLING_RATE_TRACER must be an integer for plateau upsampling')

        new_len = int(_c_raw.shape[0] * upsampling_factor)
        _t_up = _t_raw[0] + (np.arange(new_len, dtype=float) / sampling_rate_final)

        if imputation_method == 1:
            _c_up = scipy.signal.resample(_c_raw, new_len)
        else:
            half_plateau = upsampling_factor_int // 2
            plateau_idx = (np.arange(new_len) + half_plateau) // upsampling_factor_int
            plateau_idx = np.clip(plateau_idx, 0, _c_raw.shape[0] - 1)
            _c_up = _c_raw[plateau_idx]

        return interp1d(
            np.squeeze(_t_up),
            np.squeeze(_c_up),
            kind='nearest',
            bounds_error=False,
            fill_value=np.NaN,
        )(interp_time)

    return interp1d(
        np.squeeze(_t_raw),
        np.squeeze(_c_raw),
        kind='nearest',
        bounds_error=False,
        fill_value=np.NaN,
    )(interp_time)
