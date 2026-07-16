"""
spike_detection_vickers_97.py
----------------------------------
Started form a python conversion of the EDDYPRO Fortran spike detection algorithm, 
based on Vickers and Mahrt (1997) spike detection method.
Subsequently modified by Bernard Heinesch and Jonathan Bitton

"""
from typing import Tuple, Optional
import numpy as np
import matplotlib.pyplot as plt

def linear_interpolate_spikes(data: np.ndarray, is_spike: np.ndarray, error_value: Optional[float] = None) -> np.ndarray:
    """
    Replace spikes with linear interpolation using neighboring valid values.
    Started form a python conversion of the EDDYPRO Fortran spike detection algorithm, 
    based on Vickers and Mahrt (1997) spike detection method.
    Subsequently modified by Bernard Heinesch and Jonathan Bitton
    
    Parameters:
    -----------
    data : np.ndarray
        Input data array
    is_spike : np.ndarray
        Boolean array indicating spike locations
    error_value : float
        Value indicating invalid/error data
    
    Returns:
    --------
    np.ndarray
        Data array with spikes replaced by linear interpolation
    """
    data_out = data.copy()
    N = len(data)
    
    invalid_mask = np.isnan(data) | np.isinf(data)
    if error_value is not None:
        invalid_mask |= (data == error_value)

    # Find consecutive spike sequences
    spike_starts = np.where(np.diff(np.concatenate(([0], is_spike.astype(int)))) == 1)[0]
    spike_ends = np.where(np.diff(np.concatenate((is_spike.astype(int), [0]))) == -1)[0]

    for start, end in zip(spike_starts, spike_ends):
        # Find valid points before and after the spike sequence
        left_idx = start - 1
        right_idx = end + 1

        # Look for valid left point
        while left_idx >= 0 and (is_spike[left_idx] or invalid_mask[left_idx]):
            left_idx -= 1

        # Look for valid right point
        while right_idx < N and (is_spike[right_idx] or invalid_mask[right_idx]):
            right_idx += 1

        # Perform interpolation if valid points are found
        if left_idx >= 0 and right_idx < N:
            # Linear interpolation
            x = np.array([left_idx, right_idx])
            y = np.array([data[left_idx], data[right_idx]])
            x_interp = np.arange(start, end + 1)
            data_out[start:end + 1] = np.interp(x_interp, x, y)
        elif left_idx >= 0:  # Only left point available
            data_out[start:end + 1] = data[left_idx]
        elif right_idx < N:  # Only right point available
            data_out[start:end + 1] = data[right_idx]

    return data_out

# %%

def spike_detection_vickers97(data: np.ndarray,
                                   spike_mode: Optional[int] = 1,
                                   max_pass: Optional[int] = 10,
                                   avrg_len: Optional[int] = 30,
                                   ac_freq: Optional[int] = 10,
                                   spike_limit: Optional[float] = 3.5,
                                   max_consec_spikes: Optional[int] = 3,
                                   error_value: Optional[float] = None,
                                   series_label: Optional[str] = None,
                                   ctrplot: Optional[bool] = False
                                   ) -> Tuple[np.ndarray, np.ndarray, int]:
    """
    Detects and counts spikes, and replaces them by linear interpolation if requested.
    Hard-flags file for too many spikes.
    Started form a python conversion of the EDDYPRO Fortran spike detection algorithm, 
    based on Vickers and Mahrt (1997) spike detection method.
    Subsequently modified by Bernard Heinesch and Jonathan Bitton

    Parameters:
    -----------
    data : np.ndarray
        1D input data array
    spike_mode : int
        If 1: detect and flag spikes only (no interpolation; detected spikes are masked as NaN in the
        working signal for subsequent passes). If 2: detect, flag and remove spikes by linear
        interpolation.
    max_pass: int
        maximum number of passes for spikes detection
    avrg_len: int
        lenght of data, in minutes
    ac_freq: float
        Frequency (Hz)
    spike_limit : float
        Standard deviation multiplier for spike detection
    max_consec_spikes : int
        Maximum number of consecutive points that can be considered spikes
    series_label : str
        Optional label used only for diagnostic plots (ctrplot=True) to indicate which time series is
        being processed.
    ctrplot : bool
        Whether to generate and save diagnostic plots

    Returns:
    --------
    Tuple[np.ndarray, np.ndarray, int]
        - Modified data array with spikes removed (and interpolated if requested)
        - Boolean array indicating individual spike points (spike locations). The number of
          individual spike points is given by sum(is_spike).
        - Number of spike events detected (an event is a contiguous sequence of consecutive spike
          points, i.e. a run of outliers with length <= max_consec_spikes).

    Comments:
    --------
    - check step (EDDYPRO manual says half the window size. May speed up the routine).

    """

    # Parameters
    step = 100  # window advancement in samples. VM97 say 1. EP manual v7 says: 'The window moves forward
                # half its length at a time', but 100 samples are hard-coded in the version 7.0.4.
    lim_step = 0.1  # increase of inliers range

    N = len(data)

    # Calculate window length
    win_len = avrg_len // 6
    win_len = max(1, win_len)
    nn = int(win_len * ac_freq * 60)  # window length in samples
    wdw_num = (N - nn) // step + 1  # number of windows for current file

    # Initialize arrays
    is_spike = np.zeros(N, dtype=bool)
    is_spike_pass = np.zeros(N, dtype=bool)
    loc_mean = np.zeros(N)
    loc_stdev = np.zeros(N)
    data_out = data.copy()

    # Main processing loop
    passes = 0
    nspikes = 0
    adv_lim = spike_limit

    while passes < max_pass:
        passes += 1
        is_spike_pass.fill(False)
        nspikes_sng = 0
        cnt = 0  # Counter for consecutive outliers

        invalid_mask = np.isnan(data_out) | np.isinf(data_out)
        if error_value is not None:
            invalid_mask |= (data_out == error_value)

        # Process each window
        for wdw in range(wdw_num):
            # Extract window data
            start_idx = wdw * step
            window_data = data_out[start_idx:start_idx + nn]

            # Calculate window statistics
            valid_window = ~invalid_mask[start_idx:start_idx + nn]
            window_mean = np.nanmean(window_data[valid_window])
            window_std = np.nanstd(window_data[valid_window])

            # Define central points range
            imin = nn//2 - step//2 + step * wdw
            imax = nn//2 + step//2 - 1 + step * wdw

            # Assign local statistics
            loc_mean[imin:imax+1] = window_mean
            loc_stdev[imin:imax+1] = window_std

            # Handle first and last windows
            if wdw == 0:
                loc_mean[:imin] = window_mean
                loc_stdev[:imin] = window_std
            if wdw == wdw_num - 1:
                loc_mean[imax:] = window_mean
                loc_stdev[imax:] = window_std

        # Spike detection with consecutive outlier checking
        i = 0
        while i < N:
            if invalid_mask[i]:
                i += 1
                continue

            upper_limit = loc_mean[i] + adv_lim * loc_stdev[i]
            lower_limit = loc_mean[i] - adv_lim * loc_stdev[i]

            # Check if point is an outlier
            if data_out[i] > upper_limit or data_out[i] < lower_limit:
                cnt += 1
                i += 1
            else:
                # Found a valid point, check the previous sequence
                if cnt > 0 and cnt <= max_consec_spikes:
                    new_spike = False
                    # Mark the previous cnt points as spikes
                    for k in range(i-cnt, i):
                        if not is_spike[k]:
                            new_spike = True
                            nspikes_sng += 1
                            is_spike[k] = True
                            is_spike_pass[k] = True
                    if new_spike:
                        nspikes += 1
                cnt = 0  # Reset counter
                i += 1
        # Handle spike sequence reaching the end of the series
        if cnt > 0 and cnt <= max_consec_spikes:
            new_spike = False
            for k in range(N - cnt, N):
                if not is_spike[k]:
                    new_spike = True
                    nspikes_sng += 1
                    is_spike[k] = True
                    is_spike_pass[k] = True
            if new_spike:
                nspikes += 1

        # Check if another pass is needed
        if nspikes_sng == 0: 
            break
        else:
            # Replace spikes using linear interpolation
            if spike_mode == 2:
                data_out = linear_interpolate_spikes(data_out, is_spike_pass)
            elif spike_mode == 1:
                data_out[is_spike_pass] = np.nan

        # Adjust limits for next pass
        adv_lim += lim_step

    # print(f"\nDetected {np.sum(is_spike)} spikes in {passes} passes")

    # Generate diagnostic plots if requested
    if ctrplot:
        title_prefix = '' if series_label is None else f'{series_label} - '
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, 8), sharex=True)

        # Data with spikes highlighted
        ax1.plot(data, 'b-', label='Original', alpha=0.7)
        spike_points = np.where(is_spike)[0]
        if len(spike_points) > 0:
            ax1.plot(spike_points, data[spike_points], 'rx', label='Detected Spikes')
        ax1.set_title(f'{title_prefix}Original Data with Detected Spikes')
        ax1.set_ylabel('Value')
        ax1.grid(True)
        ax1.legend()

        # Final interpolated data
        ax2.plot(data_out, 'g-', label='Cleaned (Interpolated)', alpha=0.7)
        if len(spike_points) > 0:
            ax2.plot(spike_points, data[spike_points], 'rx', label='Original Spike Values')
        ax2.set_title(f'{title_prefix}Final Data with Interpolated Values')
        ax2.set_xlabel('Sample Number')
        ax2.set_ylabel('Value')
        ax2.grid(True)
        ax2.legend()

        plt.tight_layout()
        # plt.savefig('spike_detection_results.png')
        # plt.close()

    return data_out, is_spike, nspikes