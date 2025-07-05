import numpy as np
from typing import Tuple
import matplotlib.pyplot as plt

# %%


def linear_interpolate_spikes(data: np.ndarray, is_spike: np.ndarray, error_value: float) -> np.ndarray:
    """
    Replace detected spikes with linear interpolation.

    This function replaces spike values in a time series with linearly interpolated
    values using valid neighboring points. It handles both isolated spikes and
    consecutive sequences of spikes, as well as error values in the data.

    Parameters
    ----------
    data : numpy.ndarray
        Input data array containing the original time series with spikes
    is_spike : numpy.ndarray
        Boolean array of same length as data, True where spikes were detected
    error_value : float
        Special value indicating invalid or missing data points
        These points are skipped when finding valid neighbors for interpolation

    Returns
    -------
    numpy.ndarray
        Copy of input data with spikes replaced by linear interpolation

    Notes
    -----
    The interpolation strategy is:
    1. For each sequence of spikes, find valid (non-spike, non-error) points
       before and after the sequence
    2. If both points exist: perform linear interpolation
    3. If only one exists: use that value (nearest neighbor)
    4. If neither exists: spikes remain unchanged

    This implementation follows the EddyPro software's approach but is
    vectorized for better performance in Python.

    See Also
    --------
    spike_detection_vickers97 : Main spike detection algorithm

    Examples
    --------
    >>> import numpy as np
    >>> # Create sample data with spikes
    >>> data = np.array([1.0, 10.0, 1.1, np.nan, 1.2])
    >>> spikes = np.array([False, True, False, False, False])
    >>> cleaned = linear_interpolate_spikes(data, spikes, np.nan)
    >>> print(cleaned)  # [1.0, 1.05, 1.1, nan, 1.2]
    """
    data_out = data.copy()
    N = len(data)

    # Find consecutive spike sequences
    spike_starts = np.where(np.diff(np.concatenate(([0], is_spike.astype(int)))) == 1)[0]
    spike_ends = np.where(np.diff(np.concatenate((is_spike.astype(int), [0]))) == -1)[0]

    for start, end in zip(spike_starts, spike_ends):
        # Find valid points before and after the spike sequence
        left_idx = start - 1
        right_idx = end + 1

        # Look for valid left point
        while left_idx >= 0 and (is_spike[left_idx] or data[left_idx] == error_value):
            left_idx -= 1

        # Look for valid right point
        while right_idx < N and (is_spike[right_idx] or data[right_idx] == error_value):
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
                                   spike_mode: int = 1,
                                   max_pass: int = 10,
                                   avrg_len: int = 30,
                                   ac_freq: int = 10,
                                   spike_limit: float = 3.5,
                                   max_consec_spikes: int = 3,
                                   ctrplot: bool = False
                                   ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Detect and remove spikes in high-frequency eddy covariance data.

    This function implements the Vickers and Mahrt (1997) despiking algorithm for
    eddy covariance data. It uses an iterative moving window approach to identify
    outliers based on local statistics. The algorithm can either flag spikes or
    both flag and remove them through linear interpolation.

    Parameters
    ----------
    data : numpy.ndarray
        1D input data array containing high-frequency measurements
        (e.g., wind components, scalar concentrations)
    spike_mode : {1, 2}, optional
        Operation mode:
        - 1: Only detect spikes
        - 2: Detect and remove spikes via linear interpolation
        Default is 1
    max_pass : int, optional
        Maximum number of iterations for spike detection.
        Each pass may use progressively larger thresholds.
        Default is 10
    avrg_len : int, optional
        Averaging period length in minutes.
        Used to determine the window size for local statistics.
        Default is 30
    ac_freq : int, optional
        Data acquisition frequency in Hz.
        Used to calculate the number of samples in each window.
        Default is 10
    spike_limit : float, optional
        Initial threshold for spike detection in standard deviations.
        Points exceeding mean ± (spike_limit × std) are flagged.
        Default is 3.5
    max_consec_spikes : int, optional
        Maximum number of consecutive points that can be flagged as spikes.
        Longer sequences are not considered spikes.
        Default is 3
    ctrplot : bool, optional
        If True, generates diagnostic plots showing:
        - Original data with detected spikes
        - Cleaned data with interpolated values
        Default is False

    Returns
    -------
    data_out : numpy.ndarray
        If spike_mode=1: Copy of input with spikes still present
        If spike_mode=2: Data with spikes replaced by linear interpolation
    is_spike : numpy.ndarray
        Boolean array same length as input, True where spikes were detected

    Notes
    -----
    The algorithm follows these steps:
    1. Divides data into overlapping windows
    2. Calculates local mean and standard deviation
    3. Flags points exceeding threshold as potential spikes
    4. Checks for consecutive outliers
    5. Optionally interpolates across spike locations
    6. Repeats with adjusted threshold if spikes found

    The window advancement step is currently set to 100 samples, which
    differs from both the original VM97 paper (1 sample) and the EddyPro
    manual recommendation (half window size).

    See Also
    --------
    linear_interpolate_spikes : Function used to replace detected spikes

    References
    ----------
    .. [1] Vickers, D. and Mahrt, L. (1997). Quality control and flux sampling
           problems for tower and aircraft data. Journal of Atmospheric and
           Oceanic Technology, 14(3), 512-526.
    .. [2] LI-COR Biosciences (2019). EddyPro Software Instruction Manual,
           Version 7.0.4.

    Examples
    --------
    >>> # Generate sample data with artificial spikes
    >>> import numpy as np
    >>> data = np.random.normal(0, 1, 18000)  # 30 min at 10 Hz
    >>> data[1000:1002] = 10  # Add artificial spikes
    >>> cleaned, spikes = spike_detection_vickers97(
    ...     data, spike_mode=2, ctrplot=True
    ... )
    >>> print(f'Found {np.sum(spikes)} spikes')

    Author
    ------
    Written by Bernard Heinesch
    University of Liege, Gembloux Agro-Bio Tech
    """

    # Parameters
    step = 100  # window advancement in samples. VM97 say 1. EP manual v7 says: 'The window moves forward
                # half its length at a time', but 100 samples are hard-coded in the version 7.0.4.
    lim_step = 0.1  # increase of inliers range

    N = len(data)
    error_value = np.nan

    # Calculate window length
    win_len = avrg_len // 6
    win_len = max(1, win_len)
    nn = int(win_len * ac_freq * 60)  # window length in samples
    wdw_num = (N - nn) // step + 1  # number of windows for current file

    # Initialize arrays
    is_spike = np.zeros(N, dtype=bool)
    loc_mean = np.zeros(N)
    loc_stdev = np.zeros(N)
    data_out = data.copy()

    # Main processing loop
    passes = 0
    adv_lim = spike_limit

    while passes < max_pass:
        passes += 1
        nspikes = 0
        nspikes_sng = 0
        cnt = 0  # Counter for consecutive outliers

        # Process each window
        for wdw in range(wdw_num):
            # Extract window data
            start_idx = wdw * step
            window_data = data[start_idx:start_idx + nn]

            # Calculate window statistics
            valid_mask = window_data != error_value
            window_mean = np.nanmean(window_data[valid_mask])
            window_std = np.nanstd(window_data[valid_mask])

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
            if data[i] == error_value:
                i += 1
                continue

            upper_limit = loc_mean[i] + adv_lim * loc_stdev[i]
            lower_limit = loc_mean[i] - adv_lim * loc_stdev[i]

            # Check if point is an outlier
            if data[i] > upper_limit or data[i] < lower_limit:
                cnt += 1
                i += 1
            else:
                # Found a valid point, check the previous sequence
                if cnt > 0 and cnt <= max_consec_spikes:
                    # Mark the previous cnt points as spikes
                    for k in range(i-cnt, i):
                        if not is_spike[k]:
                            nspikes += 1
                            nspikes_sng += 1
                            is_spike[k] = True
                cnt = 0  # Reset counter
                i += 1

        # Replace spikes using linear interpolation
        if spike_mode == 2 and nspikes > 0:
            data_out = linear_interpolate_spikes(data, is_spike, error_value)

        # Check if another pass is needed
        if nspikes == 0:
            break

        # Adjust limits for next pass
        adv_lim += lim_step

    # print(f"\nDetected {np.sum(is_spike)} spikes in {passes} passes")

    # Generate diagnostic plots if requested
    if ctrplot:
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, 8), sharex=True)

        # Data with spikes highlighted
        ax1.plot(data, 'b-', label='Original', alpha=0.7)
        spike_points = np.where(is_spike)[0]
        if len(spike_points) > 0:
            ax1.plot(spike_points, data[spike_points], 'rx', label='Detected Spikes')
        ax1.set_title('Original Data with Detected Spikes')
        ax1.set_ylabel('Value')
        ax1.grid(True)
        ax1.legend()

        # Final interpolated data
        ax2.plot(data_out, 'g-', label='Cleaned (Interpolated)', alpha=0.7)
        if len(spike_points) > 0:
            ax2.plot(spike_points, data[spike_points], 'rx', label='Original Spike Values')
        ax2.set_title('Final Data with Interpolated Values')
        ax2.set_xlabel('Sample Number')
        ax2.set_ylabel('Value')
        ax2.grid(True)
        ax2.legend()

        plt.tight_layout()
        # plt.savefig('spike_detection_results.png')
        # plt.close()

    return data_out, is_spike
