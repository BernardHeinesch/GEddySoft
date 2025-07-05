"""Module for robust time series value interpolation.

This module provides functionality to find the closest value in a time
series to a target timestamp, with robust handling of missing values.
It implements:

- Exact timestamp matching when available
- Adaptive window median interpolation when exact match unavailable
- Progressive window size expansion until valid values found

Author
------
B. Heinesch
University of Liege, Gembloux Agro-Bio Tech
April, 2025
"""

import numpy as np
import datetime


def get_closest_value(df, target_timestamp):
    """Find closest value to target timestamp with robust handling.

    This function searches a time series for the value closest to a
    target timestamp. If an exact match is not available or contains
    NaN, it progressively expands a window around the target time
    until valid values are found, then returns their median.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame with:
        - datetime index
        - Single column of values to interpolate
        - May contain NaN values
    target_timestamp : datetime.datetime
        Target timestamp to find closest value for

    Returns
    -------
    float
        Either:
        - Exact value if timestamp match found and not NaN
        - Median of closest valid values using adaptive window

    Notes
    -----
    The search process:
    1. Try exact timestamp match first
    2. If no match or value is NaN:
       - Start with window of 10 closest timestamps
       - Expand window by 10 until valid values found
       - Return median of valid values in window
    3. Window expansion continues until either:
       - Valid values found
       - Entire series searched
    """

    closest_index = abs(df.index - target_timestamp).argmin()
    if abs(df.index - target_timestamp).min() == datetime.timedelta() and not np.isnan(df.iloc[closest_index]):
        # the var is present for this half-hour
        value = df.iloc[closest_index]  # find closest var value based on timestamp
    else:
        max_limit = len(df)
        step = 10
        found = False
        window_size = step
    
        while window_size <= max_limit:
            # Calculate distances from the target index
            valid_indices = range(len(df))
            valid_distances = [abs(i - closest_index) for i in valid_indices]
            sorted_indices = [x for _, x in sorted(zip(valid_distances, valid_indices))]
    
            # Select the closest N values
            selected_indices = sorted_indices[:window_size]
    
            # Extract values and filter out NaNs
            selected_values = df.iloc[selected_indices].values
            valid_values = selected_values[~np.isnan(selected_values)]
    
            if len(valid_values) > 0:
                value = np.median(valid_values)
                found = True
                break
            else:
                window_size += step
    
    return value
