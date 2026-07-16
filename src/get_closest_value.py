import numpy as np
import datetime


def get_closest_value(df, target_timestamp):
    """
    find the value in df that is the closest as the target timestamp

    parameters
    ----------
    df: df with timestamp as datetime index and var values as only column
    target_timestamp: the target timestamp

    returns
    -------
    value: the var value at target timestamp and if no value at that timestamp, 
           the median of the xth closest values, x starting with ten and increased 
           by steps of ten until values are found  

    Comments
    --------

    Written by B. Heinesch, April, 2025.
    University of Liege, Gembloux Agro-Bio Tech.

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
