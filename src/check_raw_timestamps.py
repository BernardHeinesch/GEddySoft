import pandas as pd

def check_raw_timestamps(time_series, threshold):
    """
    Checks if there are any gaps in a time series greater than a specified threshold,
    and if there are any overlaps in the timestamps.

    parameters
    ----------
        time_series (pd.Series): A pandas Series containing timestamps as index.
        threshold (pd.Timedelta): The maximum gap allowed between two timestamps in the time_series.

    returns
    -------
        tuple: (message, indices) where message is a string describing the issue and indices is a list
               of positions where the problem occurs. Returns (None, []) if no issues found.

    comments
    --------
    Written by B. Heinesch.
    University of Liege, Gembloux Agro-Bio Tech.
    """

    diff = time_series.diff()
    # Check for gaps
    gap_mask = diff[1:len(diff)] > threshold
    if gap_mask.any():
        gap_indices = gap_mask[gap_mask].index.tolist()
        message = f'There are gaps (max={max(diff[1:len(diff)])}) in the time series > specified threshold of {threshold} ms.'
        return message, gap_indices

    # Check for overlaps
    overlap_mask = diff[1:len(diff)] <= pd.Timedelta(milliseconds=0)
    if overlap_mask.any():
        overlap_indices = overlap_mask[overlap_mask].index.tolist()
        return "There are overlaps in the time series.", overlap_indices

    return None, []
