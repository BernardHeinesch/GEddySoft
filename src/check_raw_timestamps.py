"""Module for validating timestamp continuity in eddy covariance data.

This module provides functionality to check high-frequency eddy covariance
time series for data quality issues related to timestamps, including:

- Detection of gaps exceeding a specified threshold
- Identification of timestamp overlaps
- Reporting of problematic time periods

These checks are critical for eddy covariance processing as gaps or
overlaps in the raw data can significantly impact flux calculations.

Author
------
B. Heinesch
University of Liege, Gembloux Agro-Bio Tech
"""

import pandas as pd


def check_raw_timestamps(time_series, threshold):
    """Check for gaps and overlaps in high-frequency time series data.

    This function analyzes a time series for two types of timestamp issues:
    1. Gaps: periods where consecutive timestamps differ by more than
       the specified threshold
    2. Overlaps: periods where timestamps are not strictly increasing

    These checks are essential for eddy covariance data quality control
    as both gaps and overlaps can lead to incorrect flux calculations.

    Parameters
    ----------
    time_series : pd.Series
        Pandas Series containing timestamps as index. Should be sorted
        in ascending order.
    threshold : pd.Timedelta
        Maximum allowable gap between consecutive timestamps. Gaps
        larger than this will be flagged.

    Returns
    -------
    tuple
        Two-element tuple containing:
        - message : str or None
            Description of any issues found, or None if no issues
        - indices : list
            Indices where problems occur, empty list if no issues

    Notes
    -----
    The function performs checks in the following order:
    1. Gaps check: Identifies timestamps separated by > threshold
    2. Overlap check: Identifies non-increasing timestamps

    For gaps, the function reports the maximum gap size found.
    For overlaps, it identifies all instances of non-increasing times.

    Examples
    --------
    >>> import pandas as pd
    >>> times = pd.Series(pd.date_range('2022-01-01', periods=5, freq='100ms'))
    >>> msg, idx = check_raw_timestamps(times, pd.Timedelta('150ms'))
    >>> print(msg)  # Should be None - no issues
    None

    >>> # Create a series with a gap
    >>> times = pd.Series([pd.Timestamp('2022-01-01 00:00:00'),
    ...                   pd.Timestamp('2022-01-01 00:00:01')])
    >>> msg, idx = check_raw_timestamps(times, pd.Timedelta('100ms'))
    >>> print(msg)  # Will report gap > 100ms
    'There are gaps (max=1 s) in the time series > specified threshold of 100 ms.'
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
