"""Module for cleaning and filtering eddy covariance results data.

This module provides functionality to clean processed results
by removing entries with invalid or missing sonic anemometer
data. It handles:

- Recursive cleaning of nested data structures
- Preservation of metadata fields
- Synchronized removal across all variables
- Special handling of quality control flags

The cleaning is based on NaN values in sonic wind speed data, which
indicates periods where core measurements were unavailable or invalid.

Author
------
B. Heinesch
University of Liege, Gembloux Agro-Bio Tech
"""

import math


def remove_entries_based_on_indices(data, nan_indices, exceptions=[]):
    """Remove entries at specified indices from nested data structures.

    This helper function recursively traverses dictionaries and lists,
    removing entries at specified indices while preserving structure
    and handling exceptions.

    Parameters
    ----------
    data : dict or list
        The data structure to clean. Can be nested to any depth.
    nan_indices : list of int
        Indices of entries to remove
    exceptions : list of str, optional
        Keys in dictionaries that should not have entries removed,
        typically metadata fields

    Returns
    -------
    list or None
        If input is a list, returns filtered list
        If input is a dict, modifies in place and returns None

    Notes
    -----
    The function handles three cases:
    1. Dictionaries: Recursively process each value
    2. Lists: Filter out specified indices
    3. Exception keys: Skip processing entirely
    """
    if isinstance(data, dict):
        for key, value in data.items():
            if key in exceptions:
                continue
            if isinstance(value, list):
                data[key] = [v for i, v in enumerate(value) if i not in nan_indices]
            elif isinstance(value, dict):
                remove_entries_based_on_indices(value, nan_indices)
    elif isinstance(data, list):
        return [v for i, v in enumerate(data) if i not in nan_indices]


def clean_results(results):
    """Clean eddy covariance results by removing entries with invalid sonic data.

    This function removes data entries where sonic anemometer measurements
    were invalid or missing (indicated by NaN values in wind speed). The
    removal is synchronized across all variables to maintain data consistency.

    Parameters
    ----------
    results : dict
        GEddySoft results dictionary containing:
        - time : list
            Timestamps for each measurement period
        - MET : dict
            Meteorological measurements including:
            - wsh : list
                Wind speed measurements (used to identify invalid periods)
            - qaqc : dict
                Quality control flags
        - TRACER : dict, optional
            Tracer gas measurements and associated QC

    Returns
    -------
    dict
        Cleaned results dictionary with invalid entries removed

    Notes
    -----
    The cleaning process:
    1. Identifies invalid periods using NaN values in wind speed
    2. Removes corresponding entries from all variables
    3. Preserves structure and metadata
    4. Handles QC flags separately to maintain integrity

    The function assumes that NaN values in wind speed (wsh)
    indicate periods where sonic data was invalid or missing,
    typically due to:
    - Incomplete data files
    - Instrument malfunctions
    - Communication errors
    """


    met = results.get('MET', {})

    # Identify indices with NaN values in 'wsh'
    nan_indices = [i for i, value in enumerate(met['wsh']) if math.isnan(value)]

    # Remove corresponding entries from 'time'
    results['time'] = [time for i, time in enumerate(results['time']) if i not in nan_indices]

    # Remove entries from MET based on nan_indices
    remove_entries_based_on_indices(met, nan_indices, exceptions=['qaqc'])

    # Remove entries from MET.qaqc based on nan_indices
    remove_entries_based_on_indices(met['qaqc'], nan_indices)

    # Remove entries from each group in TRACER based on nan_indices
    tracer = results.get('TRACER', {})
    for group_key, group_value in tracer.items():
        if isinstance(group_value, dict):
            remove_entries_based_on_indices(group_value, nan_indices, exceptions=['name', 'qaqc'])
            if 'qaqc' in group_value:
                remove_entries_based_on_indices(group_value['qaqc'], nan_indices)
        elif isinstance(group_value, list):
            tracer[group_key] = remove_entries_based_on_indices(group_value, nan_indices)

    return (results)
