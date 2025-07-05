"""Module for synchronizing sonic anemometer and tracer gas analyzer files.

This module provides functionality to match sonic anemometer and tracer
gas analyzer files based on their timestamps, ensuring that flux
calculations only proceed when both measurements are available.
It handles:

- Parsing of standardized filename formats
- Half-hourly period matching
- Filtering of unpaired measurements

Author
------
B. Heinesch
University of Liege, Gembloux Agro-Bio Tech
"""

from datetime import datetime


def parse_datetime_sonic(filename):
    """Parse datetime from sonic anemometer filename.

    Parameters
    ----------
    filename : str
        Filename in format 'GHS50_YYYY_MM_DD__HH_MM_SS.hdf5'

    Returns
    -------
    datetime
        Parsed datetime object
    """
    parts = filename.split('__')  # Split by double underscore
    date_str = parts[0].split('_', 1)[1]  # Remove GHS50_ prefix and get date
    time_str = parts[1].replace('.hdf5', '')  # Get time part
    return datetime.strptime(f"{date_str}__{time_str}", "%Y_%m_%d__%H_%M_%S")


def parse_datetime_tracer(filename):
    """Parse datetime from tracer gas analyzer filename.

    Parameters
    ----------
    filename : str
        Filename in format 'K2_BE-Vie_YYYY_TOF4000_YYYY_MM_DD__HH_MM_SS.h5'

    Returns
    -------
    datetime
        Parsed datetime object
    """
    parts = filename.split('__')  # Split by double underscore
    date_str = parts[0].split('_')[-3:]  # Get last 3 parts of first section (year, month, day)
    date_str = '_'.join(date_str)  # Rejoin with underscores
    time_str = parts[1].replace('.h5', '')  # Get time part
    return datetime.strptime(f"{date_str}__{time_str}", "%Y_%m_%d__%H_%M_%S")


def map_sonic2tracer(all_sonic_files_list, all_tracer_files_list):
    """Filter sonic files to keep only those with matching tracer files.

    This function examines sonic anemometer files and keeps only those
    that have corresponding tracer gas analyzer files within the same
    half-hour period. This ensures data synchronization for flux
    calculations.

    Parameters
    ----------
    all_sonic_files_list : dict
        Dictionary of sonic files containing:
        - name : list
            Filenames
        - path : list
            Full paths to files
        - prefix : list
            File prefixes
        - date : list
            File dates
    all_tracer_files_list : dict
        Dictionary of tracer files with same structure

    Returns
    -------
    dict
        Filtered sonic files dictionary containing only entries
        with matching tracer files. Maintains same structure as
        input dictionary.

    Notes
    -----
    The matching process:
    1. For sonic files at MM=00, looks for tracer files between MM=00-29
    2. For sonic files at MM=30, looks for tracer files between MM=30-59
    3. Files must match exactly in year, month, day, and hour
    """

    # Convert tracer filenames to datetime objects
    tracer_times = [parse_datetime_tracer(name) for name in all_tracer_files_list['name']]

    # Get the indices of valid sonic files
    valid_indices = []
    for i, name in enumerate(all_sonic_files_list['name']):
        sonic_time = parse_datetime_sonic(name)

        # For each sonic file, check if there's a tracer file in the same half-hour period
        for tracer_time in tracer_times:
            # If sonic is at MM=00, look for tracer between MM=00 and MM=29
            if sonic_time.minute == 0:
                if (tracer_time.year == sonic_time.year and
                    tracer_time.month == sonic_time.month and
                    tracer_time.day == sonic_time.day and
                    tracer_time.hour == sonic_time.hour and
                    0 <= tracer_time.minute < 30):
                    valid_indices.append(i)
                    break
            # If sonic is at MM=30, look for tracer between MM=30 and MM=59
            elif sonic_time.minute == 30:
                if (tracer_time.year == sonic_time.year and
                    tracer_time.month == sonic_time.month and
                    tracer_time.day == sonic_time.day and
                    tracer_time.hour == sonic_time.hour and
                    30 <= tracer_time.minute < 60):
                    valid_indices.append(i)
                    break

    # Create filtered dictionary with same keys as input
    filtered_sonic = {}
    for key in all_sonic_files_list.keys():
        filtered_sonic[key] = [all_sonic_files_list[key][i] for i in valid_indices]

    return filtered_sonic
