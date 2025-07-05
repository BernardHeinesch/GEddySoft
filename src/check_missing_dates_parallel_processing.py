"""Module for tracking and reporting missing dates in parallel processing output.

This module provides functionality to identify and report dates that failed
to process successfully during parallel eddy covariance calculations. It
helps ensure data completeness by:

1. Checking output files against expected dates
2. Identifying missing dates for retry processing
3. Providing formatted reporting of missing dates

The module is particularly useful in parallel processing scenarios where
some days may fail due to resource constraints or data issues.

Author
------
B. Heinesch
University of Liege, Gembloux Agro-Bio Tech
"""

import os
from datetime import datetime


def print_missing_dates(missing_dates):
    """Print a formatted report of missing dates.

    This function takes a list of dates in YYYYMMDD format and prints them
    in a human-readable format (YYYY-MM-DD). It handles both valid and
    invalid date strings gracefully.

    Parameters
    ----------
    missing_dates : list of str
        List of dates in YYYYMMDD format that are missing from the output

    Notes
    -----
    The function attempts to parse each date and format it as YYYY-MM-DD.
    If parsing fails for any date (invalid format), it prints the raw string.

    Examples
    --------
    >>> print_missing_dates(['20220101', '20220102'])
    Found 2 missing dates:
      - 2022-01-01 (raw: 20220101)
      - 2022-01-02 (raw: 20220102)
    """
    if not missing_dates:
        print("No missing dates found - all files were processed successfully!")
    else:
        print(f"Found {len(missing_dates)} missing dates:")
        for date in sorted(missing_dates):
            try:
                date_obj = datetime.strptime(date, '%Y%m%d')
                formatted_date = date_obj.strftime('%Y-%m-%d')
                print(f"  - {formatted_date} (raw: {date})")
            except ValueError:
                print(f"  - {date}")


def check_missing_dates_parallel_processing(output_folder, output_prefix, unique_days):
    """Identify dates missing from parallel processing output files.

    This function compares a list of expected processing dates against
    the actual output files produced by parallel processing. It helps
    identify which dates need to be reprocessed due to failures.

    Parameters
    ----------
    output_folder : str
        Path to the directory containing output files
    output_prefix : str
        Common prefix of output filenames (e.g., 'GEddySoft_')
    unique_days : list of str
        List of expected dates in YYYY_MM_DD format to check against

    Returns
    -------
    list of str
        Dates (in YYYY_MM_DD format) that are in unique_days but not
        found in output files

    Notes
    -----
    The function performs the following steps:
    1. Lists all files in the output directory
    2. Extracts dates from filenames matching the prefix
    3. Converts dates to YYYY_MM_DD format for comparison
    4. Returns list of dates missing from outputs

    The function assumes output filenames follow the pattern:
    {prefix}{YYYYMMDD}*

    Examples
    --------
    >>> missing = check_missing_dates_parallel_processing(
    ...     'output/dir', 'GEddySoft_',
    ...     ['2022_01_01', '2022_01_02']
    ... )
    >>> print(missing)  # If 2022_01_02 is missing
    ['2022_01_02']

    See Also
    --------
    print_missing_dates : Format and print the missing dates
    """

    # Get all files in the output folder
    all_files = os.listdir(output_folder)

    # Extract dates from filenames and convert to yyyy_mm_dd format
    file_dates = []
    file_dates = []
    for filename in all_files:
        if filename.startswith(output_prefix):
            # Extract the date part (8 characters: yyyymmdd)
            date = filename[len(output_prefix):len(output_prefix)+8]
            # Convert from yyyymmdd to yyyy_mm_dd
            converted_date = f"{date[:4]}_{date[4:6]}_{date[6:]}"
            file_dates.append(converted_date)

    # Find dates in unique_days that are not in file_dates
    missing_dates = [date for date in unique_days if date not in file_dates]

    return missing_dates
