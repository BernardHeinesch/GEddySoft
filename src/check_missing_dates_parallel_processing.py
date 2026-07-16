import os
from datetime import datetime


def print_missing_dates(missing_dates):
    """Print the missing dates in a formatted way."""
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
    """
    Find dates that are missing in the output files compared to unique_days list.

    parameters
    ----------
        output_folder (str): Path to the folder containing output files
        output_prefix (str): Prefix of the output files
        unique_days (list): List of dates in 'yyyymmdd' format to check against

    returns
    -------
        list: List of dates (in yyyymmdd format) that are in unique_days but not in output files

    comments
    --------
    Written by B. Heinesch.
    University of Liege, Gembloux Agro-Bio Tech.
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
