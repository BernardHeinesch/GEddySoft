from datetime import datetime


def parse_datetime_sonic(filename):
    # Parse e.g. GHS50_2023_06_10__02_30_00.hdf5
    parts = filename.split('__')  # Split by double underscore
    date_str = parts[0].split('_', 1)[1]  # Remove GHS50_ prefix and get date
    time_str = parts[1].replace('.hdf5', '')  # Get time part
    return datetime.strptime(f"{date_str}__{time_str}", "%Y_%m_%d__%H_%M_%S")


def parse_datetime_tracer(filename):
    # Parse e.g. K2_BE-Vie_2023_TOF4000_2023_06_10__02_30_35.h5
    parts = filename.split('__')  # Split by double underscore
    date_str = parts[0].split('_')[-3:]  # Get last 3 parts of first section (year, month, day)
    date_str = '_'.join(date_str)  # Rejoin with underscores
    time_str = parts[1].replace('.h5', '')  # Get time part
    return datetime.strptime(f"{date_str}__{time_str}", "%Y_%m_%d__%H_%M_%S")


def map_sonic2tracer(all_sonic_files_list, all_tracer_files_list):
    """
    removes half-hours for which a sonic file is present but not a tracer file

    parameters
    ----------
    all_sonic_files_list: list, sonic files that were find in the folder given in the ini and in its subfolders
    all_tracer_files_list: list, tracer files that were find in the folder given in the ini and in its subfolders

    returns
    -------
    filtered_sonic: the filtered all_sonic_files_list

    comments
    --------
    Written by B. Heinesch.
    University of Liege, Gembloux Agro-Bio Tech.
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
