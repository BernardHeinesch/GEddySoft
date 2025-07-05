"""Module for date-based file listing and filtering.

This module provides functionality to list and filter files based on
embedded date strings in filenames. It supports:

- Recursive directory traversal
- Date string parsing and validation
- Time window filtering
- Sorting by date
- Flexible date formats

Typical use case: Finding all eddy covariance data files within a
specific time period.

Author
------
B. Heinesch
University of Liege, Gembloux Agro-Bio Tech
"""

from datetime import datetime
import os


def list_files_datestring(folder_path, prefixes, suffix, date_format, filter_crit='', sub_dir=True):
    """List and filter files by embedded date strings in filenames.

    This function recursively searches directories for files with specific
    prefixes and suffixes, extracts date information from filenames, and
    optionally filters them by a time window.

    Parameters
    ----------
    folder_path : str
        Root directory to start searching
    prefixes : str
        File prefix(es) to match
    suffix : str
        File suffix to match (e.g., '.dat', '.txt')
    date_format : str
        Format string for date in filename. Supported tokens:
        - yyyy: 4-digit year
        - mm: 2-digit month
        - dd: 2-digit day
        - HH: 2-digit hour (24-hour)
        - MM: 2-digit minute
        - SS: 2-digit second
        Example: 'yyyy_mm_dd_HH_MM_SS'
    filter_crit : tuple of str, optional
        (start_date, end_date) for filtering.
        Format: 'YYYY_MM_DD__HH_MM_SS'
        Default: '' (no filtering)
    sub_dir : bool, optional
        If True, search subdirectories recursively.
        Default: True

    Returns
    -------
    dict
        Dictionary with matched files info:
        - name : list of str
            Filenames
        - path : list of str
            Full paths to files
        - prefix : list of str
            Matched prefixes
        - date : list of int
            Ordinal dates (days since 01-01-0001)

    Notes
    -----
    1. Files containing 'trunk' after prefix are excluded
    2. For files with only date (no time):
       - HH, MM, SS fields can be empty
    3. Results are sorted by filename
    4. Empty lists are returned if no matches found

    Examples
    --------
    >>> files = list_files_datestring(
    ...     'data/',
    ...     'SONIC_',
    ...     '.dat',
    ...     'yyyy_mm_dd',
    ...     ('2023_01_01__00_00_00', '2023_12_31__23_59_59')
    ... )
    >>> print(len(files['name']), 'files found')
    """


    # define output dict
    files = {'name': [], 'path': [], 'prefix': [], 'date': []}

    exclude_string = "trunk"
    for path, subdirs, f in os.walk(folder_path):
        files['name'] += [each for each in f if (each.startswith(prefixes)
                                                 and each.endswith(suffix)
                                                 and exclude_string not in each[len(prefixes):])
                          ]
        files['path'] += [path for i in range(len(f)) if (f[i].startswith(prefixes) and
                                                          f[i].endswith(suffix) and
                                                          exclude_string not in f[i][len(prefixes):])
                          ]
        if not sub_dir:
            break

    # extract date from filenames
    date_num = []
    for n in range(len(files['name'])):
        date_string = files['name'][n][len(prefixes):len(prefixes)+len(date_format)]
        if date_string[date_format.find("HH"):date_format.find("HH")+2] == '' and \
           date_string[date_format.find("MM"):date_format.find("MM")+2] == '' and \
           date_string[date_format.find("SS"):date_format.find("SS")+2] == '':
            # only date is provided
            date_num += [datetime(int(date_string[date_format.find("yyyy"):date_format.find("yyyy")+4]),
                                  int(date_string[date_format.find("mm"):date_format.find("mm")+2]),
                                  int(date_string[date_format.find("dd"):date_format.find("dd")+2]),
                                  )]
        else:
            # time is also provided
            date_num += [datetime(int(date_string[date_format.find("yyyy"):date_format.find("yyyy")+4]),
                                  int(date_string[date_format.find("mm"):date_format.find("mm")+2]),
                                  int(date_string[date_format.find("dd"):date_format.find("dd")+2]),
                                  int(date_string[date_format.find("HH"):date_format.find("HH")+2]),
                                  int(date_string[date_format.find("MM"):date_format.find("MM")+2]),
                                  int(date_string[date_format.find("SS"):date_format.find("SS")+2]),
                                  )]

    # filter list of filemanes that are within the time window in filter_crit
    if len(filter_crit) != 0:
        start = datetime.strptime(filter_crit[0], "%Y_%m_%d__%H_%M_%S")
        end = datetime.strptime(filter_crit[1], "%Y_%m_%d__%H_%M_%S")
        index = [i for i in range(len(date_num)) if date_num[i] >= start and date_num[i] <= end]
        date_num = [date_num[i] for i in index]
        files['name'] = [files['name'][i] for i in index]
        files['path'] = [files['path'][i] for i in index]

    if len(files['name']) != 0:
        # sort them
        files['name'], files['path'] = (list(t) for t in zip(*sorted(zip(files['name'], files['path']))))

        # convert to ordinal
        files['date'] += [datetime.toordinal(date_num[i]) for i in range(len(files['name']))]

        # fill the other fields of the output dict
        files['prefix'] += [prefixes for i in range(len(files['name']))]

    return files
