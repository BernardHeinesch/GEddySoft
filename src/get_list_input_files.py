"""Module for locating and listing input files.

This module provides functionality to locate and list sonic anemometer
and tracer gas analyzer data files based on configuration settings.
It supports:

- Multiple input directories
- Various file formats (HDF5, GHG, SLT)
- Date-based file filtering
- Recursive subdirectory scanning for tracer files

Author
------
B. Heinesch
University of Liege, Gembloux Agro-Bio Tech
"""

import os
from list_files_datestring import list_files_datestring


def get_list_input_files(ini):
    """Locate sonic anemometer and tracer gas analyzer data files.

    This function scans specified directories for input data files from
    sonic anemometers and tracer gas analyzers, filtering by file type,
    prefix, and date range.

    Parameters
    ----------
    ini : dict
        Configuration dictionary containing:
        - files : dict
            File search parameters including:
            - sonic_files_folders : list
                Paths to search for sonic data
            - sonic_files_prefix : str
                Prefix for sonic filenames
            - tracer_files_folders : list
                Paths to search for tracer data
            - tracer_files_prefix : str
                Prefix for tracer filenames
            - tracer_files_suffix : str
                Suffix for tracer filenames
            - tracer_files_date_format : str
                Date format in tracer filenames
            - date_files_selection : str
                Date range for file filtering
        - sonic : dict
            Sonic configuration including:
            - sonic_files_type : str
                Type of sonic files ('hdf5', 'ghg', 'slt')

    Returns
    -------
    tuple
        Two-element tuple containing:
        - sonic_files_list : dict
            Dictionary of sonic files with keys:
            - name : list
                Filenames
            - path : list
                Full paths to files
            - prefix : list
                File prefixes
            - date : list
                File dates
        - tracer_files_list : dict
            Dictionary of tracer files with same structure
    """

    # get list of sonic input files in the given folder
    sonic_files_list = {'name': [], 'path': [], 'prefix': [], 'date': []}
    for n in range(len(ini['files']['sonic_files_folders'])):
        if ini['sonic']['sonic_files_type'] == 'hdf5':
            temp = list_files_datestring(ini['files']['sonic_files_folders'][n],
                                         ini['files']['sonic_files_prefix'], '.hdf5', 'yyyy_mm_dd__HH_MM_SS',
                                         ini['files']['date_files_selection'],
                                         sub_dir=False)
        elif ini['sonic']['sonic_files_type'] == 'ghg':
            temp = list_files_datestring(ini['files']['sonic_files_folders'][n],
                                         ini['files']['sonic_files_prefix'], '.ghg', 'yyyy-mm-ddTHHMMSS',
                                         ini['files']['date_files_selection'],
                                         sub_dir=False)
        elif ini['sonic']['sonic_files_type'] == 'slt':
            temp = list_files_datestring(ini['files']['sonic_files_folders'][n],
                                         ini['files']['sonic_files_prefix'], '.slt', 'xxx',
                                         ini['files']['date_files_selection'],
                                         sub_dir=False)

        for key in sonic_files_list:
            sonic_files_list[key] = sonic_files_list[key] + temp[key]

    # get list of tracer input files in the given folder (and its sub-folders)
    tracer_files_list = {'name': [], 'path': [], 'prefix': [], 'date': []}
    if ini['files']['tracer_files_folders']:
        for n in range(len(ini['files']['tracer_files_folders'])):
            temp = list_files_datestring(ini['files']['tracer_files_folders'][n],
                                         ini['files']['tracer_files_prefix'], ini['files']['tracer_files_suffix'],
                                         ini['files']['tracer_files_date_format'],
                                         ini['files']['date_files_selection'],
                                         sub_dir=True)
            for key in tracer_files_list:
                tracer_files_list[key] = tracer_files_list[key] + temp[key]
    else:
        for n in range(len(ini['files']['tracer_files'])):
            pathstr, name = os.path.split(ini['files']['tracer_files'])(n)
            tracer_files_list['name'] += name
            tracer_files_list['path'] += pathstr

    print('found ' + str(len(sonic_files_list['name']))
          + ' sonic files and ' + str(len(tracer_files_list['name'])) + ' tracer files\n')

    return sonic_files_list, tracer_files_list
