import os
from list_files_datestring import list_files_datestring


def get_list_input_files(ini):
    """
        Function to identify input raw files

        parameters
        ----------
        ini: dictionnary with all initialisation information

        returns
        -------
        sonic_files_list: list, sonic files that were find in the folder given in the ini and in its subfolders
        tracer_files_list: list, tracer files that were find in the folder given in the ini and in its subfolders

        comments
        --------
        Written by B. Heinesch.
        University of Liege, Gembloux Agro-Bio Tech.
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
        elif ini['sonic']['sonic_files_type'] == 'csv':
            temp = list_files_datestring(ini['files']['sonic_files_folders'][n],
                                         ini['files']['sonic_files_prefix'], '.csv', 'yyyymmddHHMM',
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
