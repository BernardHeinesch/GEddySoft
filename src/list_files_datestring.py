from datetime import datetime
import os


def list_files_datestring(folder_path, prefixes, suffix, date_format, filter_crit='', sub_dir=True):
    """
    lists all files in a folder and sub-folders with filenames starting with
    one of the given prefixes followed by a date string of the format given
    by date_format (e.g. "yyyy_mm_dd"), within the given time window and
    followed by a given suffix.
    Additionally sorts files according to the date string.

    parameters
    ----------
    folder_path: string
    prefixes: string
    suffix: string
    date_format: str, date format used (e.g. 'yyyy_mm_dd_HH_MM_SS')
    filter_crit: tuple of strings with the start date and the end date in a text format
    sub_dir: boolean, if True, consider also sub-folders, default is True

    returns
    -------
    files: dict with keys: 'name', 'path', 'prefix', 'date'

    comments
    --------
    Written by B. Heinesch
    University of Liege, Gembloux Agro-Bio Tech.
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

            # Check if SS exists in the date format, and set it to 00 if not found
            date_num += [datetime(int(date_string[date_format.find("yyyy"):date_format.find("yyyy")+4]),
                                  int(date_string[date_format.find("mm"):date_format.find("mm")+2]),
                                  int(date_string[date_format.find("dd"):date_format.find("dd")+2]),
                                  int(date_string[date_format.find("HH"):date_format.find("HH")+2]),
                                  int(date_string[date_format.find("MM"):date_format.find("MM")+2]),
                                  int(date_string[date_format.find("SS"):date_format.find("SS")+2] if "SS" in date_format else "00"),
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
