import numpy as np
import pandas as pd
import re
import sys


def read_metadata_files(path, OF, meteo=False, tilt=False, clock_drift=False, presc_lag=False, rh_lag=False, lpfc=False):
    """
    Read the requested metadata input files (one per call))

    parameters
    ----------
    path (string): path+name of the file
    meteo, tilt, clock_drift, presc_lag, rh_lag, lpfc (boolean or integer): used to identify the requested file

    returns
    -------
    df_xxx (pd dataframe): formatted content of the requested file

    comments
    --------
    Written by B. Heinesch.
    University of Liege, Gembloux Agro-Bio Tech.
    """

    # meteo parameters
    if meteo:

        print('using meteo file: ' + path); OF.write('using meteo file: ' + path + '\n')
        df_meteofiledata = pd.read_csv(path, header=0, names=['pressure', 'temperature', 'relative humidity'], sep=',', skiprows=1)  # read input meteo file
        df_meteofiledata.index = pd.to_datetime(df_meteofiledata.index, format='%d/%m/%Y %H:%M')  # format index as datetime
        df_meteofiledata = df_meteofiledata[~df_meteofiledata.index.duplicated()]

        return df_meteofiledata

    # tilt correction parameters
    if tilt:
        print('using tilt correction file: ' + path); OF.write('using tilt correction file: ' + path + '\n')

        # this is the file coming from the PFM preparation of eddypro
        with open(path, 'r') as file:
            lines = file.readlines()

        # Find the index of the line containing "Rotation matrices"
        rotation_matrices_index = lines.index('Rotation matrices\n') + 1

        R_tilt_PFM = {}
        sect = 1
        for i in range(rotation_matrices_index, len(lines), 4):

            # Split the text using one or more spaces as the delimiter
            elements = re.split(r'\s+', lines[i])
            # Filter out any empty strings resulting from consecutive spaces
            elements = [element for element in elements if element]

            sector = int(elements[4])

            matrix = []
            for j in range(i + 1, i + 4):
                row = list(map(float, lines[j].strip().split()))
                matrix.append(row)
            R_tilt_PFM[sector] = np.array(matrix)

            sect = sect + 1

        return R_tilt_PFM

    # clock-drift lag parameters
    if clock_drift:

        print('using clock_drift file: ' + path); OF.write('using clock_drift file: ' + path + '\n')

        # lag drift info are present and must be accounted for
        df_lag_clock_drift = pd.read_csv(path, header=0, names=['TDC-computer', 'lag drift'], sep=',')  # read input lag drift file
        df_lag_clock_drift.index = pd.to_datetime(df_lag_clock_drift.index, format='%d/%m/%Y %H:%M')  # format index as datetime
        df_lag_clock_drift = df_lag_clock_drift[~df_lag_clock_drift.index.duplicated()]

        return df_lag_clock_drift

    # prescribed time lag (clock-drift + physical)
    if presc_lag:

        if not path:
            sys.exit('LAG_DETECT_METHOD = PRESCRIBED but no lag_prescribed_filepath given')

        print('using presc_lag file: ' + path); OF.write('using presc_lag file: ' + path + '\n')

        # time lag present and must be accounted for
        df_lag_prescribed = pd.read_csv(path, header=0, names=['time lag in s'], sep=',')  # read input lag drift file
        df_lag_prescribed.index = pd.to_datetime(df_lag_prescribed.index, format='%d/%m/%Y %H:%M')  # format index as datetime
        df_lag_prescribed = df_lag_prescribed[~df_lag_prescribed.index.duplicated()]
        df_lag_prescribed = df_lag_prescribed.dropna()

        return df_lag_prescribed

    # time lag rh dependency
    if rh_lag:

        if not path:
            sys.exit('LAG_RH_DEPENDENCY = 1 but no lag_rh_dependency_filepath given')

        print('using rh_lag file: ' + path); OF.write('using rh_lag file: ' + path + '\n')

        # time lag rh dependency present and must be accounted for
        df_lag_rh_dependency = pd.read_csv(path, header=0, sep=',', index_col='RH (%)', skiprows=1)  # read input lag rh dependency file
        df_lag_rh_dependency = df_lag_rh_dependency[~df_lag_rh_dependency.index.duplicated()]
        df_lag_rh_dependency = df_lag_rh_dependency.dropna()
        
        # overwrite each column title by its mz value, rounded at the third decimal
        df_lag_rh_dependency.columns = [
            str(round(float(match.group()), 3)) if (match := re.search(r"[-+]?\d*\.\d+|\d+", col)) else col
            for col in df_lag_rh_dependency.columns
            ]

        return df_lag_rh_dependency

    # low-pass filtering correction parameters
    if lpfc == 1:

        if not path:
            sys.exit('LPFC = 1 but no COF/Massman-type lpfc_filepath given')
            
        print('using lpcf file: ' + path); OF.write('using lpcf file: ' + path + '\n')

        records = []
        stability_class = None
        with open(path, encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                # Detect section label: all,, unstable,, stable,,
                if line.startswith('all') or line.startswith('unstable') or line.startswith('stable'):
                    parts = [p.strip() for p in line.split(',')]
                    if parts[0] in ['all', 'unstable', 'stable']:
                        stability_class = parts[0]
                    continue
                # Skip header line
                if line.startswith('name,value'):
                    continue
                parts = [p.strip() for p in line.split(',')]
                if not stability_class or not parts[0]:
                    continue
                try:
                    value = float(parts[1])
                except (ValueError, IndexError):
                    continue
                records.append({
                    'stability_class': stability_class,
                    'name': parts[0],
                    'value': value
                })
        return pd.DataFrame(records)

    elif lpfc == 2:

        if not path:
            sys.exit('LPFC = 2 but no lpcf_filepath given')

        print('using lpcf file: ' + path); OF.write('using lpcf file: ' + path + '\n')

        dfs = []
        with open(path, encoding='utf-8') as f:
            lines = f.readlines()

        i = 0
        while i < len(lines):
            # Skip empty lines
            if not lines[i].strip():
                i += 1
                continue
            # Section label (unstable/stable)
            if lines[i].startswith('unstable') or lines[i].startswith('stable'):
                status = lines[i].split(',')[0].strip()
                header = lines[i].strip().split(',')
                i += 1
                data = []
                # Read until empty line or line of commas
                while i < len(lines) and lines[i].strip() and not all(x == '' for x in lines[i].strip().split(',')):
                    row = lines[i].strip().split(',')
                    if len(row) == len(header):
                        data.append(row)
                    i += 1
                # Create DataFrame for this section
                df_section = pd.DataFrame(data, columns=header)
                df_section['stability_class'] = status
                dfs.append(df_section)
            i += 1

        # Concatenate, keep only columns of interest
        df_lpfc = pd.concat(dfs, ignore_index=True)
        df_lpfc = df_lpfc[['stability_class', 'ws_max', 'CF_L']]

        # Convert to numeric where possible
        df_lpfc['ws_max'] = pd.to_numeric(df_lpfc['ws_max'], errors='coerce')
        df_lpfc['CF_L'] = pd.to_numeric(df_lpfc['CF_L'], errors='coerce')

        return df_lpfc
