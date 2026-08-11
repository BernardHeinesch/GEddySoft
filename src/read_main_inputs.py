import h5py
import numpy as np
import pandas as pd
from hdf5obj_2_nparray import hdf5obj_2_nparray
import zipfile
# import datetime
from datetime import datetime, timedelta, timezone
from copy import deepcopy
from read_GHG import read_GHG
from prepare_output import prepare_output
from reformat_date import reformat_date


def read_tracer_file(hdf5_f, ini, idx_tracers_to_process):
    """Read a PTR-TOF-MS tracer HDF5 file into arrays.

    Parameters
    ----------
    hdf5_f : h5py.File
        Open HDF5 file handle.
    ini : dict
        Parsed INI settings.
    idx_tracers_to_process : array_like of int
        Indices of tracers (m/z rows) to extract.

    Returns
    -------
    tracerdata : dict
        Dictionary containing time series and per-tracer metadata.

    Raises
    ------
    ValueError
        If ``ini['param']['TRACER_NORM'] == 1`` and required primary-ion
        information is missing (attributes FPH1/FPH2 or m/z 21.022 and 38.033
        within tolerance).
    """
    tracerdata = {}

    # convert from UNIX time format in nanoseconds to UNIX in seconds
    tracerdata['time'] = pd.Index(np.squeeze(hdf5obj_2_nparray(hdf5_f[ini['tracer']['time_column']], 'f8')) / 1e9)

    tracerdata['mz'] = np.squeeze(hdf5obj_2_nparray(hdf5_f[ini['tracer']['detected_masses_column']][idx_tracers_to_process], 'f8'))
    tracerdata['calibration'] = np.squeeze(hdf5obj_2_nparray(hdf5_f[ini['tracer']['calibration_column']][idx_tracers_to_process], 'f8'))
    tracerdata['transmission'] = np.squeeze(hdf5obj_2_nparray(hdf5_f[ini['tracer']['transmission_column']][idx_tracers_to_process], 'f8'))
    tracerdata['Xr0'] = np.squeeze(hdf5obj_2_nparray(hdf5_f[ini['tracer']['Xr0_column']][idx_tracers_to_process], 'f8'))

    tracerdata['cluster_min'] = np.squeeze(hdf5obj_2_nparray(hdf5_f[ini['tracer']['cluster_min_column']][idx_tracers_to_process], 'f8'))
    tracerdata['cluster_max'] = np.squeeze(hdf5obj_2_nparray(hdf5_f[ini['tracer']['cluster_max_column']][idx_tracers_to_process], 'f8'))

    if 'default_CC_kinetic' in hdf5_f.attrs:
        tracerdata['default_CC_kinetic'] = hdf5_f.attrs['default_CC_kinetic']
    if 'FPH1' in hdf5_f.attrs:
        tracerdata['FPH1'] = hdf5_f.attrs['FPH1']
    if 'FPH2' in hdf5_f.attrs:
        tracerdata['FPH2'] = hdf5_f.attrs['FPH2']
    
    for key in ('k_reac', 'FY', 'IF', 'sensitivity'):
        ini_key = f'{key}_column'
        if ini_key in ini['tracer']:
            column = ini['tracer'][ini_key]
            if column in hdf5_f:
                tracerdata[key] = np.squeeze(
                    hdf5obj_2_nparray(
                        hdf5_f[column][idx_tracers_to_process],
                        'f8'
                    )
                )

    signal = hdf5obj_2_nparray(hdf5_f[ini['tracer']['conc_column']][:, idx_tracers_to_process], 'f8')
    signal_primary_ions = hdf5obj_2_nparray(hdf5_f['Signal_primaryIons'][:, idx_tracers_to_process], 'f8')
    cols_primary = ~np.all(np.isnan(signal_primary_ions), axis=0)
    signal[:, cols_primary] = signal_primary_ions[:, cols_primary]
    tracerdata['conc'] = signal

    signal_prec = hdf5obj_2_nparray(hdf5_f[ini['tracer']['conc_prec_column']][:, idx_tracers_to_process], 'f8')
    signal_primary_ions_prec = hdf5obj_2_nparray(hdf5_f['Signal_primaryIons_prec'][:, idx_tracers_to_process], 'f8')
    cols_primary_prec = ~np.all(np.isnan(signal_primary_ions_prec), axis=0)
    signal_prec[:, cols_primary_prec] = signal_primary_ions_prec[:, cols_primary_prec]
    tracerdata['conc_prec'] = signal_prec

    tracerdata['conc_acc'] = hdf5obj_2_nparray(hdf5_f[ini['tracer']['conc_acc_column']][:, idx_tracers_to_process], 'f8')
    tracerdata['zero_prec'] = hdf5obj_2_nparray(hdf5_f[ini['tracer']['zero_prec_column']][:, idx_tracers_to_process], 'f8')

    if ini['param']['TRACER_NORM'] == 1:
        if ('FPH1' not in tracerdata) or ('FPH2' not in tracerdata):
            raise ValueError(
                "TRACER_NORM=1 requested but required HDF5 root attributes are missing: "
                f"FPH1 present={('FPH1' in tracerdata)}, FPH2 present={('FPH2' in tracerdata)}. "
                f"File: {getattr(hdf5_f, 'filename', 'unknown')}"
            )

        # Primary ion time-series needed for TRACER_NORM (m/z 21.022 & 38.033), read outside idx_tracers_to_process selection.
        mz_21 = 21.022
        mz_38 = 38.033
        mz_tol = 0.001
        mz_all = np.squeeze(hdf5obj_2_nparray(hdf5_f[ini['tracer']['detected_masses_column']][:], 'f8'))
        idx_21 = int(np.nanargmin(np.abs(mz_all - mz_21)))
        idx_38 = int(np.nanargmin(np.abs(mz_all - mz_38)))

        mz_found_21 = float(mz_all[idx_21])
        mz_found_38 = float(mz_all[idx_38])
        dmz_21 = float(abs(mz_found_21 - mz_21))
        dmz_38 = float(abs(mz_found_38 - mz_38))
        if dmz_21 > mz_tol:
            raise ValueError(
                "TRACER_NORM=1 requested but primary ion m/z 21.022 was not found within tolerance. "
                f"Nearest m/z={mz_found_21:.6f}, |dmz|={dmz_21:.6f} (> {mz_tol}). "
                f"File: {getattr(hdf5_f, 'filename', 'unknown')}"
            )
        if dmz_38 > mz_tol:
            raise ValueError(
                "TRACER_NORM=1 requested but primary ion m/z 38.033 was not found within tolerance. "
                f"Nearest m/z={mz_found_38:.6f}, |dmz|={dmz_38:.6f} (> {mz_tol}). "
                f"File: {getattr(hdf5_f, 'filename', 'unknown')}"
            )

        conc_21 = hdf5obj_2_nparray(hdf5_f[ini['tracer']['conc_column']][:, [idx_21]], 'f8')
        conc_21_primary = hdf5obj_2_nparray(hdf5_f['Signal_primaryIons'][:, [idx_21]], 'f8')
        if not np.all(np.isnan(conc_21_primary)):
            conc_21 = conc_21_primary
        tracerdata['conc_primary_21_022'] = np.squeeze(conc_21)

        conc_38 = hdf5obj_2_nparray(hdf5_f[ini['tracer']['conc_column']][:, [idx_38]], 'f8')
        conc_38_primary = hdf5obj_2_nparray(hdf5_f['Signal_primaryIons'][:, [idx_38]], 'f8')
        if not np.all(np.isnan(conc_38_primary)):
            conc_38 = conc_38_primary
        tracerdata['conc_primary_38_033'] = np.squeeze(conc_38)

    return tracerdata


def read_main_inputs(filepath, filename, filetype, ini, OF, idx_tracers_to_process=None, tracer_files_list=None, results=None, cov_data=None, out_len=None):
    """
    Reads the sonic, GHG or the tracer input files. The positional args are needed only for reading the tracer files.

    parameters
    ----------

    filetype: str, gives the type of input file that should be read ('sonic' or 'tracer')
    ini: dict, initialisation information
    OF: obj, logfile
    idx_tracers_to_process: np arr of int, indices of selected mz, updated only if first reading of a tracer file for the running day
    tracer_files_list: dict, list of the tracer files to be processed
    results: dict, info to be outputed in the hdf5 files, updated only if first reading of a tracer file for the running day
    cov_data: dict, covariance function (cov in function of the lag), updated only if first reading of a tracer file for the running day
    out_len: int, number of sonic files to be processed

    returns
    -------
    sonicdata: np array of floats, input data from the sonic and the IRGA
                The structure is the following:
                    np arrays of float64, size (x,y), with
                    x the number of records
                    y the variables time in UNIX format (UTC), u, v, w, T, flag
    error_code: int, error flag for the reading of the input file (0 if failed, 1 if successful)
    tracerdata: dict, input data from the PTR-TOF-MS
                The structure is the following:
                    - time: np arrays of float64, size (x,), time in UNIX format (UTC)
                    - mz: np arrays of float64, size (y,), mz values
                    - conc: np arrays of float64, size (x,y), concentrations
                    - calibration: np arrays of float64, size (y,), calibration coefficients for concentrations
                    - transmission: np arrays of float64, size (y,), transmission coefficients for concentrations
    idx_tracers_to_process: np arr of int, indices of selected mz
    tracer_file_index: int, index of the tracer file to be read in the file list
    results: dict, info to be outputed in the hdf5 files, updated only if first reading of a tracer file for the running day
    cov_data: dict, covariance function (cov in function of the lag), updated only if first reading of a tracer file for the running day

    Comments
    --------
    input file formats:
        
    wind sonic data format (hdf5) is 
        col 0: timestamp (UNIX format, in nanoseconds)
        col 1: U wind component (m s-1)
        col 2: V wind component (m s-1)
        col 3: W wind component (m s-1)
        col 4: sonic temperature (deg C)
        col 5: status adress (see HS50 manual)
        col 6: status data (see HS50 manual)
        col 5: error (see HS50 manual)

    GHG data format is the eddypro™ output
    
    tracer data format (hdf5) is (see field attributes in the file for more info)
        - Time-series datasets:
            * time: Unix timestamp array [UTC]
            * Signal: concentration data
            * concentration_accuracy, concentration_precision: Data quality metrics
            * Signal_primaryIons, Signal_primaryIons_prec: Primary ion measurements
            * zero, zero_precision: Zero calibration data
            
        - Calibration/Configuration datasets:
            * mz: Mass-to-charge ratios
            * calibration: Calibration coefficients
            * k_reac: Reaction rate coefficients
            * transmission: Transmission factors
            * FY, IF: Fragmentation/Ion data
            * cluster_min, cluster_max: Cluster bounds
            * Xr0: Reference mixing ratios
            * MeasurementsTimeInterval: Time intervals
        
        - Key attributes:
            * PAP processing version
            * Timezone info for filename (UTC+01:00)
            * Configuration parameters for:
                - Zero calibration intervals
                - Processing intervals
                - Calibration settings
                - Data averaging periods

    Written by B. Heinesch.
    University of Liege, Gembloux Agro-Bio Tech.
    """
    error_code = 0

    if filetype == 'sonic':
        # load sonic file in sonicdata np
        msg = 'sonic  file ' + filename
        print(msg); OF.write(msg + "\n")
        if ini['sonic']['sonic_files_type'] == 'hdf5':
            with h5py.File(filepath + '/' + filename, 'r') as hdf5_f:
                sonicdata = hdf5obj_2_nparray(hdf5_f['Data'], 'f8')
                sonicdata = sonicdata[:, ini['sonic']['sonic_columns']]
        elif ini['sonic']['sonic_files_type'] == 'ghg':
            if zipfile.is_zipfile(filepath + '/' + filename):
                _, sonicdata, _, _ = read_GHG(filepath + '/' + filename, 'ghg', filepath + r'/unzipped_GHG')
                sonicdata = sonicdata.to_numpy()
                sonicdata = sonicdata[:, ini['sonic']['sonic_columns'] + ini['irga']['irga_columns']]
                # convert from UNIX time format in nanoseconds to UNIX in seconds
                sonicdata[:, 0] += sonicdata[:, 1] / 1e9
                sonicdata = np.delete(sonicdata, 1, axis=1)
            else:
                e = 'ERROR on file: ' + filepath + '/' + filename + ': cannot be unzipped'
                print(e)
                error_code = 1
            sonicdata = np.array(sonicdata, dtype=float)
        elif ini['sonic']['sonic_files_type'] == 'csv':
            sonicdata = np.loadtxt(filepath + '/' + filename, delimiter=',', skiprows=1)

            sonicdata[sonicdata == -9999] = np.nan

            # Filter out negative values from both sonic_columns and irga_columns and slice the sonicdata using the valid indices
            sonic_columns = [col for col in ini['sonic']['sonic_columns'] if col >= 0]
            irga_columns = [col for col in ini['irga']['irga_columns'] if col >= 0]
            valid_columns = sonic_columns + irga_columns
            sonicdata = sonicdata[:, valid_columns]

            # add time based on filename
            date_format='yyyymmddHHMM'
            date_string = filename[len(ini['files']['sonic_files_prefix']):len(ini['files']['sonic_files_prefix'])+len(date_format)]
            date_num = datetime(int(date_string[date_format.find("yyyy"):date_format.find("yyyy")+4]),
                                          int(date_string[date_format.find("mm"):date_format.find("mm")+2]),
                                          int(date_string[date_format.find("dd"):date_format.find("dd")+2]),
                                          int(date_string[date_format.find("HH"):date_format.find("HH")+2]),
                                          int(date_string[date_format.find("MM"):date_format.find("MM")+2]),
                                          int(date_string[date_format.find("SS"):date_format.find("SS")+2] if "SS" in date_format else "00"),
                                          )
            
            # Convert date_num (UTC+1) to UTC (subtract 1 hour)
            start_datetime_utc = date_num - timedelta(hours=1)

            # Convert start_datetime_utc to Unix timestamp (in seconds)
            start_datetime_utc = start_datetime_utc.replace(tzinfo=timezone.utc)
            start_timestamp_seconds = int(start_datetime_utc.timestamp())
            start_timestamp_seconds = int(start_datetime_utc.timestamp())
            
            # Generate Unix timestamps for each data point using the SAMPLING_RATE_SONIC
            num_samples = sonicdata.shape[0]
            time_diff_seconds = 1 / ini['param']['SAMPLING_RATE_SONIC']

            unix_timestamps = np.zeros(num_samples, dtype=np.float64)
            for i in range(num_samples):
                unix_timestamps[i] = start_timestamp_seconds + i * time_diff_seconds
            
            # Pre-allocate the final sonicdata array with 3 additional columns (1 for timestamp, 2 for quality flags)
            final_sonicdata = np.full((num_samples, sonicdata.shape[1] + 3), np.nan)
            
            # Insert Unix timestamps in the first column, former sonicdata and the flag columns
            final_sonicdata[:, 0] = unix_timestamps
            final_sonicdata[:, 1:5] = sonicdata[:, 0:4]
            final_sonicdata[:, 5] = 0
            final_sonicdata[:, 6:8] = sonicdata[:, 4:6]
            final_sonicdata[:, 8] = 0

            sonicdata = final_sonicdata

        return (sonicdata, error_code)

    if filetype == 'tracer':

        # find file in tracer_file_list with the corresponding timestamp
        current_timestamp = reformat_date(filename[-25:-5], 'yyyy_mm_dd__HH_MM_SS', ini['files']['tracer_files_date_format'])
        # remove seconds in order to allow processing of tracer file that do not start exactly at the half-hour
        current_timestamp = current_timestamp[:-4]
        tracer_file_index = None

        if tracer_files_list:
            if any(current_timestamp in name for name in tracer_files_list['name']):
                tracer_file_index = next(idx for idx, name in enumerate(tracer_files_list['name']) if current_timestamp in name)
                msg = 'tracer file ' + tracer_files_list['name'][tracer_file_index]
                print(msg); OF.write(msg + "\n")
            else:
                msg = 'corresponding tracer file not found\n'
                error_code = 1
                print(msg); OF.write(msg + "\n")
                return (None, error_code, idx_tracers_to_process, None)

        if error_code == 0:  # corresponding tracer file was found
            tracer_file_path = tracer_files_list['path'][tracer_file_index] + '/' + tracer_files_list['name'][tracer_file_index]

            if 'TRACER' not in results:
                # get hdf5_nb_tracers, hdf5_mz_tracer and set idx_calibration
                with h5py.File(tracer_file_path, 'r') as hdf5_f:
                    ds = hdf5_f[ini['tracer']['detected_masses_column']]
                    hdf5_nb_tracers = ds.shape[0]
                    hdf5_mz_tracer = np.array(ds[:])

                    ds = hdf5_f[ini['tracer']['calibration_column']]
                    idx_calibration = np.where(~np.isnan(ds[:]))[0]

                if len(ini['tracer']['tracer_mz']) == 0:
                    # process all channels
                    idx_tracers_to_process = list(range(hdf5_nb_tracers))
                else:
                    # process channels having the closest mz to the list proposed in the ini
                    idx_tracers_to_process = [min(range(len(hdf5_mz_tracer)), key=lambda i: abs(hdf5_mz_tracer[i] - x)) for x in ini['tracer']['tracer_mz']]

                # append TRACER part to results dict
                results = prepare_output(ini, out_len, hdf5_nb_tracers=len(idx_tracers_to_process), hdf5_mz_tracer=hdf5_mz_tracer[idx_tracers_to_process], results=results)

                temp = {'cov': [[np.NaN] * (2 * ini['param']['LAG_OUTER_WINDOW_SIZE'] + 1)] * out_len}
                cov_data['TRACER'] = dict(zip(map(str, range(len(idx_tracers_to_process))),
                                              [deepcopy(temp) for _ in range(len(idx_tracers_to_process))]))

            # load corresponding tracer file
            with h5py.File(tracer_file_path, 'r') as hdf5_f:
                tracerdata = read_tracer_file(hdf5_f, ini, idx_tracers_to_process)

            if 'TRACER' not in results:
                return (tracerdata, error_code, idx_tracers_to_process, tracer_file_index, results, cov_data)
            else:
                return (tracerdata, error_code, idx_tracers_to_process, tracer_file_index)
