import h5py
import numpy as np
import pandas as pd
from hdf5obj_2_nparray import hdf5obj_2_nparray
import zipfile
import datetime
from copy import deepcopy
from read_GHG import read_GHG
from prepare_output import prepare_output
from reformat_date import reformat_date


def read_tracer_file(hdf5_f, ini, idx_tracers_to_process):
    tracerdata = {}

    # convert from UNIX time format in nanoseconds to UNIX in seconds
    tracerdata['time'] = pd.Index(np.squeeze(hdf5obj_2_nparray(hdf5_f[ini['tracer']['time_column']], 'f8')) / 1e9)

    tracerdata['mz'] = np.squeeze(hdf5obj_2_nparray(hdf5_f[ini['tracer']['detected_masses_column']][idx_tracers_to_process], 'f8'))
    tracerdata['calibration'] = np.squeeze(hdf5obj_2_nparray(hdf5_f[ini['tracer']['calibration_column']][idx_tracers_to_process], 'f8'))
    tracerdata['transmission'] = np.squeeze(hdf5obj_2_nparray(hdf5_f[ini['tracer']['transmission_column']][idx_tracers_to_process], 'f8'))
    tracerdata['Xr0'] = np.squeeze(hdf5obj_2_nparray(hdf5_f[ini['tracer']['Xr0_column']][idx_tracers_to_process], 'f8'))

    tracerdata['default_CC_kinetic'] = hdf5_f.attrs['default_CC_kinetic']
    tracerdata['cluster_min'] = np.squeeze(hdf5obj_2_nparray(hdf5_f[ini['tracer']['cluster_min_column']][idx_tracers_to_process], 'f8'))
    tracerdata['cluster_max'] = np.squeeze(hdf5obj_2_nparray(hdf5_f[ini['tracer']['cluster_max_column']][idx_tracers_to_process], 'f8'))
    tracerdata['k_reac'] = np.squeeze(hdf5obj_2_nparray(hdf5_f[ini['tracer']['k_reac_column']][idx_tracers_to_process], 'f8'))
    tracerdata['FY'] = np.squeeze(hdf5obj_2_nparray(hdf5_f[ini['tracer']['FY_column']][idx_tracers_to_process], 'f8'))
    tracerdata['IF'] = np.squeeze(hdf5obj_2_nparray(hdf5_f[ini['tracer']['IF_column']][idx_tracers_to_process], 'f8'))

    tracerdata['conc'] = hdf5obj_2_nparray(hdf5_f[ini['tracer']['conc_column']][:, idx_tracers_to_process], 'f8')
    tracerdata['conc_acc'] = hdf5obj_2_nparray(hdf5_f[ini['tracer']['conc_acc_column']][:, idx_tracers_to_process], 'f8')
    tracerdata['conc_prec'] = hdf5obj_2_nparray(hdf5_f[ini['tracer']['conc_prec_column']][:, idx_tracers_to_process], 'f8')
    tracerdata['zero_prec'] = hdf5obj_2_nparray(hdf5_f[ini['tracer']['zero_prec_column']][:, idx_tracers_to_process], 'f8')

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
            with h5py.File(filepath + '\\' + filename, 'r') as hdf5_f:
                sonicdata = hdf5obj_2_nparray(hdf5_f['Data'], 'f8')
                sonicdata = sonicdata[:, ini['sonic']['sonic_columns']]
        elif ini['sonic']['sonic_files_type'] == 'ghg':
            if zipfile.is_zipfile(filepath + '\\' + filename):
                _, sonicdata, _, _ = read_GHG(filepath + '\\' + filename, 'ghg', filepath + r'\unzipped_GHG')
                sonicdata = sonicdata.to_numpy()
                sonicdata = sonicdata[:, ini['sonic']['sonic_columns'] + ini['irga']['irga_columns']]
                # convert from UNIX time format in nanoseconds to UNIX in seconds
                sonicdata[:, 0] += sonicdata[:, 1] / 1e9
                sonicdata = np.delete(sonicdata, 1, axis=1)
            else:
                e = 'ERROR on file: ' + filepath + '\\' + filename + ': cannot be unzipped'
                print(e)
                error_code = 1
            sonicdata = np.array(sonicdata, dtype=float)
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
            tracer_file_path = tracer_files_list['path'][tracer_file_index] + '\\' + tracer_files_list['name'][tracer_file_index]

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
