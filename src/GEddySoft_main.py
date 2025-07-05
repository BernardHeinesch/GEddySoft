"""
GEddySoft: A Python Package for Eddy Covariance Processing

This module contains the main processing function for the GEddySoft eddy covariance
software package. It handles the complete workflow of processing raw high-frequency
data from sonic anemometers and trace gas analyzers (e.g., PTR-TOF-MS) to calculate
fluxes and quality control metrics.

Key Features
------------
* Flexible input handling for sonic and trace gas data
* Comprehensive quality control including:

  * Spike detection (Vickers & Mahrt 1997)
  * Stationarity tests (Foken & Wichura 1996)
  * Integral turbulence characteristics (Thomas & Foken 2002)
  * Instrument diagnostics

* Advanced processing options:

  * Multiple coordinate rotation methods
  * Time lag optimization
  * Spectral corrections
  * Flux uncertainty estimation

* Parallel processing support for large datasets
* HDF5-based data storage with rich metadata

Author
------
Bernard Heinesch
University of Liège, Gembloux Agro-Bio Tech

License
-------
GEddySoft is available under the terms of the MIT License.
"""

# %% Package and functions import

import numpy as np
import pandas as pd
import datetime
import os
import sys
import scipy.signal
from scipy.interpolate import interp1d
from copy import deepcopy
from print_progress import print_progress
import hdfdict  # not initially available in the anaconda distribution
import re
import h5py

# import user functions/modules
from read_GHG import read_diag_val
from read_main_inputs import read_main_inputs
from read_metadata_files import read_metadata_files
from check_raw_timestamps import check_raw_timestamps
from get_list_input_files import get_list_input_files
from prepare_output import prepare_output
from get_closest_value import get_closest_value
from logBinSpectrum import logBinSpectrum
from nanlinfit import nanlinfit
from nandetrend import nandetrend
from xcov import xcov
from find_covariance_peak import find_covariance_peak
from test_steady_state_FW96 import test_steady_state_FW96
from test_steady_state_M98 import test_steady_state_M98
from test_steady_state_D99 import test_steady_state_D99
from test_ITC import test_ITC
from compute_wind_direction import compute_wind_direction
from cospectrum import cospectrum
from inst_prob_test import inst_prob_test
from add_attributes import add_attributes
from flux_uncertainties import flux_uncertainties
from clean_results import clean_results
from correction_factor_lpf import correction_factor_lpf
from spike_detection_vickers97 import spike_detection_vickers97
from map_sonic2tracer import map_sonic2tracer
from wind_rotations import wind_rotations
from flux_unit_conversion import flux_unit_conversion

# %% read inputs and prepare the process of individual days


def GEddySoft_main(str_day, ini, log_filename):
    """
    Process eddy covariance data for a single day (when used in multithread mode)
    or for the time period given in the ini file.

    This function handles the complete processing chain for eddy covariance data,
    from raw high-frequency measurements to quality-controlled fluxes. It supports
    both single-threaded and parallel processing modes.

    Parameters
    ----------
    str_day : str
        Date string in format 'YYYYMMDD' for the day to process.
        Use 'no_multithread' for single-threaded mode.
    ini : dict
        Configuration dictionary containing all processing parameters:

        * files: Input/output paths and file patterns
        * param: Processing parameters (QC thresholds, methods, etc.)
        * run_param: Runtime parameters
        * sonic: Sonic anemometer configuration
        * irga: Gas analyzer configuration

    log_filename : str
        Path to the log file where processing messages will be written

    Returns
    -------
    unique_days : list
        List of successfully processed days

    Notes
    -----
    The processing workflow includes:

    1. Data loading and synchronization

       * Sonic anemometer data
       * Trace gas measurements
       * Auxiliary meteorological data

    2. Quality control

       * Spike detection and removal
       * Diagnostic value checking
       * Missing data handling

    3. Coordinate rotations

       * Double rotation or
       * Planar fit method

    4. Time lag optimization

       * Covariance maximization
       * RH-dependent correction

    5. Flux calculations

       * Detrending options
       * Webb-Pearman-Leuning terms

    6. Quality assessment

       * Stationarity tests
       * Integral turbulence characteristics
       * Spectral analysis

    7. Uncertainty estimation

       * Random error (Finkelstein & Sims)
       * Detection limits

    The function creates two types of output files:

    1. Results file (HDF5)

       * Fluxes and means
       * Quality flags
       * Processing parameters

    2. Covariance file (HDF5, optional)

       * Raw covariance functions
       * Time lag information

    Examples
    --------
    >>> # Process a single day
    >>> ini = read_main_inputs('config.txt')
    >>> log_file = 'processing.log'
    >>> processed = GEddySoft_main('20240101', ini, log_file)
    """

    # if multi-processing, mute the console and overwrite the file selection of the ini file
    if str_day != 'no_multithread':
        # Save the current stdout
        original_stdout = sys.stdout
        # Redirect stdout to devnull
        sys.stdout = open(os.devnull, 'w')

        ini['files']['date_files_selection'] = "('{}__00_00_00','{}__23_30_00')".format(str_day, str_day)

    # Open logfile in append mode (already opened but see no other solution when multiprocessing)
    OF = open(log_filename, 'a')

    # load metadata files if provided
    if ini['files']['meteo_filepath']:
        df_meteofiledata = read_metadata_files(ini['files']['meteo_filepath'], OF, meteo=True)

    if ini['param']['TILT_CORRECTION_MODE'] == 2:
        R_tilt_PFM = read_metadata_files(ini['files']['tilt_correction_filepath'], OF, tilt=True)

    if ini['files']['lag_clock_drift_filepath']:
        # lag drift info are present and must be accounted for
        df_lag_clock_drift = read_metadata_files(ini['files']['lag_clock_drift_filepath'], OF, clock_drift=True)

    if ini['param']['LAG_DETECT_METHOD'] == 'PRESCRIBED':
        df_lag_prescribed = read_metadata_files(ini['files']['lag_prescribed_filepath'], OF, presc_lag=True)

    if ini['param']['LAG_RH_DEPENDENCY'] == 1:
        df_lag_rh_dependency = read_metadata_files(ini['files']['lag_rh_dependency_filepath'], OF, rh_lag=True)



    if ini['param']['LPFC'] != 0:
        df_lpfc = read_metadata_files(ini['files']['lpfc_filepath'], OF, lpfc=ini['param']['LPFC'])

    print(); OF.write('\n')
    
    # check output folder
    if not os.path.exists((ini['files']['output_folder'])):
        sys.exit('output folder not found: ' + ini['files']['output_folder'])
    # create sub-folder for covariance files if not existing
    cov_folder = os.path.join(ini['files']['output_folder'], 'cov')
    if not os.path.exists(cov_folder):
        os.makedirs(cov_folder)
        print(f"Created missing folder: {cov_folder}\n"); OF.write(f"Created missing folder: {cov_folder}\n\n")
    
    # get list of sonic and tracer input files (in the given dir for sonic and also in subdir for tracer)
    all_sonic_files_list, all_tracer_files_list = get_list_input_files(ini)

    # map sonic file list to tracer file list if asked for
    if ini['run_param']['MAP_SONIC2TRACER'] and ini['run_param']['CONCENTRATION_TYPE']:
        all_sonic_files_list = map_sonic2tracer(all_sonic_files_list, all_tracer_files_list)
        print('sonic files mapped to tracer files\n')
    
    print('*************** processing data ***************' + '\n')
    OF.write('*************** processing data ***************' + '\n' + '\n')

    # get uniques yyyy_mm_dd in the sonic file list
    unique_days = list(set(list(map(lambda x: x[len(ini['files']['sonic_files_prefix']):len(ini['files']['sonic_files_prefix']) + 10], all_sonic_files_list['name']))))
    unique_days.sort()

    idx_tracers_to_process = None

# %% loop on days

    for d in range(len(unique_days)):

        process_irga_data_day = False

        # select only the given day

        # for sonic files
        condition = lambda e: unique_days[d] in e
        indices = [i for i, elem in enumerate(all_sonic_files_list['name']) if condition(elem)]
        sonic_files_list = {'name': [], 'path': [], 'prefix': [], 'date': []}
        for key in all_sonic_files_list:
            sonic_files_list[key] = [all_sonic_files_list[key][i] for i in indices]

        if ini['run_param']['CONCENTRATION_TYPE'] == 1:
            # and for tracer files
            # Extract year, month, day from the input_date
            year, month, day = unique_days[d].split('_')
            # Define a pattern to identify time-related placeholders
            time_placeholders = re.compile(r'HH|MM|SS')
            # Remove time-related placeholders from the user_format
            cleaned_format = time_placeholders.sub('', ini['files']['tracer_files_date_format'])
            # Replace placeholders for date components in the cleaned format
            formatted_date = cleaned_format.replace('yyyy', year).replace('mm', month).replace('dd', day)
            # Remove any leftover separators
            formatted_date = re.sub(r'[^a-zA-Z0-9]+$', '', formatted_date)
    
            condition = lambda e: formatted_date in e
            indices = [i for i, elem in enumerate(all_tracer_files_list['name']) if condition(elem)]
    
            tracer_files_list = {'name': [], 'path': [], 'prefix': [], 'date': []}
            for key in tracer_files_list:
                tracer_files_list[key] = [all_tracer_files_list[key][i] for i in indices]

        # prepare some variables
        out_len = len(sonic_files_list['name'])
        upsampling_factor = int(ini['param']['SAMPLING_RATE_SONIC']/ini['param']['SAMPLING_RATE_TRACER'])

        # prepare results dict structure
        results = prepare_output(ini, out_len)

        # prepare covariance data dict to be stored to file
        cov_data = {}
        temp = {'cov': [[np.NaN] * (2*ini['param']['LAG_OUTER_WINDOW_SIZE']+1)] * out_len}
        cov_data['IRGA'] = dict(zip(tuple(str(x) for x in tuple(range(0, len(ini['irga']['irga_columns'])))),
                                    [deepcopy(temp) for i in range(len(ini['irga']['irga_columns']))]))  # struct containing results for each tracer

        # prepare frequency axis for cospectrum
        f = np.arange(0, np.ceil(ini['param']['WINDOW_LENGTH']*ini['param']['SAMPLING_RATE_FINAL']/2)+1) / ini['param']['WINDOW_LENGTH']  # linear frequency
        results['freq'] = logBinSpectrum(f, np.zeros(len(f)), ini['param']['NUM_FREQ_BINS'], f[1], f[-1])[0]  # log-spaced frequency bins

# %% process sonic data

        # Loop on sonic files
        for n in range(len(sonic_files_list['name'])):

            # store current time
            proc_start_time_hh = datetime.datetime.now()

            # load sonic file in sonicdata np
            sonicdata, error_code = read_main_inputs(sonic_files_list['path'][n], sonic_files_list['name'][n], 'sonic', ini, OF)
            if error_code: continue

            # store timestamp
            # Filename is in UTC+1, start of the half-hour
            # Timestamp in the input file is in UNIX format (number of seconds since Epoch, 1/1/1970 00:00:00 UTC).
            # The instruction utcfromtimestamp converts time in the UTC timezone
            # +3600 : to convert to UTC +1
            first_ts = round(sonicdata[0, 0])
            last_ts = round(next((i for i in reversed(sonicdata[:, 0]) if i != 0), None))
            results['time'][n] = datetime.datetime.utcfromtimestamp(last_ts + 3600)

            # remove sonicdata trailing zeros (the file is initially oversized and filled with zero values)
            sonicdata = sonicdata[:len(sonicdata[sonicdata[:, 0] != 0, :]), :]

            # check completeness of sonic data
            completeness_sonicdata = sum(np.isfinite(sonicdata[:, 1])) / \
                                     (ini['param']['WINDOW_LENGTH']*ini['param']['SAMPLING_RATE_SONIC'])

            # if sonic file not complete enough, skip this file
            if completeness_sonicdata < ini['param']['COMPLETENESS_THRESHOLD']:
                print('sonic file incomplete (' + "{:.2f}".format(completeness_sonicdata) + '%), sonic  process skipped')
                continue

            # reduce sonic file to expected lenght
            # if the sonic output frequency is slightly higher that it's nominal one, "trailing" records will be discarded
            # acceptable only is these nb of records are very limited (e.g. GHS50 files have between 90000 and 90005 records instead of 90000)
            # the case of the sonic output frequency slightly lower that it's nominal one is not considered
            sonicdata = sonicdata[:ini['param']['WINDOW_LENGTH']*ini['param']['SAMPLING_RATE_SONIC'], :]

            # check validity of timestamp series, if present
            time_series = pd.Series(np.vectorize(datetime.datetime.utcfromtimestamp)(sonicdata[:, 0]))
            threshold = pd.Timedelta(milliseconds=3/ini['param']['SAMPLING_RATE_SONIC']*1000)
            msg = check_raw_timestamps(time_series, threshold)
            if msg[0]:
                print('sonic file: ' + msg[0]); OF.write('sonic file: ' + msg[0] + '\n')

            # check sonic flag and set flagged data to default value
            # if SPIKE_MODE = 2 and not more than 3 consecutive raw data with error code, then handled by the spike replacement
            if np.any(sonicdata[:, 5] != 0):
                print('Warning: detected sonic data with error code activated')
                sonicdata[sonicdata[:, 5] != 0, 1:5] = -999

            # check IRGA flag and set flagged data to NaN
            if (ini['run_param']['CONCENTRATION_TYPE'] == 0 and ini['irga']['irga_columns'][2] > 0):
                diag_count = read_diag_val(sonicdata[:, 8])
                if (diag_count.drop(columns=['AGC', 'Sync', 'Aux_input']) != 0).any(axis=1).any():
                    sonicdata[:, 6:8] = np.NaN
                    print('IRGA data flagged')

            # subsample sonic data from SAMPLING_RATE_SONIC to SAMPLING RATE FINAL if needed
            # note : InnFlux uses upfirdn routine instead of resample
            if ((ini['param']['SAMPLING_RATE_SONIC'] > ini['param']['SAMPLING_RATE_FINAL'])):

                if np.sum(np.isnan(sonicdata)):
                    print(str(np.sum(np.isnan(sonicdata))) + ' Nans have been detected in this sonic file. Resampling will not work.')  # debug purposes, see readme for nans and resample

                downsampling_factor = ini['param']['SAMPLING_RATE_SONIC']/ini['param']['SAMPLING_RATE_FINAL']
                sonicdata = scipy.signal.resample(sonicdata, int(len(sonicdata[:, 0])/downsampling_factor))

                # reconstruct the timestamp colum (used later on by interp1d).
                # also allows to correct for possible variable transmission delays between the sonic and the final
                # storage on the datalogger/computer, assuming that the sonic is sending data at a fixed frequency
                sonicdata[:, 0] = (first_ts +
                                   np.asarray([x for x in range(0, len(sonicdata))])/ini['param']['SAMPLING_RATE_FINAL'])  # first column contains timestamp

            # detect w and T spikes (and remove if requested)
            if ini['param']['SPIKE_MODE'] != 0:
                sonicdata[:, 3], is_spike = spike_detection_vickers97(
                                             sonicdata[:, 3],
                                             spike_mode = ini['param']['SPIKE_MODE'],
                                             avrg_len = int(ini['param']['WINDOW_LENGTH']/60),
                                             ac_freq = ini['param']['SAMPLING_RATE_SONIC'],
                                             spike_limit = ini['param']['SPIKE_DETECTION_THRESHOLD_W'],
                                             max_consec_spikes = ini['param']['MAX_CONSEC_SPIKES'],
                                             ctrplot = ini['run_param']['PLOT_SPIKE_DETECTION']
                                             )
                num_spikes_w = int(sum(is_spike))
            else:
                num_spikes_w = np.nan

            if ini['param']['SPIKE_MODE'] != 0:
                sonicdata[:, 4], is_spike = spike_detection_vickers97(sonicdata[:, 4], 
                                             spike_mode = ini['param']['SPIKE_MODE'],
                                             avrg_len = int(ini['param']['WINDOW_LENGTH']/60),
                                             ac_freq = ini['param']['SAMPLING_RATE_SONIC'],
                                             spike_limit = ini['param']['SPIKE_DETECTION_THRESHOLD_T'],
                                             max_consec_spikes = ini['param']['MAX_CONSEC_SPIKES'],
                                             ctrplot = ini['run_param']['PLOT_SPIKE_DETECTION']
                                             )
                num_spikes_T = int(sum(is_spike))
            else:
                num_spikes_T = np.nan

            # test on instrument malfunctions from Vitale 2020
            IPT_w = inst_prob_test(sonicdata[:, 3], False, ini['param']['SAMPLING_RATE_FINAL'], ini['run_param']['PLOT_INST_PROB_TEST'], 'W')
            IPT_T = inst_prob_test(sonicdata[:, 4], ini['param']['DETREND_TEMPERATURE_SIGNAL'], ini['param']['SAMPLING_RATE_FINAL'], ini['run_param']['PLOT_INST_PROB_TEST'], 'T_SONIC')

            # mean and standard deviation of horizontal wind speed and direction
            # note: wind direction is the angle where the wind is coming from, with north = 0 and east = 90 degrees
            mean_uvw = np.nanmean(sonicdata[:, 1:3+1], 0)
            mean_wind_speed = np.sqrt(mean_uvw[0]**2 + mean_uvw[1]**2 + mean_uvw[2]**2)
            std_wind_speed = np.std(np.sqrt(sonicdata[:, 1]**2 + sonicdata[:, 2]**2 + sonicdata[:, 3]**2))
            mean_wind_dir = compute_wind_direction(mean_uvw[0], mean_uvw[1], ini['param']['SONIC_TYPE'], ini['param']['SONIC_ORIENTATION'])

            # Compute wind directions using vectorized operations
            wind_directions = compute_wind_direction(sonicdata[:, 1],
                                                     sonicdata[:, 2],
                                                     ini['param']['SONIC_TYPE'], ini['param']['SONIC_ORIENTATION'])

            # Compute standard deviation
            std_wind_dir = np.std(wind_directions)

            # Coordinate rotations
            if (ini['param']['TILT_CORRECTION_MODE'] != 2):
                uvw_rot, rot_phi, R_tilt = wind_rotations(ini, sonicdata, mean_uvw)
            else:
                uvw_rot, rot_phi, R_tilt = wind_rotations(ini, sonicdata, mean_uvw, R_tilt_PFM, mean_wind_dir)
            mean_uvw_rot = np.nanmean(uvw_rot, axis=0)

            # create a common time vector for aligning sonic and irga/tracer data
            interp_time = sonicdata[0, 0] + np.transpose(np.arange(0, ini['param']['WINDOW_LENGTH'], 1/ini['param']['SAMPLING_RATE_FINAL']))

            # interpolate wind and scalar data (keep NaN and extrapolate using NaN)
            interp_u = interp1d(sonicdata[:, 0], uvw_rot[:, 0], kind='nearest', bounds_error=False, fill_value=np.NaN)(interp_time)
            interp_v = interp1d(sonicdata[:, 0], uvw_rot[:, 1], kind='nearest', bounds_error=False, fill_value=np.NaN)(interp_time)
            interp_w = interp1d(sonicdata[:, 0], uvw_rot[:, 2], kind='nearest', bounds_error=False, fill_value=np.NaN)(interp_time)
            interp_Tsonic = interp1d(sonicdata[:, 0], sonicdata[:, 4], kind='nearest', bounds_error=False, fill_value=np.NaN)(interp_time) + 273.15

            # check data completeness
            completeness_T = sum(np.isfinite(interp_Tsonic))/(ini['param']['WINDOW_LENGTH']*ini['param']['SAMPLING_RATE_FINAL'])

            # Reynolds decomposition of wind components (e.g. u = u_mean + u_prime)
            u_prime = interp_u - np.nanmean(interp_u)
            v_prime = interp_v - np.nanmean(interp_v)
            w_prime = interp_w - np.nanmean(interp_w)

            # standard deviation of wind
            std_uvw = [np.nanstd(u_prime, ddof=1), np.nanstd(v_prime, ddof=1), np.nanstd(w_prime, ddof=1)]

            # check for Nans and replace NaN by zero
            if np.isnan(u_prime).any():
                print("Warning: NaN values detected in u_prime. Replacing with 0.")
                u_prime = np.nan_to_num(u_prime, nan=0)
            if np.isnan(v_prime).any():
                print("Warning: NaN values detected in v_prime. Replacing with 0.")
                v_prime = np.nan_to_num(v_prime, nan=0)
            if np.isnan(w_prime).any():
                print("Warning: NaN values detected in w_prime. Replacing with 0.")
                w_prime = np.nan_to_num(w_prime, nan=0)

            # calculate covariances of wind components
            uw = np.cov(u_prime, w_prime)[0][1]
            vw = np.cov(v_prime, w_prime)[0][1]
            uv = np.cov(u_prime, v_prime)[0][1]
            uu = np.cov(u_prime, u_prime)[0][1]
            vv = np.cov(v_prime, v_prime)[0][1]
            ww = np.cov(w_prime, w_prime)[0][1]

            # calculate friction velocity, u*
            u_star = (uw**2 + vw**2)**0.25

            # calculate temperature T from virtual sonic temperature T_sonic and water vapor concentration
            # T_sonic = T*(1 + 0.32*w), were w is the volume mixing ratio of water vapor in air
            if ini['run_param']['CONCENTRATION_TYPE'] == 0:
                 interp_H2O = interp1d(sonicdata[:,0], sonicdata[:, 4], kind='nearest', bounds_error=False, fill_value=np.NaN)(interp_time)
                 completeness_H2O = sum(np.isfinite(interp_H2O))/(ini['param']['WINDOW_LENGTH']*ini['param']['SAMPLING_RATE_FINAL'])
            else:
                completeness_H2O = 0

            if ini['run_param']['CONCENTRATION_TYPE'] == 0 and completeness_H2O >= ini['param']['COMPLETENESS_THRESHOLD']:

                 if (completeness_T >= ini['param']['COMPLETENESS_THRESHOLD'] and completeness_H2O >= ini['param']['COMPLETENESS_THRESHOLD']):
                    interp_T = interp_Tsonic/(1 + 0.32*1e-3*interp_H2O) # note that H2O data is given in permille
                    T_mean = np.nanmean(interp_T)
                    dT_dt = nanlinfit(interp_T)[0]
                    dT_dt = dT_dt*ini['param']['SAMPLING_RATE_FINAL']
                    if (ini['param']['DETREND_TEMPERATURE_SIGNAL']):
                        T_prime = nandetrend(interp_T)
                    else:
                        T_prime = interp_T - T_mean
                    std_T = np.nanstd(T_prime, ddof=1)
                    T_prime = np.nan_to_num(T_prime, 0)

            else:

                # H2O concentration not available, use T_sonic
                interp_T = interp_Tsonic
                T_mean = np.nanmean(interp_T)
                dT_dt = nanlinfit(interp_T)[0]
                dT_dt = dT_dt*ini['param']['SAMPLING_RATE_FINAL']
                if (ini['param']['DETREND_TEMPERATURE_SIGNAL']):
                    T_prime = nandetrend(interp_T)
                else:
                    T_prime = interp_T - T_mean
                std_T = np.nanstd(T_prime, ddof=1)
                T_prime = np.nan_to_num(T_prime, 0)

            # get mean atmospheric pressure
            if ini['files']['meteo_filepath'] and ini['meteo']['pressure_column']:
                p_mean = get_closest_value(df_meteofiledata.iloc[:, ini['meteo']['pressure_column']-1], results['time'][n])
            else:
                p_mean = ini['param']['DEFAULT_PRESSURE']

            # get mean air temperature (in K, only for improved computation of air molar concentration)
            T_meteo = np.NaN
            if ini['files']['meteo_filepath'] and ini['meteo']['temperature_column']:
                T_meteo = get_closest_value(df_meteofiledata.iloc[:, ini['meteo']['temperature_column']-1], results['time'][n]) + 273.15
            else:
                T_meteo = np.NaN

            # get mean air relative humidity (in %, for lag-RH dependency, if any)
            if ini['files']['meteo_filepath'] and ini['meteo']['relative_humidity_column']:
                rh_meteo = get_closest_value(df_meteofiledata.iloc[:, ini['meteo']['relative_humidity_column']-1], results['time'][n])
            else:
                rh_meteo = np.NaN

            # calculate air molar concentration (mol m-3)
            if not pd.isna(T_meteo):
                air_mol_conc = 1/((T_meteo/273.15) * (1013.25/p_mean) * 22.414e-3)
            else:
                air_mol_conc = 1/((T_mean/273.15) * (1013.25/p_mean) * 22.414e-3)
                
            if (completeness_T >= ini['param']['COMPLETENESS_THRESHOLD']):

                # calculate temperature flux <w'T'>
                wT = xcov(w_prime, T_prime, [0, 0])

                # calculate spectrum for T
                [spec_T, spec_T_scaled, _] = cospectrum(ini, T_prime, T_prime, f, mean_wind_speed, plot=ini['run_param']['PLOT_COSPECTRUM'], filename=sonic_files_list['name'][n], spectrum_type='spec')
                # calculate cospectrum for wT
                [cospec_wT, cospec_wT_scaled, _] = cospectrum(ini, w_prime, T_prime, f, mean_wind_speed, plot=ini['run_param']['PLOT_COSPECTRUM'], filename=sonic_files_list['name'][n], spectrum_type='cospec')

                # steady state test for wT according to Foken and Wichura 1996 and to Mahrt 1998
                steady_state_wT_FW96 = test_steady_state_FW96(w_prime, T_prime)
                steady_state_wT_M98 = test_steady_state_M98(w_prime, T_prime)
                steady_state_wT_D99 = test_steady_state_D99(w_prime, T_prime)

            else:

                wT = np.NaN
                spec_T = [np.NaN] * ini['param']['NUM_FREQ_BINS']
                spec_T_scaled = [np.NaN] * ini['param']['NUM_FREQ_BINS']
                cospec_wT = [np.NaN] * ini['param']['NUM_FREQ_BINS']
                cospec_wT_scaled = [np.NaN] * ini['param']['NUM_FREQ_BINS']
                steady_state_wT_FW96 = np.NaN
                steady_state_wT_M98 = np.NaN
                steady_state_wT_D99 = np.NaN

            # calculate Obukhov length, L
            L = -T_mean*u_star**3/(0.4*9.81*wT)

            # calculate stability parameter, z/L
            zoL = ini['param']['SENSOR_HEIGHT']/L

            # calculate potential temperature, theta
            #           virtual potential temperature, theta_v
            # and associated wT covariance
            p_0 = 1000  # reference pressure for potential temperature [hPa]
            if ((completeness_T >= ini['param']['COMPLETENESS_THRESHOLD']) and p_mean):
                theta_mean = T_mean*(p_0/p_mean)**(0.286)
                theta_prime = T_prime*(p_0/p_mean)**(0.286)
                wtheta = xcov(w_prime, theta_prime, [0, 0])
                if (completeness_H2O >= ini['param']['COMPLETENESS_THRESHOLD']):
                    theta_v = interp_T*(p_0/p_mean)**0.286*(1 + 0.38*1e-3*interp_H2O)  # note that H2O data is given in permille
                    theta_v_mean = np.nanmean(theta_v)
                    if (ini['param']['DETREND_TEMPERATURE_SIGNAL']):
                        theta_v_prime = nandetrend(theta_v)
                    else:
                        theta_v_prime = theta_v - theta_v_mean
                    theta_v_prime = np.nan_to_num(theta_v_prime, 0)
                    wtheta_v = xcov(w_prime, theta_v_prime, [0, 0])
                else:
                    theta_v_mean = np.NaN
                    wtheta_v = np.NaN
            else:
                theta_mean = np.NaN
                theta_v_mean = np.NaN
                wtheta = np.NaN
                wtheta_v = np.NaN

            # test on developed turbulent conditions
            ITC_w, ITC_u, ITC_T = test_ITC(w_prime, u_prime, T_prime, wT, zoL, u_star, ini['param']['LATITUDE'])

            # store results for sonic data
            results['MET']['u_unrot'][n] = mean_uvw[0]
            results['MET']['v_unrot'][n] = mean_uvw[1]
            results['MET']['w_unrot'][n] = mean_uvw[2]
            results['MET']['u_rot'][n] = mean_uvw_rot[0]
            results['MET']['v_rot'][n] = mean_uvw_rot[1]
            results['MET']['w_rot'][n] = mean_uvw_rot[2]
            results['MET']['std_u'][n] = std_uvw[0]
            results['MET']['std_v'][n] = std_uvw[1]
            results['MET']['std_w'][n] = std_uvw[2]
            results['MET']['wsh'][n] = mean_wind_speed
            results['MET']['std_wsh'][n] = std_wind_speed
            results['MET']['wdir'][n] = mean_wind_dir
            results['MET']['std_wdir'][n] = std_wind_dir
            results['MET']['phi'][n] = rot_phi*180/3.1415
            results['MET']['tilt'][str(n)] = R_tilt
            results['MET']['uw'][n] = uw
            results['MET']['vw'][n] = vw
            results['MET']['uv'][n] = uv
            results['MET']['uu'][n] = uu
            results['MET']['vv'][n] = vv
            results['MET']['ww'][n] = ww
            results['MET']['ust'][n] = u_star
            results['MET']['T'][n] = T_mean
            results['MET']['std_T'][n] = std_T
            results['MET']['wT'][n] = wT
            results['MET']['L'][n] = L
            results['MET']['zoL'][n] = zoL
            results['MET']['spec_T'][n] = spec_T
            results['MET']['spec_T_scaled'][n] = spec_T_scaled
            results['MET']['cospec_wT'][n] = cospec_wT
            results['MET']['cospec_wT_scaled'][n] = cospec_wT_scaled
            results['MET']['p'][n] = p_mean
            results['MET']['rho_air_molar'][n] = air_mol_conc
            results['MET']['theta'][n] = theta_mean
            results['MET']['theta_v'][n] = theta_v_mean
            results['MET']['wtheta'][n] = wtheta
            results['MET']['wtheta_v'][n] = wtheta_v
            results['MET']['P_air'][n] = p_mean
            results['MET']['T_air'][n] = T_meteo
            results['MET']['RH_air'][n] = rh_meteo
            results['MET']['qaqc']['completeness_sonic'][n] = completeness_sonicdata
            results['MET']['qaqc']['num_spikes_w'][n] = num_spikes_w
            results['MET']['qaqc']['num_spikes_T'][n] = num_spikes_T
            results['MET']['qaqc']['SST_wT_FW96'][n] = steady_state_wT_FW96
            results['MET']['qaqc']['SST_wT_M98'][n] = steady_state_wT_M98
            results['MET']['qaqc']['SST_wT_D99'][n] = steady_state_wT_D99
            results['MET']['qaqc']['ITC_w'][n] = ITC_w
            results['MET']['qaqc']['ITC_u'][n] = ITC_u
            results['MET']['qaqc']['ITC_T'][n] = ITC_T
            results['MET']['qaqc']['IPT_w'][n] = IPT_w
            results['MET']['qaqc']['IPT_T'][n] = IPT_T

# %% process IRGA data

            # IRGA tracers
            if ini['run_param']['CONCENTRATION_TYPE'] == 0:

                process_irga_data_day = True

                lag_clock_drift_samples = 0
                # get clock-drift time lag drift if input file provided
                if ini['files']['lag_clock_drift_filepath']:
                    # time lag drift info are present and must be accounted for
                    lag_clock_drift_samples = int(-df_lag_clock_drift['TDC-computer'][abs(df_lag_clock_drift.index - results['time'][n]).argmin()] * ini['param']['SAMPLING_RATE_FINAL'])  # find closest time lag based on timestamp

                # get prescribed time lag (instrumental + physical) if input file provided
                if ini['param']['LAG_DETECT_METHOD'] == 'PRESCRIBED':
                    # get prescribed time lag (clock-drift + physical)
                    closest_index = abs(df_lag_prescribed.index - results['time'][n]).argmin()
                    if abs(df_lag_prescribed.index - results['time'][n]).min() == datetime.timedelta():
                        # a lag is present for this half-hour
                        prescribed_lag_samples = int(df_lag_prescribed['time lag in s'][closest_index] * ini['param']['SAMPLING_RATE_FINAL'])  # find closest time lag based on timestamp
                    else:                        
                        # a lag is not present for this half-hour, compute the median lag of the ten closest values
                        # the median must be computed on the physical lag

                        # Synchronize on df_lag_prescribed index
                        df_lag_clock_drift_synchronized = df_lag_clock_drift.reindex(df_lag_prescribed.index, method='nearest')

                        # Calculate indices for the non-NaN values around the gap
                        valid_indices = range(len(df_lag_prescribed))
                        valid_distances = [abs(i - closest_index) for i in valid_indices]
                        sorted_indices = [x for _, x in sorted(zip(valid_distances, valid_indices))]

                        # Select the closest 10 values for median calculation
                        selected_indices = sorted_indices[:10]

                        # Calculate the median of the selected values
                        prescribed_lag_samples = (np.median(df_lag_prescribed.iloc[selected_indices,0].values + df_lag_clock_drift_synchronized.iloc[selected_indices,0].values) - df_lag_clock_drift_synchronized.iloc[closest_index,0])* ini['param']['SAMPLING_RATE_FINAL']

                # compute flux correction factor for low-pass filtering
                if ini['param']['LPFC'] == 0:
                    cf_lpf = 1
                elif ini['param']['LPFC'] == 1:
                    cf_lpf = correction_factor_lpf(mean_wind_speed, zoL, ini, df_lpfc=df_lpfc, ctrplot=ini['run_param']['PLOT_CORRECTION_FACTOR_LPF'])
                elif ini['param']['LPFC'] == 2:
                    cf_lpf = correction_factor_lpf(mean_wind_speed, zoL, ini, df_lpfc=df_lpfc)
                results['IRGA']['cf_lpf'][n] = cf_lpf

                # loop on concentrations
                for i in range(len(ini['irga']['irga_columns'])-1):

                    lag_phys_wd_center = ini['param']['LAG_WINDOW_CENTER']

                    # replace by the rh-dependent physical lag, if the tracer is present in the rh_dependency file
                    if ini['param']['LAG_RH_DEPENDENCY'] == 1:
                        searched_mz = ini['irga']['irga_names'][i]
                        matching_cols = df_lag_rh_dependency.columns[
                            (df_lag_rh_dependency.columns.astype(float) >= searched_mz - 0.001) &
                            (df_lag_rh_dependency.columns.astype(float) <= searched_mz + 0.001)
                        ]
                        if not matching_cols.empty:
                            closest_index = abs(df_lag_rh_dependency.index - min(rh_meteo, 100)).argmin()
                            lag_phys_wd_center = df_lag_rh_dependency.loc[closest_index, matching_cols[0]]

                    # detect spikes (and remove if requested)
                    if ini['param']['SPIKE_MODE'] != 0:
                        sonicdata[:, i+6], is_spike=spike_detection_vickers97(sonicdata[:, i+6],
                                                     spike_mode=ini['param']['SPIKE_MODE'],
                                                     avrg_len=int(ini['param']['WINDOW_LENGTH']/60),
                                                     ac_freq=ini['param']['SAMPLING_RATE_FINAL'],
                                                     spike_limit=ini['param']['SPIKE_DETECTION_THRESHOLD_IRGA'],
                                                     max_consec_spikes=ini['param']['MAX_CONSEC_SPIKES'],
                                                     ctrplot=ini['run_param']['PLOT_SPIKE_DETECTION']
                                                     )
                        num_spikes = int(sum(is_spike))
                    else:
                        num_spikes = np.nan

                    # test on instrument malfunctions from Vitale 2020
                    IPT = inst_prob_test(sonicdata[:,i+6], ini['param']['DETREND_TRACER_SIGNAL'], ini['param']['SAMPLING_RATE_FINAL'], ini['run_param']['PLOT_INST_PROB_TEST'], 'c or h')
                    # IPT = [np.nan for a in range(12)]
    
                    interp_c = interp1d(sonicdata[:,0], sonicdata[:,i+6], kind='nearest', bounds_error=False, fill_value=np.NaN)(interp_time)
                    completeness_IRGA = sum(np.isfinite(interp_c))/(ini['param']['WINDOW_LENGTH']*ini['param']['SAMPLING_RATE_FINAL'])

                    if completeness_IRGA >= ini['param']['COMPLETENESS_THRESHOLD']:

                        # calculate slope of IRGA tracer signal
                        dc_dt = nanlinfit(interp_c)[0]
                        dc_dt = dc_dt*ini['param']['SAMPLING_RATE_FINAL']

                        # Reynolds decomposition / detrending
                        c_mean = np.nanmean(interp_c)
                        if (ini['param']['DETREND_TRACER_SIGNAL']):
                            c_prime = nandetrend(interp_c)
                        else:
                            c_prime = interp_c - c_mean
                        std_c = np.nanstd(c_prime, ddof=1)
                        c_prime = np.nan_to_num(c_prime, 0)

                        # calculate cross-covariance over a given lag-window
                        cov_wc = xcov(w_prime, c_prime, [lag_phys_wd_center + lag_clock_drift_samples - ini['param']['LAG_OUTER_WINDOW_SIZE'], lag_phys_wd_center + lag_clock_drift_samples + ini['param']['LAG_OUTER_WINDOW_SIZE']])

                        # find time lag
                        if ini['param']['LAG_DETECT_METHOD'] == 'CONSTANT':
                            lag_samples = int(lag_phys_wd_center + lag_clock_drift_samples)
                        elif ini['param']['LAG_DETECT_METHOD'] == 'MAX':
                            lag_samples = find_covariance_peak(cov_wc, ini['param']['LAG_OUTER_WINDOW_SIZE'], ini['param']['LAG_INNER_WINDOW_SIZE'], lag_phys_wd_center + lag_clock_drift_samples, ini['param']['LAG_COVPEAK_FILTER_LENGTH'], 0, ini['run_param']['PLOT_FIND_COVARIANCE_PEAK'], sonic_files_list['name'][n])  # lag in samples
                        elif ini['param']['LAG_DETECT_METHOD'] == 'MAX_WITH_DEFAULT':
                            lag_samples = find_covariance_peak(cov_wc, ini['param']['LAG_OUTER_WINDOW_SIZE'], ini['param']['LAG_INNER_WINDOW_SIZE'], lag_phys_wd_center + lag_clock_drift_samples, ini['param']['LAG_COVPEAK_FILTER_LENGTH'], 0, ini['run_param']['PLOT_FIND_COVARIANCE_PEAK'], sonic_files_list['name'][n])  # lag in samples
                            # if lag was found at the extremum of the search window, replace by the default
                            if (lag_samples == lag_phys_wd_center + lag_clock_drift_samples - ini['param']['LAG_INNER_WINDOW_SIZE'] + 1 or 
                                lag_samples == lag_phys_wd_center + lag_clock_drift_samples + ini['param']['LAG_INNER_WINDOW_SIZE']):
                                lag_samples = lag_phys_wd_center + lag_clock_drift_samples
                        elif ini['param']['LAG_DETECT_METHOD'] == 'PRESCRIBED':
                            lag_samples = int(prescribed_lag_samples)

                        lagtime_individual = lag_samples/ini['param']['SAMPLING_RATE_FINAL']  # seconds

                        # calculate flux in umol/m^2/s or mmol/m^2/s
                        flux_individual = (cov_wc[ini['param']['LAG_OUTER_WINDOW_SIZE'] - (ini['param']['LAG_WINDOW_CENTER'] + lag_clock_drift_samples) + lag_samples])

                        # align c_prime for time-lag and trim the circular part, corrected time-series being needed for further computations
                        if lag_samples > 0:
                            # Shift c_prime forward (to the right), remove front (wrapped)
                            c_prime_trim = np.roll(c_prime, lag_samples)[lag_samples:]
                            # Remove first lag_samples from w_prime to align it in time
                            w_prime_trim = w_prime[lag_samples:]
                        elif lag_samples < 0:
                            lag_abs = abs(lag_samples)
                            # Shift c_prime backward (to the left), remove end (wrapped)
                            c_prime_trim = np.roll(c_prime, lag_samples)[:-lag_abs]
                            # Remove last lag_abs from w_prime to align
                            w_prime_trim = w_prime[:-lag_abs]

                        # calculate spectrum for c and cospectrum for wc
                        [spec_c, spec_c_scaled, _] = cospectrum(ini, c_prime_trim, c_prime_trim, f, mean_wind_speed, spectrum_type='spec')
                        [cospec_wc, cospec_wc_scaled, _] = cospectrum(ini, w_prime_trim, c_prime_trim, f, mean_wind_speed, spectrum_type='cospec')

                        # steady state test according to Foken and Wichura 1996 and to Mahrt 1998
                        steady_state_wc_FW96 = test_steady_state_FW96(w_prime_trim, c_prime_trim)
                        steady_state_wc_M98 = test_steady_state_M98(w_prime_trim, c_prime_trim)
                        steady_state_wc_D99 = test_steady_state_D99(w_prime_trim, c_prime_trim)

                        # compute flux uncertainties
                        [flux_noise_mean, flux_noise_std, flux_noise_rmse, random_error_FS, random_flux, random_error_noise] = \
                            flux_uncertainties(ini, w_prime_trim, c_prime_trim)

                        # convert IRGA flux from kinematic units to desired units
                        flux_individual, flux_noise_mean, flux_noise_std, flux_noise_rmse, random_error_FS, random_flux, random_error_noise = \
                            flux_unit_conversion(flux_individual, flux_noise_mean, flux_noise_std, flux_noise_rmse, random_error_FS, random_flux, random_error_noise, 
                                                 'to umol m-2 s-1', air_mol_conc)

                        # store results for individual IRGA tracer
                        results['IRGA'][str(i)]['conc_mean'][n] = c_mean
                        results['IRGA'][str(i)]['conc_std'][n] = std_c
                        results['IRGA'][str(i)]['lagtime'][n] = lagtime_individual
                        results['IRGA'][str(i)]['lagtime_clock_drift'][n] = lag_clock_drift_samples/ini['param']['SAMPLING_RATE_FINAL']
                        results['IRGA'][str(i)]['flux'][n] = flux_individual
                        results['IRGA'][str(i)]['spec'][n] = spec_c
                        results['IRGA'][str(i)]['spec_scaled'][n] = spec_c_scaled
                        results['IRGA'][str(i)]['cospec'][n] = cospec_wc
                        results['IRGA'][str(i)]['cospec_scaled'][n] = cospec_wc_scaled
                        results['IRGA'][str(i)]['qaqc']['completeness_IRGA'][n] = completeness_IRGA
                        results['IRGA'][str(i)]['qaqc']['num_spikes'][n] = num_spikes
                        results['IRGA'][str(i)]['qaqc']['SST_FW96'][n] = steady_state_wc_FW96
                        results['IRGA'][str(i)]['qaqc']['SST_M98'][n] = steady_state_wc_M98
                        results['IRGA'][str(i)]['qaqc']['SST_D99'][n] = steady_state_wc_D99
                        results['IRGA'][str(i)]['qaqc']['IPT'][n] = IPT
                        results['IRGA'][str(i)]['qaqc']['flux_SNR'][n] = abs(flux_individual - flux_noise_mean)/flux_noise_std
                        results['IRGA'][str(i)]['qaqc']['flux_noise_mean'][n] = flux_noise_mean
                        results['IRGA'][str(i)]['qaqc']['flux_noise_std'][n] = flux_noise_std
                        results['IRGA'][str(i)]['qaqc']['flux_noise_rmse'][n] = flux_noise_rmse
                        results['IRGA'][str(i)]['qaqc']['random_error_FS'][n] = random_error_FS
                        results['IRGA'][str(i)]['qaqc']['random_flux'][n] = random_flux
                        results['IRGA'][str(i)]['qaqc']['random_error_noise'][n] = random_error_noise
                        cov_data['IRGA'][str(i)]['cov'][n] = cov_wc


# %% process TRACER data

            if ini['run_param']['CONCENTRATION_TYPE'] == 1:
                # load corresponding tracer file
                tracerdata, error_code, idx_tracers_to_process, tracer_file_index, *rest = read_main_inputs(sonic_files_list['path'][n], sonic_files_list['name'][n], 'tracer', ini, OF, idx_tracers_to_process, tracer_files_list, results, cov_data, out_len)
                if rest:
                    results = rest[0]
                    cov_data = rest[1]
    
                # process tracer data
                if not(error_code):  # corresponding tracer file was found
                
                    # compute completeness of tracer data
                    completeness_tracerdata = sum(np.isfinite(tracerdata['conc']), 1)/(ini['param']['WINDOW_LENGTH']*ini['param']['SAMPLING_RATE_FINAL'])
    
                    # check validity of timestamp series, if present
                    time_series = pd.Series(np.vectorize(datetime.datetime.utcfromtimestamp)(tracerdata['time']))
                    threshold = pd.Timedelta(milliseconds=2/ini['param']['SAMPLING_RATE_TRACER']*1000)
                    msg = check_raw_timestamps(time_series, threshold)
                    if msg[0]:
                        print('tracer file: ' + msg[0]); OF.write('tracer file: ' + msg[0] + '\n')
    
                    # # debug, check for time differences between sonic and tracer timestamps
                    # if completeness_tracerdata[0] > ini['param']['COMPLETENESS_THRESHOLD']:
                    #     soniclen=len(sonicdata[:,0])
                    #     tracerlen=len(tracerdata['time'])
                    #     timedif = sonicdata[0:min(soniclen,tracerlen),0] - tracerdata['time'][0:min(soniclen,tracerlen)]
                    #     if abs(max(timedif)) > 0.1:
                    #         import matplotlib.pyplot as plt
                    #         plt.plot(timedif)
                    #         print('timedif detected for file : ' + sonic_files_list['name'][n] )  # + '\n'

                    # upsample data to final sampling rate if necessary
                    if (ini['param']['SAMPLING_RATE_TRACER'] < ini['param']['SAMPLING_RATE_FINAL']):
                        first_ts = tracerdata['time'][0, 0]

                        upsampling_factor = ini['param']['SAMPLING_RATE_FINAL']/ini['param']['SAMPLING_RATE_TRACER']
                        tracerdata['conc'] = scipy.signal.resample(tracerdata['conc'], int(tracerdata['conc'].shape[0]*upsampling_factor))

                        # reconstruct the timestamp colum (used later on by interp1d). 
                        tracerdata['time'][:] = (first_ts +
                                                 np.asarray([x for x in range(0, len(tracerdata['time']))])/ini['param']['SAMPLING_RATE_FINAL'])

                    lag_clock_drift_samples = 0
                    # get clock-drift time lag drift if input file provided
                    if ini['files']['lag_clock_drift_filepath']:
                        # time lag drift info are present and must be accounted for
                        lag_clock_drift_samples = int(-df_lag_clock_drift['TDC-computer'][abs(df_lag_clock_drift.index - results['time'][n]).argmin()] * ini['param']['SAMPLING_RATE_FINAL'])  # find closest time lag based on timestamp

                    # get prescribed time lag (instrumental + physical) if input file provided
                    if ini['param']['LAG_DETECT_METHOD'] == 'PRESCRIBED':

                        closest_index = abs(df_lag_prescribed.index - results['time'][n]).argmin()
                        if abs(df_lag_prescribed.index - results['time'][n]).min() == datetime.timedelta():
                            # a lag is present for this half-hour
                            prescribed_lag_samples = int(df_lag_prescribed['time lag in s'][closest_index] * ini['param']['SAMPLING_RATE_FINAL'])  # find closest time lag based on timestamp
                        else:                 
                            # a lag is not present for this half-hour, compute the median lag of the ten closest values
                            # the median must be computed on the physical lag

                            # Synchronize on df_lag_prescribed index
                            df_lag_clock_drift_synchronized = df_lag_clock_drift.reindex(df_lag_prescribed.index, method='nearest')

                            # Calculate indices for the non-NaN values around the gap
                            valid_indices = range(len(df_lag_prescribed))
                            valid_distances = [abs(i - closest_index) for i in valid_indices]
                            sorted_indices = [x for _, x in sorted(zip(valid_distances, valid_indices))]

                            # Select the closest 10 values for median calculation
                            selected_indices = sorted_indices[:10]

                            # Calculate the median of the selected values
                            prescribed_lag_samples = (np.median(df_lag_prescribed.iloc[selected_indices,0].values + df_lag_clock_drift_synchronized.iloc[selected_indices,0].values) - df_lag_clock_drift_synchronized.iloc[closest_index,0])* ini['param']['SAMPLING_RATE_FINAL']

                    # compute flux correction factor for low-pass filtering
                    if ini['param']['LPFC'] == 0:
                        cf_lpf = 1
                    elif ini['param']['LPFC'] == 1:
                        cf_lpf = correction_factor_lpf(mean_wind_speed, zoL, ini, df_lpfc=df_lpfc, ctrplot=ini['run_param']['PLOT_CORRECTION_FACTOR_LPF'])
                    elif ini['param']['LPFC'] == 2:
                        cf_lpf = correction_factor_lpf(mean_wind_speed, zoL, ini, df_lpfc=df_lpfc)
                    results['TRACER']['cf_lpf'][n] = cf_lpf

                    # save config parameters
                    results['TRACER']['default_CC_kinetic'][n] = tracerdata['default_CC_kinetic']

                    # iterate over tracers
                    msg_tracer = ''
                    for i in range(tracerdata['conc'].shape[1]):

                        # check if whole values are NaNs and skip process if the case
                        if np.all(np.isnan(tracerdata['conc'][:, i])):
                            msg_tracer = msg_tracer + f"Warning: full NaN values found in tracerdata[{i}]\n"
                            OF.write((f"\nWarning: full NaN values found in tracerdata[{i}]"))
                            continue
    
                        # check for NaNs, send a warning if present and skip process if the case
                        if np.any(np.isnan(tracerdata['conc'][:, i])):
                            msg_tracer = msg_tracer + f"Warning: NaN values found in tracerdata[{i}]\n"
                            OF.write((f"Warning: NaN values found in tracerdata[{i}]\n"))
                            continue

                        # check completeness of tracer data
                        if(completeness_tracerdata[i] < ini['param']['COMPLETENESS_THRESHOLD']):
                            if i == 0:
                                msg_tracer = msg_tracer + 'tracer file incomplete (' + str(int(completeness_tracerdata[i] * 100)) + '%), tracer process skipped\n'
                                continue

                        lag_phys_wd_center = ini['param']['LAG_WINDOW_CENTER']

                        # replace by the rh-dependent physical lag, if the tracer is present in the rh_dependency file
                        if ini['param']['LAG_RH_DEPENDENCY'] == 1:
                            searched_mz = round(float(tracerdata['mz'][i]), 3)
                            matching_cols = df_lag_rh_dependency.columns[
                                (df_lag_rh_dependency.columns.astype(float) >= searched_mz - 0.001) &
                                (df_lag_rh_dependency.columns.astype(float) <= searched_mz + 0.001)
                            ]
                            if not matching_cols.empty:
                                closest_index = abs(df_lag_rh_dependency.index - min(rh_meteo, 100)).argmin()
                                lag_phys_wd_center = df_lag_rh_dependency.loc[closest_index, matching_cols[0]]
  
                        # statistics on tracer signal and zero
                        c_mean = np.nanmean(tracerdata['conc'][:, i])
                        c_std = np.nanstd(tracerdata['conc'][:, i])
                        c_Q5 = np.nanpercentile(tracerdata['conc'][:, i], 5)
                        c_Q95 = np.nanpercentile(tracerdata['conc'][:, i], 95)
                        c_acc = np.nanmean(abs(tracerdata['conc_acc'][:, i]/tracerdata['conc'][:, i]))
                        c_prec = np.sqrt(np.nansum(tracerdata['conc_prec'][:, i]**2) / np.count_nonzero(~np.isnan(tracerdata['conc_prec'][:, i])))
                        c_LOD = 3*np.sqrt(np.nansum(tracerdata['zero_prec'][:, i]**2) / np.count_nonzero(~np.isnan(tracerdata['zero_prec'][:, i])))

                        # detect spikes (and remove if requested)
                        if ini['param']['SPIKE_MODE'] != 0:
                            tracerdata['conc'][:, i], is_spike = spike_detection_vickers97(tracerdata['conc'][:, i], 
                                                  spike_mode = ini['param']['SPIKE_MODE'],
                                                  avrg_len = int(ini['param']['WINDOW_LENGTH']/60),
                                                  ac_freq = ini['param']['SAMPLING_RATE_FINAL'],
                                                  spike_limit = ini['param']['SPIKE_DETECTION_THRESHOLD_TRACER'],
                                                  max_consec_spikes = ini['param']['MAX_CONSEC_SPIKES'],
                                                  ctrplot = ini['run_param']['PLOT_SPIKE_DETECTION']
                                                  )
                            num_spikes = int(sum(is_spike))
                        else:
                            num_spikes = np.nan

                        # test on instrument malfunctions from Vitale 2020
                        if ini['run_param']['PROCESS_IPT_TRACERS']:
                            IPT = inst_prob_test(tracerdata['conc'][:, i], ini['param']['DETREND_TRACER_SIGNAL'], ini['param']['SAMPLING_RATE_FINAL'], ini['run_param']['PLOT_INST_PROB_TEST'], 'tracer')
                        else:
                            IPT = [np.nan for a in range(12)]

                        # caclulate slope of tracer signal
                        dc_dt = nanlinfit(tracerdata['conc'][:, i])[0]
                        dc_dt = dc_dt*ini['param']['SAMPLING_RATE_TRACER']

                        # interpolate concentration data (keep NaNs and extrapolate using NaNs)
                        interp_c = interp1d(np.squeeze(tracerdata['time'][:]), np.squeeze(tracerdata['conc'][:, i]), kind='nearest', bounds_error=False, fill_value=np.NaN)(interp_time)

                        # Reynolds decomposition / detrending
                        if ini['param']['DETREND_TRACER_SIGNAL']:
                            c_prime = nandetrend(interp_c)
                        else:
                            c_prime = interp_c - c_mean

                        # trim w_prime and c_prime for NaNs, keeping the alignment
                        first, last = np.argmax(~np.isnan(w_prime)), len(w_prime) - np.argmax(~np.isnan(w_prime)[::-1]) - 1
                        w_prime_trim = w_prime[first:last + 1]; c_prime_trim = c_prime[first:last + 1]
                        first, last = np.argmax(~np.isnan(c_prime_trim)), len(c_prime_trim) - np.argmax(~np.isnan(c_prime_trim)[::-1]) - 1
                        w_prime_trim = w_prime_trim[first:last + 1]; c_prime_trim = c_prime_trim[first:last + 1]

                        # calculate cross-covariance over a given lag-window
                        cov_wc = xcov(w_prime_trim, c_prime_trim, [lag_phys_wd_center + lag_clock_drift_samples - ini['param']['LAG_OUTER_WINDOW_SIZE'], lag_phys_wd_center + lag_clock_drift_samples + ini['param']['LAG_OUTER_WINDOW_SIZE']])

                        # find time lag
                        if ini['param']['LAG_DETECT_METHOD'] == 'CONSTANT':
                            lag_samples = int(lag_phys_wd_center + lag_clock_drift_samples)
                        elif ini['param']['LAG_DETECT_METHOD'] == 'MAX':
                            lag_samples = find_covariance_peak(cov_wc, ini['param']['LAG_OUTER_WINDOW_SIZE'], ini['param']['LAG_INNER_WINDOW_SIZE'], lag_phys_wd_center + lag_clock_drift_samples, ini['param']['LAG_COVPEAK_FILTER_LENGTH'], 0, ini['run_param']['PLOT_FIND_COVARIANCE_PEAK'], tracer_files_list['name'][tracer_file_index])  # lag in samples
                        elif ini['param']['LAG_DETECT_METHOD'] == 'MAX_WITH_DEFAULT':
                            lag_samples = find_covariance_peak(cov_wc, ini['param']['LAG_OUTER_WINDOW_SIZE'], ini['param']['LAG_INNER_WINDOW_SIZE'], lag_phys_wd_center + lag_clock_drift_samples, ini['param']['LAG_COVPEAK_FILTER_LENGTH'], 0, ini['run_param']['PLOT_FIND_COVARIANCE_PEAK'], tracer_files_list['name'][tracer_file_index])  # lag in samples
                            # if lag was found at the extremum of the search window, replace by the default
                            if (lag_samples == lag_phys_wd_center + lag_clock_drift_samples - ini['param']['LAG_INNER_WINDOW_SIZE'] + 1 or 
                                lag_samples == lag_phys_wd_center + lag_clock_drift_samples + ini['param']['LAG_INNER_WINDOW_SIZE']):
                                lag_samples = lag_phys_wd_center + lag_clock_drift_samples
                        elif ini['param']['LAG_DETECT_METHOD'] == 'PRESCRIBED':
                            lag_samples = int(prescribed_lag_samples)

                        lagtime_individual = lag_samples/ini['param']['SAMPLING_RATE_FINAL']  # seconds

                        # calculate flux in ug/m^2/s
                        flux_individual = cov_wc[ini['param']['LAG_OUTER_WINDOW_SIZE'] - (lag_phys_wd_center + lag_clock_drift_samples) + lag_samples]

                        # correct flux for low-pass filtering
                        flux_individual = flux_individual * cf_lpf

                        # align c_prime for time-lag and trim the circular part, corrected time-series being needed for further computations
                        if lag_samples > 0:
                            # Shift c_prime forward (to the right), remove front (wrapped)
                            c_prime_trim = np.roll(c_prime_trim, lag_samples)[lag_samples:]
                            # Remove first lag_samples from w_prime to align it in time
                            w_prime_trim = w_prime_trim[lag_samples:]
                        elif lag_samples < 0:
                            lag_abs = abs(lag_samples)
                            # Shift c_prime backward (to the left), remove end (wrapped)
                            c_prime_trim = np.roll(c_prime_trim, lag_samples)[:-lag_abs]
                            # Remove last lag_abs from w_prime to align
                            w_prime_trim = w_prime_trim[:-lag_abs]

                        # calculate spectrum for c and cospectrum for wc
                        [spec_c, spec_c_scaled, _] = cospectrum(ini, c_prime_trim, c_prime_trim, f, mean_wind_speed, spectrum_type='spec')
                        [cospec_wc, cospec_wc_scaled, _] = cospectrum(ini, w_prime_trim, c_prime_trim, f, mean_wind_speed, spectrum_type='cospec')

                        # steady state test according to Foken and Wichura 1996 and to Mahrt 1998
                        steady_state_wc_FW96 = test_steady_state_FW96(w_prime_trim, c_prime_trim)
                        steady_state_wc_M98 = test_steady_state_M98(w_prime_trim, c_prime_trim)
                        steady_state_wc_D99 = test_steady_state_D99(w_prime_trim, c_prime_trim)

                        # compute flux uncertainties
                        [flux_noise_mean, flux_noise_std, flux_noise_rmse, random_error_FS, random_flux, random_error_noise] = \
                            flux_uncertainties(ini, w_prime_trim, c_prime_trim)

                        # convert tracer flux and associated variables from kinematic units to desired units
                        flux_individual, flux_noise_mean, flux_noise_std, flux_noise_rmse, random_error_FS, random_flux, random_error_noise = \
                            flux_unit_conversion(flux_individual, flux_noise_mean, flux_noise_std, flux_noise_rmse, random_error_FS, random_flux, random_error_noise, 
                                                 'to ug m-2 s-1', air_mol_conc, tracerdata['mz'][i] - 1)

                        # store results for individual tracer
                        results['TRACER'][str(i)]['conc_mean'][n] = c_mean
                        results['TRACER'][str(i)]['conc_std'][n] = c_std
                        results['TRACER'][str(i)]['conc_Q5'][n] = c_Q5
                        results['TRACER'][str(i)]['conc_Q95'][n] = c_Q95
                        results['TRACER'][str(i)]['conc_acc'][n] = c_acc
                        results['TRACER'][str(i)]['conc_prec'][n] = c_prec
                        results['TRACER'][str(i)]['conc_LOD'][n] = c_LOD
                        results['TRACER'][str(i)]['lagtime'][n] = lagtime_individual
                        results['TRACER'][str(i)]['lagtime_clock_drift'][n] = lag_clock_drift_samples/ini['param']['SAMPLING_RATE_FINAL']
                        results['TRACER'][str(i)]['flux'][n] = flux_individual
                        results['TRACER'][str(i)]['spec'][n] = spec_c
                        results['TRACER'][str(i)]['spec_scaled'][n] = spec_c_scaled
                        results['TRACER'][str(i)]['cospec'][n] = cospec_wc
                        results['TRACER'][str(i)]['cospec_scaled'][n] = cospec_wc_scaled
                        results['TRACER'][str(i)]['qaqc']['completeness_TRACER'][n] = completeness_tracerdata[i]
                        results['TRACER'][str(i)]['qaqc']['num_spikes'][n] = num_spikes
                        results['TRACER'][str(i)]['qaqc']['SST_FW96'][n] = steady_state_wc_FW96
                        results['TRACER'][str(i)]['qaqc']['SST_M98'][n] = steady_state_wc_M98
                        results['TRACER'][str(i)]['qaqc']['SST_D99'][n] = steady_state_wc_D99
                        results['TRACER'][str(i)]['qaqc']['IPT'][n] = IPT
                        results['TRACER'][str(i)]['qaqc']['flux_SNR'][n] = abs(flux_individual - flux_noise_mean)/flux_noise_std
                        results['TRACER'][str(i)]['qaqc']['flux_noise_mean'][n] = flux_noise_mean
                        results['TRACER'][str(i)]['qaqc']['flux_noise_std'][n] = flux_noise_std
                        results['TRACER'][str(i)]['qaqc']['flux_noise_rmse'][n] = flux_noise_rmse
                        results['TRACER'][str(i)]['qaqc']['random_error_FS'][n] = random_error_FS
                        results['TRACER'][str(i)]['qaqc']['random_flux'][n] = random_flux
                        results['TRACER'][str(i)]['qaqc']['random_error_noise'][n] = random_error_noise
                        results['TRACER'][str(i)]['calibration'][n] = tracerdata['calibration'][i]
                        results['TRACER'][str(i)]['transmission'][n] = tracerdata['transmission'][i]
                        results['TRACER'][str(i)]['Xr0'][n] = tracerdata['Xr0'][i]
                        results['TRACER'][str(i)]['cluster_min'][n] = tracerdata['cluster_min'][i]
                        results['TRACER'][str(i)]['cluster_max'][n] = tracerdata['cluster_max'][i]
                        results['TRACER'][str(i)]['k_reac'][n] = tracerdata['k_reac'][i]
                        results['TRACER'][str(i)]['FY'][n] = tracerdata['FY'][i]
                        results['TRACER'][str(i)]['IF'][n] = tracerdata['IF'][i]

                        cov_data['TRACER'][str(i)]['cov'][n] = cov_wc

                        print_progress(i+1, tracerdata['conc'].shape[1], 'tracers processed', flush_=True, end_='\r', width=10)

                    # display time needed
                    proc_end_time_hh = datetime.datetime.now()
                    msg = 'process time: ' + str((proc_end_time_hh - proc_start_time_hh))
                    print('\n' + msg); OF.write("\n" + msg + "\n\n")

                    # display messages
                    print(msg_tracer)

        # %% clean results for unprocessed sonic files
        clean_results(results)

# %% output data

        # save covariance functions to file
        # if (tracer_files_list or ini['irga']['irga_columns']) and ini['run_param']['WRITE_COV_OUTPUTS']:
        if ini['run_param']['WRITE_COV_OUTPUTS']:
            ts = results['time'][0]
            cov_filename = ini['files']['output_files_prefix'] + 'cov_%d%02d%02d.hdf5' % (ts.year, ts.month, ts.day)
            hdfdict.dump(cov_data, ini['files']['output_folder'] + '\\cov\\' + cov_filename, mode='w')

            # add attributes
            with h5py.File(ini['files']['output_folder'] + '\\cov\\' + cov_filename, 'r+') as hdf5_f:
                if ini['run_param']['CONCENTRATION_TYPE'] == 0:
                    # loop on tracers
                    for i in range(len(hdf5_f["IRGA"])):
                        hdf5_f['IRGA'][str(i)]['cov'].attrs['description'] = 'covariance function (cov vs lag)'; hdf5_f['IRGA'][str(i)]['cov'].attrs['units'] = 'ppbv ncps-1 m s-1 and lag in samples'
                elif ini['run_param']['CONCENTRATION_TYPE'] == 1:
                    # loop on tracers
                    for i in range(len(hdf5_f["TRACER"])):
                        hdf5_f['TRACER'][str(i)]['cov'].attrs['description'] = 'covariance function (cov vs lag)'; hdf5_f['TRACER'][str(i)]['cov'].attrs['units'] = 'ppbv ncps-1 m s-1 and lag in samples'

        # save results to file
        date = datetime.datetime.now()
        results['file_creation_time'] = '%d%02d%02d%02d%02d' % (date.year, date.month, date.day, date.hour, date.minute)
        ts = results['time'][0]
        res_filename = ini['files']['output_files_prefix'] + '%d%02d%02d.hdf5' % (ts.year, ts.month, ts.day)
        results['time'] = [date_obj.strftime('%Y-%m-%d %H-%M-%S') for date_obj in results['time']]
        hdfdict.dump(results, ini['files']['output_folder']+'\\' + res_filename, mode='w')

        # add attributes "name" and "units"
        if  'TRACER' in results:  # TRACER part has been created in the result dict (tracer data for that day)
            add_attributes(ini['files']['output_folder'], res_filename, process_irga_data_day, len(ini['irga']['irga_columns']), True, tracerdata)
        else:
            add_attributes(ini['files']['output_folder'], res_filename, process_irga_data_day, len(ini['irga']['irga_columns']))

    return unique_days
