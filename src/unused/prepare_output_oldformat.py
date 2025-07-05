import numpy as np
from copy import deepcopy


def prepare_output(ini, out_len, hdf5_nb_tracers=False, hdf5_mz_tracer=False, results=False):
    """
    Preparation of the results dict.
    If the positional args are not provided, it will create the sonic/IRGA part of the results dict.
    If the positional args are provided, it will create the TRACER part of the results dict.

    parameters
    ----------
    ini: dict, initialisation information
    out_len: int, number of sonic files to be processed
    hdf5_nb_tracers:
    hdf5_mz_tracer:
    results: dictionnary that will be filled with all info to be outputed in the hdf5 files

    returns
    -------
    results: dictionnary that will be filled with all info to be outputed in the hdf5 files

    comments
    --------
    Written by B. Heinesch, 2 January, 2023.
    University of Liege, Gembloux Agro-Bio Tech.

    """

    if hdf5_nb_tracers is False:

        results = {}

        results['procedure_version'] = ini['version']['procedure_version']             # semantic versioning
        results['procedure_version_date'] = ini['version']['procedure_version_date']

        results['file_creation_time'] = ''

        # store parameters
        results['param'] = {
            'SONIC_ORIENTATION': ini['param']['SONIC_ORIENTATION'],                   # degree
            'SENSOR_HEIGHT': ini['param']['SENSOR_HEIGHT'],                           # meters above roughness height
            'WINDOW_LENGTH': ini['param']['WINDOW_LENGTH'],                           # samples
            'SAMPLING_RATE_SONIC': ini['param']['SAMPLING_RATE_SONIC'],               # per second sonic anemometer sampling rate
            'SAMPLING_RATE_TRACER': ini['param']['SAMPLING_RATE_TRACER'],             # per second gas analyzer sampling rate
            'DETREND_TEMPERATURE_SIGNAL': ini['param']['DETREND_TEMPERATURE_SIGNAL'], # if 1, linear detrending is applied to temperature signal
            'DETREND_TRACER_SIGNAL': ini['param']['DETREND_TRACER_SIGNAL'],           # if 1, linear detrending is applied to tracer signal
            'LAG_DETECT_METHOD': ini['param']['LAG_DETECT_METHOD'],                   # time-lag detection method to be applied; choose between: 'CONST', 'MAX', 'MAX_WITH_DEFAULT'. 
                                                                                           # If CONST, then the center of the [LAG_SEARCH_MIN,LAG_SEARCH_MAX] will be used
            'LAG_PRESCRIBED': ini['files']['lag_prescribed_filepath'],                # If a mz proposed, only the lag on this (closest) mz will be estimated, and used for all mz. 
                                                                                           # Insert 0 if this option should not be used
            'LAG_TRACER_WINDOW_CENTER': ini['param']['LAG_TRACER_WINDOW_CENTER'],     # expected physical lag (not including possible clock drift)
            'LAG_OUTER_WINDOW_SIZE': ini['param']['LAG_OUTER_WINDOW_SIZE'],           # time lag plotting window around the expected lag (in samples)
            'LAG_INNER_WINDOW_SIZE': ini['param']['LAG_INNER_WINDOW_SIZE'],           # time lag search window around the expected lag (in samples)
            'LAG_COVPEAK_FILTER_LENGTH': ini['param']['LAG_COVPEAK_FILTER_LENGTH'],   # samples length of filter applied to covariance peak prior to lag-time determination
            'LAG_RH_DEPENDENCY': ini['param']['LAG_RH_DEPENDENCY'],                   # imposed time-lag, overruling, for provided mz, the MAX/MAX_WITH_DEFAULT window center
                                                                                      # or the PRESCRIBED lag (0 if applied, 1 if not)

            'COF_LPF': ini['param']['COF_LPF'],                                       # half-power cut-off frequency of the total transfer function (Hz)

            'NUM_FREQ_BINS': ini['param']['NUM_FREQ_BINS'],                           # number of logarithmically spaced frequency bins for (co)spectra
            'TILT_CORRECTION_MODE': ini['param']['TILT_CORRECTION_MODE'],             # 0: no tilt correction, 1: double rotation, 2: angle-dependent planar fit
            'APPLY_WPL_CORRECTION': ini['param']['APPLY_WPL_CORRECTION'],             # if 1, WPL correction is applied
            'COMPLETENESS_THRESHOLD': ini['param']['COMPLETENESS_THRESHOLD'],         # threshold value for completeness of input data (0.0 - 1.0)
            'DEFAULT_PRESSURE': ini['param']['DEFAULT_PRESSURE'],                     # hPa / mbar fixed pressure used if no pressure is given in the input files
            'SPIKE_MODE': ini['param']['SPIKE_MODE'],                                 # 0 = no spikes handling, 1 = spike detection, 2 = spikes detection and replacement
            'SPIKE_DETECTION_THRESHOLD_W': ini['param']['SPIKE_DETECTION_THRESHOLD_W'],             # standard deviations threshold for detecting w spikes in the input data
            'SPIKE_DETECTION_THRESHOLD_T': ini['param']['SPIKE_DETECTION_THRESHOLD_T'],             # standard deviations threshold for detecting T spikes in the input data
            'SPIKE_DETECTION_THRESHOLD_IRGA': ini['param']['SPIKE_DETECTION_THRESHOLD_IRGA'],       # standard deviations threshold for detecting IRGA spikes in the input data
            'SPIKE_DETECTION_THRESHOLD_TRACER': ini['param']['SPIKE_DETECTION_THRESHOLD_TRACER'],   # standard deviations threshold for detecting TRACER spikes in the input data
            'MAX_CONSEC_SPIKES': ini['param']['MAX_CONSEC_SPIKES'],                   # maximum number of consecutive outliers that can be considered as spikes
            'MEAS_ID': ini['param']['MEAS_ID']                                        # measurement ID
            }

        results['time'] = [np.NaN] * out_len
        results['freq'] = []                        # frequency axis of (co)spectra [1/s]
        results['MET'] = {
                'uvw': [[np.NaN]*3]*out_len,        # mean wind speed vector (u, v, w), with u component pointing into direction of mean wind [m/s]
                'std_uvw': [[np.NaN]*3]*out_len,    # standard deviation of wind speed components u, v, w [m/s]
                'wsh': [np.NaN] * out_len,          # horizontal wind speed [m/s]
                'std_wsh': [np.NaN] * out_len,      # horizontal wind speed standard deviation [m/s]
                'wdir': [np.NaN] * out_len,         # horizontal wind direction [deg]
                'std_wdir': [np.NaN] * out_len,     # horizontal wind direction standard deviation[deg]
                'phi': [np.NaN] * out_len,          # second rotation angle [deg]
                'tilt': dict(zip(tuple(str(x) for x in tuple(range(0,out_len))),
                                 [deepcopy(np.array([])) for i in range(out_len)])),
                'uw': [np.NaN] * out_len,           # covariance of along-wind and vertical wind component, <u'w'> [m^2/s^2]
                'vw': [np.NaN] * out_len,           # covariance of cross-wind and vertical wind component, <v'w'> [m^2/s^2]
                'uv': [np.NaN] * out_len,           # covariance of along-wind and cross-wind component, <u'v'> [m^2/s^2]
                'uu': [np.NaN] * out_len,           # auto-covariance of along-wind component, <u'u'> [m^2/s^2]
                'vv': [np.NaN] * out_len,           # auto-covariance of cross-wind component, <v'v'> [m^2/s^2]
                'ww': [np.NaN] * out_len,           # auto-covariance of vertical wind component, <w'w'> [m^2/s^2]
                'ust': [np.NaN] * out_len,          # friction velocity, u* [m/s]
                'T': [np.NaN] * out_len,            # sonic temperature [K]
                'std_T': [np.NaN] * out_len,        # standard deviation of temperature [K]
                'wT': [np.NaN] * out_len,           # temperature flux, <w'T'> [K*m/s]
                'L': [np.NaN] * out_len,            # Obukhov out_length [m]
                'zoL': [np.NaN] * out_len,          # stability parameter, z/L, dimensionless
                'spec_T': [[np.NaN] * ini['param']['NUM_FREQ_BINS']] * out_len, # spectrum for T, Sp(T',f)
                'spec_T_scaled': [[np.NaN] * ini['param']['NUM_FREQ_BINS']] * out_len, # scaled spectrum for T, f*Sp(T',f)/TT
                'cospec_wT': [[np.NaN] * ini['param']['NUM_FREQ_BINS']] * out_len, # cospectrum for wT, Co(w',c',f)
                'cospec_wT_scaled': [[np.NaN] * ini['param']['NUM_FREQ_BINS']] * out_len, # scaled cospectrum for wT, f*Co(w',T',f)/wT
                'p': [np.NaN] * out_len,            # pressure [hPa]
                'rho_air_molar': [np.NaN] * out_len, # dry molar air density [mol m-3]
                'theta': [np.NaN] * out_len,        # potential temperature [K]
                'theta_v': [np.NaN] * out_len,      # virtual potential temperature [K]
                'wtheta': [np.NaN] * out_len,       # potential temperature (heat) flux, <w'theta'> [K*m/s]
                'wtheta_v': [np.NaN] * out_len,     # vitual potential temperature (buoyancy) flux, <w'theta_v'> [K*m/s]
                'P_air': [np.NaN] * out_len,        # air pressure [mbar]
                'T_air': [np.NaN] * out_len,        # air temperature [mbar]
                'RH_air': [np.NaN] * out_len,       # air relative humidity [%]
                'qaqc':    {
                    'completeness': [np.NaN] * out_len,         # fraction of sonic data used in this interval
                    'num_spikes_w': [np.NaN] * out_len,         # number of w spikes detected in this interval
                    'num_spikes_T': [np.NaN] * out_len,         # number of T spikes detected in this interval
                    'SST_wT_FW96': [np.NaN] * out_len,          # steady state test for wT after Foken & Wichura 1996
                    'SST_wT_M98': [np.NaN] * out_len,           # steady state test for wT after Mahrt 1998
                    'SST_wT_D99': [np.NaN] * out_len,           # steady state test for wT after Dutaur 1999
                    'ITC_w': [np.NaN] * out_len,                # relative model deviation of integral turbulence characteristics test for w
                    'ITC_u': [np.NaN] * out_len,                # relative model deviation of integral turbulence characteristics test for u
                    'ITC_T': [np.NaN] * out_len,                # relative model deviation of integral turbulence characteristics test for T
                    'IPT_w': [[np.NaN]*12]*out_len,             # instrumental malfunctions tests for w, from Vitale 2020
                    'IPT_T': [[np.NaN]*12]*out_len,             # instrumental malfunctions tests for T, from Vitale 2020
                            }
                          }

        if ini['irga']['irga_columns']:
            IRGA = {}
            IRGA['compound'] = []                               # tracer name
            IRGA['conc_mean'] = [np.NaN] * out_len          # mean concentration [ppbv, if IRGA concentration provided in ppbv]
            IRGA['conc_std'] = [np.NaN] * out_len           # standard deviation of concentration
            IRGA['lagtime'] = [np.NaN] * out_len            # lag time determined from individual averaging interval [s]
            IRGA['lagtime_clock_drift'] = [np.NaN] * out_len  # clock drift lag time provided as input [s]
            IRGA['flux'] = [np.NaN] * out_len               # flux from individual lagtime [umol/m^2/s, if IRGA concentration provided in ppbv]
            IRGA['spec'] = [[np.NaN] * ini['param']['NUM_FREQ_BINS']] * out_len  # spectrum for c
            IRGA['spec_scaled'] = [[np.NaN] * ini['param']['NUM_FREQ_BINS']] * out_len  # scaled spectrum for c
            IRGA['cospec'] = [[np.NaN] * ini['param']['NUM_FREQ_BINS']] * out_len  # cospectrum for wc
            IRGA['cospec_scaled'] = [[np.NaN] * ini['param']['NUM_FREQ_BINS']] * out_len  # scaled cospectrum for wc
            IRGA['qaqc'] = {
                'completeness': [np.NaN] * out_len,        # fraction of tracer data used in this interval
                'num_spikes': [np.NaN] * out_len,          # number of spikes detected in this interval
                'SST_FW96': [np.NaN] * out_len,            # steady state test for tracer after Foken & Wichura 1996
                'SST_M98': [np.NaN] * out_len,             # steady state test for tracer after Mahrt 1998
                'SST_D99': [np.NaN] * out_len,             # steady state test for tracer after Dutaur 1999
                'IPT': [[np.NaN]*12]*out_len,              # instrumental malfunctions tests for irga, from Vitale 2020
                'flux_SNR': [np.NaN] * out_len,            # flux signal to noise ratio
                'flux_noise_std': [np.NaN] * out_len,      # standard deviation of flux noise far off integral time scale
                'flux_noise_mean': [np.NaN] * out_len,     # mean flux noise far off integral time scale
                'flux_noise_rmse': [np.NaN] * out_len,     # RMSE of flux noise far off integral time scale
                'random_error_FS': [np.NaN] * out_len,     # random error as described by Finkelstein and Sims 2001
                'random_error_noise': [np.NaN] * out_len,  # random error noise estimated according to Mauder 2013
                'random_flux': [np.NaN] * out_len,         # random flux level estimated by random shuffle criterium by Billesbach 2011
                          }

            results['IRGA'] = dict(zip(tuple(str(x) for x in tuple(range(0,len(ini['irga']['irga_columns'])-1))),
                                       [deepcopy(IRGA) for i in range(len(ini['irga']['irga_columns'])-1)]))  # struct containing results for each tracer

            results['IRGA']['cf_lpf'] = [np.NaN] * out_len

            del IRGA

            for n in range(len(ini['irga']['irga_columns'])-1):
                results['IRGA'][str(n)]['name'] = ini['irga']['irga_names'][n]

    else:
        TRACER = {}
        TRACER['mz'] = []                               # protonated mass to charge ratio
        TRACER['conc_mean'] = [np.NaN] * out_len        # mean tracer concentration [ppb]
        TRACER['conc_std'] = [np.NaN] * out_len         # standard deviation of tracer concentration [ppb]
        TRACER['lagtime'] = [np.NaN] * out_len          # lag time (lag time (physical + clock_drift)) determined from individual averaging interval [s]
        TRACER['lagtime_clock_drift'] = [np.NaN] * out_len    # clock_drift lag time provided as input [s]
        TRACER['flux'] = [np.NaN] * out_len             # flux from individual lagtime [ug/m^2/s]
        TRACER['spec'] = [[np.NaN] * ini['param']['NUM_FREQ_BINS']] * out_len           # spectrum for c
        TRACER['spec_scaled'] = [[np.NaN] * ini['param']['NUM_FREQ_BINS']] * out_len    # scaled spectrum for c
        TRACER['cospec'] = [[np.NaN] * ini['param']['NUM_FREQ_BINS']] * out_len         # cospectrum for wc
        TRACER['cospec_scaled'] = [[np.NaN] * ini['param']['NUM_FREQ_BINS']] * out_len  # scaled cospectrum for wc
        TRACER['qaqc'] = {
            'completeness': [np.NaN] * out_len,        # fraction of tracer data used in this interval
            'num_spikes': [np.NaN] * out_len,          # number of spikes detected in this interval
            'SST_FW96': [np.NaN] * out_len,            # steady state test for tracer after Foken & Wichura 1996
            'SST_M98': [np.NaN] * out_len,             # steady state test for tracer after Mahrt 1998
            'SST_D99': [np.NaN] * out_len,             # steady state test for tracer after Dutaur 1999
            'IPT': [[np.NaN]*12]*out_len,              # instrumental malfunctions tests for tracer, from Vitale 2020
            'flux_SNR': [np.NaN] * out_len,            # flux signal to noise ratio
            'flux_noise_std': [np.NaN] * out_len,      # standard deviation of flux noise far off integral time scale
            'flux_noise_mean': [np.NaN] * out_len,     # mean flux noise far off integral time scale
            'flux_noise_rmse': [np.NaN] * out_len,     # RMSE of flux noise far off integral time scale
            'random_error_FS': [np.NaN] * out_len,     # random error as described by Finkelstein and Sims 2001
            'random_error_noise': [np.NaN] * out_len,  # random error noise estimated according to Mauder 2013
            'random_flux': [np.NaN] * out_len,         # random flux level estimated by random shuffle criterium (Billesbach 2011)
                      }
        TRACER['calcoef_norm'] = [np.NaN] * out_len     # transmission-corrected calibration coefficient
        TRACER['transmission'] = [np.NaN] * out_len     # transmission

        results['TRACER'] = dict(zip(tuple(str(x) for x in tuple(range(0, hdf5_nb_tracers))),
                                    [deepcopy(TRACER) for i in range(hdf5_nb_tracers)]))  # struct containing results for each tracer

        del TRACER

        # fill mz values
        for n in range(hdf5_nb_tracers):
            results['TRACER'][str(n)]['mz'] = str(hdf5_mz_tracer[n])

        results['TRACER']['cf_lpf'] = [np.NaN] * out_len
    return results
