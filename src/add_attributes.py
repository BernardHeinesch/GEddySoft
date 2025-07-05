"""Module for adding standardized metadata to GEddySoft HDF5 output files.

This module handles the addition of descriptive attributes to HDF5 datasets
produced by GEddySoft eddy covariance processing. It ensures that all output
files have consistent and well-documented metadata including:

- Variable descriptions
- Physical units
- Quality control flags
- Processing parameters
- Instrument-specific metadata

The metadata follows community standards and best practices for eddy
covariance data, making the outputs compatible with standard analysis
tools and workflows.

Author
------
B. Heinesch
University of Liege, Gembloux Agro-Bio Tech
"""

import h5py


def add_attributes(filepath, filename, process_irga_data_day, n_irga, process_tracer_data_day=False, tracerdata=False):
    """Add standardized metadata attributes to HDF5 output files.

    This function adds descriptive attributes to each dataset in the HDF5
    output file, including variable descriptions and physical units. It
    handles metadata for different data types:

    - Basic file information (creation time, frequencies)
    - Meteorological variables (wind, temperature, pressure)
    - IRGA measurements (if present)
    - Tracer measurements (if present)
    - Quality control metrics
    - Processing parameters

    Parameters
    ----------
    filepath : str
        Directory path containing the HDF5 output file
    filename : str
        Name of the HDF5 output file
    process_irga_data_day : bool
        Flag indicating whether IRGA data was processed
    n_irga : int
        Number of concentration variables from IRGA
    process_tracer_data_day : bool, optional
        Flag indicating whether tracer data was processed
    tracerdata : bool or dict, optional
        If dict, contains tracer-specific metadata

    Returns
    -------
    None
        Function modifies the HDF5 file in place

    Notes
    -----
    The function adds two key attributes to each dataset:
    - 'description': Human-readable description of the variable
    - 'units': Physical units in standard notation

    For quality control flags, additional metadata about test
    thresholds and criteria is included.

    The metadata structure follows the hierarchy:
    - Root level: Basic file information
    - /MET: Meteorological measurements
    - /IRGA: Gas analyzer data (if present)
    - /TRACER: Tracer measurements (if present)

    See Also
    --------
    h5py.File : HDF5 file handling in Python
    """


    with h5py.File(filepath + '\\' + filename, 'r+') as hdf5_f:

        hdf5_f['file_creation_time'].attrs['description'] = 'file_creation_time'; hdf5_f['file_creation_time'].attrs['units'] = 'local computer time'
        hdf5_f['freq'].attrs['description'] = 'frequency axis of co-spectra'; hdf5_f['freq'].attrs['units'] = 's-1'
        hdf5_f['param'].attrs['description'] = 'input parameters'; hdf5_f['param'].attrs['units'] = '-'
        hdf5_f['time'].attrs['description'] = 'timestamp of the end of each averaging interval'; hdf5_f['time'].attrs['units'] = 'yyyy-mm-dd hh-mm-ss'

        hdf5_f['MET']['L'].attrs['description'] = 'Obukhov length'; hdf5_f['MET']['L'].attrs['units'] = 'm'
        hdf5_f['MET']['T'].attrs['description'] = 'temperature'; hdf5_f['MET']['T'].attrs['units'] = 'K'
        hdf5_f['MET']['cospec_wT'].attrs['description'] = 'co-spectrum for wT'; hdf5_f['MET']['cospec_wT'].attrs['units'] = 'Co(w_prime,T_prime,f)'
        hdf5_f['MET']['cospec_wT_scaled'].attrs['description'] = 'scaled co-spectrum for wT'; hdf5_f['MET']['cospec_wT_scaled'].attrs['units'] = 'f*Co(w_prime,T_prime,f)/cov(w_prime,T_prime)'
        hdf5_f['MET']['p'].attrs['description'] = 'air pressure'; hdf5_f['MET']['p'].attrs['units'] = 'hP'
        hdf5_f['MET']['qaqc'].attrs['description'] = 'quality controls'; hdf5_f['MET']['qaqc'].attrs['units'] = '-'
        hdf5_f['MET']['rho_air_molar'].attrs['description'] = 'air molar concentration'; hdf5_f['MET']['rho_air_molar'].attrs['units'] = 'mol m-3'
        hdf5_f['MET']['std_T'].attrs['description'] = 'standard deviation of mean temperature'; hdf5_f['MET']['std_T'].attrs['units'] = 'K'
        hdf5_f['MET']['std_u'].attrs['description'] = 'standard deviation of rotated u wind speed components'; hdf5_f['MET']['std_u'].attrs['units'] = 'm s-1'
        hdf5_f['MET']['std_v'].attrs['description'] = 'standard deviation of rotated v wind speed components'; hdf5_f['MET']['std_v'].attrs['units'] = 'm s-1'
        hdf5_f['MET']['std_w'].attrs['description'] = 'standard deviation of rotated w wind speed components'; hdf5_f['MET']['std_w'].attrs['units'] = 'm s-1'
        hdf5_f['MET']['theta'].attrs['description'] = 'potential temperature'; hdf5_f['MET']['theta'].attrs['units'] = 'K'
        hdf5_f['MET']['theta_v'].attrs['description'] = 'virtual potential temperature'; hdf5_f['MET']['theta_v'].attrs['units'] = 'K'
        hdf5_f['MET']['tilt'].attrs['description'] = 'tilt correction matrix applied to wind vectors'; hdf5_f['MET']['tilt'].attrs['units'] = '-'
        hdf5_f['MET']['phi'].attrs['description'] = 'second rotation angle'; hdf5_f['MET']['phi'].attrs['units'] = 'deg'
        hdf5_f['MET']['ust'].attrs['description'] = 'friction velocity'; hdf5_f['MET']['ust'].attrs['units'] = 'm s-1'
        hdf5_f['MET']['u_unrot'].attrs['description'] = 'unrotated mean u wind speed component'; hdf5_f['MET']['u_unrot'].attrs['units'] = 'm s-1'
        hdf5_f['MET']['v_unrot'].attrs['description'] = 'unrotated mean v wind speed component'; hdf5_f['MET']['v_unrot'].attrs['units'] = 'm s-1'
        hdf5_f['MET']['w_unrot'].attrs['description'] = 'unrotated mean w wind speed component'; hdf5_f['MET']['w_unrot'].attrs['units'] = 'm s-1'
        hdf5_f['MET']['u_rot'].attrs['description'] = 'rotated mean u wind speed component'; hdf5_f['MET']['u_rot'].attrs['units'] = 'm s-1'
        hdf5_f['MET']['v_rot'].attrs['description'] = 'rotated mean v wind speed component'; hdf5_f['MET']['v_rot'].attrs['units'] = 'm s-1'
        hdf5_f['MET']['w_rot'].attrs['description'] = 'rotated mean w wind speed component'; hdf5_f['MET']['w_rot'].attrs['units'] = 'm s-1'
        hdf5_f['MET']['uu'].attrs['description'] = 'auto-covariance of along-wind component'; hdf5_f['MET']['uu'].attrs['units'] = 'm2 s-2'
        hdf5_f['MET']['uv'].attrs['description'] = 'covariance of along-wind and crosswind component'; hdf5_f['MET']['uv'].attrs['units'] = 'm2 s-2'
        hdf5_f['MET']['uw'].attrs['description'] = 'covariance of along-wind and vertical wind component'; hdf5_f['MET']['uw'].attrs['units'] = 'm2 s-2'
        hdf5_f['MET']['vv'].attrs['description'] = 'auto-covariance of crosswind component'; hdf5_f['MET']['vv'].attrs['units'] = 'm2 s-2'
        hdf5_f['MET']['vw'].attrs['description'] = 'covariance of crosswind and vertical wind component'; hdf5_f['MET']['vw'].attrs['units'] = 'm2 s-2'
        hdf5_f['MET']['ww'].attrs['description'] = 'auto-covariance of vertical wind component'; hdf5_f['MET']['ww'].attrs['units'] = 'm2 s-2'
        hdf5_f['MET']['wT'].attrs['description'] = 'temperature flux'; hdf5_f['MET']['wT'].attrs['units'] = 'm s-1 K'
        hdf5_f['MET']['wtheta'].attrs['description'] = 'potential temperature (heat) flux'; hdf5_f['MET']['wtheta'].attrs['units'] = 'm s-1 K'
        hdf5_f['MET']['wtheta_v'].attrs['description'] = 'virtual potential temperature (buoyancy) flux'; hdf5_f['MET']['wtheta_v'].attrs['units'] = 'm s-1 K'
        hdf5_f['MET']['P_air'].attrs['description'] = 'air pressure'; hdf5_f['MET']['P_air'].attrs['units'] = 'hPa'
        hdf5_f['MET']['T_air'].attrs['description'] = 'air temperature'; hdf5_f['MET']['T_air'].attrs['units'] = 'deg C'
        hdf5_f['MET']['RH_air'].attrs['description'] = 'air relative humidity'; hdf5_f['MET']['RH_air'].attrs['units'] = '%'
        hdf5_f['MET']['wdir'].attrs['description'] = 'horizontal wind direction'; hdf5_f['MET']['wdir'].attrs['units'] = 'deg'
        hdf5_f['MET']['std_wdir'].attrs['description'] = 'standard deviation of horizontal wind direction'; hdf5_f['MET']['std_wdir'].attrs['units'] = 'deg'
        hdf5_f['MET']['wsh'].attrs['description'] = 'horizontal wind speed'; hdf5_f['MET']['wsh'].attrs['units'] = 'm s-1'
        hdf5_f['MET']['std_wsh'].attrs['description'] = 'standard deviation of horizontal wind speed'; hdf5_f['MET']['std_wsh'].attrs['units'] = 'm s-1'
        hdf5_f['MET']['zoL'].attrs['description'] = 'stability parameter'; hdf5_f['MET']['zoL'].attrs['units'] = '-'
        hdf5_f['MET']['spec_T'].attrs['description'] = 'spectrum for T'; hdf5_f['MET']['spec_T'].attrs['units'] = 'Sp(T_prime,f)'
        hdf5_f['MET']['spec_T_scaled'].attrs['description'] = 'scaled spectrum for T'; hdf5_f['MET']['spec_T_scaled'].attrs['units'] = 'f*Sp(T_prime,f)/var(T_prime)'

        hdf5_f['MET']['qaqc']['IPT_T'].attrs['description'] = 'instrumental problems tests for T (S_VM97, K_VM97, KID0, KID, HF5, HF10, HD5, HD10, AL1, DDI, DIP, OOR)'; hdf5_f['MET']['qaqc']['IPT_T'].attrs['units'] = '-'
        hdf5_f['MET']['qaqc']['IPT_w'].attrs['description'] = 'instrumental problems tests for w (S_VM97, K_VM97, KID0, KID, HF5, HF10, HD5, HD10, AL1, DDI, DIP, OOR)'; hdf5_f['MET']['qaqc']['IPT_w'].attrs['units'] = '-'
        hdf5_f['MET']['qaqc']['ITC_T'].attrs['description'] = 'relative model deviation of integral turbulence characteristics test for T'; hdf5_f['MET']['qaqc']['ITC_T'].attrs['units'] = '%'
        hdf5_f['MET']['qaqc']['ITC_u'].attrs['description'] = 'relative model deviation of integral turbulence characteristics test for u'; hdf5_f['MET']['qaqc']['ITC_u'].attrs['units'] = '%'
        hdf5_f['MET']['qaqc']['ITC_w'].attrs['description'] = 'relative model deviation of integral turbulence characteristics test for w'; hdf5_f['MET']['qaqc']['ITC_w'].attrs['units'] = '%'
        hdf5_f['MET']['qaqc']['SST_wT_FW96'].attrs['description'] = 'relative deviation in steady-state test for wT according to Foken and Wichura 1996'; hdf5_f['MET']['qaqc']['SST_wT_FW96'].attrs['units'] = '-'
        hdf5_f['MET']['qaqc']['SST_wT_M98'].attrs['description'] = 'relative deviation in steady-state test for wT according to Mahrt 1998'; hdf5_f['MET']['qaqc']['SST_wT_M98'].attrs['units'] = '-'
        hdf5_f['MET']['qaqc']['SST_wT_D99'].attrs['description'] = 'relative deviation in steady-state test for wT according to Dutaur 1999'; hdf5_f['MET']['qaqc']['SST_wT_D99'].attrs['units'] = '-'
        hdf5_f['MET']['qaqc']['completeness_sonic'].attrs['description'] = 'fraction of sonic data used in this averaging interval'; hdf5_f['MET']['qaqc']['completeness_sonic'].attrs['units'] = '%'
        hdf5_f['MET']['qaqc']['num_spikes_w'].attrs['description'] = 'number of spikes detected on w raw data'; hdf5_f['MET']['qaqc']['num_spikes_w'].attrs['units'] = '-'
        hdf5_f['MET']['qaqc']['num_spikes_T'].attrs['description'] = 'number of spikes detected on T raw data'; hdf5_f['MET']['qaqc']['num_spikes_T'].attrs['units'] = '-'

        if process_irga_data_day:

            hdf5_f['IRGA']['cf_lpf'].attrs['description'] = 'flux correction factor for low-pass filtering'; hdf5_f['IRGA']['cf_lpf'].attrs['units'] = '-'

            for i in range(n_irga-1):
                hdf5_f['IRGA'][str(i)]['qaqc']['flux_SNR'].attrs['description'] = 'flux signal-to-noise ratio'; hdf5_f['IRGA'][str(i)]['qaqc']['flux_SNR'].attrs['units'] = '-'
                hdf5_f['IRGA'][str(i)]['cospec'].attrs['description'] = 'co-spectrum for w and tracer concentration'; hdf5_f['IRGA'][str(i)]['cospec'].attrs['units'] = 'Co(w_prime,T_prime,f)'
                hdf5_f['IRGA'][str(i)]['cospec_scaled'].attrs['description'] = 'scaled co-spectrum for w and tracer concentration'; hdf5_f['IRGA'][str(i)]['cospec_scaled'].attrs['units'] = 'f*Co(w_prime,T_prime,f)/cov(w_prime,T_prime)'
                hdf5_f['IRGA'][str(i)]['lagtime_clock_drift'].attrs['description'] = 'clock-drift lag'; hdf5_f['IRGA'][str(i)]['lagtime_clock_drift'].attrs['units'] = 's'
                hdf5_f['IRGA'][str(i)]['lagtime'].attrs['description'] = 'lag time (physical + clock-drift)'; hdf5_f['IRGA'][str(i)]['lagtime'].attrs['units'] = 's'
                hdf5_f['IRGA'][str(i)]['compound'].attrs['description'] = 'mz measured value'; hdf5_f['IRGA'][str(i)]['compound'].attrs['units'] = '-'
                hdf5_f['IRGA'][str(i)]['qaqc'].attrs['description'] = 'quality analysis/quality control tests'; hdf5_f['IRGA'][str(i)]['qaqc'].attrs['units'] = '-'

                hdf5_f['IRGA'][str(i)]['qaqc']['IPT'].attrs['description'] = 'instrumental problems tests for IRGA (S_VM97, K_VM97, KID0, KID, HF5, HF10, HD5, HD10, AL1, DDI, DIP, OOR)'; hdf5_f['IRGA'][str(i)]['qaqc']['IPT'].attrs['units'] = '-'
                hdf5_f['IRGA'][str(i)]['qaqc']['SST_FW96'].attrs['description'] = 'relative deviation in steady-state test for wT according to Foken and Wichura 1996'; hdf5_f['IRGA'][str(i)]['qaqc']['SST_FW96'].attrs['units'] = '%'
                hdf5_f['IRGA'][str(i)]['qaqc']['SST_M98'].attrs['description'] = 'relative deviation in steady-state test for wT according to Mahrt 1998'; hdf5_f['IRGA'][str(i)]['qaqc']['SST_M98'].attrs['units'] = '-'
                hdf5_f['IRGA'][str(i)]['qaqc']['SST_D99'].attrs['description'] = 'relative deviation in steady-state test for wT according to Dutaur 1999'; hdf5_f['IRGA'][str(i)]['qaqc']['SST_D99'].attrs['units'] = '-'
                hdf5_f['IRGA'][str(i)]['qaqc']['completeness_IRGA'].attrs['description'] = 'fraction of IRGA data used in this averaging interval'; hdf5_f['IRGA'][str(i)]['qaqc']['completeness_IRGA'].attrs['units'] = '%'
                hdf5_f['IRGA'][str(i)]['qaqc']['num_spikes'].attrs['description'] = 'number of spikes on raw data'; hdf5_f['IRGA'][str(i)]['qaqc']['num_spikes'].attrs['units'] = '-'
                
                if i == 0:  # CO2
                    hdf5_f['IRGA'][str(i)]['conc_mean'].attrs['description'] = 'mean concentration'; hdf5_f['IRGA'][str(i)]['conc_mean'].attrs['units'] = 'ppm'
                    hdf5_f['IRGA'][str(i)]['conc_std'].attrs['description'] = 'standard deviation of (detrended) tracer concentration'; hdf5_f['IRGA'][str(i)]['conc_std'].attrs['units'] = 'ppm'
                    hdf5_f['IRGA'][str(i)]['flux'].attrs['description'] = 'flux'; hdf5_f['IRGA'][str(i)]['flux'].attrs['units'] = 'umol m-2 s-1'
                    hdf5_f['IRGA'][str(i)]['qaqc']['flux_noise_mean'].attrs['description'] = 'mean flux noise far off the integral timescale'; hdf5_f['IRGA'][str(i)]['qaqc']['flux_noise_mean'].attrs['units'] = 'umol m-2 s-1'
                    hdf5_f['IRGA'][str(i)]['qaqc']['flux_noise_rmse'].attrs['description'] = 'RMSE of flux noise far off the integral timescale'; hdf5_f['IRGA'][str(i)]['qaqc']['flux_noise_rmse'].attrs['units'] = 'umol m-2 s-1'
                    hdf5_f['IRGA'][str(i)]['qaqc']['flux_noise_std'].attrs['description'] = 'standard deviation of flux noise far off the integral timescale'; hdf5_f['IRGA'][str(i)]['qaqc']['flux_noise_std'].attrs['units'] = 'umol m-2 s-1'
                    hdf5_f['IRGA'][str(i)]['qaqc']['random_error_FS'].attrs['description'] = 'random error as described by Finkelstein and Sims (2001);'; hdf5_f['IRGA'][str(i)]['qaqc']['random_error_FS'].attrs['units'] = 'umol m-2 s-1'
                    hdf5_f['IRGA'][str(i)]['qaqc']['random_error_noise'].attrs['description'] = 'random error noise estimated according to Mauder (2013)'; hdf5_f['IRGA'][str(i)]['qaqc']['random_error_noise'].attrs['units'] = 'umol m-2 s-1'
                    hdf5_f['IRGA'][str(i)]['qaqc']['random_flux'].attrs['description'] = 'random flux level estimated by random shuffle criteria (Billesbach, 2011)'; hdf5_f['IRGA'][str(i)]['qaqc']['random_flux'].attrs['units'] = 'umol m-2 s-1'

                elif i == 1:  # H2O
                    hdf5_f['IRGA'][str(i)]['conc_mean'].attrs['description'] = 'mean concentration'; hdf5_f['IRGA'][str(i)]['conc_mean'].attrs['units'] = 'ppt'
                    hdf5_f['IRGA'][str(i)]['conc_std'].attrs['description'] = 'standard deviation of (detrended) tracer concentration'; hdf5_f['IRGA'][str(i)]['conc_std'].attrs['units'] = 'ppt'
                    hdf5_f['IRGA'][str(i)]['flux'].attrs['description'] = 'flux'; hdf5_f['IRGA'][str(i)]['flux'].attrs['units'] = 'mmol m-2 s-1'
                    hdf5_f['IRGA'][str(i)]['qaqc']['flux_noise_mean'].attrs['description'] = 'mean flux noise far off the integral timescale'; hdf5_f['IRGA'][str(i)]['qaqc']['flux_noise_mean'].attrs['units'] = 'mmol m-2 s-1'
                    hdf5_f['IRGA'][str(i)]['qaqc']['flux_noise_rmse'].attrs['description'] = 'RMSE of flux noise far off the integral timescale'; hdf5_f['IRGA'][str(i)]['qaqc']['flux_noise_rmse'].attrs['units'] = 'mmol m-2 s-1'
                    hdf5_f['IRGA'][str(i)]['qaqc']['flux_noise_std'].attrs['description'] = 'standard deviation of flux noise far off the integral timescale'; hdf5_f['IRGA'][str(i)]['qaqc']['flux_noise_std'].attrs['units'] = 'mmol m-2 s-1'
                    hdf5_f['IRGA'][str(i)]['qaqc']['random_error_FS'].attrs['description'] = 'random error as described by Finkelstein and Sims (2001);'; hdf5_f['IRGA'][str(i)]['qaqc']['random_error_FS'].attrs['units'] = 'mmol m-2 s-1'
                    hdf5_f['IRGA'][str(i)]['qaqc']['random_error_noise'].attrs['description'] = 'random error noise estimated according to Mauder (2013)'; hdf5_f['IRGA'][str(i)]['qaqc']['random_error_noise'].attrs['units'] = 'mmol m-2 s-1'
                    hdf5_f['IRGA'][str(i)]['qaqc']['random_flux'].attrs['description'] = 'random flux level estimated by random shuffle criteria (Billesbach, 2011)'; hdf5_f['IRGA'][str(i)]['qaqc']['random_flux'].attrs['units'] = 'mmol m-2 s-1'
                    

        if process_tracer_data_day:
            
            hdf5_f['TRACER']['cf_lpf'].attrs['description'] = 'flux correction factor for low-pass filtering'; hdf5_f['TRACER']['cf_lpf'].attrs['units'] = '-'
            hdf5_f['TRACER']['default_CC_kinetic'].attrs['description'] = 'default_CC_kinetic'; hdf5_f['TRACER']['default_CC_kinetic'].attrs['units'] = '[trncps ppbv^-1]'

            
            # loop on tracers
            for key in hdf5_f["TRACER"].keys():
                if key.isdigit():
                    hdf5_f['TRACER'][key]['calibration'].attrs['description'] = 'concentration calibration coefficent'; hdf5_f['TRACER'][key]['calibration'].attrs['units'] = 'ppbv ncps-1'

                    hdf5_f['TRACER'][key]['conc_mean'].attrs['description'] = 'mean (detrended) concentration'; hdf5_f['TRACER'][key]['conc_mean'].attrs['units'] = 'ppbv'
                    hdf5_f['TRACER'][key]['conc_std'].attrs['description'] = 'standard deviation of (detrended) concentration'; hdf5_f['TRACER'][key]['conc_std'].attrs['units'] = 'ppbv'
                    hdf5_f['TRACER'][key]['conc_acc'].attrs['description'] = 'mean accuracy of (detrended) concentration'; hdf5_f['TRACER'][key]['conc_acc'].attrs['units'] = 'ppbv'
                    hdf5_f['TRACER'][key]['conc_prec'].attrs['description'] = 'mean precision of (detrended) concentration ((sum(conc_prec^2)/n)^0.5)'; hdf5_f['TRACER'][key]['conc_prec'].attrs['units'] = 'ppbv'
                    hdf5_f['TRACER'][key]['conc_LOD'].attrs['description'] = 'limit of detection of (detrended) tracer concentration (3*(sum(zero_prec^2)/n)^0.5)'; hdf5_f['TRACER'][key]['conc_LOD'].attrs['units'] = 'ppbv'
                    hdf5_f['TRACER'][key]['conc_Q5'].attrs['description'] = 'quantile 5 of (detrended) concentration'; hdf5_f['TRACER'][key]['conc_Q5'].attrs['units'] = 'ppbv'
                    hdf5_f['TRACER'][key]['conc_Q95'].attrs['description'] = 'quantile 95 of (detrended) concentration'; hdf5_f['TRACER'][key]['conc_Q95'].attrs['units'] = 'ppbv'
                        
                    hdf5_f['TRACER'][key]['cospec'].attrs['description'] = 'co-spectrum for w and tracer concentration'; hdf5_f['TRACER'][key]['cospec'].attrs['units'] = 'Co(w_prime,T_prime,f)'
                    hdf5_f['TRACER'][key]['cospec_scaled'].attrs['description'] = 'scaled co-spectrum for w and tracer concentration'; hdf5_f['TRACER'][key]['cospec_scaled'].attrs['units'] = 'f*Co(w_prime,T_prime,f)/cov(w_prime,T_prime)'
                    hdf5_f['TRACER'][key]['spec'].attrs['description'] = 'spectrum for tracer concentration'; hdf5_f['TRACER'][key]['spec'].attrs['units'] = 'Sp(T_prime,f)'
                    hdf5_f['TRACER'][key]['spec_scaled'].attrs['description'] = 'scaled spectrum for tracer concentration'; hdf5_f['TRACER'][key]['spec_scaled'].attrs['units'] = 'f*Sp(T_prime,f)/var(T_prime)'
                    hdf5_f['TRACER'][key]['flux'].attrs['description'] = 'flux (corrected for HF losses, see cf_lpf)'; hdf5_f['TRACER'][key]['flux'].attrs['units'] = 'ug m-2 s-1'
                    hdf5_f['TRACER'][key]['lagtime_clock_drift'].attrs['description'] = 'clock-drift lag'; hdf5_f['TRACER'][key]['lagtime_clock_drift'].attrs['units'] = 's'
                    hdf5_f['TRACER'][key]['lagtime'].attrs['description'] = 'lag time (physical + clock-drift)'; hdf5_f['TRACER'][key]['lagtime'].attrs['units'] = 's'
                    hdf5_f['TRACER'][key]['mz'].attrs['description'] = 'mass-to-charge ratio'; hdf5_f['TRACER'][key]['mz'].attrs['units'] = '-'
                    hdf5_f['TRACER'][key]['qaqc'].attrs['description'] = 'quality analysis/quality control tests'; hdf5_f['TRACER'][key]['qaqc'].attrs['units'] = '-'
                    hdf5_f['TRACER'][key]['transmission'].attrs['description'] = 'Transmission coefficient'; hdf5_f['TRACER'][key]['transmission'].attrs['units'] = 'Transmission relative to Tr_21'
                    hdf5_f['TRACER'][key]['Xr0'].attrs['description'] = 'Xr0'; hdf5_f['TRACER'][key]['Xr0'].attrs['units'] = '-'
                    hdf5_f['TRACER'][key]['cluster_min'].attrs['description'] = 'Cluster min'; hdf5_f['TRACER'][key]['cluster_min'].attrs['units'] = 'atomic mass unit'
                    hdf5_f['TRACER'][key]['cluster_max'].attrs['description'] = 'Cluster max'; hdf5_f['TRACER'][key]['cluster_max'].attrs['units'] = 'atomic mass unit'
                    hdf5_f['TRACER'][key]['k_reac'].attrs['description'] = 'Ion/Molecule reaction rate constant'; hdf5_f['TRACER'][key]['k_reac'].attrs['units'] = '1.e-9 cm3 molecule-1 s-1'
                    hdf5_f['TRACER'][key]['FY'].attrs['description'] = 'Fragmentation yield to correct signal due to fragmentation in the drift tube'; hdf5_f['TRACER'][key]['FY'].attrs['units'] = '-'
                    hdf5_f['TRACER'][key]['IF'].attrs['description'] = 'Isotopic factor to correct isotopic ratio'; hdf5_f['TRACER'][key]['IF'].attrs['units'] = '-'
    
                    hdf5_f['TRACER'][key]['qaqc']['IPT'].attrs['description'] = 'instrumental problems tests for TRACER (S_VM97, K_VM97, KID0, KID, HF5, HF10, HD5, HD10, AL1, DDI, DIP, OOR)'; hdf5_f['TRACER'][key]['qaqc']['IPT'].attrs['units'] = '-'
                    hdf5_f['TRACER'][key]['qaqc']['SST_FW96'].attrs['description'] = 'relative deviation in steady-state test for wc according to Foken and Wichura 1996'; hdf5_f['TRACER'][key]['qaqc']['SST_FW96'].attrs['units'] = '%'
                    hdf5_f['TRACER'][key]['qaqc']['SST_M98'].attrs['description'] = 'relative deviation in steady-state test for wc according to Mahrt 1998'; hdf5_f['TRACER'][key]['qaqc']['SST_M98'].attrs['units'] = '-'
                    hdf5_f['TRACER'][key]['qaqc']['SST_D99'].attrs['description'] = 'relative deviation in steady-state test for wc according to Dutaur 1999'; hdf5_f['TRACER'][key]['qaqc']['SST_D99'].attrs['units'] = '-'
                    hdf5_f['TRACER'][key]['qaqc']['completeness_TRACER'].attrs['description'] = 'fraction of tracer data used in this averaging interval'; hdf5_f['TRACER'][key]['qaqc']['completeness_TRACER'].attrs['units'] = '%'
                    hdf5_f['TRACER'][key]['qaqc']['num_spikes'].attrs['description'] = 'number of spikes on raw data'; hdf5_f['TRACER'][key]['qaqc']['num_spikes'].attrs['units'] = '-'
                    hdf5_f['TRACER'][key]['qaqc']['flux_SNR'].attrs['description'] = 'flux signal-to-noise ratio'; hdf5_f['TRACER'][key]['qaqc']['flux_SNR'].attrs['units'] = '-'
                    hdf5_f['TRACER'][key]['qaqc']['flux_noise_mean'].attrs['description'] = 'mean flux noise far off the integral timescale'; hdf5_f['TRACER'][key]['qaqc']['flux_noise_mean'].attrs['units'] = 'ug m-2 s-1'
                    hdf5_f['TRACER'][key]['qaqc']['flux_noise_rmse'].attrs['description'] = 'RMSE of flux noise far off the integral timescale'; hdf5_f['TRACER'][key]['qaqc']['flux_noise_rmse'].attrs['units'] = 'ug m-2 s-1'
                    hdf5_f['TRACER'][key]['qaqc']['flux_noise_std'].attrs['description'] = 'standard deviation of flux noise far off the integral timescale'; hdf5_f['TRACER'][key]['qaqc']['flux_noise_std'].attrs['units'] = 'ug m-2 s-1'
                    hdf5_f['TRACER'][key]['qaqc']['num_spikes'].attrs['description'] = 'number of spikes (tbd)'; hdf5_f['TRACER'][key]['qaqc']['num_spikes'].attrs['units'] = '-'
                    hdf5_f['TRACER'][key]['qaqc']['random_error_FS'].attrs['description'] = 'random error as described by Finkelstein and Sims (2001);'; hdf5_f['TRACER'][key]['qaqc']['random_error_FS'].attrs['units'] = 'ug m-2 s-1'
                    hdf5_f['TRACER'][key]['qaqc']['random_error_noise'].attrs['description'] = 'random error noise estimated according to Mauder (2013)'; hdf5_f['TRACER'][key]['qaqc']['random_error_noise'].attrs['units'] = 'ug m-2 s-1'
                    hdf5_f['TRACER'][key]['qaqc']['random_flux'].attrs['description'] = 'random flux level estimated by random shuffle criteria (Billesbach, 2011)'; hdf5_f['TRACER'][key]['qaqc']['random_flux'].attrs['units'] = 'ug m-2 s-1'
            
    return None
