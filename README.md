# GEddySoft

Project to compute turbulent fluxes for IRGA or PTR-TOF-MS

## Name
GEddySoft (Gembloux Eddy-covariance Software)

## Documentation
The HTML documentation is available online via GitHub Pages:
[https://bernardheinesch.github.io/GEddySoft/](https://bernardheinesch.github.io/GEddySoft/)

## Description
Initially inspired by InnFLUX, which is an an eddy covariance software package developed
at the University of Innsbruck, Austria by Marcus Striednig et al.
(innFLUX - an open-source code for conventional and disjunct eddy covariance
analysis of trace gas measurements: an urban test case,
Atmos. Meas. Tech., 13, 1447–1465, https://doi.org/10.5194/amt-13-1447-2020,
2020.)

Translated to python by Bernard Heinesch (University of Liege, Belgium), starting in March 2022,
using Python 3.11.7 64-bit | Qt 5.9.7 | PyQt5 5.9.2 | Windows 10,
starting from innFLUX v1.1.0 2021-02-23, downloaded on github at https://git.uibk.ac.at/acinn/apc/innflux.

Improvements include 
- more flexibility in input file formats and imperfect content
- addition of QAQC like spike detection and replacement, IPT tests, etc)
- addition of some flux uncertainty assesments
- addition of several options for timelag estimation
- addition of a spectral correction for high-frequency losses
- addition of some plotting routines
- addition of a multithreading option for improved performances

Description of the different subroutines/features:

1. GEddySoft_start.py: starting routine, reading the ini and allowing the possibility of parallel processing on days 

2. GEddySoft_main.py: main routine

3. a ini file is gathering all user-defined choices. 
   The only user setting, besides this ini, is the ini file name to be provided in GEddySoft_start. 

4. mandatory input data files are
    - sonic/IRGA GHG format data read by the read_main_inputs.py routine
    or
    - sonic data and TRACER data (meaning VOC data from the TOF) read also by the read_main_inputs.py routine
      If TRACER data are not provided, sonic data can be processed but not the reverse.
    GEddySoft cannot be run on both the IRGA and TRACER. You need to do separate runs.
    A limited flexibility in the input format is allowed through paramtrization in the .ini file

   optional metadata input files can be activated in the .ini:
    - meteo (air pressure, air temperature) for a more precise computation of air molar concentration
    - tilt_correction for an application of the Planar-Fit rotation method
    - lag clock-drift to correct for a computer/datalogger clock drift in tha acquisition of TRACER data
    - lag prescribed for the use of a user-defined lag for each time interval
    - lag-air relative humidity relation to correct for lag-RH dependency
    - low-pass filter correction

## Performances

Typical running time for a SONIC file at 50 Hz and TRACER file at 10Hz with 80 ions,
all CPU-consuming options being activated, is 
- 15 s/file (so per half-hour) in the multiprocessing mode, on a computer 
  11th Gen Intel(R) Core(TM) i7-1165G7 @ 2.80GHz with 8 cores 
- about 3 times more in the non-multiprocessing mode

## Visuals

## Installation
The codes must be run on a Python environment.
If using Anaconda, make sure you have installed the following libraries on your Python environment:
* hdfdict
* diptest


## Contributing

For contributions, please consider working on a separate branch.

Once the user wants to contribute to the main branch, please create a merge request indicating @Bernard Heinesch as reviewer.


## Versioning
    v1.0 2022-03-xx
    ---------------
    Conversion of InnFLUX from MATLAB to Python. No modification brought to the functionnalities.
    This version was not running nor optimized.
    
    v2.0 2022-06-xx
    ---------------
    First running version on BE-Vie data
    
    v3.0 2023-0x-xx
    ---------------
    - /

    v3.1 2025-05-16
    ---------------
    Corrected mistakes:
    - Mistake in Billebash random uncertainty corrected 
    - Mistake in flux conversion to ug corrected (the protonated mz value was used instead of the non-protonated one) 
    - Mistake in air temperature choice corrected: sonic temperature was always taken for improved molar air density 
      computation, instead of 'meteo' air temperature when available 
    - Mistake in call to inst_prob_test corrected: argument detrend was set arbitrarily  

    Improved existing features: 
    - improved spikes detection/replacement (detect_num_spikes replaced by test_spike_detection_vickers97)
    - Subsampling of sonicdata from SAMPLING_RATE_SONIC to SAMPLING_RATE_FINAL is now made through
      the resample function instead of taking one sample over x
    - Prescribed_lag strategy adapted (in case no value available, median of the 10 closest values)
    - centralized flux unit conversion 

    New functionnalities:
    - Option implemented to avoid processing sonic files if tracer file not present 
    - Multithreading implemented, one process per day if option activated in the ini
    - Corrections for low-pass filtering implemented with two options available:
      one using a cut-off frequency and the co-spectral approach with a Massman fitted reference co-spectrum 
      one using a correction factor vs wind speed table, separately for unstable and stable atmospheric conditions
    - RH-dependent expected lag implemented
    - Clean_results def implemented, for cleaning unuseful lines from the final output
    - Option implemented to use external meteodata (air pressure and air temperature)  
      for a more precise computation of molar air concentration and therefore fluxes 
    - Creation of the read_main_inputs def and addition of more flexibility on input format for tracer
    - Improvements on the GEddySoft_parameters routine and replaced by an ini text file using configparser
    - IRGA processing implemented
    - Check on the consistency of the timestamp in the raw input data files implemented
    - Dutaur 1999 stationarity test implemented
    - Out of range test (OOR) implemented, for U, V, W and T_SONIC
    - Spikes detection (num_spikes) implemented for W, IRGA and TRACER
    - Second rotation angle added to the outputs
    - File selection adapted to handle also the trunk-space EC and a MEAS_ID added to discriminate 
      top and trunk output files
    - Planar Fit Method option implemented (starting from the EddyPro PFM file)
    - Optimization of the compute_wind_direction routine through vectorization
    - Disjunct option removed everywhere
    
    v4.0 2025-06-27
    ---------------
    Corrected mistakes:
    - 27 June 2025: correction of the error introduced in v4.0 when implementing time-series wrapping after rolling
    - 28 June 2025: correction for rh dependent lag : it was never found due to string-type search of mz. Now tolerance of 0.001 for matching
    - 29 June 2025: correction for get_closest_index : in case of missing/nan/empty value, added the search in a window of increasing size 
                    (not only 10 values), until non nan/empty values are found
                    
    v4.1 2026-07-16
    ---------------
    - Added configurable UTC offset through UTC_OFFSET parameter in ini settings
    - Simplifications brought to date handling in read_main_inputs
    - Added one accepted input file format (csv from QCL)
    - Added IMPUTATION_METHOD option (0/1/2) for TRACER concentration alignment on the SAMPLING_RATE_FINAL grid. Method 0 maps only measured TRACER timestamps (no imputation), method 1 uses FFT resampling, method 2 uses plateau (zero-order hold) upsampling.
    - Refactored lag/flux computation: CONST/PRESCRIBED flux uses xcov at a single lag; MAX and MAX_WITH_DEFAULT share a common code path. Def 'compute_time_lag also created, to improve modularity.

    Corrected mistakes:
    - 16 August 2025: correction on conc_prec and conc_LOD computation. Division by number of samples must be out of the sqrt.
    - w_prime_trim was used instead of w_prime in the IRGA part when implementing time-series wrapping after rolling. Idem for c_prime_trim.
    - 23 December 2025: added correction for low-pass filtering on flux uncertainties for TRACER : flux uncertainties were not corrected while the flux was.
    - 26 December 2025: added correction for prescribed lag which was corrected for clock-drift altough the input file is physical + clock drift lag
    - 26 December 2025: make time-zone aware all the timestamps for the metadata input files
    - 28 May 2026: fixed spike_detection_vickers97 multi-pass despiking (now updates the working signal between passes) and improved invalid-data handling (NaN/Inf via error_value)
    - 29 May 2026: wind_rotations: replaced arctan by arctan2 for yaw/pitch (correct quadrant + improved stability); no impact expected on key outputs (mean wind speed/direction and scalar fluxes), except rare edge cases when u_mean≈0.

    v4.2 2026-11-12
    ---------------

    Improved existing features:
    - GEddySoft_start (parallel processing): replaced the multiprocessing.Pool backend by joblib.Parallel (loky backend), which provides simpler code, live progress report, worker stdout/stderr forwarding to the console, and more robust worker management. joblib is now a required dependency.
    - GEddySoft_main (T_prime): add warning when NaN values are detected in T_prime before replacing by 0.
    - compatibility: replaced np.NaN by np.nan across the code, integration routine naming: replaced trapz by scipy.integrate.trapezoid, replaced numpy.string_ by numpy.bytes_ in hdfdict, suggestions by Kim Hoang.

    New functionalities:
    - GEddySoft_start: can now be launched from command line with an optional ini filepath argument (argparse). If no argument is provided, the example ini path in this repository is used.
    - Added standalone daily HDF5-to-NetCDF concatenation utility (concat_hdf5.py).

    Corrected mistakes:
    - compute_time_lag(MAX_WITH_DEFAULT): corrected edge fallback conditions (off-by-one) to match the lag range returned by find_covariance_peak. Bug detected by Pedro Coimbra.
    - xcov: lag=0 now uses the full arrays (no dropped sample) and negative lags now trim the wrapped samples correctly. Bug detected by Pedro Coimbra. Measured effect at a fixed lag (~-110 samples in an 18000-sample window): median 0.000 %, 90th percentile ~3 %, worst 27 % — negligible on strong fluxes, large on weak ones, and it also feeds the peak search.
    - GEddySoft_main (wT): 'wT' computation was using 'xcov(w_prime, T_prime,[0, 0])' which drops the last sample. Now with np.cov, same estimator as other sonic covariances, and without dropping the last sample. That one-sample split is the whole of a measured 0.067 % disagreement on wT while every neighbouring sonic quantity reproduces to ~1e-13 %. Bug detected by Pedro Coimbra.
    - GEddySoft_main (sonic despiking): the two sonic `spike_detection_vickers97` calls (w and T_sonic) passed `ac_freq = SAMPLING_RATE_SONIC` instead of SAMPLING_RATE_FINAL (after resampling). Bug detected by Pedro Coimbra.
    - GEddySoft_main (TRACER conc_std ddof): tracer conc_std was using ddof=0 instead of ddof=1. Worth `sqrt(N/(N-1)) - 1` = 0.0028 % at N = 17900. Bug detected by Pedro Coimbra.
    - GEddySoft_main (PTR-TOF-MS neutral mass convention): for conversion of TRACER fluxes to ug m-2 s-1 (neutrals), molar mass is now derived as `mz - 1.00728` (proton mass) instead of `mz - 1`. Suggestion by Pedro Coimbra.
    - nanlinfit: non-finite (NaN/Inf) samples were deleted before building the fit x-axis (`t = 0..N_valid-1`), so the slope was computed per valid sample instead of per real sample index — overestimating dT_dt/dc_dt when multiplied by the sampling rate, and inconsistent with nandetrend which evaluates the fit on the original indices. Now the fit is performed on the original sample indices (isfinite mask). Bug detected by Matthias Cuntz.
    - GEddySoft_main (air_mol_conc): flux unit conversion used the moist-air molar concentration p/(R*T) instead of the dry-air concentration (p-ea)/(R*T) — appropriate for dry-air mole fractions (QCL/PTR-TOF-MS convention). All converted fluxes were overestimated by ea/p (~1 % typical, up to ~4 % hot & humid). Water vapour pressure ea is now estimated with a fallback chain: (1) meteo-file relative humidity via new function esat (Magnus formula), (2) measured H2O mixing ratio, (3) ea=0 dry-air assumption (identical to previous behaviour). The fixed standard molar-volume constant (`22.414e-3` m3 mol-1) previously used in the conversion is removed as a consequence. Bug detected by Matthias Cuntz.
    - read_GHG (read_diag_val, LI-7200 diagnostics): 16-bit flag conversion `'000' + format(x, 'b')` gave wrong bit alignment whenever the value had fewer than 13 significant bits. Now `format(x, '016b')`. Bug detected by Matthias Cuntz.