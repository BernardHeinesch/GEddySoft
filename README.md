# GEddySoft

Project to compute turbulent fluxes for IRGA or PTR-TOF-MS

## Name
GEddySoft (Gembloux Eddy-covariance Software)

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
                    
    v4.0.1 2025-08-16
    -----------------
    Corrected mistakes:
    - 16 August 2025: correction on conc_prec and conc_LOD computation. Division by number of samples must be out of the sqrt.
