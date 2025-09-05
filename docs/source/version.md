# Version History

## v4.0.2 (2025-08-16)

### Improvements
- Added configurable UTC offset through UTC_OFFSET parameter in ini settings
- Simplifications brought to date handling in read_main_inputs

### Corrected Mistakes
- w_prime_trim was used instead of w_prime in the IRGA part when implementing time-series wrapping after rolling. Idem for c_prime_trim.

## v4.0.1 (2025-08-16)

### Corrected Mistakes
- correction on conc_prec and conc_LOD computation. Division by number of samples is now out of the sqrt

## v4.0 (2025-06-27)

### Corrected Mistakes
- Correction of the error introduced in v3.1 when implementing time-series wrapping after rolling
- Correction for RH dependent lag - fixed string-type search of mz with 0.001 tolerance for matching

### Improvements
- Improved get_closest_index to search in windows of increasing size for missing/nan/empty values

## v3.1 (2025-05-16)

### Corrected Mistakes
- Fixed Billebash random uncertainty calculation
- Fixed flux conversion to μg (now using non-protonated mz value)
- Fixed air temperature selection for improved molar air density computation
- Fixed inst_prob_test argument detrend setting

### Improvements
- Enhanced spike detection/replacement (using test_spike_detection_vickers97)
- Improved sonic data subsampling using resample function
- Adapted prescribed_lag strategy for missing values
- Centralized flux unit conversion

### New Features
- Option to skip sonic files without tracer file
- Multithreading support (one process per day)
- Low-pass filtering corrections with two options:
  - Cut-off frequency with Massman fitted reference co-spectrum
  - Wind speed based correction factor for stable/unstable conditions
- RH-dependent expected lag
- Clean_results function for output cleanup
- External meteodata support for precise molar air concentration
- More flexible input format for tracer data
- Replaced parameters routine with configparser-based INI file
- IRGA processing support
- Input timestamp consistency checks
- Additional quality tests:
  - Dutaur 1999 stationarity test
  - Out of range test for U, V, W and T_SONIC
  - Spike detection for W, IRGA and TRACER
- Second rotation angle in outputs
- Support for trunk-space EC with MEAS_ID
- Planar Fit Method from EddyPro PFM file
- Optimized wind direction computation
- Removed disjunct option

## v2.0 (2022-06)
First running version on BE-Vie data

## v1.0 (2022-03)
Initial conversion of InnFLUX from MATLAB to Python. No functionality modifications.
This version was not running nor optimized.
