# Changelog

All notable changes to GEddySoft will be documented in this file.

## [4.0.0] - 2025-06-27
### Fixed
- Correction of error introduced when implementing time-series wrapping after rolling (27 June 2025)
- Fixed RH dependent lag detection - added tolerance of 0.001 for matching mz (28 June 2025)
- Improved get_closest_index to search in windows of increasing size for missing/nan/empty values (29 June 2025)

## [3.1.0] - 2025-05-16
### Fixed
- Billebash random uncertainty calculation
- Flux conversion to ug (now using non-protonated mz value)
- Air temperature selection (now correctly using 'meteo' air temperature when available)
- Argument 'detrend' in inst_prob_test call

### Improved
- Enhanced spike detection/replacement (test_spike_detection_vickers97)
- Better subsampling using resample function instead of simple decimation
- Improved prescribed_lag strategy for missing values
- Centralized flux unit conversion

### Added
- Option to skip sonic file processing when tracer file absent
- Multithreading support (one process per day)
- Low-pass filtering corrections with two options:
  - Co-spectral approach with Massman fitted reference
  - Wind-speed based correction factors for stable/unstable conditions
- RH-dependent expected lag
- Clean_results function for output cleanup
- External meteodata support for precise molar air concentration
- Flexible input format for tracer data
- IRGA processing
- Timestamp consistency checking
- Dutaur 1999 stationarity test
- Out of range (OOR) test for U, V, W and T_SONIC
- Spike detection for W, IRGA and TRACER
- Second rotation angle in outputs
- Trunk-space EC support with MEAS_ID
- Planar Fit Method
- Optimized wind direction computation

### Removed
- Disjunct option

## [3.0.0] - 2023
- No changes documented

## [2.0.0] - 2022-06
- First running version on BE-Vie data

## [1.0.0] - 2022-03
- Initial conversion of InnFLUX from MATLAB to Python
- Basic functionality port, not optimized
