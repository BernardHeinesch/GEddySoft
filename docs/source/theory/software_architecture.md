# Software Architecture

GEddySoft follows a modular architecture designed for efficient eddy covariance data processing. The software is organized into core processing modules and utility functions, with a clear data flow between components.

## Processing Flow

```
GEddySoft_start
    |
    v
GEddySoft_main
    |
    |-> Input Processing
    |     |- read_main_inputs
    |     |- read_metadata_files
    |     |- get_list_input_files
    |          |
    |          v
    |-> Raw Data Loading
    |     |- read_sonic -> check_raw_timestamps
    |     |- read_GHG -> map_sonic2tracer
    |          |
    |          v
    |-> Quality Control
    |     |- spike_detection_vickers97
    |     |- inst_prob_test
    |     |- test_steady_state_FW96
    |     |- test_ITC
    |          |
    |          v
    |-> Flux Calculation
    |     |- wind_rotations
    |     |- xcov
    |     |- find_covariance_peak
    |     |- flux_unit_conversion
    |          |
    |          v
    |-> Analysis & Corrections
    |     |- cospectrum
    |     |- correction_factor_lpf
    |     |- flux_uncertainties
    |          |
    |          v
    |-> Output Generation
          |- prepare_output
          |- clean_results
          |- add_attributes
```

## Module Categories

### Core Processing Modules

These modules implement essential eddy covariance algorithms and scientific calculations:

1. **Main Processing**
   - `GEddySoft_main`: Central processing workflow
   - `flux_uncertainties`: Uncertainty calculations
   - `cospectrum`: Spectral analysis

2. **Quality Control**
   - `spike_detection_vickers97`: Despiking algorithm
   - `test_steady_state_FW96`: Stationarity test
   - `test_steady_state_M98`: Alternative stationarity test
   - `test_steady_state_D99`: Alternative stationarity test
   - `test_ITC`: Integral turbulence characteristics

3. **Flux Processing**
   - `flux_unit_conversion`: Unit conversions
   - `wind_rotations`: Coordinate rotations
   - `correction_factor_lpf`: Low-pass filtering corrections

### Utility Functions

These modules provide supporting functionality:

1. **Data Handling**
   - `clean_results`: Result filtering
   - `get_closest_value`: Value interpolation
   - `nandetrend`: NaN-aware detrending
   - `nanlinfit`: NaN-aware linear fitting

2. **File Operations**
   - `get_ini`: Configuration parsing
   - `get_list_input_files`: File discovery
   - `read_GHG`: Gas analyzer data reading
   - `read_main_inputs`: Input file parsing

3. **Data Transformations**
   - `compute_wind_direction`: Wind vector processing
   - `correction_factor_lpf`: Spectral corrections
   - `logBinSpectrum`: Spectral binning

4. **Metadata Handling**
   - `add_attributes`: HDF5 attribute management
   - `read_metadata_files`: Metadata parsing

5. **Support Functions**
   - `prepare_output`: Output formatting
   - `print_progress`: Progress reporting
   - `reformat_date`: Date string handling

## Data Flow

1. **Input Stage**
   - Configuration loading
   - File discovery
   - Raw data reading

2. **Preprocessing Stage**
   - Timestamp validation
   - Data synchronization
   - Quality control

3. **Processing Stage**
   - Coordinate rotations
   - Covariance calculation
   - Spectral corrections
   - Uncertainty estimation

4. **Output Stage**
   - Result cleaning
   - Metadata addition
   - HDF5 file generation

## Performance Considerations

- Parallel processing for multiple days
- Memory-efficient array operations
- Robust error handling and logging
- Progress tracking for long operations
