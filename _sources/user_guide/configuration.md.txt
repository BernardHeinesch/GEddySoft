# Configuration

GEddySoft is configured through an INI file that controls all aspects of the processing. This guide focuses on the most important settings you'll need to adjust.

## INI File Structure

The INI file is organized into sections:

```ini
[version]
# Version information

[files]
# Input/output file paths and settings

[sonic]
# Sonic anemometer configuration

[irga]
# IRGA data configuration

[tracer]
# Tracer data configuration

[meteo]
# Meteorological data configuration

[param]
# Processing parameters

[run_param]
# Runtime configuration
```

## Key Settings

### Concentration Type
In the `[run_param]` section, `CONCENTRATION_TYPE` determines which type of concentration data to process:

```ini
[run_param]
# Choice of concentration files to be processed
# 0 for GHG (IRGA)
# 1 for VOCs (TRACER)
CONCENTRATION_TYPE = 0
```

### Multiprocessing
In the `[run_param]` section, `MULTIPROCESSING` controls parallel processing:

```ini
[run_param]
# Run in multiprocessing mode, one process per day of data
# 0 for conventional processing
# 1 for multiprocessing
MULTIPROCESSING = 0
```

**Important Notes about Multiprocessing:**
- Multiprocessing runs one thread per day of data
- Console output is muted when multiprocessing is enabled
- **Recommended Workflow:**
  1. First test with `MULTIPROCESSING = 0` to verify everything works
  2. Once confirmed working, stop the process
  3. Enable multiprocessing by setting to `1` for faster processing of multiple days
  4. Re-run the processing

## Other Important Settings

In addition to the basic settings above, there are several other important configuration options:

Optional metadata input files can be activated in the .ini:
- meteo (air pressure, air temperature) for a more precise computation of air molar concentration
- tilt_correction for an application of the Planar-Fit rotation method
- lag clock-drift to correct for a computer/datalogger clock drift in tha acquisition of TRACER data
- lag prescribed for the use of a user-defined lag for each time interval
- lag-air relative humidity relation to correct for lag-RH dependency
- low-pass filter correction

### Time Window
```ini
[files]
# (start date, end date), use format yyyy_mm_dd__hh_mm_ss
# Only files within this time range will be processed
# Leave empty if no selection to be applied
date_files_selection = ('2023_06_12__00_00_00','2023_06_14__23_30_00')
```

For a complete list of all available settings and their descriptions, refer to the example INI files provided in the `input/metadata` directory.
