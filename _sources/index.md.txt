# Welcome to GEddySoft Documentation

GEddySoft (Gembloux Eddy-covariance Software) is a Python package for computing turbulent fluxes from IRGA or PTR-TOF-MS measurements.

The package is approximately 5.8 GB in size, with the vast majority (5.7 GB) consisting of example datasets provided in the `examples` directory. These example files are high-frequency eddy covariance measurements that can be used to test and validate the software's functionality.

Initially inspired by InnFLUX (a MATLAB-based eddy covariance software from the University of Innsbruck, v1.1.0 2021-02-23, https://git.uibk.ac.at/acinn/apc/innflux), GEddySoft began as a Python translation of the original code and has since been enhanced with:
- more flexible input file formats and content handling
- improved quality control (e.g. spike detection and replacement, instrumental problem tests)
- two additional stationarity tests (Mauder and Dutour)
- an additional flux uncertainty assessments (Lenschow metric)
- two additional time lag estimation options (handling of a clock-drift lag on one time series and handling of a humidity dependance)
- improved spectral corrections (computation of spectra in addition to co-spectra and possibility of application of a correction factor)
- multithreading support
- some plotting routines (instrumental tests, lag search, spike detection, (co)spectrum and correction factor)

## Contents

```{toctree}
:maxdepth: 2
:caption: User Guide

user_guide/installation
user_guide/quickstart
user_guide/configuration
user_guide/performance
user_guide/running_example
user_guide/input_files
user_guide/processing
user_guide/outputs
```

```{toctree}
:maxdepth: 2
:caption: Project Info

version
contributing
citation
```

```{toctree}
:maxdepth: 2
:caption: Implementation Details

theory/software_architecture
```

```{toctree}
:caption: API Reference
:maxdepth: 2

api/index
```
