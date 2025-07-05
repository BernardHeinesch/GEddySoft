# Welcome to GEddySoft Documentation

GEddySoft (Gembloux Eddy-covariance Software) is a Python package for computing turbulent fluxes from IRGA or PTR-TOF-MS measurements.

Initially inspired by InnFLUX (a MATLAB-based eddy covariance software from the University of Innsbruck, v1.1.0 2021-02-23, https://git.uibk.ac.at/acinn/apc/innflux), GEddySoft began as a Python translation of the original code and has since been enhanced with:
- More flexible input file formats and content handling
- Improved quality control (spike detection, IPT tests, etc.)
- Flux uncertainty assessments
- Multiple timelag estimation options (e.g. smoothing of the covariance function, possibility of handling a clock-drift lag on one time series)
- Spectral correction for high-frequency losses
- Multithreading support
- Plotting routines

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
