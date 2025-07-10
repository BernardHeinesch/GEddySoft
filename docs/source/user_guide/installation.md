# Installation

## Requirements

GEddySoft requires Python 3.11.7 or higher and the following dependencies:

### Core Dependencies
All core dependencies are listed in requirements.txt and will be automatically installed:
- numpy: For numerical computations
- pandas: For data handling and time series operations
- scipy: For signal processing and statistical functions
- matplotlib: For plotting routines
- h5py: For HDF5 file operations
- hdfdict: For HDF5 dictionary operations
- diptest: For statistical tests

### Documentation Dependencies
These are only needed if you want to build the documentation locally:
- sphinx
- myst-parser
- sphinx-rtd-theme

All dependencies with their version requirements are listed in requirements.txt:

```{include} ../../../requirements.txt
:code: text
```

## Installation Steps

1. Clone the repository:
   ```bash
   git clone https://github.com/BernardHeinesch/GEddySoft.git
   cd GEddySoft
   ```

2. Create and activate a virtual environment (recommended):
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```

3. Install required packages:
   ```bash
   pip install -r requirements.txt
   ```

4. Install GEddySoft in development mode:
   ```bash
   pip install -e .
   ```

## Verification

To verify the installation:

1. Check Python version:
   ```bash
   python --version  # Should be 3.11.7 or higher
   ```

2. Try importing key dependencies:
   ```python
   import numpy
   import pandas
   import h5py
   import hdfdict
   import diptest
   ```

## Troubleshooting

If you encounter any issues during installation:

1. Ensure you have Python 3.11.7 or higher
2. Check that all required dependencies are installed
3. Verify your system meets the minimum requirements
4. Contact the maintainer if problems persist
