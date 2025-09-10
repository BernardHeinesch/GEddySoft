# Installation

## Requirements

GEddySoft requires Python 3.7 or higher and the following dependencies:

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

## Environment Setup

GEddySoft requires Python 3.7 or higher. We recommend using a dedicated Python environment to avoid conflicts with other packages. You can choose either of these approaches:

### Using Anaconda (Recommended for Scientific Users)

Anaconda provides a user-friendly way to manage Python environments and includes many scientific packages by default.

1. Install Anaconda from [https://www.anaconda.com/download](https://www.anaconda.com/download)
2. Create a new environment for GEddySoft:
   ```bash
   conda create -n geddysoft python=3.11  # or any version >=3.7
   conda activate geddysoft
   ```

### Using Python venv

Python's built-in `venv` module is a lighter alternative if you prefer not to install Anaconda.

1. Ensure you have Python >=3.7 installed
2. Create and activate a virtual environment:
   ```bash
   # On Windows
   python -m venv geddysoft-env
   .\geddysoft-env\Scripts\activate

   # On Linux/MacOS
   python3 -m venv geddysoft-env
   source geddysoft-env/bin/activate
   ```

## Installation Steps

Once you have set up and activated your Python environment (see Environment Setup above):

1. Clone the repository:
   ```bash
   git clone https://github.com/BernardHeinesch/GEddySoft.git
   cd GEddySoft
   ```

2. Install required packages:
   ```bash
   pip install -r requirements.txt
   ```
   Watch the console output during installation to ensure there are no error messages.

3. Install GEddySoft:
   ```bash
   pip install .
   ```

   Note: For developers who plan to modify the source code, use the development mode instead:
   ```bash
   pip install -e .
   ```

## Verification

To verify the installation:

1. Check Python version:
   ```bash
   python --version  # Should be 3.7 or higher
   ```

## Troubleshooting

If you encounter any issues during installation:

1. Ensure you have Python 3.7 or higher
2. Check that all required dependencies are installed
3. Verify your system meets the minimum requirements
4. Contact the maintainer if problems persist
