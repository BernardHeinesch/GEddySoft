# Installation

## Requirements

GEddySoft requires Python 3.11.7 or higher and the following packages:
- hdfdict
- diptest
- Other dependencies listed in requirements.txt:

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
