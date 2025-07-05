# Quick Start

Get started with GEddySoft in a few simple steps.

## Basic Usage

1. Set up your input files according to the expected format (see Input Files section)

2. Configure the INI file with your processing parameters (see Configuration section)

3. Run the processing:
   - Open `GEddySoft_start.py` in your Python editor
   - Locate the following line near the top of the file:
     ```python
     ini = 'path/to/your/config.ini'
     ```
   - Replace the path with the full path to your INI file. Example INI files can be found in the `input/metadata` directory (e.g., `GEddySoft_parameters_IRGA.ini` for GHG data). This is the only modification you need to make to the code.
   - Run `GEddySoft_start.py` using Python:
     ```bash
     python GEddySoft_start.py
     ```
   - The script will read your INI file, process the data according to your configuration, and generate output files in the specified directory.

For more detailed information about each step, refer to the relevant sections in the User Guide.
