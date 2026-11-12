# Quick Start

Get started with GEddySoft in a few simple steps.

## Basic Usage

1. Use example data:
- Download example data from [https://doi.org/10.18758/iypn6fnm](https://doi.org/10.18758/iypn6fnm) (see 'Example Dataset' page for some info on this test dataset) 
- unzip the downloaded file in a folder of your choice (e.g. \examples\data)

2. Configure the INI file with your processing parameters (see Configuration section)

3. Run the processing:
   - Open `src/GEddySoft_start.py` in your Python editor
   - Locate the following line near the top of the file:
     ```python
     ini = 'path/to/your/config.ini'
     ```
   - Replace this path with the full path to your INI file. Example INI files can be found in the `examples/metadata` directory (e.g., `GEddySoft_parameters_VOC_top.ini` for VOC data). This is the only modification you need to make to the code.
   - adapt the 'sonic_files_folders' and 'tracer_files_folders' fields of the ini to point to the folders where you downloaded the data 
   - Run `GEddySoft_start.py` using Python from the `src` directory.
   - The script will read your INI file, process the data according to your configuration, and generate output files in the specified directory.

For more detailed information about each step, refer to the Implementation Details and API Reference sections.
