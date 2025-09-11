# Example Dataset

A test dataset is provided with GEddySoft, containing 5 days of measurements from the Vielsalm forest ICOS station (BE-Vie) collected in 2023. The dataset can be downloaded from [https://doi.org/10.18758/iypn6fnm](https://doi.org/10.18758/iypn6fnm) and includes:

## Available Data Files

- Sonic anemometer data: `GHS50_yyyy_mm_dd__HH_MM_SS.hdf5`, 277 Mb
- PTR-TOF-MS measurements: `K8_BE-Vie_yyyy_TOF4000_yyyy_mm_dd_HH_MM_SS.h5`, 1100 Mb
- GHG data from SmartFlux: `yyyy-mm-ddTHHMMSS_OTV.ghg`, 709 Mb

## Processing Examples

You can test GEddySoft processing capabilities using:

1. **BVOC Flux Processing**
   - Use the configuration file: `GEddySoft_parameters_VOC_top.ini`
   - This will process the PTR-TOF-MS data (`K8_BE-Vie_*.h5` containing signals for 20 masses) and sonic data (`GHS50_*.hdf5`) for BVOC flux calculations
   - Typical running time: 14 minutes (on a i7 processor, without the multiprocessing option activated in the ini)

2. **CO2/H2O Flux Processing**
   - Use the configuration file: `GEddySoft_parameters_IRGA.ini`
   - This will process the SmartFlux GHG data (`*_OTV.ghg`) for CO2 and H2O flux calculations
   - Typical running time: 9 minutes (on a i7 processor, without the multiprocessing option activated in the ini)
