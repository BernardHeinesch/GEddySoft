# Running on Example Dataset

A test dataset is provided with GEddySoft, containing 3 days of measurements from the Vielsalm forest ICOS station (BE-Vie) collected in 2023. The dataset includes:

## Available Data Files

- Sonic anemometer data: `GHS50_yyyy_mm_dd__HH_MM_SS.hdf5`
- PTR-TOF-MS measurements: `K8_BE-Vie_yyyy_TOF4000_yyyy_mm_dd_HH_MM_SS.h5`
- GHG data from SmartFlux: `yyyy-mm-ddTHHMMSS_OTV.ghg`

## Processing Examples

You can test GEddySoft processing capabilities using:

1. **BVOC Flux Processing**
   - Use the configuration file: `GEddySoft_parameters_VOC_top.ini`
   - This will process the PTR-TOF-MS data (`K8_BE-Vie_*.h5`) and sonic data (`GHS50_*.hdf5`) for BVOC flux calculations

2. **CO2/H2O Flux Processing**
   - Use the configuration file: `GEddySoft_parameters_IRGA`
   - This will process the SmartFlux GHG data (`*_OTV.ghg`) for CO2 and H2O flux calculations
