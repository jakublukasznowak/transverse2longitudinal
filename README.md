# longi_vs_trans

This code was written by Jakub L. Nowak for the analysis of the transverse-to-longitudinal ratio of 2nd order turbulent velocity statistics (structure functions and power spectra) derived from aircraft measurements.

There are 4 field experiments and 3 aircrafts sampling 2 marine boundary layer regimes included:
- EUREC4A - ATR42 aircraft - shallow trade-wind convection
- RICO - C130 aircraft - shallow trade-wind convection
- VOCALS-REx - C130 aircraft - subtropical stratocumulus
- POST - Twin Otter aircraft - subtropical stratocumulus

The analysis requires access to the datasets which can be downloaded from the public data repositories:
- [Lothon, M. & Brilouet, P. (2020). SAFIRE ATR42: Turbulence Data 25 Hz. Aeris. doi.org/10.25326/128](https://observations.ipsl.fr/aeris/eurec4a-data/AIRCRAFT/ATR/SAFIRE-TURB/PROCESSED/) (here L3 v1.9 is used)
- [NCAR/NSF C-130 Navigation, State Parameter, and Microphysics HRT Data - 25 Hz. Version 1.0. doi.org/10.5065/D64J0CDM](https://data.eol.ucar.edu/dataset/87.049)
- [NCAR/NSF C-130 Navigation, State Parameter, and Microphysics HRT (25 sps) Data. doi.org/10.5065/D69K48JK](https://data.eol.ucar.edu/dataset/89.002)
- [NSF/NCAR C130 Radar, Lidar and Radiometer Integrated Dataset. doi.org/10.26023/8KEJ-BQNG-W808](https://data.eol.ucar.edu/dataset/89.159)
- [UC Irvine 40-hz Probes - netCDF format. doi.org/10.26023/KP56-KFJS-VC07](https://data.eol.ucar.edu/dataset/111.033)
and Table 1 from Carman, at al. 2012: Observational constraints on entrainment and the entrainment interface layer in stratocumulus(https://doi.org/10.5194/acp-12-11135-2012) which can be found here in 'aux_data'.


The code was developed in MATLAB R2019b. The functionality in other versions of this environment was not tested.

Two external packages are used:
- YAML 1.1 parser and emitter for MATLAB by Martin Koch from https://www.mathworks.com/matlabcentral/fileexchange/106765-yaml
- boxplotGroup by Adam Danz from https://www.mathworks.com/matlabcentral/fileexchange/74437-boxplotgroup

