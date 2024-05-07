The repository contains MATLAB code developed for the purpose of the analysis presented in the manuscript

*The ratio of transverse to longitudinal turbulent velocity statistics for aircraft measurements*

by Jakub L. Nowak, Marie Lothon, Donald H. Lenschow, and Szymon P. Malinowski


## Introduction

The classical theory of homogeneous isotropic turbulence predicts the ratio of transverse to longitudinal structure functions or power spectra equal to 4/3 in the inertial subrange. For the typical turbulence cascade in the inertial subrange, it also predicts a power law scaling with an exponent of +2/3 and -5/3 for the structure functions and the power spectra, respectively.

We estimate those ratios and exponents from in-situ high-rate turbulence measurements collected by three research aircraft during four field experiments in two regimes of the marine atmospheric boundary layer: shallow trade-wind convection and subtropical stratocumulus. The results were derived by fitting power law formulas to the structure functions and power spectra computed separately for the three components of the turbulent wind velocity measured in horizontal flight segments.


## Data

The measurements were performed during four field experiments:
- EUREC4A (Elucidating the role of cloud–circulation coupling in climate) in Jan - Feb 2020 in trade-wind cumulus regime in northwestern Atlantic,
- RICO (Rain in Cumulus Over Ocean) in Nov 2004 - Jan 2005 in trade-wind cumulus regime in northwestern Atlantic,
- VOCALS-REx (Variability of the American Monsoon Systems Ocean-Cloud-Atmosphere-Land Study Regional Experiment) in Oct-Nov 2008 in subtropical stratocumulus regime in southeastern Pacific,
- POST (Physics of the Stratocumulus Top) in Jul-Aug 2008 in subtropical stratocumulus regime in northeastern Pacific.

The measurements were obtained with three research aircraft:
- SAFIRE (the French facility for airborne research) ATR42 during EUREC4A,
- NSF/NCAR (National Science Foundation - National Center for Atmospheric Research) C130 during RICO and VOCALS-REx,
- NPS CIRPAS (Naval Postgraduate School - Center for Interdisciplinary Remotely-Piloted Aircraft Studies) Twin Otter during POST.

The relavant datasets can be downloaded from the public repositories:
- [Lothon, M. & Brilouet, P. (2020). SAFIRE ATR42: Turbulence Data 25 Hz.](https://observations.ipsl.fr/aeris/eurec4a-data/AIRCRAFT/ATR/SAFIRE-TURB/PROCESSED/) doi:10.25326/128,
- [NCAR/NSF C-130 Navigation, State Parameter, and Microphysics HRT Data - 25 Hz. Version 1.0.](https://data.eol.ucar.edu/dataset/87.049) doi:10.5065/D64J0CDM,
- [NCAR/NSF C-130 Navigation, State Parameter, and Microphysics HRT (25 sps) Data.](https://data.eol.ucar.edu/dataset/89.002) doi:10.5065/D69K48JK,
- [NSF/NCAR C130 Radar, Lidar and Radiometer Integrated Dataset.](https://data.eol.ucar.edu/dataset/89.159) doi:10.26023/8KEJ-BQNG-W808,
- [UC Irvine 40-hz Probes - netCDF format.](https://data.eol.ucar.edu/dataset/111.033) doi:10.26023/KP56-KFJS-VC07,
- Table 1 from [Carman, at al. 2012](https://doi.org/10.5194/acp-12-11135-2012) can be also found here in `aux_data`.


## Code

The code was developed in MATLAB R2019b. The functionality in other versions of this environment was not tested.

The script `main.m` performs the analysis described throughout the manuscript. The script `sensitivity.m` performs an additional sensitivity study on the choice of the fitting range described in the appendix.

Two external packages are used:
- [YAML 1.1 parser and emitter for MATLAB](https://www.mathworks.com/matlabcentral/fileexchange/106765-yaml) by Martin Koch
- [boxplotGroup](https://www.mathworks.com/matlabcentral/fileexchange/74437-boxplotgroup) by Adam Danz


