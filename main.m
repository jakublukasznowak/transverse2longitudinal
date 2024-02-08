
%% Introduction

% This code was written by Jakub L. Nowak for the analysis of the
% transverse-to-longitudinal ratio of 2nd order turbulent velocity
% statistics (structure functions and power spectra) derived from aircraft
% measurements. There are 4 field experiments and 3 aircrafts sampling 2 marine
% boundary layer regimes included:
% - EUREC4A - ATR42 aircraft - shallow trade-wind convection
% - RICO - C130 aircraft - shallow trade-wind convection
% - VOCALS-REx - C130 aircraft - subtropical stratocumulus
% - POST - Twin Otter aircraft - subtropical stratocumulus
%
% The analysis requires access to the datasets which can be downloaded from
% the public data repositories.
%
% MYPROJECTPATH is the location where you downloaded the codes
% MYDATAPATH is the path where you downloaded the datasets
%
%
% MYDATAPATH/ATR-EUREC4A/TURBLENCE
%
% Lothon, M. & Brilouet, P. (2020). SAFIRE ATR42: Turbulence Data 25 Hz. Aeris.
% doi.org/10.25326/128
% https://observations.ipsl.fr/aeris/eurec4a-data/AIRCRAFT/ATR/SAFIRE-TURB/PROCESSED/
% 
% In this code 'longlegs' L3 v1.9 is used.
%
%
% MYDATAPATH/C130-RICO/TURBULENCE
% 
% NCAR/NSF C-130 Navigation, State Parameter, and Microphysics HRT Data - 25 Hz. Version 1.0
% doi.org/10.5065/D64J0CDM
% https://data.eol.ucar.edu/dataset/87.049
% 
%
% MYDATAPATH/C130-VOCALS/TURBULENCE
% 
% NCAR/NSF C-130 Navigation, State Parameter, and Microphysics HRT (25 sps) Data
% doi.org/10.5065/D69K48JK
% https://data.eol.ucar.edu/dataset/89.002
%
% MYDATAPATH/C130-VOCALS/LIDAR
% 
% NSF/NCAR C130 Radar, Lidar and Radiometer Integrated Dataset
% doi.org/10.26023/8KEJ-BQNG-W808 
% https://data.eol.ucar.edu/dataset/89.159
%
%
% MYDATAPATH/TO-POST/TURBULENCE
%
% UC Irvine 40-hz Probes - netCDF format
% doi.org/10.26023/KP56-KFJS-VC07 
% https://data.eol.ucar.edu/dataset/111.033
%
% MYDATAPATH/TO-POST/cloud_tops.txt
%
% Table 1 from Carman, J. K., Rossiter, D. L., Khelif, D., Jonsson, H. H.,
% Faloona, I. C., and Chuang, P. Y.: Observational constraints on entrainment
% and the entrainment interface layer in stratocumulus, Atmos. Chem. Phys.,
% 12, 11135–11152, https://doi.org/10.5194/acp-12-11135-2012, 2012. 
% (as tab-delimited text file). Can be copied from the AUX_DATA folder in
% the code repository.
%
%
%
% The code was developed in MATLAB R2019b. The functionality in other
% versions of this environment was not tested.
%
% Two external packages are used
%       (1) YAML 1.1 parser and emitter for MATLAB by Martin Koch
%           from https://www.mathworks.com/matlabcentral/fileexchange/106765-yaml
%       (2) boxplotGroup by Adam Danz
%           from https://www.mathworks.com/matlabcentral/fileexchange/74437-boxplotgroup



%% Select experiment
% Uncomment one of the four options below.

% plane = 'ATR-EUREC4A';
% plane = 'C130-RICO';
% plane = 'C130-VOCALS-REx'; 
% plane = 'TO-POST';     



%% Load datasets

% Prepare paths

addpath(genpath(myprojectpath))
datapath = [mydatapath,filesep,plane];

plotpath = [myprojectpath,filesep,'figures'];
if ~isfolder(plotpath), mkdir(plotpath), end
plotpath = [plotpath,filesep,plane];


% Load experiment-specific datasets
% + specify fitting ranges for structure functions and power spectra
% + select a segment to plot examples of those statistics
% + choose whether to distinguish between segments flown along/across mean wind in plots

if strcmp(plane,'ATR-EUREC4A')
    
    sfc_fit_range = [8 40];
    psd_fit_range = [16 80];
    ex_s = ["RF12","R2B"];
    ifdirs = true;
    
    % Load L3 v1.9 'longlegs' data:
    % - segmentation timestamps from yaml files
    % - means/moments: variables ALT, MEAN_TAS, MEAN_THDG, MEAN_WDIR
    % - turb. fluctuations: variables UX_DET, VY_DET, W_DET
    % for flights RF09-RF19
    % and levels: cloud-base, top-subcloud, mid-subcloud, near-surface
    [TURB,MOM,levels] = load_eureca_all(datapath);

elseif strcmp(plane,'C130-RICO')
    
    sfc_fit_range = [8 40];
    psd_fit_range = [16 80];
    ex_s = ["RF06","SC01"];
    ifdirs = false;
    
    % Load turbulence data, variables: UXC, VYC, WIC, GGALTC, TASX, THDG, PSXC
    % Apply segmentation algorithm to mark horizontal segments
    % Cut segments out of turbulence data
    % Compute mean values
    % Classify segments into levels and select: cloud-layer, cloud-base, sub-cloud, near-surface
    [TURB,MOM,levels] = load_rico_all(datapath);

elseif strcmp(plane,'C130-VOCALS-REx')
    
    sfc_fit_range = [8 40];
    psd_fit_range = [16 80];
    ex_s = ["RF09","C6"];
    ifdirs = false;
    
    % Load turbulence data, variables: UXC, VYC, WIC, ALTX, TASX, THDG
    % Load segment timestamps from lidar dataset, rename them
    %    and select levels: in-cloud, cloud-base, sub-cloud
    % Cut segments out of turbulence data
    % Compute mean values
    [TURB,MOM,levels] = load_vocals_all(datapath);
    
elseif strcmp(plane,'TO-POST')
    
    sfc_fit_range = [4 40];
    psd_fit_range = [8 80];
    ex_s = ["RF12","CB01"];
    ifdirs = true;
    
    % Load turbulence data, variables: WX, WY, WZ, RADALT, TAS, GTRK
    % Load cloud bas/top height table
    % Apply segmentation algorithm to mark horizontal segments
    % Cut segments out of turbulence data
    % Rotate eastward/northward WX, WY to longitudinal/lateral UX, VY
    % Compute mean values
    % Classify segments into levels and select: cloud-top,cloud-base, sub-cloud
    [TURB,MOM,levels] = load_post_all(datapath);

end



%% Compute structure functions and power spectra

% Settings

sfc_method = "logmean";
sfc_fit_points = 5;

psd_method = "logmean";
psd_fit_points = 5;
psd_win_length = 1000; % m
psd_win_overlap = 500; % m

vars = {'UX','VY','W'};


% Compute statistics and fit power laws

Nvar = numel(vars);
Nseg = size(MOM,1);
fprintf('Number of segments: %d\n',Nseg)
disp('Compute structure functions and power spectra ...')

for i_v = 1:Nvar
    var = vars{i_v}; fprintf('%2s',var)
    
    for i_s = 1:Nseg
        fprintf(' %d',i_s)
        dr = MOM.dr(i_s);

        [MOM.(['off_sfc_',var])(i_s),MOM.(['slp_sfc_',var])(i_s),es] = ...
            fit_sfc( detrend(TURB(i_s).(var)), dr, sfc_fit_range, ...
            'Method',sfc_method, 'FitPoints',sfc_fit_points );
        
        MOM.(['e_off_sfc_',var])(i_s) = es.O;
        MOM.(['e_slp_sfc_',var])(i_s) = es.slp;
        
        [MOM.(['off_psd_',var])(i_s),MOM.(['slp_psd_',var])(i_s),ep] = ...
            fit_psd( detrend(TURB(i_s).(var)), dr, psd_fit_range, ...
            'Method',psd_method, 'FitPoints',psd_fit_points, ...
            'WindowLength',floor(psd_win_length/dr), 'WindowOverlap',floor(psd_win_overlap/dr) );
        
        MOM.(['e_off_psd_',var])(i_s) = ep.O;
        MOM.(['e_slp_psd_',var])(i_s) = ep.slp;
    end
    
    MOM.(['slp_psd_',var]) = - MOM.(['slp_psd_',var]);
    
    fprintf('\n')
end


% Calculate ratios

MOM.ar_sfc_VU = MOM.off_sfc_VY./MOM.off_sfc_UX;
MOM.ar_sfc_WU = MOM.off_sfc_W ./MOM.off_sfc_UX;
MOM.ar_sfc_WV = MOM.off_sfc_W ./MOM.off_sfc_VY;

MOM.ar_psd_VU = MOM.off_psd_VY./MOM.off_psd_UX;
MOM.ar_psd_WU = MOM.off_psd_W ./MOM.off_psd_UX;
MOM.ar_psd_WV = MOM.off_psd_W ./MOM.off_psd_VY;

% Error propagation for ratios

MOM.e_ar_sfc_VU = MOM.ar_sfc_VU .* sqrt( (MOM.e_off_sfc_VY./MOM.off_sfc_VY).^2 + (MOM.e_off_sfc_UX./MOM.off_sfc_UX).^2 );
MOM.e_ar_sfc_WU = MOM.ar_sfc_WU .* sqrt( (MOM.e_off_sfc_W ./MOM.off_sfc_W ).^2 + (MOM.e_off_sfc_UX./MOM.off_sfc_UX).^2 );
MOM.e_ar_sfc_WV = MOM.ar_sfc_WV .* sqrt( (MOM.e_off_sfc_W ./MOM.off_sfc_W ).^2 + (MOM.e_off_sfc_VY./MOM.off_sfc_VY).^2 );

MOM.e_ar_psd_VU = MOM.ar_psd_VU .* sqrt( (MOM.e_off_psd_VY./MOM.off_psd_VY).^2 + (MOM.e_off_psd_UX./MOM.off_psd_UX).^2 );
MOM.e_ar_psd_WU = MOM.ar_psd_WU .* sqrt( (MOM.e_off_psd_W ./MOM.off_psd_W ).^2 + (MOM.e_off_psd_UX./MOM.off_psd_UX).^2 );
MOM.e_ar_psd_WV = MOM.ar_psd_WV .* sqrt( (MOM.e_off_psd_W ./MOM.off_psd_W ).^2 + (MOM.e_off_psd_VY./MOM.off_psd_VY).^2 );



%% Compute integral length scale

disp('Compute integral length scale ...')

var = 'W'; fprintf('%2s','L')
for i_s = 1:Nseg
    fprintf(' %d',i_s)
    MOM.int_scale(i_s) = int_ls_short(detrend(TURB(i_s).(var)),'Method','e-decay')*MOM.dr(i_s);
end
fprintf('\n')



%% PLOTS

plot_all



%% Print out summary of results for each level

% Segment info: number, average altitude / length / integral scale
print_table(MOM,{'alt','int_scale','length_km'},1,0)

% Transverse-to-longitudinal ratios and their uncertainties
print_table(MOM,{'ar_sfc_VU','ar_psd_VU','ar_sfc_WU','ar_psd_WU'})
print_table(MOM,{'e_ar_sfc_VU','e_ar_psd_VU','e_ar_sfc_WU','e_ar_psd_WU'},0,2,["median","std"])

% Scaling exponents and their uncertainties
print_table(MOM,{'slp_sfc_UX','slp_sfc_VY','slp_sfc_W',...
                 'slp_psd_UX','slp_psd_VY','slp_psd_W'})
print_table(MOM,{'e_slp_sfc_UX','e_slp_sfc_VY','e_slp_sfc_W',...
                 'e_slp_psd_UX','e_slp_psd_VY','e_slp_psd_W'},0,2,["median","std"])
