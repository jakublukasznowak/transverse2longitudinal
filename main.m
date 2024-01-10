
% Select plane: ATR, C130 or TO

plane = 'ATR';


% Prepare paths

projectpath = pwd;

addpath(genpath(projectpath))

datapath = [mydatapath,filesep,plane];

plotpath = [projectpath,filesep,'figures'];
if ~isfolder(plotpath)
    mkdir(plotpath)
end


% MYDATAPATH is the path where you downloaded the datasets:
%
% MYDATAPATH/ATR/TURBLENCE
%
% Lothon, M. & Brilouet, P. (2020). SAFIRE ATR42: Turbulence Data 25 Hz. Aeris.
% doi.org/10.25326/128
% https://observations.ipsl.fr/aeris/eurec4a-data/AIRCRAFT/ATR/SAFIRE-TURB/PROCESSED/
% 
% In this code 'longlegs' L3 v1.9 is used.
%
%
% MYDATAPATH/C130/TURBULENCE
% 
% NCAR/NSF C-130 Navigation, State Parameter, and Microphysics HRT (25 sps) Data
% doi.org/10.5065/D69K48JK
% https://data.eol.ucar.edu/dataset/89.002
%
% MYDATAPATH/C130/LIDAR
% 
% NSF/NCAR C130 Radar, Lidar and Radiometer Integrated Dataset
% doi.org/10.26023/8KEJ-BQNG-W808 
% https://data.eol.ucar.edu/dataset/89.159
%
%
% MYDATAPATH/TO/TURBULENCE
%
% UC Irvine 40-hz Probes - netCDF format
% doi.org/10.26023/KP56-KFJS-VC07 
% https://data.eol.ucar.edu/dataset/111.033
%
% MYDATAPATH/TO/cloud_tops.txt
%
% Table 1 from Carman, J. K., Rossiter, D. L., Khelif, D., Jonsson, H. H.,
% Faloona, I. C., and Chuang, P. Y.: Observational constraints on entrainment
% and the entrainment interface layer in stratocumulus, Atmos. Chem. Phys.,
% 12, 11135–11152, https://doi.org/10.5194/acp-12-11135-2012, 2012. 
% (as tab-delimited text file).



%% Load datasets


if strcmp(plane,'ATR')
    
    fit_range = [16 80];
    
    % List of levels
    levels  = {'cloud-base','top-subcloud','mid-subcloud','near-surface'};

    % List of flights
    flight_ids = num2cell(num2str((9:19)','RF%02d'),2); % RF09 - RF19

    % List of variables from turbulent moments dataset
    mom_vars = {'alt','MEAN_WSPD','MEAN_WDIR','MEAN_TAS','MEAN_THDG'};

    % List of variables from turbulent fluctuations dataset
    turb_vars = {'W','W_DET';
                 'UX','UX_DET';
                 'VY','VY_DET'};


    % Read data files

    SEG = load_atr_seg(datapath,'v1.9','longlegs');                          % Flight segmentation
    [MOM,mom_info] = load_atr_mom(datapath,'L3','v1.9','longlegs',mom_vars); % Mean values and moments
    
    SEG = SEG(ismember(SEG.flight,flight_ids) & ismember(SEG.level,levels),:); % Select segments
    MOM = join(SEG,MOM,'Keys',{'start','end'});

    TURB = load_atr_turb(MOM,datapath,'L3','v1.9',turb_vars(:,2),turb_vars(:,1)); % Load signals for the selected segments only

    
elseif strcmp(plane,'C130')
    
    fit_range = [16 80];
    
    % List of levels
    levels  = {'in-cloud','cloud-base','sub-cloud'};
    
    % List of variables from turbulence dataset
    turb_vars = {'time','Time';
                 'ALT','ALTX';
                 'TAS','TASX';
                 'THDG','THDG';
                 'U','UIC';
                 'V','VIC';
                 'W','WIC';
                 'UX','UXC';
                 'VY','VYC'};
    
             
    % Read data files
    
    [DATA,turb_info] = read_turb([datapath,filesep,'TURBULENCE'],turb_vars(:,2),turb_vars(:,1));
    [SEG,seg_info] = load_c130_seg(datapath); % Segmentation from auxiliary dataset
    
    
    % Process
    
    SEG = SEG(ismember(SEG.level,levels),:); % Select segments
    TURB = calc_turb(SEG,DATA);              % Cut signals to the segments
    MOM = calc_mom(TURB);                    % Calculate mean segment values

    
elseif strcmp(plane,'TO')
    
    fit_range = [8 80];
    
    % List of levels
    levels  = {'cloud-top','cloud-base','sub-cloud','near-surface'};
    
    % List of variables from turbulence dataset
    turb_vars = {'time','Time';
                 'ALT','RADALT';
                 'TAS','TAS';
                 'THDG','GTRK'; 
                 'U','WX';
                 'V','WY';
                 'W','WZ'};
   
    
    % Read data files
    
    [DATA,turb_info] = read_turb([datapath,filesep,'TURBULENCE'],turb_vars(:,2),turb_vars(:,1));
    CT = readtable([datapath,filesep,'cloud_tops.txt']);
    
    
    % Process
    
    SEG = calc_seg(DATA);       % Algorithmic detection of horizontal segments
    TURB = calc_turb(SEG,DATA); % Cut signals to the segments
    TURB = uv2uxvy(TURB);       % Wind rotation from U,V to UX,VY
    MOM = calc_mom(TURB);       % Calculate mean segment values
    
    % Append cloud heights to MOM
    for i_s = 1:size(MOM,1)
        ind_f = find(strcmp(CT.flight,MOM.flight(i_s)));
        MOM.cloud_base(i_s) = CT.cloud_base(ind_f);
        MOM.cloud_top(i_s)  = CT.cloud_top(ind_f);
        MOM.cloud_top_std(i_s) = CT.cloud_top_std(ind_f);
    end
    MOM.cloud_mid = mean([MOM.cloud_base MOM.cloud_top],2);
    
    % Classify segments according to height 
    MOM.level = repmat("",size(MOM,1),1);
    MOM.level(MOM.alt<60) = "near-surface";
    MOM.level(MOM.alt>=60 & MOM.alt<MOM.cloud_base) = "sub-cloud";
    MOM.level(MOM.alt>=MOM.cloud_base & MOM.alt<MOM.cloud_mid) = "cloud-base";
    MOM.level(MOM.alt>=MOM.cloud_mid & MOM.alt<MOM.cloud_top+MOM.cloud_top_std) = "cloud-top";
    MOM.level(MOM.alt>=MOM.cloud_top+MOM.cloud_top_std) = "free-troposphere";
    MOM = movevars(MOM,{'flight','level','alt','length',...
        'cloud_base','cloud_top','cloud_top_std'},'Before',1);
    
    % Select segments
    ind_s = ismember(MOM.level,levels) & MOM.length>=10e3;
    TURB = TURB(ind_s,:);
    MOM  = MOM(ind_s,:);

end

MOM.dr = MOM.MEAN_TAS./[TURB.fsamp]';

clear SEG DATA


% Plot overview of the segments

plot_seg_overview(MOM,levels);
title(plane)


%%

% Constants

B_L = 2.0; B_T = 2.6;
C_L = 0.5; C_T = 0.66;

% Settings

sfc_method = "logmean";
sfc_fit_points = 6;

psd_method = "logmean";
psd_fit_points = 6;
psd_win_length = 1000; % m
psd_win_overlap = 500; % m

vars = {'W','UX','VY'};
B = [B_T B_L B_T];
C = [C_T C_L C_T];

Nvar = numel(vars);


% Example

ex_s = 32;

i_s = ex_s;
dr = MOM.dr(i_s);

for i_v = 1:Nvar
    var = vars{i_v};
    
    edr_sfc( detrend(TURB(i_s).(var)), dr,fit_range,B(i_v),'Method',sfc_method,...
        'FitPoints',sfc_fit_points,'Plot',true,'PlotRange',[dr 1000] );
    
    title(join([plane,MOM.flight(i_s),round(MOM.alt(i_s)),'m',var]))
    print(gcf,join([[plotpath,filesep,'ex'],plane,'sfc',var,string(i_s)],'_'),'-dpng','-r300')
end

for i_v = 1:Nvar
    var = vars{i_v};
    
    edr_psd( detrend(TURB(i_s).(var)), dr,fit_range,B(i_v),'Method',psd_method,...
        'FitPoints',psd_fit_points,'Plot',true,'PlotRange',[2*dr 1000],...
        'WindowLength',floor(psd_win_length/dr),'WindowOverlap',floor(psd_win_overlap/dr) );
    
    title(join([plane,MOM.flight(i_s),round(MOM.alt(i_s)),'m',var]))
    print(gcf,join([[plotpath,filesep,'ex'],plane,'psd',var,string(i_s)],'_'),'-dpng','-r300')
end


%% Calculate dissipation

Nseg = size(MOM,1);

for i_s = 1:Nseg
    dr = MOM.dr(i_s);
    
    for i_v = 1:Nvar
        var = vars{i_v};

        [MOM.(['edr_sfc_',var])(i_s),MOM.(['slp_sfc_',var])(i_s)] = edr_sfc( detrend(TURB(i_s).(var)),...
            dr,fit_range,B(i_v),'Method',sfc_method,'FitPoints',sfc_fit_points );
        
        [MOM.(['edr_psd_',var])(i_s),MOM.(['slp_psd_',var])(i_s)] = edr_psd( detrend(TURB(i_s).(var)),...
            dr,fit_range,C(i_v),'Method',psd_method,'FitPoints',psd_fit_points,...
            'WindowLength',floor(psd_win_length/dr),'WindowOverlap',floor(psd_win_overlap/dr) );
    end
end

