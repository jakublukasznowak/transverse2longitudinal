%% Select campaign

% plane = 'ATR-EUREC4A'; ifdirs = true;
% plane = 'C130-RICO';   ifdirs = false;
% plane = 'C130-VOCALS-REx'; ifdirs = false;
% plane = 'TO-POST';     ifdirs = true;



%% Prepare paths

% MYPROJECTPATH is the path where you downloaded the codes
%
%
% MYDATAPATH is the path where you downloaded the datasets:
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
% MYDATAPATH/ATR-EUREC4A/TURBLENCE
%
% Table 3 from Bony et al. 2022: EUREC4A observations from the SAFIRE ATR42 aircraft,
% Earth Syst. Sci. Data, 14, 2021–2064, doi.org/10.5194/essd-14-2021-2022, 2022.
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
% (as tab-delimited text file).


addpath(genpath(myprojectpath))

datapath = [mydatapath,filesep,plane];

plotpath = [myprojectpath,filesep,'figures_new'];
if ~isfolder(plotpath)
    mkdir(plotpath)
end
plotpath = [plotpath,filesep,plane];



%% Load datasets

if strcmp(plane,'ATR-EUREC4A')
    
    sfc_fit_range = [8 40];
    psd_fit_range = [16 80];
    ex_s = ["RF12","R2B"];
    
    [TURB,MOM,levels] = load_eureca_all(datapath);

elseif strcmp(plane,'C130-RICO')
    
    sfc_fit_range = [8 40];
    psd_fit_range = [16 80];
    ex_s = ["RF06","SC01"];
    
    [TURB,MOM,levels] = load_rico_all(datapath);

elseif strcmp(plane,'C130-VOCALS-REx')
    
    sfc_fit_range = [8 40];
    psd_fit_range = [16 80];
    ex_s = ["RF09","C6"];
    
    [TURB,MOM,levels] = load_vocals_all(datapath);
    
elseif strcmp(plane,'TO-POST')
    
    sfc_fit_range = [4 40];
    psd_fit_range = [8 80];
    ex_s = ["RF12","CB01"];
    
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

vars = {'W','UX','VY'};


% Compute and fit

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


% Ratios

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




%% Calculate dissipation

% % Constants
% 
% B_L = 2.0; B_T = 2.6;
% C_L = 0.5; C_T = 0.66;
% 
% 
% % Settings
% 
% sfc_method = "logmean";
% sfc_fit_points = 5;
% 
% psd_method = "logmean";
% psd_fit_points = 5;
% psd_win_length = 1000; % m
% psd_win_overlap = 500; % m
% 
% vars = {'W','UX','VY'};
% B = [B_T B_L B_T];
% C = [C_T C_L C_T];
% 
% 
% % Compute
% 
% disp('Compute dissipation rate ...')
% 
% Nvar = numel(vars);
% Nseg = size(MOM,1);
% E = struct([]);
% 
% for i_v = 1:Nvar
%     var = vars{i_v}; fprintf('%2s',var)
%     
%     for i_s = 1:Nseg
%         fprintf(' %d',i_s)
%         dr = MOM.dr(i_s);
% 
%         [MOM.(['edr_sfc_',var])(i_s),MOM.(['slp_sfc_',var])(i_s),es] = edr_sfc( detrend(TURB(i_s).(var)),...
%             dr,sfc_fit_range,B(i_v),'Method',sfc_method,'FitPoints',sfc_fit_points );
% %         print(f,join([[myprojectpath,filesep,'figures',filesep,'sfc',filesep,plane],...
% %             MOM.flight(i_s),MOM.name(i_s),var],'_'),'-dpng','-r300')
%         
%         [MOM.(['edr_psd_',var])(i_s),MOM.(['slp_psd_',var])(i_s),ep] = edr_psd( detrend(TURB(i_s).(var)),...
%             dr,psd_fit_range,C(i_v),'Method',psd_method,'FitPoints',psd_fit_points,...
%             'WindowLength',floor(psd_win_length/dr),'WindowOverlap',floor(psd_win_overlap/dr) );
% %         print(f,join([[myprojectpath,filesep,'figures',filesep,'psd',filesep,plane],...
% %             MOM.flight(i_s),MOM.name(i_s),var],'_'),'-dpng','-r300')
%         
%         E(1).(['sfc_',var])(i_s) = es;
%         E(1).(['psd_',var])(i_s) = ep;
%     end
%     
%     E.(['sfc_',var]) = struct2table(E.(['sfc_',var]));
%     E.(['psd_',var]) = struct2table(E.(['psd_',var]));
%     
%     MOM.(['slp_psd_',var]) = - MOM.(['slp_psd_',var]);
%     
%     fprintf('\n')
% end



%% Dependent parameters

% % sfc and psd prefactors
% 
% for i_v = 1:Nvar
%     var = vars{i_v};
%     
%     MOM.(['off_sfc_',var]) = B(i_v)*MOM.(['edr_sfc_',var]).^(2/3);
%     MOM.(['off_psd_',var]) = C(i_v)*MOM.(['edr_psd_',var]).^(2/3);
%     
%     MOM.(['e_off_sfc_',var]) = MOM.(['off_sfc_',var]) .* E.(['sfc_',var]).offsetFixed;
%     MOM.(['e_off_psd_',var]) = MOM.(['off_psd_',var]) .* E.(['psd_',var]).offsetFixed;
%     
%     MOM.(['er_off_sfc_',var]) = E.(['sfc_',var]).offsetFixed;
%     MOM.(['er_off_psd_',var]) = E.(['psd_',var]).offsetFixed;
%     
%     MOM.(['e_slp_sfc_',var]) = E.(['sfc_',var]).slopeFree;
%     MOM.(['e_slp_psd_',var]) = E.(['psd_',var]).slopeFree;
% end
% 
% 
% % Anisotropy
% 
% MOM.ar_sfc_VU = MOM.off_sfc_VY./MOM.off_sfc_UX;
% MOM.ar_sfc_WU = MOM.off_sfc_W ./MOM.off_sfc_UX;
% MOM.ar_sfc_WV = MOM.off_sfc_W ./MOM.off_sfc_VY;
% 
% MOM.ar_psd_VU = MOM.off_psd_VY./MOM.off_psd_UX;
% MOM.ar_psd_WU = MOM.off_psd_W ./MOM.off_psd_UX;
% MOM.ar_psd_WV = MOM.off_psd_W ./MOM.off_psd_VY;
% 
% MOM.e_ar_sfc_VU = MOM.ar_sfc_VU .* sqrt( (MOM.e_off_sfc_VY./MOM.off_sfc_VY).^2 + (MOM.e_off_sfc_UX./MOM.off_sfc_UX).^2 );
% MOM.e_ar_sfc_WU = MOM.ar_sfc_WU .* sqrt( (MOM.e_off_sfc_W ./MOM.off_sfc_W ).^2 + (MOM.e_off_sfc_UX./MOM.off_sfc_UX).^2 );
% MOM.e_ar_sfc_WV = MOM.ar_sfc_WV .* sqrt( (MOM.e_off_sfc_W ./MOM.off_sfc_W ).^2 + (MOM.e_off_sfc_VY./MOM.off_sfc_VY).^2 );
% 
% MOM.e_ar_psd_VU = MOM.ar_psd_VU .* sqrt( (MOM.e_off_psd_VY./MOM.off_psd_VY).^2 + (MOM.e_off_psd_UX./MOM.off_psd_UX).^2 );
% MOM.e_ar_psd_WU = MOM.ar_psd_WU .* sqrt( (MOM.e_off_psd_W ./MOM.off_psd_W ).^2 + (MOM.e_off_psd_UX./MOM.off_psd_UX).^2 );
% MOM.e_ar_psd_WV = MOM.ar_psd_WV .* sqrt( (MOM.e_off_psd_W ./MOM.off_psd_W ).^2 + (MOM.e_off_psd_VY./MOM.off_psd_VY).^2 );
% 
% 
% % Dissipation rates after reversal of longi/trans
% 
% MOM.edr_sfc_UY = MOM.edr_sfc_UX * (B_L/B_T).^(3/2);
% MOM.edr_sfc_VX = MOM.edr_sfc_VY * (B_T/B_L).^(3/2);
% MOM.edr_psd_UY = MOM.edr_psd_UX * (C_L/C_T).^(3/2);
% MOM.edr_psd_VX = MOM.edr_psd_VY * (C_T/C_L).^(3/2);



%% Integral length scale

disp('Compute integral length scale ...')

var = 'W'; fprintf('%2s','L')
for i_s = 1:Nseg
    fprintf(' %d',i_s)
    MOM.int_scale(i_s) = int_ls_short(detrend(TURB(i_s).(var)),'Method','e-decay')*MOM.dr(i_s);
end
fprintf('\n')



%% PLOTS

plot_all


%% Summmary of segments

print_table(MOM,{'alt','int_scale','length_km'},1,0)

print_table(MOM,{'ar_sfc_VU','ar_psd_VU','ar_sfc_WU','ar_psd_WU'})
print_table(MOM,{'e_ar_sfc_VU','e_ar_psd_VU','e_ar_sfc_WU','e_ar_psd_WU'},0,2,["median","std"])

print_table(MOM,{'slp_sfc_UX','slp_sfc_VY','slp_sfc_W',...
                 'slp_psd_UX','slp_psd_VY','slp_psd_W'})
print_table(MOM,{'e_slp_sfc_UX','e_slp_sfc_VY','e_slp_sfc_W',...
                 'e_slp_psd_UX','e_slp_psd_VY','e_slp_psd_W'},0,2,["median","std"])
