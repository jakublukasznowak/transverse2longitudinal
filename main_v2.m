

%% Introduction

% The code requires access to the datasets which can be downloaded from
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



%% Settings

% List of planes/experiments
planes = {'ATR-EUREC4A','C130-RICO','C130-VOCALS-REx','TO-POST'};

% Velocity variables (suffixes)
vars_calc = {'UX','VY','W'};


% Max scale
r_max = 400; % m


% Fit range: bottom limit (= factor*dr)
sfc_bot_factors = [2 2 2 3]; % list of factors for each plane/experiment
psd_bot_factors = [4 4 4 6];

% Fit range: upper limit (= factor*L)
sfc_upp_factor = 1.0; % one number for all planes/experiments
psd_upp_factor = 2.0; 

% Fit points
sfc_fit_points = 5;
psd_fit_points = 5;

sfc_min_fit_points = 3;
psd_min_fit_points = 3;

% Power spectrum Welch window size
psd_win_length = 1000; % m
psd_win_overlap = 500; % m


% Scale-by-scale averaging points
sfc_calc_points = 10;
psd_calc_points = 16;

% Scale-by-scale interpolation points
sfc_interp_points = 10;
psd_interp_points = 16;

% Interpolation variables
vars_interp = {'UX','VY','W','ar_VU','ar_WU'};

% Normalized scale range r/L
r_L_range = [0.01 3];


% Output file
outputfile = [myprojectpath,filesep,'results.mat'];


% Turn off warnings
warning('off','LOGMEAN:EmptyBins')
warning('off','FIT_PSD:InvalidFitRange')

% Add to path
addpath(genpath(myprojectpath))



%% Computations

Npl = numel(planes);

MOM_vec = cell(Npl,1); TURB_vec = cell(Npl,1);
SFC_vec = cell(Npl,1); rawSFC_vec = cell(Npl,1); avSFC_vec = cell(Npl,1);
PSD_vec = cell(Npl,1); rawPSD_vec = cell(Npl,1); avPSD_vec = cell(Npl,1);


% Normalized grid for scale-by-scale analysis
r_L_sfc = exp(linspace( log(r_L_range(1)), log(r_L_range(2)), sfc_interp_points))';
r_L_psd = exp(linspace( log(r_L_range(1)), log(r_L_range(2)), psd_interp_points))';


for i_p = 1:Npl    
    plane = planes{i_p};
    fprintf('%s\n',plane)
    
    
    % Load datasets
    
    datapath = [mydatapath,filesep,plane];
    
    if strcmp(plane,'ATR-EUREC4A')
        
        % Load L3 v1.9 'longlegs' data:
        % - segmentation timestamps from yaml files
        % - means/moments: variables ALT, MEAN_TAS, MEAN_THDG, MEAN_WDIR
        % - turb. fluctuations: variables UX_DET, VY_DET, W_DET
        % for flights RF09-RF19
        % and levels: cloud-base, top-subcloud, mid-subcloud, near-surface
        [TURB,MOM] = load_eureca_all(datapath);
        
    elseif strcmp(plane,'C130-RICO')
        
        % Load turbulence data, variables: UXC, VYC, WIC, GGALTC, TASX, THDG, PSXC
        % Apply segmentation algorithm to mark horizontal segments
        % Cut segments out of turbulence data
        % Compute mean values
        % Classify segments into levels and select: cloud-layer, cloud-base, sub-cloud, near-surface
        [TURB,MOM] = load_rico_all(datapath);
        
    elseif strcmp(plane,'C130-VOCALS-REx')
        
        % Load turbulence data, variables: UXC, VYC, WIC, ALTX, TASX, THDG
        % Load segment timestamps from lidar dataset, rename them
        %    and select levels: in-cloud, cloud-base, sub-cloud
        % Cut segments out of turbulence data
        % Compute mean values
        [TURB,MOM] = load_vocals_all(datapath);
        
    elseif strcmp(plane,'TO-POST')
        
        % Load turbulence data, variables: WX, WY, WZ, RADALT, TAS, GTRK
        % Load cloud bas/top height table
        % Apply segmentation algorithm to mark horizontal segments
        % Cut segments out of turbulence data
        % Rotate eastward/northward WX, WY to longitudinal/lateral UX, VY
        % Compute mean values
        % Classify segments into levels and select: cloud-top, cloud-base, sub-cloud, near-surface
        [TURB,MOM] = load_post_all(datapath);
    end
    
    Nseg = size(MOM,1);
    fprintf('Number of segments: %d\n',Nseg)
    
    
    % Integral length scale
    
    disp('Integral length scales ...')
    
    MOM.int_scale = nan(Nseg,1);
    for i_s = 1:Nseg
        MOM.int_scale(i_s) = int_ls_short(detrend(TURB(i_s).W),...
            'Method','e-decay','MaxLag',floor(10e3/MOM.dr(i_s))) * MOM.dr(i_s);
    end

    
    % Structure functions
    
    fprintf('%s','Structure functions ... ')
    
    SFC = cell(Nseg,1); rawSFC = cell(Nseg,1);
    
    for i_v = 1:numel(vars_calc)
        var = vars_calc{i_v}; fprintf('%3s',var)
        
        MOM{:,["off","slp","e_off","e_slp","R2","N"]+"_sfc_"+var} = nan;
        
        for i_s = 1:Nseg
            dr = MOM.dr(i_s);
            sfc_calc_range = [dr r_max];
            sfc_fit_range = [dr*sfc_bot_factors(i_p) MOM.int_scale(i_s)*sfc_upp_factor];
            
            [rawSFC{i_s}.(var),rawSFC{i_s}.r] = calc_raw_sfc( detrend(TURB(i_s).(var)), ...
                dr, sfc_calc_range);
            
            [SFC{i_s}.r, SFC{i_s}.(var)] = logmean( rawSFC{i_s}.r, rawSFC{i_s}.(var), ...
                sfc_calc_points, sfc_calc_range );
            
%             [rv_fit, sfc_fit] = logmean( rawSFC{i_s}.r, rawSFC{i_s}.(var), ...
%                 sfc_fit_points, sfc_fit_range );
%             [MOM.(['off_sfc_',var])(i_s),MOM.(['slp_sfc_',var])(i_s),es] = ...
%                 fit_uni( rv_fit, sfc_fit, 2/3 );
            [MOM.(['off_sfc_',var])(i_s),MOM.(['slp_sfc_',var])(i_s),es] = ...
                fit_sfc( detrend(TURB(i_s).(var)), dr, sfc_fit_range, ...
                'Method','logmean', 'FitPoints',sfc_fit_points );
            
            MOM{i_s,["e_off","e_slp","R2","N"]+"_sfc_"+var} = [es.O es.slp es.R2 es.N];
            if MOM{i_s,"N_sfc_"+var} < sfc_min_fit_points
                MOM{i_s,["off","slp","e_off","e_slp","R2"]+"_sfc_"+var} = nan;
            end
        end
        
    end
   
    SFC = vertcat(SFC{:}); rawSFC = vertcat(rawSFC{:});
    
    fprintf('\n')
    
    
    % Power spectra
    
    fprintf('%s','Power spectra ... ')
    
    PSD = cell(Nseg,1); rawPSD = cell(Nseg,1);
    
    for i_v = 1:numel(vars_calc)
        var = vars_calc{i_v}; fprintf('%3s',var)
        
        MOM{:,["off","slp","e_off","e_slp","R2","N"]+"_psd_"+var} = nan;
        
        for i_s = 1:Nseg
            dr = MOM.dr(i_s);
            psd_calc_range = [2*dr r_max];
            psd_fit_range = [dr*psd_bot_factors(i_p) MOM.int_scale(i_s)*psd_upp_factor];
            
            [rawPSD{i_s}.(var),rawPSD{i_s}.k] = calc_raw_psd( detrend(TURB(i_s).(var)), ...
               dr, 'WindowLength',psd_win_length, 'WindowOverlap',psd_win_overlap);
            rawPSD{i_s}.r = 2*pi./rawPSD{i_s}.k;
            
            [PSD{i_s}.k, PSD{i_s}.(var)] = logmean( rawPSD{i_s}.k, rawPSD{i_s}.(var), ...
                psd_calc_points, 2*pi./psd_calc_range([2 1]) );
            PSD{i_s}.r = 2*pi./PSD{i_s}.k;
            
%             [kv_fit, psd_fit] = logmean( rawPSD{i_s}.k, rawPSD{i_s}.(var), ...
%                 psd_fit_points, 2*pi./psd_fit_range([2 1]) );
%             [MOM.(['off_psd_',var])(i_s),MOM.(['slp_psd_',var])(i_s),es] = ...
%                 fit_uni( kv_fit, psd_fit, -5/3 );
            [MOM.(['off_psd_',var])(i_s),MOM.(['slp_psd_',var])(i_s),es] = ...
                fit_psd( detrend(TURB(i_s).(var)), dr, psd_fit_range, ...
                'Method','logmean', 'FitPoints',psd_fit_points, ...
                'WindowLength',floor(psd_win_length/dr), 'WindowOverlap',floor(psd_win_overlap/dr) );

            MOM{i_s,["e_off","e_slp","R2","N"]+"_psd_"+var} = [es.O es.slp es.R2 es.N];
            if MOM{i_s,"N_psd_"+var} < psd_min_fit_points
                MOM{i_s,["off","slp","e_off","e_slp","R2"]+"_psd_"+var} = nan;
            end
        end
        
        MOM{:,["slp","R2"]+"_psd_"+var} = -MOM{:,["slp","R2"]+"_psd_"+var};
        
    end
    PSD = vertcat(PSD{:}); rawPSD = vertcat(rawPSD{:});
    
    fprintf('\n')
    
    
    % Ratios
    
    disp('Transverse-to-longitudinal ratios ...')
    
    for i_s = 1:Nseg
        SFC(i_s).ar_VU = SFC(i_s).VY./SFC(i_s).UX;
        PSD(i_s).ar_VU = PSD(i_s).VY./PSD(i_s).UX;
        SFC(i_s).ar_WU = SFC(i_s).W./SFC(i_s).UX;
        PSD(i_s).ar_WU = PSD(i_s).W./PSD(i_s).UX;
    end
    
    MOM.ar_sfc_VU = MOM.off_sfc_VY./MOM.off_sfc_UX;
    MOM.ar_sfc_WU = MOM.off_sfc_W ./MOM.off_sfc_UX;

    MOM.ar_psd_VU = MOM.off_psd_VY./MOM.off_psd_UX;
    MOM.ar_psd_WU = MOM.off_psd_W ./MOM.off_psd_UX;
    
    MOM.e_ar_sfc_VU = MOM.ar_sfc_VU .* sqrt( (MOM.e_off_sfc_VY./MOM.off_sfc_VY).^2 + (MOM.e_off_sfc_UX./MOM.off_sfc_UX).^2 );
    MOM.e_ar_sfc_WU = MOM.ar_sfc_WU .* sqrt( (MOM.e_off_sfc_W ./MOM.off_sfc_W ).^2 + (MOM.e_off_sfc_UX./MOM.off_sfc_UX).^2 );
    
    MOM.e_ar_psd_VU = MOM.ar_psd_VU .* sqrt( (MOM.e_off_psd_VY./MOM.off_psd_VY).^2 + (MOM.e_off_psd_UX./MOM.off_psd_UX).^2 );
    MOM.e_ar_psd_WU = MOM.ar_psd_WU .* sqrt( (MOM.e_off_psd_W ./MOM.off_psd_W ).^2 + (MOM.e_off_psd_UX./MOM.off_psd_UX).^2 );
    
    
    % Scale-by-scale level composites

    disp('Scale-by-scale level composites ...')
    
    levels = sortrows(groupsummary(MOM,{'level'},{'mean'},{'alt'}),'mean_alt','descend').level';
    Nlvl = numel(levels);
    
    avSFC = cell(Nlvl,1); avPSD = cell(Nlvl,1);
    
    % Interpolate to metric grid
    
%     sfc_interp_range = [ max(cellfun(@min,{SFC(:).r})) min(cellfun(@max,{SFC(:).r})) ];
%     psd_interp_range = [ max(cellfun(@min,{PSD(:).r})) min(cellfun(@max,{PSD(:).r})) ];
%     
%     r_sfc = exp(linspace( log(sfc_interp_range(1)), log(sfc_interp_range(2)), sfc_interp_points))';
%     r_psd = exp(linspace( log(psd_interp_range(1)), log(psd_interp_range(2)), psd_interp_points))';
%     
%     for i_s = 1:Nseg
%         for i_v = 1:numel(vars_interp)
%             var = vars_interp{i_v};     
%             SFC(i_s).([var,'_i']) = interp1( SFC(i_s).r, SFC(i_s).(var), r_sfc, 'linear' );
%             PSD(i_s).([var,'_i']) = interp1( PSD(i_s).r, PSD(i_s).(var), r_psd, 'linear' );
%         end
%         SFC(i_s).r_i = r_sfc;
%         PSD(i_s).r_i = r_psd;
%     end
    
    % Interpolate to normalized grid
    
    for i_s = 1:Nseg
        L = MOM.int_scale(i_s);
        
        for i_v = 1:numel(vars_interp)
            var = vars_interp{i_v};         
            SFC(i_s).([var,'_iL']) = interp1( SFC(i_s).r/L, SFC(i_s).(var), r_L_sfc, 'linear' );
            PSD(i_s).([var,'_iL']) = interp1( PSD(i_s).r/L, PSD(i_s).(var), r_L_psd, 'linear' );           
        end
        
        SFC(i_s).r_iL = r_L_sfc;
        PSD(i_s).r_iL = r_L_psd;
    end
    
    % Level averaging
    
    ff = fieldnames(SFC);
    vars_avg = ff(endsWith(ff,{'_i','_iL'}));
    
    for i_l = 1:Nlvl
        ind_l = (MOM.level==levels(i_l));
        avSFC{i_l}.level = levels(i_l);
        avPSD{i_l}.level = levels(i_l);
        
        for i_v = 1:numel(vars_avg)
            var = vars_avg{i_v};
            
            avSFC{i_l}.(var) = mean( horzcat(SFC(ind_l).(var)), 2 ,'omitnan');
            avPSD{i_l}.(var) = mean( horzcat(PSD(ind_l).(var)), 2, 'omitnan');
            avSFC{i_l}.([var,'_std']) = std( horzcat(SFC(ind_l).(var)), 0, 2 ,'omitnan');
            avPSD{i_l}.([var,'_std']) = std( horzcat(PSD(ind_l).(var)), 0, 2, 'omitnan');
        end
    end
    
    avSFC = vertcat(avSFC{:}); avPSD = vertcat(avPSD{:});
    
    
    % Store results
    
    MOM_vec{i_p} = MOM;
%     TURB_vec{i_p} = TURB;
    SFC_vec{i_p} = SFC; 
    PSD_vec{i_p} = PSD; 
    rawSFC_vec{i_p} = rawSFC;
    rawPSD_vec{i_p} = rawPSD;
    avSFC_vec{i_p} = avSFC;
    avPSD_vec{i_p} = avPSD;
    
    clear MOM TURB SFC PSD rawSFC rawPSD avSFC avPSD
    
    fprintf('\n')
    
end



%% Tables

disp('Summary tables ...')

for i_p = 1:Npl    
    plane = planes{i_p};
    fprintf('%s\n',plane)
    
    MOM = MOM_vec{i_p};
    
    print_table(MOM,{'level'},{'alt','int_scale','length_km'},true,true,0)

    print_table(MOM,{'level'},{'ar_sfc_VU','ar_psd_VU','ar_sfc_WU','ar_psd_WU'})
    
    print_table(MOM,{'level'},{'slp_sfc_UX','slp_sfc_VY','slp_sfc_W',...
                     'slp_psd_UX','slp_psd_VY','slp_psd_W'})
    
%     CNT = MOM(MOM.N_sfc_UX<sfc_min_fit_points,:);
%     if ~isempty(CNT)
%         print_table(CNT,[],true,0)
%     end
    
    fprintf('\n')
end



%% Save results

save(outputfile,'planes','vars_calc','vars_interp','r_max','r_L_range',...
    'sfc_bot_factors','sfc_upp_factor','sfc_fit_points','sfc_min_fit_points','sfc_calc_points','sfc_interp_points',...
    'psd_bot_factors','psd_upp_factor','psd_fit_points','psd_min_fit_points','psd_calc_points','psd_interp_points',...
    'psd_win_length','psd_win_overlap',...
    'MOM_vec','TURB_vec',...
    'SFC_vec','rawSFC_vec','avSFC_vec',...
    'PSD_vec','rawPSD_vec','avPSD_vec')
