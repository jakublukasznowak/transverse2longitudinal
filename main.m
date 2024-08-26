
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

% Fit method
sfc_method = "logmean";
psd_method = "logmean";

% Power spectrum Welch window size
psd_win_length = 1000; % m
psd_win_overlap = 500; % m


% Velocity variables (suffixes)
vars = {'UX','VY','W'};


% Turn off fit range warnings
warning('off','LOGMEAN:EmptyBins')
warning('off','FIT_PSD:InvalidFitRange')
warning('off','backtrace')


% Example segments to plot sfc/psd
examples = {["RF12","R2B"],["RF06","SC01"],["RF09","C6"],["RF12","CB01"]}; % [flight, name] for each plane/experiment


addpath(genpath(myprojectpath))
plotpath = [myprojectpath,filesep,'figures'];
if ~isfolder(plotpath), mkdir(plotpath), end



%% Computations

Npl = numel(planes);
Nvar = numel(vars);

MOM_vec = cell(Npl,1);
TURB_vec = cell(Npl,1);


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
    
    disp('Integral length scale ...')
    
    MOM.int_scale = nan(Nseg,1);
    for i_s = 1:Nseg
        MOM.int_scale(i_s) = int_ls_short(detrend(TURB(i_s).W),...
            'Method','e-decay','MaxLag',floor(10e3/MOM.dr(i_s))) * MOM.dr(i_s);
    end

    
    % Structure functions
    
    fprintf('%s','Structure functions ... ')
    
    MOM.sfc_fit_range = [MOM.dr*sfc_bot_factors(i_p) MOM.int_scale*sfc_upp_factor];
    
    for i_v = 1:Nvar
        var = vars{i_v}; fprintf('%3s',var)
        
        MOM{:,["off","slp","e_off","e_slp","R2","N"]+"_sfc_"+var} = nan;

        for i_s = 1:Nseg
            dr = MOM.dr(i_s);
            
            try
                [MOM.(['off_sfc_',var])(i_s),MOM.(['slp_sfc_',var])(i_s),es] = ...
                    fit_sfc( detrend(TURB(i_s).(var)), dr, MOM.sfc_fit_range(i_s,:), ...
                    'Method',sfc_method, 'FitPoints',sfc_fit_points );
                MOM{i_s,["e_off","e_slp","R2","N"]+"_sfc_"+var} = [es.O es.slp es.R2 es.N];
            
            catch ME
                if strcmp(ME.identifier,'FIT_SFC:TooFewFitPoints')
                    fprintf('\n')
                    warning(ME.identifier,['In seg %d: ',getReport(ME,'basic')],i_s)
                    MOM{i_s,"N_sfc_"+var} = 0;
                else
                    throw(ME)
                end
                
            end
            
            if MOM{i_s,"N_sfc_"+var} < sfc_min_fit_points
                MOM{i_s,["off","slp","e_off","e_slp","R2"]+"_sfc_"+var} = nan;
            end
        end
        
    end
    fprintf('\n')
    
    
    % Power spectra
    
    fprintf('%s','Power spectra ... ')
    
    MOM.psd_fit_range = [MOM.dr*psd_bot_factors(i_p) MOM.int_scale*psd_upp_factor];
        
    for i_v = 1:Nvar
        var = vars{i_v}; fprintf('%3s',var)
        
        MOM{:,["off","slp","e_off","e_slp","R2","N"]+"_psd_"+var} = nan;

        for i_s = 1:Nseg
            dr = MOM.dr(i_s);
            
            try
                [MOM.(['off_psd_',var])(i_s),MOM.(['slp_psd_',var])(i_s),es] = ...
                    fit_psd( detrend(TURB(i_s).(var)), dr, MOM.psd_fit_range(i_s,:), ...
                    'Method',psd_method, 'FitPoints',psd_fit_points, ...
                    'WindowLength',floor(psd_win_length/dr), 'WindowOverlap',floor(psd_win_overlap/dr) );
                MOM{i_s,["e_off","e_slp","R2","N"]+"_psd_"+var} = [es.O es.slp es.R2 es.N];
            
            catch ME
                if strcmp(ME.identifier,'FIT_PSD:TooFewFitPoints')
                    fprintf('\n')
                    warning(ME.identifier,['In seg %d: ',getReport(ME,'basic')],i_s)
                    MOM{i_s,"N_psd_"+var} = 0;
                else
                    throw(ME)
                end
            end

            if MOM{i_s,"N_psd_"+var} < psd_min_fit_points
                MOM{i_s,["off","slp","e_off","e_slp","R2"]+"_psd_"+var} = nan;
            end
        end

        MOM{:,["slp","R2"]+"_psd_"+var} = -MOM{:,["slp","R2"]+"_psd_"+var};
 
    end
    fprintf('\n')
    
    
    % Transverse-to-longitudinal ratios
    
    disp('Transverse-to-longitudinal ratios ...')

    MOM.ar_sfc_VU = MOM.off_sfc_VY./MOM.off_sfc_UX;
    MOM.ar_sfc_WU = MOM.off_sfc_W ./MOM.off_sfc_UX;

    MOM.ar_psd_VU = MOM.off_psd_VY./MOM.off_psd_UX;
    MOM.ar_psd_WU = MOM.off_psd_W ./MOM.off_psd_UX;
    
    MOM.e_ar_sfc_VU = MOM.ar_sfc_VU .* sqrt( (MOM.e_off_sfc_VY./MOM.off_sfc_VY).^2 + (MOM.e_off_sfc_UX./MOM.off_sfc_UX).^2 );
    MOM.e_ar_sfc_WU = MOM.ar_sfc_WU .* sqrt( (MOM.e_off_sfc_W ./MOM.off_sfc_W ).^2 + (MOM.e_off_sfc_UX./MOM.off_sfc_UX).^2 );
    
    MOM.e_ar_psd_VU = MOM.ar_psd_VU .* sqrt( (MOM.e_off_psd_VY./MOM.off_psd_VY).^2 + (MOM.e_off_psd_UX./MOM.off_psd_UX).^2 );
    MOM.e_ar_psd_WU = MOM.ar_psd_WU .* sqrt( (MOM.e_off_psd_W ./MOM.off_psd_W ).^2 + (MOM.e_off_psd_UX./MOM.off_psd_UX).^2 );
    
    
    % Store results in memory
    
    MOM_vec{i_p} = MOM;
    TURB_vec{i_p} = TURB;
    
    clear MOM TURB
    
    fprintf('\n')
end


    
%% Tables

disp('Summary tables ...')

for i_p = 1:Npl    
    plane = planes{i_p};
    fprintf('%s\n',plane)
    
    MOM = MOM_vec{i_p};
    
    print_table(MOM,{'alt','int_scale','length_km'},true,0)

    print_table(MOM,{'ar_sfc_VU','ar_psd_VU','ar_sfc_WU','ar_psd_WU'})
    
    print_table(MOM,{'slp_sfc_UX','slp_sfc_VY','slp_sfc_W',...
                     'slp_psd_UX','slp_psd_VY','slp_psd_W'})
    
%     CNT = MOM(MOM.N_sfc_UX<sfc_min_fit_points,:);
%     if ~isempty(CNT)
%         print_table(CNT,[],true,0)
%     end
    
    fprintf('\n')
end

%%

% vars = {'ar_sfc_VU','ar_psd_VU','ar_sfc_WU','ar_psd_WU',...
%         'slp_sfc_UX','slp_sfc_VY','slp_sfc_W',...
%         'slp_psd_UX','slp_psd_VY','slp_psd_W'};
% 
% for i_p = 1:Npl
%     plane = planes{i_p};
%     
%     if ismember(plane,{'ATR-EUREC4A','TO-POST'})
%         fprintf('%s\n',plane)
%         
%         MOM = MOM_vec{i_p};
%         
%         DIRSUM = groupsummary(MOM,{'level','dir2'},{'std'},vars);
%         for i_v = 1:numel(vars)
%             var = vars{i_v};
%             DIRSUM{:,['std_',var]} = DIRSUM{:,['std_',var]}./sqrt(DIRSUM.GroupCount);
%         end
%         DIRSUM
%         
%         levels = sortrows(groupsummary(MOM,{'level'},{'mean'},{'alt'}),'mean_alt','descend').level';
%         for i_l = 1:numel(levels)
%             DIRSUM.level_id(DIRSUM.level==levels(i_l)) = i_l;
%         end
%         
%         [fig,ax] = plot_xy_uni(DIRSUM,{'std_ar_psd_VU'},{'std_ar_sfc_VU'},'level_id','','dir2',true,{'cross1'},10);
%         legend(levels,'Location','northwest','Interpreter','latex')
%         xlabel('$P_v/P_u$','Interpreter','latex')
%         ylabel('$D_v/D_u$','Interpreter','latex')
%         title(plane)
% %         print(fig,[plotpath_res,filesep,plane,'_ar_uv'],'-dpng','-r300')
%     end
% end



%% Plot examples

disp('Plot example segments ...')

plotpath_ex = [plotpath,filesep,'examples'];
if ~isfolder(plotpath_ex), mkdir(plotpath_ex), end


for i_p = 1:Npl    
    plane = planes{i_p};
    fprintf('%s\n',plane)
    
    if ~isfolder([plotpath_ex,filesep,plane]), mkdir([plotpath_ex,filesep,plane]), end
    
    MOM = MOM_vec{i_p};
    TURB = TURB_vec{i_p};
    
    if ~isempty(examples)
        ind_s = find( MOM.flight==examples{i_p}(1) & MOM.name==examples{i_p}(2) );
        plotpath_plane = [plotpath_ex,filesep,plane,'_'];
        ylim_sfc = [1e-2 0.3];
        ylim_psd = [1e-4 1];
    else % if there is no example list, plot all segments
        ind_s = 1:Nseg;
        plotpath_plane = [plotpath_ex,filesep,plane,filesep];
        ylim_sfc = [-inf inf];
        ylim_psd = [-inf inf];
    end
    
    
    for ii_s = 1:numel(ind_s)
        i_s = ind_s(ii_s);
        dr = MOM.dr(i_s);
    
        
        % Structure functions

        for i_v = 1:Nvar
            var = vars{i_v};

            fit_sfc( detrend(TURB(i_s).(var)), dr, MOM.sfc_fit_range(i_s,:), ...
                'Method',sfc_method, 'FitPoints',sfc_fit_points, ...
                'Plot',true, 'PlotXLim',[dr 400], 'PlotYLim',ylim_sfc);

            ylabel(['$D_',lower(var(1)),'\,[\mathrm{m^2\,s^{-2}}]$'],'Interpreter','latex')
            title(join([plane,MOM.flight(i_s),MOM.name(i_s),round(MOM.alt(i_s)),'m']))
            print(gcf,join([[plotpath_plane,'sfc'],MOM.level(i_s),var,string(i_s)],'_'),'-dpng','-r300')
        end


        % Power spectra

        for i_v = 1:Nvar
            var = vars{i_v};

            fit_psd( detrend(TURB(i_s).(var)), dr, MOM.psd_fit_range(i_s,:), ...
                'Method',psd_method, 'FitPoints',psd_fit_points, ...
                'WindowLength',floor(psd_win_length/dr), 'WindowOverlap',floor(psd_win_overlap/dr), ...
                'Plot',true, 'PlotXLim',[2*dr 400], 'PlotYLim',ylim_psd);       

            ylabel(['$P_',lower(var(1)),'\,[\mathrm{m^2\,s^{-2}\,rad^{-1}}]$'],'Interpreter','latex')
            title(join([plane,MOM.flight(i_s),MOM.name(i_s),round(MOM.alt(i_s)),'m']))
            print(gcf,join([[plotpath_plane,'psd'],MOM.level(i_s),var,string(i_s)],'_'),'-dpng','-r300')
        end
        
        close all
    end
    
    close all
end



%% Plot results

disp('Plot results ...')

plotpath_res = [plotpath,filesep,'main'];
if ~isfolder(plotpath_res), mkdir(plotpath_res), end

% Ax limits
ratio_lim = [0 1.6];
p_lim = [0.4 2.5];
s_lim = [0 1.3];


for i_p = 1:Npl    
    plane = planes{i_p};
    fprintf('%s\n',plane)
    
    MOM = MOM_vec{i_p};
    
    if ismember(plane,{'ATR-EUREC4A','TO-POST'})
        dirvar = 'dir2';
    else
        dirvar = '';
    end
    
    levels = sortrows(groupsummary(MOM,{'level'},{'mean'},{'alt'}),'mean_alt','descend').level';
    for i_l = 1:numel(levels)
        MOM.level_id(MOM.level==levels(i_l)) = i_l;
    end
    
    
    % (Pv/Pu,Dv/Du)
    
    [fig,ax] = plot_xy_uni(MOM,{'ar_psd_VU'},{'ar_sfc_VU'},'level_id','',dirvar,true,...
        {'ver3/4','ver4/3','hor3/4','hor4/3','cross1'},[],'XLim',ratio_lim,'YLim',ratio_lim);
    plot(ax,4/3,4/3,'d','Color',"#77AC30",'MarkerFaceColor',"#77AC30",'MarkerSize',12)
    legend(horzcat(levels,{'HIT'}),'Location','northwest','Interpreter','latex')
    xlabel('$P_v/P_u$','Interpreter','latex')
    ylabel('$D_v/D_u$','Interpreter','latex')
    title(plane)
    print(fig,[plotpath_res,filesep,plane,'_ar_uv'],'-dpng','-r300')
    
    
    % (Pw/Pu,Dw/Du)
    
    [fig,ax] = plot_xy_uni(MOM,{'ar_psd_WU'},{'ar_sfc_WU'},'level_id','',dirvar,true,...
        {'ver3/4','ver4/3','hor3/4','hor4/3','cross1'},[],'XLim',ratio_lim,'YLim',ratio_lim);
    plot(ax,4/3,4/3,'d','Color',"#77AC30",'MarkerFaceColor',"#77AC30",'MarkerSize',12)
    legend(horzcat(levels,{'HIT'}),'Location','northwest','Interpreter','latex')
    xlabel('$P_w/P_u$','Interpreter','latex')
    ylabel('$D_w/D_u$','Interpreter','latex')
    title(plane)
    print(fig,[plotpath_res,filesep,plane,'_ar_uw'],'-dpng','-r300')
    
    
    % (p,s)

    [fig,ax] = plot_xy_uni(MOM,{'slp_psd_UX','slp_psd_VY','slp_psd_W'},...
        {'slp_sfc_UX','slp_sfc_VY','slp_sfc_W'},'yN','level_id',dirvar,false,...
        {'hor2/3','ver5/3'},[],'XLim',p_lim,'YLim',s_lim);
    plot(ax,5/3,2/3,'d','Color',"#77AC30",'MarkerFaceColor',"#77AC30",'MarkerSize',12)
    legend(horzcat({'u','v','w'},levels,{'K41'}),'Location','northwest','Interpreter','latex')
    xlabel('$p$','Interpreter','latex')
    ylabel('$s$','Interpreter','latex')
    title(plane)
    print(fig,[plotpath_res,filesep,plane,'_slp'],'-dpng','-r300')
    
end



%% Plot uncertainties in box-whisker

% Ax limits
e_ratio_lim = [0 0.3];
e_s_lim = [0 0.1];
e_p_lim = [0 0.3];


for i_p = 1:Npl    
    plane = planes{i_p};
    fprintf('%s\n',plane)
    
    MOM = MOM_vec{i_p};

    levels = sortrows(groupsummary(MOM,{'level'},{'mean'},{'alt'}),'mean_alt','descend').level';
    for i_l = 1:numel(levels)
        MOM.level_id(MOM.level==levels(i_l)) = i_l;
    end

    
    h = plot_whisker(MOM,{'e_ar_sfc_VU','e_ar_psd_VU','e_ar_sfc_WU','e_ar_psd_WU'},...
        levels,0,'PrimaryLabels',{'$D_v/D_u$','$P_v/P_u$','$D_w/D_u$','$P_w/P_u$'},'DataLim',e_ratio_lim);
    hold on
    h.axis.YLim = e_ratio_lim;
    ylabel('Uncertainty','Interpreter','latex')
    title(plane)
    print(h.figure,[plotpath_res,filesep,plane,'_e_wsk_ar'],'-dpng','-r300')
    
    h = plot_whisker(MOM,{'e_slp_sfc_UX','e_slp_sfc_VY','e_slp_sfc_W'},...
        levels,0,'PrimaryLabels',{'$s_u$','$s_v$','$s_w$'},'DataLim',e_s_lim);
    hold on
    h.axis.YLim = e_s_lim;
    ylabel('Uncertainty','Interpreter','latex')
    title(plane)
    print(h.figure,[plotpath_res,filesep,plane,'_e_wsk_slp_sfc'],'-dpng','-r300')
    
    h = plot_whisker(MOM,{'e_slp_psd_UX','e_slp_psd_VY','e_slp_psd_W'},...
        levels,0,'PrimaryLabels',{'$p_u$','$p_v$','$p_w$'},'DataLim',e_p_lim);
    hold on
    h.axis.YLim = e_p_lim;
    ylabel('Uncertainty','Interpreter','latex')
    title(plane)
    print(h.figure,[plotpath_res,filesep,plane,'_e_wsk_slp_psd'],'-dpng','-r300')

end