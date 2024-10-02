
%% Settings

% List of planes/experiments
planes = {'ATR-EUREC4A','C130-RICO','C130-VOCALS-REx','TO-POST'};

% Velocity variables (suffixes)
vars_calc = {'UX','VY','W'};


% Fit range: bottom limit (= factor*dr)
sfc_bot_factors = [2 2 2 3]; % list of factors for each plane/experiment
psd_bot_factors = [4 4 4 6];

% Fit range: upper limits (= factor*L)
sfc_upp_factors = [0.6 0.7 0.8 1.0 1.2 1.4]; % list of factors to iterate
psd_upp_factors = sfc_upp_factors*2;

% Fit points
sfc_fit_points = 5;
psd_fit_points = 5;

sfc_min_fit_points = 3;
psd_min_fit_points = 3;

% Power spectrum Welch window size
psd_win_length = 1000; % m
psd_win_overlap = 500; % m


% Output file
outputfile = [myprojectpath,filesep,'sensitivity.mat'];


% Turn off warnings
warning('off','LOGMEAN:EmptyBins')
warning('off','FIT_SFC:InvalidFitRange')
warning('off','FIT_PSD:InvalidFitRange')
warning('off','backtrace')

% Add to path
addpath(genpath(myprojectpath))



%% Computations

Npl = numel(planes);
Nfc = numel(sfc_upp_factors);

MOM_matrix = cell(Npl,Nfc);


for i_p = 1:Npl    
    plane = planes{i_p};
    fprintf('%s\n',plane)
    
    
    % Load datasets
    
    datapath = [mydatapath,filesep,plane];
    
    if strcmp(plane,'ATR-EUREC4A')
        [TURB,MOM] = load_eureca_all(datapath);
    elseif strcmp(plane,'C130-RICO')
        [TURB,MOM] = load_rico_all(datapath);
    elseif strcmp(plane,'C130-VOCALS-REx')
        [TURB,MOM] = load_vocals_all(datapath);
    elseif strcmp(plane,'TO-POST')
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
    
    
    % Iterations
    
    for i_f = 1:Nfc
        testname = sprintf('sfc %.1fL psd %.1fL',sfc_upp_factors(i_f),psd_upp_factors(i_f));
        fprintf('%s',testname)
        
        MOM.plane(:) = string(plane);
        MOM.sfc_bot_factor(:) = sfc_bot_factors(i_p);
        MOM.psd_bot_factor(:) = psd_bot_factors(i_p);
        
        MOM.testname(:) = string(testname);
        MOM.sfc_upp_factor(:) = sfc_upp_factors(i_f);
        MOM.psd_upp_factor(:) = psd_upp_factors(i_f);
        
        
        % Structure functions
        
        fprintf('%s','  SFC ... ')
        
        for i_v = 1:numel(vars_calc)
            var = vars_calc{i_v}; fprintf('%3s',var)
            
            MOM{:,["off","slp","e_off","e_slp","R2","N"]+"_sfc_"+var} = nan;
            
            for i_s = 1:Nseg
                dr = MOM.dr(i_s);
                sfc_fit_range = [dr*sfc_bot_factors(i_p) MOM.int_scale(i_s)*sfc_upp_factors(i_f)];
                
                try
                    [MOM.(['off_sfc_',var])(i_s),MOM.(['slp_sfc_',var])(i_s),es] = ...
                        fit_sfc( detrend(TURB(i_s).(var)), dr, sfc_fit_range, ...
                        'Method','logmean', 'FitPoints',sfc_fit_points );
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
        
        
        % Power spectra
        
        fprintf('%s','  PSD ... ')
        
        for i_v = 1:numel(vars_calc)
            var = vars_calc{i_v}; fprintf('%3s',var)
            
            MOM{:,["off","slp","e_off","e_slp","R2","N"]+"_psd_"+var} = nan;
            
            for i_s = 1:Nseg
                dr = MOM.dr(i_s);
                psd_fit_range = [dr*psd_bot_factors(i_p) MOM.int_scale(i_s)*psd_upp_factors(i_f)];
                
                try
                    [MOM.(['off_psd_',var])(i_s),MOM.(['slp_psd_',var])(i_s),es] = ...
                        fit_psd( detrend(TURB(i_s).(var)), dr, psd_fit_range, ...
                        'Method','logmean', 'FitPoints',psd_fit_points, ...
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
        
        
        MOM.R2_sfc_tot = MOM.R2_sfc_UX .* MOM.R2_sfc_VY .* MOM.R2_sfc_W;
        MOM.R2_psd_tot = MOM.R2_psd_UX .* MOM.R2_psd_VY .* MOM.R2_psd_W;
        MOM.R2_tot_tot = MOM.R2_sfc_tot .* MOM.R2_psd_tot;
  
        
        % Transverse-to-longitudinal ratios
        
        MOM.ar_sfc_VU = MOM.off_sfc_VY./MOM.off_sfc_UX;
        MOM.ar_sfc_WU = MOM.off_sfc_W ./MOM.off_sfc_UX;

        MOM.ar_psd_VU = MOM.off_psd_VY./MOM.off_psd_UX;
        MOM.ar_psd_WU = MOM.off_psd_W ./MOM.off_psd_UX;
        
        MOM.e_ar_sfc_VU = MOM.ar_sfc_VU .* sqrt( (MOM.e_off_sfc_VY./MOM.off_sfc_VY).^2 + (MOM.e_off_sfc_UX./MOM.off_sfc_UX).^2 );
        MOM.e_ar_sfc_WU = MOM.ar_sfc_WU .* sqrt( (MOM.e_off_sfc_W ./MOM.off_sfc_W ).^2 + (MOM.e_off_sfc_UX./MOM.off_sfc_UX).^2 );
       
        MOM.e_ar_psd_VU = MOM.ar_psd_VU .* sqrt( (MOM.e_off_psd_VY./MOM.off_psd_VY).^2 + (MOM.e_off_psd_UX./MOM.off_psd_UX).^2 );
        MOM.e_ar_psd_WU = MOM.ar_psd_WU .* sqrt( (MOM.e_off_psd_W ./MOM.off_psd_W ).^2 + (MOM.e_off_psd_UX./MOM.off_psd_UX).^2 );
        
        
        % Store results in memory
        
        MOM_matrix{i_p,i_f} = MOM;
        
        fprintf('\n')
    end
    
    clear MOM TURB
    fprintf('\n')
    
end



%% Summary

sum_vars = {'ar_sfc_VU','ar_psd_VU','ar_sfc_WU','ar_psd_WU',...
    'slp_sfc_UX','slp_sfc_VY','slp_sfc_W',...
    'slp_psd_UX','slp_psd_VY','slp_psd_W',...
    'R2_sfc_tot','R2_psd_tot','R2_tot_tot',...
    'N_sfc_UX','N_psd_UX'};

grp_vars = {'testname','sfc_upp_factor','psd_upp_factor','plane'};


temp =  cellfun(@(x) x(:,horzcat(grp_vars,sum_vars)),MOM_matrix(:),'UniformOutput',false);
MOM = vertcat(temp{:});

SUM = groupsummary(MOM,grp_vars,{'mean'},sum_vars);



%% Save results

save(outputfile,'planes','vars_calc',...
    'sfc_bot_factors','sfc_upp_factors','sfc_fit_points','sfc_min_fit_points',...
    'psd_bot_factors','psd_upp_factors','psd_fit_points','psd_min_fit_points',...
    'psd_win_length','psd_win_overlap',...
    'MOM_matrix','MOM','SUM')

