
%% Settings

% List of planes/experiments
planes = {'ATR-EUREC4A','C130-RICO','C130-VOCALS-REx','TO-POST'};

% Velocity variables (suffixes)
vars_calc = {'UX','VY','W'};

% Max scale
r_max = 400;

% Logmean averaging points
sfc_calc_points = 10;
psd_calc_points = 16;

% Power spectrum Welch window size
psd_win_length  = 1000; % m
psd_win_overlap = 500; % m


addpath(genpath(myprojectpath))



%% Computations

Npl = numel(planes);

MOM_vec = cell(Npl,1); TURB_vec = cell(Npl,1);
SFC_vec = cell(Npl,1); rawSFC_vec = cell(Npl,1);
PSD_vec = cell(Npl,1); rawPSD_vec = cell(Npl,1);


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
        
        for i_s = 1:Nseg    
            sfc_calc_range = [MOM.dr(i_s) r_max];
            
            [rawSFC{i_s}.(var),rawSFC{i_s}.r] = calc_raw_sfc( detrend(TURB(i_s).(var)), ...
                MOM.dr(i_s), [MOM.dr(i_s) r_max]);
            
            [SFC{i_s}.r, SFC{i_s}.(var)] = logmean( rawSFC{i_s}.r, rawSFC{i_s}.(var), ...
                sfc_calc_points, sfc_calc_range );
        end
        
    end
   
    SFC = vertcat(SFC{:}); rawSFC = vertcat(rawSFC{:});
    
    fprintf('\n')
    
    
    % Power spectra
    
    fprintf('%s','Power spectra ... ')
    
    PSD = cell(Nseg,1); rawPSD = cell(Nseg,1);
    
    for i_v = 1:numel(vars_calc)
        var = vars_calc{i_v}; fprintf('%3s',var)
        
        for i_s = 1:Nseg
            psd_calc_range = [2*MOM.dr(i_s) r_max];
            
            [rawPSD{i_s}.(var),rawPSD{i_s}.k] = calc_raw_psd( detrend(TURB(i_s).(var)), ...
                MOM.dr(i_s), 'WindowLength',psd_win_length, 'WindowOverlap',psd_win_overlap);
            rawPSD{i_s}.r = 2*pi./rawPSD{i_s}.k;
            
            [PSD{i_s}.k, PSD{i_s}.(var)] = logmean( rawPSD{i_s}.k, rawPSD{i_s}.(var), ...
                psd_calc_points, 2*pi./psd_calc_range([2 1]) );
            PSD{i_s}.r = 2*pi./PSD{i_s}.k;
        end
        
    end
    PSD = vertcat(PSD{:}); rawPSD = vertcat(rawPSD{:});
    
    fprintf('\n')
    
    
    % Ratios
    
    for i_s = 1:Nseg
        SFC(i_s).ar_VU = SFC(i_s).VY./SFC(i_s).UX;
        PSD(i_s).ar_VU = PSD(i_s).VY./PSD(i_s).UX;
        SFC(i_s).ar_WU = SFC(i_s).W./SFC(i_s).UX;
        PSD(i_s).ar_WU = PSD(i_s).W./PSD(i_s).UX;
    end
    
    
    % Save results
    
    MOM_vec{i_p} = MOM;
%     TURB_vec{i_p} = TURB;
    SFC_vec{i_p} = SFC; 
    PSD_vec{i_p} = PSD; 
    rawSFC_vec{i_p} = rawSFC;
    rawPSD_vec{i_p} = rawPSD;
    
    clear MOM TURB SFC PSD rawSFC rawPSD
    
    fprintf('\n')
    
end



%% Save/load

save([myprojectpath,filesep,'scale_by_scale.mat'],'planes','r_max',...
    'sfc_calc_points','psd_calc_points',...
    'MOM_vec','SFC_vec','PSD_vec','rawSFC_vec','rawPSD_vec')
