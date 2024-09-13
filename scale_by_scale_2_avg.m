
%% Settings

% Interpolation variables
vars_interp = {'UX','VY','W','ar_VU','ar_WU'};

% Interpolation points
sfc_interp_points = 10;
psd_interp_points = 16;

% Normalized scale range r/L
r_L_range = [0.01 3];



%% Load data

addpath(genpath(myprojectpath))

load([myprojectpath,filesep,'scale_by_scale.mat'])



%% Interpolation to common grids

disp('Interpolate to common scale grids ...')


% Define normalized grid
r_L_sfc = exp(linspace( log(r_L_range(1)), log(r_L_range(2)), sfc_interp_points))';
r_L_psd = exp(linspace( log(r_L_range(1)), log(r_L_range(2)), psd_interp_points))';


Npl = numel(planes);
for i_p = 1:Npl    
    
    MOM = MOM_vec{i_p};
    SFC = SFC_vec{i_p};
    PSD = PSD_vec{i_p};
    Nseg = size(MOM,1);
    
    
    % Define metric grid
    
    sfc_interp_range = [ max(cellfun(@min,{SFC(:).r})) min(cellfun(@max,{SFC(:).r})) ];
    psd_interp_range = [ max(cellfun(@min,{PSD(:).r})) min(cellfun(@max,{PSD(:).r})) ];
    
    r_sfc = exp(linspace( log(sfc_interp_range(1)), log(sfc_interp_range(2)), sfc_interp_points))';
    r_psd = exp(linspace( log(psd_interp_range(1)), log(psd_interp_range(2)), psd_interp_points))';
    
    
    % Interpolate to metric grid
    
    for i_s = 1:Nseg
        
        for i_v = 1:numel(vars_interp)
            var = vars_interp{i_v};     
            SFC(i_s).([var,'_i']) = interp1( SFC(i_s).r, SFC(i_s).(var), r_sfc, 'linear' );
            PSD(i_s).([var,'_i']) = interp1( PSD(i_s).r, PSD(i_s).(var), r_psd, 'linear' );
        end
        
        SFC(i_s).r_i = r_sfc;
        PSD(i_s).r_i = r_psd;
    end
            
    
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
    
    
    % Save
    
    SFC_vec{i_p} = SFC;
    PSD_vec{i_p} = PSD;
    
    clear MOM SFC PSD
    
end



%% Level averaging

disp('Average at levels ...')

avSFC_vec = cell(Npl,1);
avPSD_vec = cell(Npl,1);


for i_p = 1:Npl    
    
    MOM = MOM_vec{i_p};
    SFC = SFC_vec{i_p};
    PSD = PSD_vec{i_p}; 
    Nseg = size(MOM,1);
    
    
    ff = fieldnames(SFC);
    vars_avg = ff(endsWith(ff,{'_i','_iL'}));
    
    
    levels = sortrows(groupsummary(MOM,{'level'},{'mean'},{'alt'}),'mean_alt','descend').level';
    Nlvl = numel(levels);
    
    
    avSFC = cell(Nlvl,1);
    avPSD = cell(Nlvl,1);
    
    for i_l = 1:Nlvl
        ind_l = (MOM.level==levels(i_l));
        avSFC{i_l}.level = levels(i_l);
        avPSD{i_l}.level = levels(i_l);
        
        for i_v = 1:numel(vars_avg)
            var = vars_avg{i_v};
            
            avSFC{i_l}.(var) = mean( horzcat(SFC(ind_l).(var)), 2 ,'omitnan');
            avPSD{i_l}.(var) = mean( horzcat(PSD(ind_l).(var)), 2, 'omitnan');
        end
    end
    avSFC = vertcat(avSFC{:});
    avPSD = vertcat(avPSD{:});
        
    
    % Save
    
    avSFC_vec{i_p} = avSFC;
    avPSD_vec{i_p} = avPSD;
    
    clear MOM SFC PSD avSFC avPSD
    
end


