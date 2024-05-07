
%% Settings

% List of planes/experiments
planes = {'ATR-EUREC4A','C130-RICO','C130-VOCALS-REx','TO-POST'};


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


addpath(genpath(myprojectpath))
plotpath_base = [myprojectpath,filesep,'figures',filesep,'sensitivity'];
if ~isfolder(plotpath_base), mkdir(plotpath_base), end



%% Computations

Npl = numel(planes);
Nvar = numel(vars);
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
        
        
        % Structure functions
        
        fprintf('%s','SFC ... ')

        MOM.sfc_fit_range = [MOM.dr*sfc_bot_factors(i_p) MOM.int_scale*sfc_upp_factors(i_f)];
        
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
        
        
        % Power spectra
        
        fprintf('%s','PSD ... ')
        
        MOM.psd_fit_range = [MOM.dr*psd_bot_factors(i_p) MOM.int_scale*psd_upp_factors(i_f)];
        
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



%% Save/load

save('sensitivity.mat','MOM_matrix','planes',...
    'psd_upp_factors','sfc_upp_factors','psd_bot_factors','sfc_bot_factors',...
    'plotpath_base')

% addpath(genpath(myprojectpath))
% load('sensitivity.mat')



%% Plot individual results

% Ax limits
ratio_lim = [0 2];
p_lim = [0 2.5];
s_lim = [0 1.5];


Npl = size(MOM_matrix,1);
Nfc = size(MOM_matrix,2);


for i_p = 1:Npl 
    plane = planes{i_p};
    
    if ismember(plane,{'ATR-EUREC4A','TO-POST'})
        dirvar = 'dir2';
    else
        dirvar = '';
    end
    
    
    for i_f = 1:Nfc
        testname = sprintf('sfc %.1fL psd %.1fL',sfc_upp_factors(i_f),psd_upp_factors(i_f));
        
        plotpath = [plotpath_base,filesep,replace(testname,{'.',' ',',','/'},{'p','_','',''})];
        if ~isfolder(plotpath), mkdir(plotpath), end
        
        MOM = MOM_matrix{i_p,i_f};

        levels = sortrows(groupsummary(MOM,{'level'},{'mean'},{'alt'}),'mean_alt','descend').level';
        for i_l = 1:numel(levels)
            MOM.level_id(MOM.level==levels(i_l)) = i_l;
        end
               

        % (Pv/Pu,Dv/Du)
        
        [fig,ax] = plot_xy_uni(MOM,{'ar_psd_VU'},{'ar_sfc_VU'},'level_id','',dirvar,true,...
            {'ver3/4','ver4/3','hor3/4','hor4/3','cross1'},[],'XLim',ratio_lim,'YLim',ratio_lim);
        plot(ax,4/3,4/3,'d','Color',"#77AC30",'MarkerFaceColor',"#77AC30",'MarkerSize',12)
        legend(horzcat(levels,{'HIT'}),'Location','northwest','Interpreter','latex')
        xlabel(sprintf('$P_v/P_u$ in range [%.0f$\\Delta r$, %.1f$L$]',psd_bot_factors(i_p),psd_upp_factors(i_f)),'Interpreter','latex')
        ylabel(sprintf('$D_v/D_u$ in range [%.0f$\\Delta r$, %.1f$L$]',sfc_bot_factors(i_p),sfc_upp_factors(i_f)),'Interpreter','latex')
        title(plane)
        print(fig,[plotpath,filesep,plane,'_ar_uv'],'-dpng','-r300')

        
        % (Pw/Pu,Dw/Du)
        
        [fig,ax] = plot_xy_uni(MOM,{'ar_psd_WU'},{'ar_sfc_WU'},'level_id','',dirvar,true,...
            {'ver3/4','ver4/3','hor3/4','hor4/3','cross1'},[],'XLim',ratio_lim,'YLim',ratio_lim);
        plot(ax,4/3,4/3,'d','Color',"#77AC30",'MarkerFaceColor',"#77AC30",'MarkerSize',12)
        legend(horzcat(levels,{'HIT'}),'Location','northwest','Interpreter','latex')
        xlabel(sprintf('$P_w/P_u$ in range [%.0f$\\Delta r$, %.1f$L$]',psd_bot_factors(i_p),psd_upp_factors(i_f)),'Interpreter','latex')
        ylabel(sprintf('$D_w/D_u$ in range [%.0f$\\Delta r$, %.1f$L$]',sfc_bot_factors(i_p),sfc_upp_factors(i_f)),'Interpreter','latex')
        title(plane)
        print(fig,[plotpath,filesep,plane,'_ar_uw'],'-dpng','-r300')

        
        % (p,s)
        
        [fig,ax] = plot_xy_uni(MOM,{'slp_psd_UX','slp_psd_VY','slp_psd_W'},...
            {'slp_sfc_UX','slp_sfc_VY','slp_sfc_W'},'yN','level_id',dirvar,false,...
            {'hor2/3','ver5/3'},[],'XLim',p_lim,'YLim',s_lim);
        plot(ax,5/3,2/3,'d','Color',"#77AC30",'MarkerFaceColor',"#77AC30",'MarkerSize',12)
        legend(horzcat({'u','v','w'},levels,{'K41'}),'Location','northwest','Interpreter','latex')
        xlabel(sprintf('$p$ in range [%.0f$\\Delta r$, %.1f$L$]',psd_bot_factors(i_p),psd_upp_factors(i_f)),'Interpreter','latex')
        ylabel(sprintf('$s$ in range [%.0f$\\Delta r$, %.1f$L$]',sfc_bot_factors(i_p),sfc_upp_factors(i_f)),'Interpreter','latex')
        title(plane)
        print(fig,[plotpath,filesep,plane,'_slp'],'-dpng','-r300')
        
      
        clear MOM levels     
    end
    
    close all
    
end


% Printout a piece of a jam script

jamplots = {'ar_uv','ar_uw','slp'};
for i_j = 1:numel(jamplots)
    fprintf('pdfjam --nup %dx%d --papersize ''{%dcm,%dcm}''',Npl,Nfc,Npl*16,Nfc*12)
    for i_f = 1:Nfc
        testname = sprintf('sfc %.1fL psd %.1fL',sfc_upp_factors(i_f),psd_upp_factors(i_f));
        fprintf(' ./%s/*%s.png',replace(testname,{'.',' ',',','/'},{'p','_','',''}),jamplots{i_j})
    end
    fprintf(' --outfile ./jam/%s.pdf\n',jamplots{i_j})
end
fprintf('pdfjam --papersize ''{%dcm,%dcm}''',Npl*16,Nfc*12)
fprintf(' ./jam/%s.pdf',jamplots{:})
fprintf(' --outfile ./jam/scatter.pdf\n')

        

%% Summary

sum_vars = {'ar_sfc_VU','ar_psd_VU','ar_sfc_WU','ar_psd_WU',...
    'slp_sfc_UX','slp_sfc_VY','slp_sfc_W',...
    'slp_psd_UX','slp_psd_VY','slp_psd_W',...
    'R2_sfc_UX','R2_sfc_VY','R2_sfc_W',...
    'R2_psd_UX','R2_psd_VY','R2_psd_W',...
    'R2_sfc_tot','R2_psd_tot','R2_tot_tot',...
    'N_sfc_UX','N_psd_UX'};


% Propagate settings

Npl = size(MOM_matrix,1);
Nfc = size(MOM_matrix,2);

for i_p = 1:Npl 
    plane = planes{i_p};
      
    for i_f = 1:Nfc
        testname = sprintf('sfc %.1fL psd %.1fL',sfc_upp_factors(i_f),psd_upp_factors(i_f));
        
        MOM_matrix{i_p,i_f}.plane(:) = string(plane);
        MOM_matrix{i_p,i_f}.sfc_bot_factor(:) = sfc_bot_factors(i_p);
        MOM_matrix{i_p,i_f}.psd_bot_factor(:) = psd_bot_factors(i_p);
        
        MOM_matrix{i_p,i_f}.testname(:) = string(testname);
        MOM_matrix{i_p,i_f}.sfc_upp_factor(:) = sfc_upp_factors(i_f);
        MOM_matrix{i_p,i_f}.psd_upp_factor(:) = psd_upp_factors(i_f);
        
        MOM_matrix{i_p,i_f}.R2_sfc_tot = MOM_matrix{i_p,i_f}.R2_sfc_UX .* MOM_matrix{i_p,i_f}.R2_sfc_VY .* MOM_matrix{i_p,i_f}.R2_sfc_W;
        MOM_matrix{i_p,i_f}.R2_psd_tot = MOM_matrix{i_p,i_f}.R2_psd_UX .* MOM_matrix{i_p,i_f}.R2_psd_VY .* MOM_matrix{i_p,i_f}.R2_psd_W;
        MOM_matrix{i_p,i_f}.R2_tot_tot = MOM_matrix{i_p,i_f}.R2_sfc_tot .* MOM_matrix{i_p,i_f}.R2_psd_tot;
    end
    
end


% Generate group summary table

keep_vars = horzcat({'testname','sfc_upp_factor','psd_upp_factor',...
    'plane','flight','name','level','alt'},sum_vars);

temp =  cellfun(@(x) x(:,keep_vars),MOM_matrix(:),'UniformOutput',false);
MOM = vertcat(temp{:});

MOM.rel_upp_factor(:) = "shifted";
MOM.rel_upp_factor(MOM.sfc_upp_factor==MOM.psd_upp_factor) = "equal";


CNT = groupcounts(MOM,{'testname','sfc_upp_factor','psd_upp_factor',...
    'plane','N_sfc_UX'});

SUM2 = groupsummary(MOM,{'rel_upp_factor','testname','sfc_upp_factor','plane'},...
    {'mean'},sum_vars);



%% Summary plots

% Ax limits
ratio_lim = [0 1.6];
p_lim = [0 2.5];
s_lim = [0 1.3];

rel_upp_factor = "shifted";
SUM = SUM2(SUM2.rel_upp_factor==rel_upp_factor,:);

% rel_upp_factor = "both";
% SUM = SUM2;

plotpath = [plotpath_base,filesep,char(rel_upp_factor),filesep];
if ~isfolder(plotpath), mkdir(plotpath), end


mks = 10;

tests = unique(SUM.testname)';
planes = unique(SUM.plane)';


% Ratios

fig = plot_xy_uni(SUM,{'mean_ar_psd_VU'},{'mean_ar_sfc_VU'},...
    'testname','plane',[],true,{'cross1','ver3/4','hor3/4','ver4/3','hor4/3'},mks,...
    'XLim',ratio_lim,'YLim',ratio_lim);
xlabel('$P_v/P_u$','Interpreter','latex')
ylabel('$D_v/D_u$','Interpreter','latex')
print(fig,[plotpath,'ar_vu'],'-dpng','-r300')

fig = plot_xy_uni(SUM,{'mean_ar_psd_WU'},{'mean_ar_sfc_WU'},...
    'testname','plane',[],true,{'cross1','ver3/4','hor3/4','hor3/4','ver4/3','hor4/3'},mks,...
    'XLim',ratio_lim,'YLim',ratio_lim);
xlabel('$P_w/P_u$','Interpreter','latex')
ylabel('$D_w/D_u$','Interpreter','latex')
print(fig,[plotpath,'ar_wu'],'-dpng','-r300')

[fig,ax] = plot_xy_uni(SUM,{'mean_ar_psd_VU'},{'mean_ar_sfc_VU'},...
    'testname','plane',[],true,{'cross1','ver3/4','hor3/4','hor3/4','ver4/3','hor4/3'},mks);
ax.Position = ax.Position + [-10 0 0 0];
legend( horzcat(tests,planes),'Position',[0.5 0.5 0 0],'Interpreter','latex')
print(fig,[plotpath,'legend'],'-dpng','-r300')


% Exponents

fig = plot_xy_uni(SUM,{'mean_slp_psd_UX'},{'mean_slp_sfc_UX'},...
    'testname','plane',[],false,{'ver5/3','hor2/3'},mks,'XLim',p_lim,'YLim',s_lim);
xlabel('$p_u$','Interpreter','latex')
ylabel('$s_u$','Interpreter','latex')
print(fig,[plotpath,'slp_u'],'-dpng','-r300')

fig = plot_xy_uni(SUM,{'mean_slp_psd_VY'},{'mean_slp_sfc_VY'},...
    'testname','plane',[],false,{'ver5/3','hor2/3'},mks,'XLim',p_lim,'YLim',s_lim);
xlabel('$p_v$','Interpreter','latex')
ylabel('$s_v$','Interpreter','latex')
print(fig,[plotpath,'slp_v'],'-dpng','-r300')

fig = plot_xy_uni(SUM,{'mean_slp_psd_W'},{'mean_slp_sfc_W'},...
    'testname','plane',[],false,{'ver5/3','hor2/3'},mks,'XLim',p_lim,'YLim',s_lim);
xlabel('$p_w$','Interpreter','latex')
ylabel('$s_w$','Interpreter','latex')
print(fig,[plotpath,'slp_w'],'-dpng','-r300')


% R2

fig = plot_xy_uni(SUM,{'mean_R2_psd_UX'},{'mean_R2_sfc_UX'},...
    'testname','plane',[],true,{'cross1'},mks);
xlabel('$R_{Pu}^2$','Interpreter','latex')
ylabel('$R_{Du}^2$','Interpreter','latex')
print(fig,[plotpath,'R2_u'],'-dpng','-r300')

fig = plot_xy_uni(SUM,{'mean_R2_psd_VY'},{'mean_R2_sfc_VY'},...
    'testname','plane',[],true,{'cross1'},mks);
xlabel('$R_{Pv}^2$','Interpreter','latex')
ylabel('$R_{Dv}^2$','Interpreter','latex')
print(fig,[plotpath,'R2_v'],'-dpng','-r300')

fig = plot_xy_uni(SUM,{'mean_R2_psd_W'},{'mean_R2_sfc_W'},...
    'testname','plane',[],true,{'cross1'},mks);
xlabel('$R_{Pw}^2$','Interpreter','latex')
ylabel('$R_{Dw}^2$','Interpreter','latex')
print(fig,[plotpath,'R2_w'],'-dpng','-r300')

fig = plot_xy_uni(SUM,{'mean_R2_psd_tot'},{'mean_R2_sfc_tot'},...
    'testname','plane',[],true,{'cross1'},mks);
xlabel('$R_{Pu}^2R_{Pv}^2R_{Pw}^2$','Interpreter','latex')
ylabel('$R_{Du}^2R_{Dv}^2R_{Dw}^2$','Interpreter','latex')
print(fig,[plotpath,'R2_tot'],'-dpng','-r300')


[fig,~,co] = fig16x12;
for i_p = 1:Npl
    ind = SUM.plane==planes{i_p};
    plot(SUM.sfc_upp_factor(ind),SUM.mean_R2_tot_tot(ind),'-o','Color',co(i_p,:))
end
% for i_p = 1:Npl
%     ind = SUM.plane==planes{i_p};
%     plot(SUM.sfc_upp_factor(ind),SUM.mean_R2_sfc_tot(ind),'--','Color',co(i_p,:),'HandleVisibility','off')
%     plot(SUM.sfc_upp_factor(ind),SUM.mean_R2_psd_tot(ind),':','Color',co(i_p,:),'HandleVisibility','off')
% end
legend(planes,'Location','best')
xlabel('$F$ sfc','Interpreter','latex')
ylabel('$R_{Du}^2R_{Dv}^2R_{Dw}^2R_{Pu}^2R_{Pv}^2R_{Pw}^2$','Interpreter','latex')
print(fig,[plotpath,'R2_tot_tot'],'-dpng','-r300')


% Rejected segments

sfc_min_fit_points = 3;
psd_min_fit_points = 3;

CNTrej = CNT;
CNTrej.GroupCount(CNTrej.N_sfc_UX>=sfc_min_fit_points) = 0;
CNTrej = groupsummary(CNTrej,{'sfc_upp_factor','plane'},{'sum'},{'GroupCount'});

fig = fig16x12;
for i_p = 1:Npl
    ind = CNTrej.plane==planes{i_p};
    plot(CNTrej.sfc_upp_factor(ind),CNTrej.sum_GroupCount(ind),'-o')
end
legend(planes,'Location','best')
xlabel('$F$ sfc','Interpreter','latex')
ylabel('# rejected segments')
print(fig,[plotpath,'rej_tot'],'-dpng','-r300')


