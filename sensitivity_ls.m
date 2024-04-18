
planes = {'ATR-EUREC4A','C130-RICO','C130-VOCALS-REx','TO-POST'};
% planes = {'TO-POST'};

% psd_bottom = [16 16 16 8];
% sfc_bottom = psd_bottom/2;

psd_upp_factor = [0.8 1 1.2 1.4 2 3];
sfc_upp_factor = [psd_upp_factor/2 psd_upp_factor];
psd_upp_factor = [psd_upp_factor psd_upp_factor];



sfc_method = "logmean";
sfc_fit_points = 5;

psd_method = "logmean";
psd_fit_points = 5;
psd_win_length = 1000; % m
psd_win_overlap = 500; % m

warning('off','LOGMEAN:EmptyBins')
warning('off','FIT_PSD:InvalidFitRange')
warning('off','backtrace')

addpath(genpath(myprojectpath))
plotpathB = [myprojectpath,filesep,'figures',filesep,'sensitivity_ls'];
if ~isfolder(plotpathB), mkdir(plotpathB), end



%% Iterate computations

Npl = numel(planes);
Nex = numel(psd_upp_factor);
MOM_matrix = cell(Npl,Nex);


for i_p = 1:Npl    
    plane = planes{i_p};
    datapath = [mydatapath,filesep,plane];
    fprintf('%s\n',plane)
    
    
    sfc_bot_factor = 2;
    psd_bot_factor = 4;
    
    if strcmp(plane,'ATR-EUREC4A')
        dirvar = 'dir2';
        [TURB,MOM] = load_eureca_all(datapath);
    elseif strcmp(plane,'C130-RICO')
        dirvar = '';
        [TURB,MOM] = load_rico_all(datapath);
    elseif strcmp(plane,'C130-VOCALS-REx')
        dirvar = '';
        [TURB,MOM] = load_vocals_all(datapath);
    elseif strcmp(plane,'TO-POST')
        dirvar = 'dir2';
        [TURB,MOM] = load_post_all(datapath);
        sfc_bot_factor = 3;
        psd_bot_factor = 6;
    end
    
    levels = sortrows(groupsummary(MOM,{'level'},{'mean'},{'alt'}),'mean_alt','descend').level';
    for i_l = 1:numel(levels)
        MOM.level_id(MOM.level==levels(i_l)) = i_l;
    end
    
    
    Nseg = size(MOM,1);
    
    fprintf('%2s','L ')
    MOM.int_scale = nan(Nseg,1);
    for i_s = 1:Nseg
%         fprintf(' %d',i_s)
        MOM.int_scale(i_s) = int_ls_short(detrend(TURB(i_s).W),...
            'Method','e-decay','dr',MOM.dr(i_s),'MaxLag',floor(10e3/MOM.dr(i_s)));
    end
    fprintf('\n')

    
    for i_e = 1:Nex
        ex = sprintf('sfc %.1fL psd %.1fL',sfc_upp_factor(i_e),psd_upp_factor(i_e));
        fprintf('%s\n',ex)

        MOM.sfc_fit_range = [MOM.dr*sfc_bot_factor MOM.int_scale*sfc_upp_factor(i_e)];
        MOM.psd_fit_range = [MOM.dr*psd_bot_factor MOM.int_scale*psd_upp_factor(i_e)];
        
        vars = {'UX','VY','W'};
        Nvar = numel(vars);
        
        for i_v = 1:Nvar
            var = vars{i_v}; fprintf('%2s',var)
            
            for i_s = 1:Nseg
%                 fprintf(' %d',i_s)
                dr = MOM.dr(i_s);
                
                try
                    [MOM.(['off_sfc_',var])(i_s),MOM.(['slp_sfc_',var])(i_s),es] = ...
                        fit_sfc( detrend(TURB(i_s).(var)), dr, MOM.sfc_fit_range(i_s,:), ...
                        'Method',sfc_method, 'FitPoints',sfc_fit_points );
                    MOM.(['e_off_sfc_',var])(i_s) = es.O;
                    MOM.(['e_slp_sfc_',var])(i_s) = es.slp;
                    MOM.sfc_fit_points(i_s) = es.N;
                catch ME
                    if strcmp(ME.identifier,'FIT_SFC:TooFewFitPoints')
                        fprintf('\n')
                        warning(ME.identifier,['In seg %d: ',getReport(ME,'basic'),' Assign NaN.'],i_s)
                        MOM.(['off_sfc_',var])(i_s) = nan;
                        MOM.(['slp_sfc_',var])(i_s) = nan; 
                        MOM.(['e_off_sfc_',var])(i_s) = nan;
                        MOM.(['e_slp_sfc_',var])(i_s) = nan;
                        MOM.sfc_fit_points(i_s) = 0;
                    else
                        throw(ME)
                    end
                end
                
                try
                    [MOM.(['off_psd_',var])(i_s),MOM.(['slp_psd_',var])(i_s),ep] = ...
                        fit_psd( detrend(TURB(i_s).(var)), dr, MOM.psd_fit_range(i_s,:), ...
                        'Method',psd_method, 'FitPoints',psd_fit_points, ...
                        'WindowLength',floor(psd_win_length/dr), 'WindowOverlap',floor(psd_win_overlap/dr) );
                    MOM.(['e_off_psd_',var])(i_s) = ep.O;
                    MOM.(['e_slp_psd_',var])(i_s) = ep.slp;
                    MOM.psd_fit_points(i_s) = ep.N;
                catch ME
                    if strcmp(ME.identifier,'FIT_PSD:TooFewFitPoints')
                        fprintf('\n')
                        warning(ME.identifier,['In seg %d: ',getReport(ME,'basic'),' Assign NaN.'],i_s)
                        MOM.(['off_psd_',var])(i_s) = nan;
                        MOM.(['slp_psd_',var])(i_s) = nan; 
                        MOM.(['e_off_psd_',var])(i_s) = nan;
                        MOM.(['e_slp_psd_',var])(i_s) = nan;
                        MOM.psd_fit_points(i_s) = 0;
                    else
                        throw(ME)
                    end
                end
                
                if MOM.sfc_fit_points(i_s) < 3
                    MOM.(['off_sfc_',var])(i_s) = nan;
                    MOM.(['slp_sfc_',var])(i_s) = nan; 
                    MOM.(['e_off_sfc_',var])(i_s) = nan;
                    MOM.(['e_slp_sfc_',var])(i_s) = nan;
                end
                
                if MOM.psd_fit_points(i_s) < 3
                    MOM.(['off_psd_',var])(i_s) = nan;
                    MOM.(['slp_psd_',var])(i_s) = nan; 
                    MOM.(['e_off_psd_',var])(i_s) = nan;
                    MOM.(['e_slp_psd_',var])(i_s) = nan;
                end
                    
            end
 
            MOM.(['slp_psd_',var]) = - MOM.(['slp_psd_',var]);
            fprintf('\n')
        end
        
        MOM.ar_sfc_VU = MOM.off_sfc_VY./MOM.off_sfc_UX;
        MOM.ar_sfc_WU = MOM.off_sfc_W ./MOM.off_sfc_UX;
        MOM.ar_sfc_WV = MOM.off_sfc_W ./MOM.off_sfc_VY;

        MOM.ar_psd_VU = MOM.off_psd_VY./MOM.off_psd_UX;
        MOM.ar_psd_WU = MOM.off_psd_W ./MOM.off_psd_UX;
        MOM.ar_psd_WV = MOM.off_psd_W ./MOM.off_psd_VY;
        
        MOM.e_ar_sfc_VU = MOM.ar_sfc_VU .* sqrt( (MOM.e_off_sfc_VY./MOM.off_sfc_VY).^2 + (MOM.e_off_sfc_UX./MOM.off_sfc_UX).^2 );
        MOM.e_ar_sfc_WU = MOM.ar_sfc_WU .* sqrt( (MOM.e_off_sfc_W ./MOM.off_sfc_W ).^2 + (MOM.e_off_sfc_UX./MOM.off_sfc_UX).^2 );
        MOM.e_ar_sfc_WV = MOM.ar_sfc_WV .* sqrt( (MOM.e_off_sfc_W ./MOM.off_sfc_W ).^2 + (MOM.e_off_sfc_VY./MOM.off_sfc_VY).^2 );

        MOM.e_ar_psd_VU = MOM.ar_psd_VU .* sqrt( (MOM.e_off_psd_VY./MOM.off_psd_VY).^2 + (MOM.e_off_psd_UX./MOM.off_psd_UX).^2 );
        MOM.e_ar_psd_WU = MOM.ar_psd_WU .* sqrt( (MOM.e_off_psd_W ./MOM.off_psd_W ).^2 + (MOM.e_off_psd_UX./MOM.off_psd_UX).^2 );
        MOM.e_ar_psd_WV = MOM.ar_psd_WV .* sqrt( (MOM.e_off_psd_W ./MOM.off_psd_W ).^2 + (MOM.e_off_psd_VY./MOM.off_psd_VY).^2 );
        
        
        MOM.plane(:) = string(plane);
        MOM.ex(:) = string(ex);
        MOM_matrix{i_p,i_e} = MOM;
        
        
        plotpathF = [plotpathB,filesep,replace(ex,{'.',' ',',','/'},{'p','_','',''})];
        if ~isfolder(plotpathF), mkdir(plotpathF), end
        plotpath = [plotpathF,filesep,plane];
        
        % (Pv/Pu,Dv/Du) scatter
        [fig,ax] = plot_xy_uni(MOM,{'ar_psd_VU'},{'ar_sfc_VU'},'level_id','',dirvar,true,...
            {'ver3/4','ver4/3','hor3/4','hor4/3','cross1'},[],'XLim',[0 1.5],'YLim',[0 1.5]);
        plot(ax,4/3,4/3,'d','Color',"#77AC30",'MarkerFaceColor',"#77AC30",'MarkerSize',12)
        legend(horzcat(levels,{'HIT'}),'Location','northwest','Interpreter','latex')
        xlabel(sprintf('$P_v/P_u$ in range [%.0f$\\Delta r$, %.1f$L$]',psd_bot_factor,psd_upp_factor(i_e)),'Interpreter','latex')
        ylabel(sprintf('$D_v/D_u$ in range [%.0f$\\Delta r$, %.1f$L$]',sfc_bot_factor,sfc_upp_factor(i_e)),'Interpreter','latex')
        title(plane)
        print(fig,[plotpath,'_ar_uv'],'-dpng','-r300')

        % (Pw/Pu,Dw/Du) scatter
        [fig,ax] = plot_xy_uni(MOM,{'ar_psd_WU'},{'ar_sfc_WU'},'level_id','',dirvar,true,...
            {'ver3/4','ver4/3','hor3/4','hor4/3','cross1'},[],'XLim',[0 1.5],'YLim',[0 1.5]);
        plot(ax,4/3,4/3,'d','Color',"#77AC30",'MarkerFaceColor',"#77AC30",'MarkerSize',12)
        legend(horzcat(levels,{'HIT'}),'Location','northwest','Interpreter','latex')
        xlabel(sprintf('$P_w/P_u$ in range [%.0f$\\Delta r$, %.1f$L$]',psd_bot_factor,psd_upp_factor(i_e)),'Interpreter','latex')
        ylabel(sprintf('$D_w/D_u$ in range [%.0f$\\Delta r$, %.1f$L$]',sfc_bot_factor,sfc_upp_factor(i_e)),'Interpreter','latex')
        title(plane)
        print(fig,[plotpath,'_ar_uw'],'-dpng','-r300')

        % (p,s) scatter
        [fig,ax] = plot_xy_uni(MOM,{'slp_psd_UX','slp_psd_VY','slp_psd_W'},...
            {'slp_sfc_UX','slp_sfc_VY','slp_sfc_W'},'yN','level_id',dirvar,false,...
            {'hor2/3','ver5/3'},[],'XLim',[0 2.5],'YLim',[0 1.5]);
        plot(ax,5/3,2/3,'d','Color',"#77AC30",'MarkerFaceColor',"#77AC30",'MarkerSize',12)
        legend(horzcat({'u','v','w'},levels,{'K41'}),'Location','northwest','Interpreter','latex')
        xlabel(sprintf('$p$ in range [%.0f$\\Delta r$, %.1f$L$]',psd_bot_factor,psd_upp_factor(i_e)),'Interpreter','latex')
        ylabel(sprintf('$s$ in range [%.0f$\\Delta r$, %.1f$L$]',sfc_bot_factor,sfc_upp_factor(i_e)),'Interpreter','latex')
        title(plane)
        print(fig,[plotpath,'_slp'],'-dpng','-r300')
        
    end
    
%     close all
    clear MOM TURB levels
        
end


% save('sensitivity_ls.mat','planes','psd_factor','sfc_factor','MOM_matrix')



%% Summary

Nex = numel(psd_upp_factor);
Npl = numel(planes);


sum_vars = {'ar_sfc_VU','ar_psd_VU','ar_sfc_WU','ar_psd_WU',...
    'slp_sfc_UX','slp_sfc_VY','slp_sfc_W',...
    'slp_psd_UX','slp_psd_VY','slp_psd_W',...
    'sfc_fit_points','psd_fit_points'};

keep_vars = horzcat({'ex','plane','flight','name','level','alt'},sum_vars);


temp =  cellfun(@(x) x(:,keep_vars),MOM_matrix(:),'UniformOutput',false);
MOM = vertcat(temp{:});

CNT = groupcounts(MOM,{'ex','plane','sfc_fit_points'});
CNT3 = CNT(CNT.sfc_fit_points<3,:);

MOM = MOM(1:size(MOM,1)/2,:);
% MOM = MOM(1+size(MOM,1)/2:end,:);
SUM = groupsummary(MOM,{'ex','plane'},{'nummissing','mean','std'},sum_vars);



%% Comparison plots

plotpath = [plotpathB,filesep,'comparison1',filesep];
if ~isfolder(plotpath), mkdir(plotpath), end

mks = 10;

experiments = unique(SUM.ex)';
planes = unique(SUM.plane)';


fig = plot_xy_uni(SUM,{'mean_ar_psd_VU'},{'mean_ar_sfc_VU'},...
    'ex','plane',[],true,{'cross1','ver3/4','hor3/4','ver4/3','hor4/3'},mks,...
    'XLim',[0 1.5],'YLim',[0 1.5]);
xlabel('$P_v/P_u$','Interpreter','latex')
ylabel('$D_v/D_u$','Interpreter','latex')
print(fig,[plotpath,'ar_vu'],'-dpng','-r300')

fig = plot_xy_uni(SUM,{'mean_ar_psd_WU'},{'mean_ar_sfc_WU'},...
    'ex','plane','yN',true,{'cross1','ver3/4','hor3/4','hor3/4','ver4/3','hor4/3'},mks,...
    'XLim',[0 1.5],'YLim',[0 1.5]);
xlabel('$P_w/P_u$','Interpreter','latex')
ylabel('$D_w/D_u$','Interpreter','latex')
print(fig,[plotpath,'ar_wu'],'-dpng','-r300')

[fig,ax] = plot_xy_uni(SUM,{'mean_ar_psd_VU'},{'mean_ar_sfc_VU'},...
    'ex','plane','yN',true,{'cross1','ver3/4','hor3/4','hor3/4','ver4/3','hor4/3'},mks);
ax.Position = ax.Position + [-10 0 0 0];
legend( horzcat(experiments,planes),'Position',[0.5 0.5 0 0],'Interpreter','latex')
print(fig,[plotpath,'legend'],'-dpng','-r300')


fig = plot_xy_uni(SUM,{'mean_slp_psd_UX'},{'mean_slp_sfc_UX'},...
    'ex','plane',[],false,{'ver5/3','hor2/3'},mks,'XLim',[0 2.5],'YLim',[0 1.5]);
xlabel('$p_u$','Interpreter','latex')
ylabel('$s_u$','Interpreter','latex')
print(fig,[plotpath,'slp_u'],'-dpng','-r300')

fig = plot_xy_uni(SUM,{'mean_slp_psd_VY'},{'mean_slp_sfc_VY'},...
    'ex','plane',[],false,{'ver5/3','hor2/3'},mks,'XLim',[0 2.5],'YLim',[0 1.5]);
xlabel('$p_v$','Interpreter','latex')
ylabel('$s_v$','Interpreter','latex')
print(fig,[plotpath,'slp_v'],'-dpng','-r300')

fig = plot_xy_uni(SUM,{'mean_slp_psd_W'},{'mean_slp_sfc_W'},...
    'ex','plane',[],false,{'ver5/3','hor2/3'},mks,'XLim',[0 2.5],'YLim',[0 1.5]);
xlabel('$p_w$','Interpreter','latex')
ylabel('$s_w$','Interpreter','latex')
print(fig,[plotpath,'slp_w'],'-dpng','-r300')