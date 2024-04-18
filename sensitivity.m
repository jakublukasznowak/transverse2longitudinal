

experiments = {'reference','shift-down','shift-up','half-down','half-up'};

planes = {'ATR-EUREC4A','C130-RICO','C130-VOCALS-REx','TO-POST'};

sfc_basic = {[8 40],[8 40],[8 40],[4 20]};
psd_basic = {[16 80],[16 80],[16 80],[8 40]};


Nex = numel(experiments);
Npl = numel(planes);

MOM_matrix = cell(Nex,Npl);
sfc_matrix = cell(Nex,Npl);
psd_matrix = cell(Nex,Npl);

for i_e = 1:Nex
    ex = experiments{i_e};
    
    for i_p = 1:Npl    
        plane = planes{i_p};
        
        if strcmp(ex,'reference')
            sfc_fit_range = sfc_basic{i_p};
            psd_fit_range = psd_basic{i_p};
        elseif strcmp(ex,'shift-down')
            sfc_fit_range = sfc_basic{i_p}/2;
            psd_fit_range = psd_basic{i_p}/2;
        elseif strcmp(ex,'shift-up')
            sfc_fit_range = sfc_basic{i_p}*2;
            psd_fit_range = psd_basic{i_p}*2;
        elseif strcmp(ex,'half-down')
%             sfc_fit_range = [sfc_basic{i_p}(1) sqrt(prod(sfc_basic{i_p}))];
%             psd_fit_range = [psd_basic{i_p}(1) sqrt(prod(psd_basic{i_p}))];
            sfc_fit_range = [sfc_basic{i_p}(1) mean(sfc_basic{i_p})];
            psd_fit_range = [psd_basic{i_p}(1) mean(psd_basic{i_p})];
        elseif strcmp(ex,'half-up')
%             sfc_fit_range = [sqrt(prod(sfc_basic{i_p})) sfc_basic{i_p}(2)];
%             psd_fit_range = [sqrt(prod(psd_basic{i_p})) psd_basic{i_p}(2)];
            sfc_fit_range = [mean(sfc_basic{i_p}) sfc_basic{i_p}(2)];
            psd_fit_range = [mean(psd_basic{i_p}) psd_basic{i_p}(2)];
        end
        
        plotpath = [myprojectpath,filesep,'figures',filesep,'sensitivity',filesep,ex];
        
        main
        
        MOM.plane(:) = string(plane);
        MOM.ex(:) = string(ex);
        
        sfc_matrix{i_e,i_p} = sfc_fit_range;
        psd_matrix{i_e,i_p} = psd_fit_range;
        MOM_matrix{i_e,i_p} = MOM;
        
        close all
        clear MOM TURB levels

    end
end

save('sensitivity','experiments','planes','sfc_basic','psd_basic',...
    'MOM_matrix','sfc_matrix','psd_matrix')



%% Summary

Nex = numel(experiments);
Npl = numel(planes);

% Propagate reference ex to every record

prop_vars = {'ar_sfc_VU','ar_psd_VU','ar_sfc_WU','ar_psd_WU',...
    'slp_sfc_UX','slp_sfc_VY','slp_sfc_W',...
    'slp_psd_UX','slp_psd_VY','slp_psd_W'};

for i_e = 1:Nex
    for i_p = 1:Npl
        for i_v = 1:numel(prop_vars)
            MOM_matrix{i_e,i_p}.(['ref_',prop_vars{i_v}]) = MOM_matrix{1,i_p}.(prop_vars{i_v});
        end
    end
end


% Keep selected fields and concatenate between experiments and planes

keep_vars = horzcat({'ex','plane','flight','name','level','alt'},...
    prop_vars,strcat('ref_',prop_vars));

temp =  cellfun(@(x) x(:,keep_vars),MOM_matrix(:),'UniformOutput',false);
MOM = vertcat(temp{:});


% Group summary and stats

sum_vars = horzcat(prop_vars,strcat('ref_',prop_vars));

SUM = groupsummary(MOM,{'ex','plane'},{'mean','std'},sum_vars);

SUM = SUM([find(SUM.ex=='reference');find(SUM.ex~='reference')],:);
experiments = unique(SUM.ex,'stable')';
planes = unique(SUM.plane,'stable')';



%% Comparison plots vs reference

addpath(genpath(myprojectpath))

plotpath = [myprojectpath,filesep,'figures',filesep,'sensitivity',filesep,'comparison',filesep];


fig = plot_xy_uni(SUM,{'mean_ref_ar_sfc_VU','mean_ref_ar_psd_VU'},{'mean_ar_sfc_VU','mean_ar_psd_VU'},...
    'ex','plane','yN',true,{'cross1','ver3/4','hor3/4','ver4/3','hor4/3'},'XLim',[0 1.5],'YLim',[0 1.5]);
xlabel('$D_v/D_u, P_v/P_u$ reference','Interpreter','latex')
ylabel('$D_v/D_u, P_v/P_u$','Interpreter','latex')
print(fig,[plotpath,'ref_ar_vu'],'-dpng','-r300')

fig = plot_xy_uni(SUM,{'mean_ref_ar_sfc_WU','mean_ref_ar_psd_WU'},{'mean_ar_sfc_WU','mean_ar_psd_WU'},...
    'ex','plane','yN',true,{'cross1','ver3/4','hor3/4','hor3/4','ver4/3','hor4/3'},'XLim',[0 1.5],'YLim',[0 1.5]);
xlabel('$D_w/D_u, P_w/P_u$ reference','Interpreter','latex')
ylabel('$D_w/D_u, P_w/P_u$','Interpreter','latex')
print(fig,[plotpath,'ref_ar_wu'],'-dpng','-r300')

[fig,ax] = plot_xy_uni(SUM,{'mean_ref_ar_sfc_WU','mean_ref_ar_psd_WU'},{'mean_ar_sfc_WU','mean_ar_psd_WU'},...
    'ex','plane','yN',true,{'cross1','ver3/4','hor3/4','hor3/4','ver4/3','hor4/3'});
ax.Position = ax.Position + [-10 0 0 0];
legend( horzcat(experiments,planes,{'$D$ ratio','$P$ ratio'}),...
    'Position',[0.5 0.5 0 0],'Interpreter','latex')
print(fig,[plotpath,'ref_legend'],'-dpng','-r300')


fig = plot_xy_uni(SUM,{'mean_ref_slp_sfc_UX'},{'mean_slp_sfc_UX'},...
    'ex','plane',[],true,{'cross1','ver2/3','hor2/3'},'XLim',[0 1.5],'YLim',[0 1.5]);
xlabel('$s_u$ reference','Interpreter','latex')
ylabel('$s_u$','Interpreter','latex')
print(fig,[plotpath,'ref_slp_sfc_u'],'-dpng','-r300')

fig = plot_xy_uni(SUM,{'mean_ref_slp_sfc_VY'},{'mean_slp_sfc_VY'},...
    'ex','plane',[],true,{'cross1','ver2/3','hor2/3'},'XLim',[0 1.5],'YLim',[0 1.5]);
xlabel('$s_v$ reference','Interpreter','latex')
ylabel('$s_v$','Interpreter','latex')
print(fig,[plotpath,'ref_slp_sfc_v'],'-dpng','-r300')

fig = plot_xy_uni(SUM,{'mean_ref_slp_sfc_W'},{'mean_slp_sfc_W'},...
    'ex','plane',[],true,{'cross1','ver2/3','hor2/3'},'XLim',[0 1.5],'YLim',[0 1.5]);
xlabel('$s_w$ reference','Interpreter','latex')
ylabel('$s_w$','Interpreter','latex')
print(fig,[plotpath,'ref_slp_sfc_w'],'-dpng','-r300')

fig = plot_xy_uni(SUM,{'mean_ref_slp_psd_UX'},{'mean_slp_psd_UX'},...
    'ex','plane',[],true,{'cross1','ver5/3','hor5/3'},'XLim',[0 2.5],'YLim',[0 2.5]);
xlabel('$p_u$ reference','Interpreter','latex')
ylabel('$p_u$','Interpreter','latex')
print(fig,[plotpath,'ref_slp_psd_u'],'-dpng','-r300')

fig = plot_xy_uni(SUM,{'mean_ref_slp_psd_VY'},{'mean_slp_psd_VY'},...
    'ex','plane',[],true,{'cross1','ver5/3','hor5/3'},'XLim',[0 2.5],'YLim',[0 2.5]);
xlabel('$p_v$ reference','Interpreter','latex')
ylabel('$p_v$','Interpreter','latex')
print(fig,[plotpath,'ref_slp_psd_v'],'-dpng','-r300')

fig = plot_xy_uni(SUM,{'mean_ref_slp_psd_W'},{'mean_slp_psd_W'},...
    'ex','plane',[],true,{'cross1','ver5/3','hor5/3'},'XLim',[0 2.5],'YLim',[0 2.5]);
xlabel('$p_w$ reference','Interpreter','latex')
ylabel('$p_w$','Interpreter','latex')
print(fig,[plotpath,'ref_slp_psd_w'],'-dpng','-r300')



%% Comparison plots absolute

fig = plot_xy_uni(SUM,{'mean_ar_psd_VU'},{'mean_ar_sfc_VU'},...
    'ex','plane',[],true,{'cross1','ver3/4','hor3/4','ver4/3','hor4/3'},'XLim',[0 1.5],'YLim',[0 1.5]);
xlabel('$P_v/P_u$','Interpreter','latex')
ylabel('$D_v/D_u$','Interpreter','latex')
print(fig,[plotpath,'ar_vu'],'-dpng','-r300')

fig = plot_xy_uni(SUM,{'mean_ar_psd_WU'},{'mean_ar_sfc_WU'},...
    'ex','plane','yN',true,{'cross1','ver3/4','hor3/4','hor3/4','ver4/3','hor4/3'},'XLim',[0 1.5],'YLim',[0 1.5]);
xlabel('$P_w/P_u$','Interpreter','latex')
ylabel('$D_w/D_u$','Interpreter','latex')
print(fig,[plotpath,'ar_wu'],'-dpng','-r300')

[fig,ax] = plot_xy_uni(SUM,{'mean_ar_psd_VU'},{'mean_ar_sfc_VU'},...
    'ex','plane','yN',true,{'cross1','ver3/4','hor3/4','hor3/4','ver4/3','hor4/3'});
ax.Position = ax.Position + [-10 0 0 0];
legend( horzcat(experiments,planes),'Position',[0.5 0.5 0 0],'Interpreter','latex')
print(fig,[plotpath,'legend'],'-dpng','-r300')


fig = plot_xy_uni(SUM,{'mean_slp_psd_UX'},{'mean_slp_sfc_UX'},...
    'ex','plane',[],false,{'ver5/3','hor2/3'},'XLim',[0 2.5],'YLim',[0 1.5]);
xlabel('$p_u$','Interpreter','latex')
ylabel('$s_u$','Interpreter','latex')
print(fig,[plotpath,'slp_u'],'-dpng','-r300')

fig = plot_xy_uni(SUM,{'mean_slp_psd_VY'},{'mean_slp_sfc_VY'},...
    'ex','plane',[],false,{'ver5/3','hor2/3'},'XLim',[0 2.5],'YLim',[0 1.5]);
xlabel('$p_v$','Interpreter','latex')
ylabel('$s_v$','Interpreter','latex')
print(fig,[plotpath,'slp_v'],'-dpng','-r300')

fig = plot_xy_uni(SUM,{'mean_slp_psd_W'},{'mean_slp_sfc_W'},...
    'ex','plane',[],false,{'ver5/3','hor2/3'},'XLim',[0 2.5],'YLim',[0 1.5]);
xlabel('$p_w$','Interpreter','latex')
ylabel('$s_w$','Interpreter','latex')
print(fig,[plotpath,'slp_w'],'-dpng','-r300')