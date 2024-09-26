


% Axes limits
lim_sum.ar = [0.6 1.4];
lim_sum.s  = [0.3 1.1];
lim_sum.p  = [1.0 2.2];
lim_bulk.ar= [0 2];
lim_bulk.s = [0 1.5];
lim_bulk.p = [0 2.5];


% Plot style
mks = 10;
lw = 2;

% Plot format
pfrm = '-dpng';

% Plot resolution
pres = '-r300';


% Output path
addpath(genpath(myprojectpath))
plotpath = [myprojectpath,filesep,'figures',filesep,'sensitivity'];
if ~isfolder(plotpath), mkdir(plotpath), end

% Load results
load([myprojectpath,filesep,'sensitivity.mat'])
Npl = size(MOM_matrix,1);
Nfc = size(MOM_matrix,2);



%% Temporary


sum_vars = {'ar_sfc_VU','ar_psd_VU','ar_sfc_WU','ar_psd_WU',...
    'slp_sfc_UX','slp_sfc_VY','slp_sfc_W',...
    'slp_psd_UX','slp_psd_VY','slp_psd_W',...
    'R2_sfc_tot','R2_psd_tot','R2_tot_tot',...
    'N_sfc_UX','N_psd_UX'};

grp_vars = {'testname','sfc_upp_factor','psd_upp_factor','plane'};


% Propagate settings



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

temp =  cellfun(@(x) x(:,horzcat(grp_vars,sum_vars)),MOM_matrix(:),'UniformOutput',false);
MOM = vertcat(temp{:});

SUM = groupsummary(MOM,grp_vars,{'mean'},sum_vars);

CNT = groupcounts(MOM,horzcat(grp_vars,{'N_sfc_UX'}));



%% Summary plots

plotpath_sum = [plotpath,filesep,'summary'];
if ~isfolder(plotpath_sum), mkdir(plotpath_sum), end


tests = unique(SUM.testname)';
planes = unique(SUM.plane)';


% Ratios

fig = plot_xy_uni(SUM,{'mean_ar_psd_VU'},{'mean_ar_sfc_VU'},...
    'testname','plane',[],true,{'cross1','ver4/3','hor4/3'},mks,...
    'XLim',lim_sum.ar,'YLim',lim_sum.ar);
xlabel('$P_v/P_u$','Interpreter','latex')
ylabel('$D_v/D_u$','Interpreter','latex')
print(fig,[plotpath_sum,filesep,'ar_vu'],pfrm,pres)

fig = plot_xy_uni(SUM,{'mean_ar_psd_WU'},{'mean_ar_sfc_WU'},...
    'testname','plane',[],true,{'cross1','ver4/3','hor4/3'},mks,...
    'XLim',lim_sum.ar,'YLim',lim_sum.ar);
xlabel('$P_w/P_u$','Interpreter','latex')
ylabel('$D_w/D_u$','Interpreter','latex')
print(fig,[plotpath_sum,filesep,'ar_wu'],pfrm,pres)


% Legend

[fig,ax] = plot_xy_uni(SUM,{'mean_ar_psd_VU'},{'mean_ar_sfc_VU'},...
    'testname','plane',[],true,{'cross1','ver4/3','hor4/3'},mks);
ax.Position = ax.Position + [-10 0 0 0];
legend( horzcat(tests,planes),'Position',[0.5 0.5 0 0],'Interpreter','latex','FontSize',18)
print(fig,[plotpath_sum,filesep,'legend'],pfrm,pres)


% Exponents

fig = plot_xy_uni(SUM,{'mean_slp_psd_UX'},{'mean_slp_sfc_UX'},...
    'testname','plane',[],false,{'ver5/3','hor2/3'},mks,'XLim',lim_sum.p,'YLim',lim_sum.s);
xlabel('$p_u$','Interpreter','latex')
ylabel('$s_u$','Interpreter','latex')
print(fig,[plotpath_sum,filesep,'slp_u'],pfrm,pres)

fig = plot_xy_uni(SUM,{'mean_slp_psd_VY'},{'mean_slp_sfc_VY'},...
    'testname','plane',[],false,{'ver5/3','hor2/3'},mks,'XLim',lim_sum.p,'YLim',lim_sum.s);
xlabel('$p_v$','Interpreter','latex')
ylabel('$s_v$','Interpreter','latex')
print(fig,[plotpath_sum,filesep,'slp_v'],pfrm,pres)

fig = plot_xy_uni(SUM,{'mean_slp_psd_W'},{'mean_slp_sfc_W'},...
    'testname','plane',[],false,{'ver5/3','hor2/3'},mks,'XLim',lim_sum.p,'YLim',lim_sum.s);
xlabel('$p_w$','Interpreter','latex')
ylabel('$s_w$','Interpreter','latex')
print(fig,[plotpath_sum,filesep,'slp_w'],pfrm,pres)


% R2

% fig = plot_xy_uni(SUM,{'mean_R2_psd_tot'},{'mean_R2_sfc_tot'},...
%     'testname','plane',[],true,{'cross1'},mks);
% xlabel('$R_{Pu}^2R_{Pv}^2R_{Pw}^2$','Interpreter','latex')
% ylabel('$R_{Du}^2R_{Dv}^2R_{Dw}^2$','Interpreter','latex')
% print(fig,[plotpath_sum,filesep,'R2_tot'],'-dpng','-r300')
% 
% [fig,~,co] = fig16x12;
% for i_p = 1:Npl
%     ind = SUM.plane==planes{i_p};
%     plot(SUM.sfc_upp_factor(ind),SUM.mean_R2_tot_tot(ind),'-o','Color',co(i_p,:))
% end
% legend(planes,'Location','best')
% xlabel('$F$ sfc','Interpreter','latex')
% ylabel('$R_{Du}^2R_{Dv}^2R_{Dw}^2R_{Pu}^2R_{Pv}^2R_{Pw}^2$','Interpreter','latex')
% print(fig,[plotpath_sum,filesep,'R2_tot_tot'],pfrm,pres)


% Rejected segments

% CNTrej = CNT;
% CNTrej.GroupCount(CNTrej.N_sfc_UX>=sfc_min_fit_points) = 0;
% CNTrej = groupsummary(CNTrej,{'sfc_upp_factor','plane'},{'sum'},{'GroupCount'});
% 
% fig = fig16x12;
% for i_p = 1:Npl
%     ind = CNTrej.plane==planes{i_p};
%     plot(CNTrej.sfc_upp_factor(ind),CNTrej.sum_GroupCount(ind),'-o')
% end
% legend(planes,'Location','best')
% xlabel('$F$ sfc','Interpreter','latex')
% ylabel('# rejected segments')
% print(fig,[plotpath_sum,filesep,'rej_tot'],pfrm,pres)



%% Bulk ratios for each fitting range


for i_p = 1:Npl 
    plane = planes{i_p};
    
    if ismember(plane,{'ATR-EUREC4A','TO-POST'})
        dirvar = 'dir2';
    else
        dirvar = '';
    end
    
    
    for i_f = 1:Nfc
        MOM = MOM_matrix{i_p,i_f};
        testname = char(MOM.testname(1));
        
        plotpath_test = [plotpath,filesep,replace(testname,{'.',' ',',','/'},{'p','_','',''})];
        if ~isfolder(plotpath_test), mkdir(plotpath_test), end
        
        levels = sortrows(groupsummary(MOM,{'level'},{'mean'},{'alt'}),'mean_alt','descend').level';
        for i_l = 1:numel(levels)
            MOM.level_id(MOM.level==levels(i_l)) = i_l;
        end
               

        % (Pv/Pu,Dv/Du)
        
        [fig,ax] = plot_xy_uni(MOM,{'ar_psd_VU'},{'ar_sfc_VU'},'level_id','',dirvar,true,...
            {'ver4/3','hor4/3','cross1'},[],'XLim',lim_bulk.ar,'YLim',lim_bulk.ar);
        plot(ax,4/3,4/3,'d','Color',"#77AC30",'MarkerFaceColor',"#77AC30",'MarkerSize',12)
        legend(horzcat(levels,{'HIT'}),'Location','northwest','Interpreter','latex')
        xlabel(sprintf('$P_v/P_u$ in range [%.0f$\\Delta r$, %.1f$L$]',psd_bot_factors(i_p),psd_upp_factors(i_f)),'Interpreter','latex')
        ylabel(sprintf('$D_v/D_u$ in range [%.0f$\\Delta r$, %.1f$L$]',sfc_bot_factors(i_p),sfc_upp_factors(i_f)),'Interpreter','latex')
        title(plane)
        print(fig,[plotpath_test,filesep,plane,'_ar_uv'],pfrm,pres)

        
        % (Pw/Pu,Dw/Du)
        
        [fig,ax] = plot_xy_uni(MOM,{'ar_psd_WU'},{'ar_sfc_WU'},'level_id','',dirvar,true,...
            {'ver4/3','hor4/3','cross1'},[],'XLim',lim_bulk.ar,'YLim',lim_bulk.ar);
        plot(ax,4/3,4/3,'d','Color',"#77AC30",'MarkerFaceColor',"#77AC30",'MarkerSize',12)
        legend(horzcat(levels,{'HIT'}),'Location','southeast','Interpreter','latex')
        xlabel(sprintf('$P_w/P_u$ in range [%.0f$\\Delta r$, %.1f$L$]',psd_bot_factors(i_p),psd_upp_factors(i_f)),'Interpreter','latex')
        ylabel(sprintf('$D_w/D_u$ in range [%.0f$\\Delta r$, %.1f$L$]',sfc_bot_factors(i_p),sfc_upp_factors(i_f)),'Interpreter','latex')
        title(plane)
        print(fig,[plotpath_test,filesep,plane,'_ar_uw'],pfrm,pres)

        
        % (p,s)
        
        [fig,ax] = plot_xy_uni(MOM,{'slp_psd_UX','slp_psd_VY','slp_psd_W'},...
            {'slp_sfc_UX','slp_sfc_VY','slp_sfc_W'},'yN','level_id',dirvar,false,...
            {'hor2/3','ver5/3'},[],'XLim',lim_bulk.p,'YLim',lim_bulk.s);
        plot(ax,5/3,2/3,'d','Color',"#77AC30",'MarkerFaceColor',"#77AC30",'MarkerSize',12)
        legend(horzcat({'u','v','w'},levels,{'K41'}),'Location','northwest','Interpreter','latex')
        xlabel(sprintf('$p$ in range [%.0f$\\Delta r$, %.1f$L$]',psd_bot_factors(i_p),psd_upp_factors(i_f)),'Interpreter','latex')
        ylabel(sprintf('$s$ in range [%.0f$\\Delta r$, %.1f$L$]',sfc_bot_factors(i_p),sfc_upp_factors(i_f)),'Interpreter','latex')
        title(plane)
        print(fig,[plotpath_test,filesep,plane,'_slp'],pfrm,pres)
        
      
        clear MOM levels     
    end
    
    close all
    
end


% Printout a piece of a jam script

% jamplots = {'ar_uv','ar_uw','slp'};
% for i_j = 1:numel(jamplots)
%     fprintf('pdfjam --nup %dx%d --papersize ''{%dcm,%dcm}''',Npl,Nfc,Npl*16,Nfc*12)
%     for i_f = 1:Nfc
%         testname = sprintf('sfc %.1fL psd %.1fL',sfc_upp_factors(i_f),psd_upp_factors(i_f));
%         fprintf(' ./%s/*%s.png',replace(testname,{'.',' ',',','/'},{'p','_','',''}),jamplots{i_j})
%     end
%     fprintf(' --outfile ./jam/%s.pdf\n',jamplots{i_j})
% end
% fprintf('pdfjam --papersize ''{%dcm,%dcm}''',Npl*16,Nfc*12)
% fprintf(' ./jam/%s.pdf',jamplots{:})
% fprintf(' --outfile ./jam/scatter.pdf\n')
