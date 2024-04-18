
if ~exist('ex','var')
    ex = '';
end
    


%% Overview of the segments

% plot_seg_overview(MOM,levels,ifdirs,false);
% title(plane)
% print(gcf,[plotpath,'_overview'],'-dpng','-r300')


%% Examples: structure functions

% i_s = find( MOM.flight==ex_s(1) & MOM.name==ex_s(2) );
% dr = MOM.dr(i_s);
% 
% for i_v = 1:Nvar
%     var = vars{i_v};
%     
%     fit_sfc( detrend(TURB(i_s).(var)), dr, sfc_fit_range, ...
%         'Method',sfc_method, 'FitPoints',sfc_fit_points, ...
%         'Plot',true, 'PlotXLim',[dr 400], 'PlotYLim',[1e-2 0.3]);
%     
%     ylabel(['$D_',lower(var(1)),'\,[\mathrm{m^2\,s^{-2}}]$'],'Interpreter','latex')
%     title(join([plane,MOM.flight(i_s),MOM.name(i_s),round(MOM.alt(i_s)),'m',ex]))
%     print(gcf,join([plotpath,'ex','sfc',var,string(i_s)],'_'),'-dpng','-r300')
% end


%% Examples: power spectra

% for i_v = 1:Nvar
%     var = vars{i_v};
%     
%     fit_psd( detrend(TURB(i_s).(var)), dr, psd_fit_range, ...
%         'Method',psd_method, 'FitPoints',psd_fit_points, ...
%         'WindowLength',floor(psd_win_length/dr), 'WindowOverlap',floor(psd_win_overlap/dr), ...
%         'Plot',true, 'PlotXLim',[2*dr 400], 'PlotYLim',[1e-4 1]);
%     
%     ylabel(['$P_',lower(var(1)),'\,[\mathrm{m^2\,s^{-2}\,rad^{-1}}]$'],'Interpreter','latex')
%     title(join([plane,MOM.flight(i_s),MOM.name(i_s),round(MOM.alt(i_s)),'m',ex]))
%     print(gcf,join([plotpath,'ex','psd',var,string(i_s)],'_'),'-dpng','-r300')
% end


%

[fig,ax] = plot_xy(MOM,{'ar_psd_VU'},{'ar_sfc_VU'},levels,'levels',ifdirs,1,...
    {'ver3/4','ver4/3','hor3/4','hor4/3','cross1'},'XLim',[0 1.5],'YLim',[0 1.5]);
plot(ax,4/3,4/3,'d','Color',"#77AC30",'MarkerFaceColor',"#77AC30",'MarkerSize',12)
legend(cat(2,levels,{'HIT'}),'Location','northwest','Interpreter','latex')
xlabel('$P_v/P_u$','Interpreter','latex')
ylabel('$D_v/D_u$','Interpreter','latex')
title([plane,' ',ex])
print(fig,[plotpath,'_ar_uv'],'-dpng','-r300')


%% (Pw/Pu,Dw/Du) scatter

[fig,ax] = plot_xy(MOM,{'ar_psd_WU'},{'ar_sfc_WU'},levels,'levels',ifdirs,1,...
    {'ver3/4','ver4/3','hor3/4','hor4/3','cross1'},'XLim',[0 1.5],'YLim',[0 1.5]);
plot(ax,4/3,4/3,'d','Color',"#77AC30",'MarkerFaceColor',"#77AC30",'MarkerSize',12)
legend(cat(2,levels,{'HIT'}),'Location','northwest','Interpreter','latex')
xlabel('$P_w/P_u$','Interpreter','latex')
ylabel('$D_w/D_u$','Interpreter','latex')
title([plane,' ',ex])
print(fig,[plotpath,'_ar_uw'],'-dpng','-r300')


%% (p,s) scatter

[fig,ax] = plot_xy(MOM,{'slp_psd_UX','slp_psd_VY','slp_psd_W'},...
    {'slp_sfc_UX','slp_sfc_VY','slp_sfc_W'},levels,'vars',ifdirs,0,{'hor2/3','ver5/3'},...
    'XLim',[0 2.5],'YLim',[0 1.5]);
plot(ax,5/3,2/3,'d','Color',"#77AC30",'MarkerFaceColor',"#77AC30",'MarkerSize',12)
legend(cat(2,{'u','v','w'},levels,{'K41'}),'Location','northwest','Interpreter','latex')
xlabel('$p$','Interpreter','latex')
ylabel('$s$','Interpreter','latex')
title([plane,' ',ex])
print(fig,[plotpath,'_slp'],'-dpng','-r300')


%% (Pv/Pu,Dv/Du,Pw/Pu,Dw/Du) whisker

% h = plot_whisker(MOM,{'ar_sfc_VU','ar_psd_VU','ar_sfc_WU','ar_psd_WU'},...
%     levels,0,'PrimaryLabels',{'$D_v/D_u$','$P_v/P_u$','$D_w/D_u$','$P_w/P_u$'},'DataLim',[0 1.5]);
% hold on
% h.axis.YLim = [0 1.5];
% xlim = h.axis.XLim;
% plot(h.axis,xlim,4/3*[1 1],'Color','black','LineWidth',1.5,'LineStyle','--')
% plot(h.axis,xlim,3/4*[1 1],'Color','black','LineWidth',1.5,'LineStyle','--')
% h.axis.XLim = xlim;
% title([plane,' ',ex])
% print(h.figure,[plotpath,'_ar_wsk'],'-dpng','-r300')


%% (p,s) whisker

% h = plot_whisker(MOM,{'slp_sfc_UX','slp_sfc_VY','slp_sfc_W'},...
%     levels,0,'PrimaryLabels',{'u','v','w'},'DataLim',[0 1.5]);
% hold on
% h.axis.YLim = [0 1.5];
% xlim = h.axis.XLim;
% plot(h.axis,xlim,2/3*[1 1],'Color','black','LineWidth',1.5,'LineStyle','--')
% h.axis.XLim = xlim;
% ylabel('$s$','Interpreter','latex')
% title([plane,' ',ex])
% print(h.figure,[plotpath,'_slp_sfc_wsk'],'-dpng','-r300')
% 
% h = plot_whisker(MOM,{'slp_psd_UX','slp_psd_VY','slp_psd_W'},...
%     levels,0,'PrimaryLabels',{'u','v','w'},'DataLim',[0 2.5]);
% hold on
% h.axis.YLim = [0 2.5];
% xlim = h.axis.XLim;
% plot(h.axis,xlim,5/3*[1 1],'Color','black','LineWidth',1.5,'LineStyle','--')
% h.axis.XLim = xlim;
% ylabel('$p$','Interpreter','latex')
% title([plane,' ',ex])
% print(h.figure,[plotpath,'_slp_psd_wsk'],'-dpng','-r300')


%% Uncertainties whisker

% h = plot_whisker(MOM,{'e_ar_sfc_VU','e_ar_psd_VU','e_ar_sfc_WU','e_ar_psd_WU'},...
%     levels,0,'PrimaryLabels',{'$D_v/D_u$','$P_v/P_u$','$D_w/D_u$','$P_w/P_u$'},'DataLim',[0 0.2]);
% hold on
% h.axis.YLim = [0 0.2];
% ylabel('Uncertainty','Interpreter','latex')
% title([plane,' ',ex])
% print(h.figure,[plotpath,'_e_wsk_ar'],'-dpng','-r300')
% 
% h = plot_whisker(MOM,{'e_slp_sfc_UX','e_slp_sfc_VY','e_slp_sfc_W'},...
%     levels,0,'PrimaryLabels',{'$s_u$','$s_v$','$s_w$'},'DataLim',[0 0.1]);
% hold on
% h.axis.YLim = [0 0.1];
% ylabel('Uncertainty','Interpreter','latex')
% title([plane,' ',ex])
% print(h.figure,[plotpath,'_e_wsk_slp_sfc'],'-dpng','-r300')
% 
% h = plot_whisker(MOM,{'e_slp_psd_UX','e_slp_psd_VY','e_slp_psd_W'},...
%     levels,0,'PrimaryLabels',{'$p_u$','$p_v$','$p_w$'},'DataLim',[0 0.3]);
% hold on
% h.axis.YLim = [0 0.3];
% ylabel('Uncertainty','Interpreter','latex')
% title([plane,' ',ex])
% print(h.figure,[plotpath,'_e_wsk_slp_psd'],'-dpng','-r300')


%% Uncertainties scatter

% [fig,ax] = plot_xy(MOM,{'off_sfc_UX','off_sfc_VY','off_sfc_W'},{'e_off_sfc_UX','e_off_sfc_VY','e_off_sfc_W'},levels,'vars',ifdirs,0);
% legend(cat(2,{'ux','vy','w'},levels),'Location','best')
% xlabel('$O\,\textrm{sfc}$','Interpreter','latex')
% ylabel('Uncertainty','Interpreter','latex')
% title(plane)
% print(fig,[plotpath,'_e_off_sfc_xy'],'-dpng','-r300')
% 
% [fig,ax] = plot_xy(MOM,{'off_psd_UX','off_psd_VY','off_psd_W'},{'e_off_psd_UX','e_off_psd_VY','e_off_psd_W'},levels,'vars',ifdirs,0);
% legend(cat(2,{'ux','vy','w'},levels),'Location','best')
% xlabel('$O\,\textrm{psd}$','Interpreter','latex')
% ylabel('Uncertainty','Interpreter','latex')
% title(plane)
% print(fig,[plotpath,'_e_off_psd_xy'],'-dpng','-r300')
% 
% 
% [fig,ax] = plot_xy(MOM,{'slp_sfc_UX','slp_sfc_VY','slp_sfc_W'},{'e_slp_sfc_UX','e_slp_sfc_VY','e_slp_sfc_W'},levels,'vars',ifdirs,0);
% legend(cat(2,{'ux','vy','w'},levels),'Location','best')
% xlabel('$s\,\textrm{sfc}$','Interpreter','latex')
% ylabel('Uncertainty','Interpreter','latex')')
% title(plane)
% print(fig,[plotpath,'_e_slp_sfc_xy'],'-dpng','-r300')
% 
% [fig,ax] = plot_xy(MOM,{'slp_psd_UX','slp_psd_VY','slp_psd_W'},{'e_slp_psd_UX','e_slp_psd_VY','e_slp_psd_W'},levels,'vars',ifdirs,0);
% legend(cat(2,{'ux','vy','w'},levels),'Location','best')
% xlabel('$s\,\textrm{psd}$','Interpreter','latex')
% ylabel('Uncertainty','Interpreter','latex'))
% title(plane)
% print(fig,[plotpath,'_e_slp_psd_xy'],'-dpng','-r300')
% 
% 
% [fig,ax] = plot_xy(MOM,{'ar_sfc_VU','ar_sfc_WU','ar_psd_VU','ar_psd_WU'},{'e_ar_sfc_VU','e_ar_sfc_WU','e_ar_psd_VU','e_ar_psd_WU'},levels,'vars',ifdirs,0);
% legend(cat(2,{'Dv/Du','Dw/Du','Pv/Pu','Pw/Pu'},levels),'Location','best')
% xlabel('Ratio')
% ylabel('Uncertainty','Interpreter','latex')
% title(plane)
% print(fig,[plotpath,'_e_ar_xy'],'-dpng','-r300')
