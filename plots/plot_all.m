

%% Overview of the segments

plot_seg_overview(MOM,levels,ifdirs,false);
title(plane)
print(gcf,[plotpath,'_overview'],'-dpng','-r300')


%% D examples

i_s = find( MOM.flight==ex_s(1) & MOM.name==ex_s(2) );
dr = MOM.dr(i_s);

for i_v = 1:Nvar
    var = vars{i_v};
    
    edr_sfc( detrend(TURB(i_s).(var)), dr,sfc_fit_range,B(i_v),'Method',sfc_method,...
        'FitPoints',sfc_fit_points,'Plot',true,'PlotXLim',[dr 400],'PlotYLim',[1e-2 0.3]);
    
    ylabel(['$D_',lower(var(1)),'\,[\mathrm{m^2\,s^{-2}}]$'],'Interpreter','latex')
    title(join([plane,MOM.flight(i_s),MOM.name(i_s),round(MOM.alt(i_s)),'m']))
    print(gcf,join([plotpath,'ex','sfc',var,string(i_s)],'_'),'-dpng','-r300')
end


%% P examples

for i_v = 1:Nvar
    var = vars{i_v};
    
    edr_psd( detrend(TURB(i_s).(var)), dr,psd_fit_range,C(i_v),'Method',psd_method,...
        'FitPoints',psd_fit_points,'Plot',true,'PlotXLim',[2*dr 400],'PlotYLim',[1e-4 1],...
        'WindowLength',floor(psd_win_length/dr),'WindowOverlap',floor(psd_win_overlap/dr) );
    
    ylabel(['$P_',lower(var(1)),'\,[\mathrm{m^2\,s^{-2}\,rad^{-1}}]$'],'Interpreter','latex')
    title(join([plane,MOM.flight(i_s),MOM.name(i_s),round(MOM.alt(i_s)),'m']))
    print(gcf,join([plotpath,'ex','psd',var,string(i_s)],'_'),'-dpng','-r300')
end


%% (Pv/Pu,Dv/Du) scatter

[fig,ax] = plot_xy(MOM,{'ar_psd_VU'},{'ar_sfc_VU'},levels,'levels',ifdirs,1,...
    {'ver3/4','ver4/3','hor3/4','hor4/3'},'XLim',[0 1.5],'YLim',[0 1.5]);
plot(ax,4/3,4/3,'d','Color',"#77AC30",'MarkerFaceColor',"#77AC30",'MarkerSize',12)
legend(cat(2,levels,{'K41'}),'Location','northwest','Interpreter','latex')
xlabel('$P_v/P_u$','Interpreter','latex')
ylabel('$D_v/D_u$','Interpreter','latex')
title(plane)
print(fig,[plotpath,'_ar_uv'],'-dpng','-r300')


%% (Pw/Pu,Dw/Du) scatter

[fig,ax] = plot_xy(MOM,{'ar_psd_WU'},{'ar_sfc_WU'},levels,'levels',ifdirs,1,...
    {'ver3/4','ver4/3','hor3/4','hor4/3'},'XLim',[0 1.5],'YLim',[0 1.5]);
plot(ax,4/3,4/3,'d','Color',"#77AC30",'MarkerFaceColor',"#77AC30",'MarkerSize',12)
legend(cat(2,levels,{'K41'}),'Location','northwest','Interpreter','latex')
xlabel('$P_w/P_u$','Interpreter','latex')
ylabel('$D_w/D_u$','Interpreter','latex')
title(plane)
print(fig,[plotpath,'_ar_uw'],'-dpng','-r300')


%% (p,s) scatter

[fig,ax] = plot_xy(MOM,{'slp_psd_UX','slp_psd_VY','slp_psd_W'},...
    {'slp_sfc_UX','slp_sfc_VY','slp_sfc_W'},levels,'vars',ifdirs,0,{'hor2/3','ver5/3'},...
    'XLim',[0 2.5],'YLim',[0 1.5]);
plot(ax,5/3,2/3,'d','Color',"#77AC30",'MarkerFaceColor',"#77AC30",'MarkerSize',12)
legend(cat(2,{'u','v','w'},levels,{'K41'}),'Location','northwest','Interpreter','latex')
xlabel('$p$','Interpreter','latex')
ylabel('$s$','Interpreter','latex')
title(plane)
print(fig,[plotpath,'_slp'],'-dpng','-r300')


%% (Pv/Pu,Dv/Du,Pw/Pu,Dw/Du) whisker

h = plot_whisker(MOM,{'ar_sfc_VU','ar_psd_VU','ar_sfc_WU','ar_psd_WU'},...
    levels,0,'PrimaryLabels',{'$D_v/D_u$','$P_v/P_u$','$D_w/D_u$','$P_w/P_u$'},'DataLim',[0 1.5]);
hold on
h.axis.YLim = [0 1.5];
xlim = h.axis.XLim;
plot(h.axis,xlim,4/3*[1 1],'Color','black','LineWidth',1.5,'LineStyle','--')
plot(h.axis,xlim,3/4*[1 1],'Color','black','LineWidth',1.5,'LineStyle','--')
h.axis.XLim = xlim;
title(plane)
print(h.figure,[plotpath,'_ar_wsk'],'-dpng','-r300')


%% (p,s) whisker

h = plot_whisker(MOM,{'slp_sfc_UX','slp_sfc_VY','slp_sfc_W'},...
    levels,0,'PrimaryLabels',{'u','v','w'},'DataLim',[0 1.5]);
hold on
h.axis.YLim = [0 1.5];
xlim = h.axis.XLim;
plot(h.axis,xlim,2/3*[1 1],'Color','black','LineWidth',1.5,'LineStyle','--')
h.axis.XLim = xlim;
ylabel('$s$','Interpreter','latex')
title(plane)
print(h.figure,[plotpath,'_slp_sfc_wsk'],'-dpng','-r300')

h = plot_whisker(MOM,{'slp_psd_UX','slp_psd_VY','slp_psd_W'},...
    levels,0,'PrimaryLabels',{'u','v','w'},'DataLim',[0 2.5]);
hold on
h.axis.YLim = [0 2.5];
xlim = h.axis.XLim;
plot(h.axis,xlim,5/3*[1 1],'Color','black','LineWidth',1.5,'LineStyle','--')
h.axis.XLim = xlim;
ylabel('$p$','Interpreter','latex')
title(plane)
print(h.figure,[plotpath,'_slp_psd_wsk'],'-dpng','-r300')


%% Uncertainties whisker

h = plot_whisker(MOM,{'e_ar_sfc_VU','e_ar_psd_VU','e_ar_sfc_WU','e_ar_psd_WU'},...
    levels,0,'PrimaryLabels',{'$D_v/D_u$','$P_v/P_u$','$D_w/D_u$','$P_w/P_u$'},'DataLim',[0 0.2]);
hold on
h.axis.YLim = [0 0.2];
ylabel('Uncertainty','Interpreter','latex')
title(plane)
print(h.figure,[plotpath,'_e_wsk_ar'],'-dpng','-r300')

h = plot_whisker(MOM,{'e_slp_sfc_UX','e_slp_sfc_VY','e_slp_sfc_W'},...
    levels,0,'PrimaryLabels',{'$s_u$','$s_v$','$s_w$'},'DataLim',[0 0.1]);
hold on
h.axis.YLim = [0 0.1];
ylabel('Uncertainty','Interpreter','latex')
title(plane)
print(h.figure,[plotpath,'_e_wsk_slp_sfc'],'-dpng','-r300')

h = plot_whisker(MOM,{'e_slp_psd_UX','e_slp_psd_VY','e_slp_psd_W'},...
    levels,0,'PrimaryLabels',{'$p_u$','$p_v$','$p_w$'},'DataLim',[0 0.3]);
hold on
h.axis.YLim = [0 0.3];
ylabel('Uncertainty','Interpreter','latex')
title(plane)
print(h.figure,[plotpath,'_e_wsk_slp_psd'],'-dpng','-r300')
















%% (u,v) offset

% [fig,ax] = plot_xy(MOM,{'off_sfc_UX'},{'off_sfc_VY'},levels,'levels',ifdirs,1,{'cross3/4','cross1','cross4/3'});
% legend(cat(2,levels),'Location','northwest')
% xlabel('$B_L\epsilon_u^{2/3}\,\textrm{sfc}$','Interpreter','latex')
% ylabel('$B_T\epsilon_v^{2/3}\,\textrm{sfc}$','Interpreter','latex')
% title(plane)
% print(fig,[plotpath,'_off_uv_sfc'],'-dpng','-r300')
% % ax.XLim = ax.XLim(1) + [0 0.5*diff(ax.XLim)];
% % ax.YLim = ax.YLim(1) + [0 0.5*diff(ax.YLim)];
% % print(fig,[plotpath,'_uv23_sfc_zoom'],'-dpng','-r300')
% 
% [fig,ax] = plot_xy(MOM,{'off_psd_UX'},{'off_psd_VY'},levels,'levels',ifdirs,1,{'cross3/4','cross1','cross4/3'});
% legend(cat(2,levels),'Location','northwest')
% xlabel('$C_L\epsilon_u^{2/3}\,\textrm{psd}$','Interpreter','latex')
% ylabel('$C_T\epsilon_v^{2/3}\,\textrm{psd}$','Interpreter','latex')
% title(plane)
% print(fig,[plotpath,'_off_uv_psd'],'-dpng','-r300')
% % ax.XLim = ax.XLim(1) + [0 0.5*diff(ax.XLim)];
% % ax.YLim = ax.YLim(1) + [0 0.5*diff(ax.YLim)];
% % print(fig,[plotpath,'_uv23_psd_zoom'],'-dpng','-r300')


%% (u,w) offset

% [fig,ax] = plot_xy(MOM,{'off_sfc_UX'},{'off_sfc_W'},levels,'levels',ifdirs,1,{'cross3/4','cross1','cross4/3'});
% legend(cat(2,levels),'Location','northwest')
% xlabel('$B_L\epsilon_u^{2/3}\,\textrm{sfc}$','Interpreter','latex')
% ylabel('$B_T\epsilon_w^{2/3}\,\textrm{sfc}$','Interpreter','latex')
% title(plane)
% print(fig,[plotpath,'_off_uw_sfc'],'-dpng','-r300')
% % ax.XLim = ax.XLim(1) + [0 0.5*diff(ax.XLim)];
% % ax.YLim = ax.YLim(1) + [0 0.5*diff(ax.YLim)];
% % print(fig,[plotpath,'_uw23_sfc_zoom'],'-dpng','-r300')
% 
% [fig,ax] = plot_xy(MOM,{'off_psd_UX'},{'off_psd_W'},levels,'levels',ifdirs,1,{'cross3/4','cross1','cross4/3'});
% legend(cat(2,levels),'Location','northwest')
% xlabel('$C_L\epsilon_u^{2/3}\,\textrm{psd}$','Interpreter','latex')
% ylabel('$C_T\epsilon_w^{2/3}\,\textrm{psd}$','Interpreter','latex')
% title(plane)
% print(fig,[plotpath,'_off_uw_psd'],'-dpng','-r300')
% % ax.XLim = ax.XLim(1) + [0 0.5*diff(ax.XLim)];
% % ax.YLim = ax.YLim(1) + [0 0.5*diff(ax.YLim)];
% % print(fig,[plotpath,'_uw23_psd_zoom'],'-dpng','-r300')


%% (w/u,w/v) offset

% [fig,ax] = plot_xy(MOM,{'ar_sfc_WU'},{'ar_sfc_WV'},levels,'levels',ifdirs,1,...
%     {'cross3/4','cross4/3','ver3/4','ver4/3','hor3/4','hor4/3'});
% legend(cat(2,levels),'Location','southeast')
% xlabel('$B_T\epsilon_w^{2/3}/(B_L\epsilon_u^{2/3})\,\textrm{sfc}$','Interpreter','latex')
% ylabel('$B_T\epsilon_w^{2/3}/(B_T\epsilon_v^{2/3})\,\textrm{sfc}$','Interpreter','latex')
% title(plane)
% print(fig,[plotpath,'_ar_sfc'],'-dpng','-r300')
% % ax.XLim = [0.5 2]; ax.YLim = [0.5 2];
% % print(fig,[plotpath,'_ar23_sfc_zoom'],'-dpng','-r300')
% 
% [fig,ax] = plot_xy(MOM,{'ar_psd_WU'},{'ar_psd_WV'},levels,'levels',ifdirs,1,...
%     {'cross3/4','cross4/3','ver3/4','ver4/3','hor3/4','hor4/3'});
% legend(cat(2,levels),'Location','southeast')
% xlabel('$C_T\epsilon_w^{2/3}/(C_L\epsilon_u^{2/3})\,\textrm{psd}$','Interpreter','latex')
% ylabel('$C_T\epsilon_w^{2/3}/(C_T\epsilon_v^{2/3})\,\textrm{psd}$','Interpreter','latex')
% title(plane)
% print(fig,[plotpath,'_ar_psd'],'-dpng','-r300')
% % ax.XLim = [0.4 1.6]; ax.YLim = [0.4 1.6];
% % print(fig,[plotpath,'_ar23_psd_zoom'],'-dpng','-r300')





%% Dissipation rate

% [fig,ax] = plot_xy(MOM,{'edr_sfc_UX'},{'edr_sfc_VY'},levels,'levels',ifdirs,1,...
%     {'cross(3/4)^3','cross1','cross(4/3)^3'});
% legend(cat(2,levels),'Location','best')
% xlabel('$\epsilon_u\,[\mathrm{m^2s^{-3}}]\,\textrm{sfc}$','Interpreter','latex')
% ylabel('$\epsilon_v\,[\mathrm{m^2s^{-3}}]\,\textrm{sfc}$','Interpreter','latex')
% title(plane)
% print(fig,[plotpath,'_edr_uv_sfc'],'-dpng','-r300')
% 
% [fig,ax] = plot_xy(MOM,{'edr_psd_UX'},{'edr_psd_VY'},levels,'levels',ifdirs,1,...
%     {'cross(3/4)^3','cross1','cross(4/3)^3'});
% legend(cat(2,levels),'Location','best')
% xlabel('$\epsilon_u\,[\mathrm{m^2s^{-3}}]\,\textrm{psd}$','Interpreter','latex')
% ylabel('$\epsilon_v\,[\mathrm{m^2s^{-3}}]\,\textrm{psd}$','Interpreter','latex')
% title(plane)
% print(fig,[plotpath,'_edr_uv_psd'],'-dpng','-r300')
% 
% [fig,ax] = plot_xy(MOM,{'edr_sfc_UX'},{'edr_sfc_W'},levels,'levels',ifdirs,1,...
%     {'cross(3/4)^3','cross1','cross(4/3)^3'});
% legend(cat(2,levels),'Location','best')
% xlabel('$\epsilon_u\,[\mathrm{m^2s^{-3}}]\,\textrm{sfc}$','Interpreter','latex')
% ylabel('$\epsilon_w\,[\mathrm{m^2s^{-3}}]\,\textrm{sfc}$','Interpreter','latex')
% title(plane)
% print(fig,[plotpath,'_edr_uw_sfc'],'-dpng','-r300')
% 
% [fig,ax] = plot_xy(MOM,{'edr_psd_UX'},{'edr_psd_W'},levels,'levels',ifdirs,1,...
%     {'cross(3/4)^3','cross1','cross(4/3)^3'});
% legend(cat(2,levels),'Location','best')
% xlabel('$\epsilon_u\,[\mathrm{m^2s^{-3}}]\,\textrm{psd}$','Interpreter','latex')
% ylabel('$\epsilon_w\,[\mathrm{m^2s^{-3}}]\,\textrm{psd}$','Interpreter','latex')
% title(plane)
% print(fig,[plotpath,'_edr_uw_psd'],'-dpng','-r300')



%% integral length scale stats

% h = plot_whisker(MOM,{'ls_W'},levels,0,'PrimaryLabels',{''});
% ylabel('$L\,[\mathrm{m}]$','Interpreter','latex')
% title(plane)
% print(h.figure,[plotpath,'_ls_wsk'],'-dpng','-r300')


%% Length

% MOM.length_km = MOM.length/1000;
% 
% h = plot_whisker(MOM,{'length_km'},levels,0,'PrimaryLabels',{''});
% ylabel('Segment length $[\mathrm{km}]$','Interpreter','latex')
% title(plane)
% print(h.figure,[plotpath,'_length_wsk'],'-dpng','-r300')


%% Altitude

% MOM.vstop = MOM.top - MOM.alt;
% 
% h = plot_whisker(MOM,{'alt','vstop'},levels,0,'PrimaryLabels',{'alt','vstop'});
% ylabel('Altitude $[\mathrm{m}]$','Interpreter','latex')
% hold on
% xlim = h.axis.XLim;
% plot(h.axis,xlim,prctile(MOM.top,50)*[1 1],'Color','red','LineWidth',1)
% plot(h.axis,xlim,prctile(MOM.top,25)*[1 1],'Color','red','LineWidth',1,'LineStyle','--')
% plot(h.axis,xlim,prctile(MOM.top,75)*[1 1],'Color','red','LineWidth',1,'LineStyle','--')
% % plot(h.axis,xlim,min(MOM.top)*[1 1],'Color','red','LineWidth',1,'LineStyle',':')
% % plot(h.axis,xlim,max(MOM.top)*[1 1],'Color','red','LineWidth',1,'LineStyle',':')
% h.axis.XLim = xlim;
% h.axis.View = [0 90];
% title(plane)
% print(h.figure,[plotpath,'_alt_wsk'],'-dpng','-r300')


%% Uncertainties scatter

% [fig,ax] = plot_xy(MOM,{'off_sfc_UX','off_sfc_VY','off_sfc_W'},{'e_off_sfc_UX','e_off_sfc_VY','e_off_sfc_W'},levels,'vars',ifdirs,0);
% legend(cat(2,{'ux','vy','w'},levels),'Location','best')
% xlabel('$O\,\textrm{sfc}$','Interpreter','latex')
% ylabel('$\textrm{Uncertainty}$','Interpreter','latex')
% title(plane)
% print(fig,[plotpath,'_e_off_sfc_xy'],'-dpng','-r300')
% 
% [fig,ax] = plot_xy(MOM,{'off_psd_UX','off_psd_VY','off_psd_W'},{'e_off_psd_UX','e_off_psd_VY','e_off_psd_W'},levels,'vars',ifdirs,0);
% legend(cat(2,{'ux','vy','w'},levels),'Location','best')
% xlabel('$O\,\textrm{psd}$','Interpreter','latex')
% ylabel('$\textrm{Uncertainty}$','Interpreter','latex')
% title(plane)
% print(fig,[plotpath,'_e_off_psd_xy'],'-dpng','-r300')
% 
% 
% [fig,ax] = plot_xy(MOM,{'slp_sfc_UX','slp_sfc_VY','slp_sfc_W'},{'e_slp_sfc_UX','e_slp_sfc_VY','e_slp_sfc_W'},levels,'vars',ifdirs,0);
% legend(cat(2,{'ux','vy','w'},levels),'Location','best')
% xlabel('$s\,\textrm{sfc}$','Interpreter','latex')
% ylabel('$\textrm{Uncertainty}$','Interpreter','latex')
% title(plane)
% print(fig,[plotpath,'_e_slp_sfc_xy'],'-dpng','-r300')
% 
% [fig,ax] = plot_xy(MOM,{'slp_psd_UX','slp_psd_VY','slp_psd_W'},{'e_slp_psd_UX','e_slp_psd_VY','e_slp_psd_W'},levels,'vars',ifdirs,0);
% legend(cat(2,{'ux','vy','w'},levels),'Location','best')
% xlabel('$s\,\textrm{psd}$','Interpreter','latex')
% ylabel('$\textrm{Uncertainty}$','Interpreter','latex')
% title(plane)
% print(fig,[plotpath,'_e_slp_psd_xy'],'-dpng','-r300')
% 
% 
% [fig,ax] = plot_xy(MOM,{'ar_sfc_VU','ar_sfc_WU','ar_psd_VU','ar_psd_WU'},{'e_ar_sfc_VU','e_ar_sfc_WU','e_ar_psd_VU','e_ar_psd_WU'},levels,'vars',ifdirs,0);
% legend(cat(2,{'Dv/Du','Dw/Du','Pv/Pu','Pw/Pu'},levels),'Location','best')
% xlabel('Value')
% ylabel('Uncertainty')
% title(plane)
% print(fig,[plotpath,'_e_ar_xy'],'-dpng','-r300')

% [fig,ax] = plot_xy(MOM,{'ar_sfc_VU','ar_sfc_WU'},{'e_ar_sfc_VU','e_ar_sfc_WU'},levels,'vars',ifdirs,0);
% legend(cat(2,{'v/u','w/u'},levels),'Location','best')
% xlabel('$D_v/D_u,\,D_w/D_u$','Interpreter','latex')
% ylabel('$\textrm{Uncertainty}$','Interpreter','latex')
% title(plane)
% print(fig,[plotpath,'_e_xy_ar_sfc'],'-dpng','-r300')
% 
% [fig,ax] = plot_xy(MOM,{'ar_psd_VU','ar_psd_WU'},{'e_ar_psd_VU','e_ar_psd_WU'},levels,'vars',ifdirs,0);
% legend(cat(2,{'v/u','w/u'},levels),'Location','best')
% xlabel('$P_v/P_u,\,P_w/P_u$','Interpreter','latex')
% ylabel('$\textrm{Uncertainty}$','Interpreter','latex')
% title(plane)
% print(fig,[plotpath,'_e_xy_ar_psd'],'-dpng','-r300')




%% Uncertainties whisker


% h = plot_whisker(MOM,{'er_off_sfc_UX','er_off_sfc_VY','er_off_sfc_W','er_off_psd_UX','er_off_psd_VY','er_off_psd_W'},...
%     levels,0,'PrimaryLabels',{'ux sfc','vy sfc','w sfc','ux psd','vy psd','w psd'});
% ylabel('$\Delta O/O$','Interpreter','latex')
% title(plane)
% print(h.figure,[plotpath,'_e_off_wsk'],'-dpng','-r300')

% h = plot_whisker(MOM,{'er_off_sfc_UX','er_off_sfc_VY','er_off_sfc_W'},levels,0,'PrimaryLabels',{'ux','vy','w'});
% ylabel('$\Delta O/O\,\textrm{sfc}$','Interpreter','latex')
% title(plane)
% print(h.figure,[plotpath,'_e_wsk_off_sfc'],'-dpng','-r300')
% 
% h = plot_whisker(MOM,{'er_off_psd_UX','er_off_psd_VY','er_off_psd_W'},levels,0,'PrimaryLabels',{'ux','vy','w'});
% ylabel('$\Delta O/O\,\textrm{psd}$','Interpreter','latex')
% title(plane)
% print(h.figure,[plotpath,'_e_wsk_off_psd'],'-dpng','-r300')


% h = plot_whisker(MOM,{'e_slp_sfc_UX','e_slp_sfc_VY','e_slp_sfc_W','e_slp_psd_UX','e_slp_psd_VY','e_slp_psd_W'},...
%     levels,0,'PrimaryLabels',{'ux sfc','vy sfc','w sfc','ux psd','vy psd','w psd'});
% ylabel('$\Delta s$','Interpreter','latex')
% title(plane)
% print(h.figure,[plotpath,'_e_slp_wsk'],'-dpng','-r300')

% h = plot_whisker(MOM,{'e_slp_sfc_UX','e_slp_sfc_VY','e_slp_sfc_W'},levels,0,'PrimaryLabels',{'ux','vy','w'});
% ylabel('$\Delta s\,\textrm{sfc}$','Interpreter','latex')
% title(plane)
% print(h.figure,[plotpath,'_e_wsk_slp_sfc'],'-dpng','-r300')
% 
% h = plot_whisker(MOM,{'e_slp_psd_UX','e_slp_psd_VY','e_slp_psd_W'},levels,0,'PrimaryLabels',{'ux','vy','w'});
% ylabel('$\Delta s\,\textrm{psd}$','Interpreter','latex')
% title(plane)
% print(h.figure,[plotpath,'_e_wsk_slp_psd'],'-dpng','-r300')


% h = plot_whisker(MOM,{'e_ar_sfc_VU','e_ar_sfc_WU','e_ar_psd_VU','e_ar_psd_WU'},...
%     levels,0,'PrimaryLabels',{'Dv/Du','Dw/Du','Pv/Pu','Pw/Pu'});
% ylabel('Uncertainty')
% title(plane)
% print(h.figure,[plotpath,'_e_ar_wsk'],'-dpng','-r300')