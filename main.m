

%% Select campaign

% plane = 'ATR-EUREC4A';
% plane = 'C130-VOCALS';
% plane = 'C130-RICO';
% plane = 'TO-POST';



%% Prepare paths

% MYPROJECTPATH is the path where you downloaded the codes
%
%
% MYDATAPATH is the path where you downloaded the datasets:
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
% MYDATAPATH/ATR-EUREC4A/TURBLENCE
%
% Table 3 from Bony et al. 2022: EUREC4A observations from the SAFIRE ATR42 aircraft,
% Earth Syst. Sci. Data, 14, 2021–2064, doi.org/10.5194/essd-14-2021-2022, 2022.
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
% (as tab-delimited text file).


addpath(genpath(myprojectpath))

datapath = [mydatapath,filesep,plane];

plotpath = [myprojectpath,filesep,'figures'];
if ~isfolder(plotpath)
    mkdir(plotpath)
end
plotpath = [plotpath,filesep,plane];



%% Load datasets

if strcmp(plane,'ATR-EUREC4A')
    
    fit_range = [16 80];
    ex_s = ["RF12","R2B"];
    
    [TURB,MOM,levels] = load_eureca_all(datapath);

elseif strcmp(plane,'C130-VOCALS')
    
    fit_range = [16 80];
    ex_s = ["RF09","C6"];
    
    [TURB,MOM,levels] = load_vocals_all(datapath);

elseif strcmp(plane,'C130-RICO')
    
    fit_range = [16 80];
    ex_s = ["RF06","SC01"];
    
    [TURB,MOM,levels] = load_rico_all(datapath);
    
elseif strcmp(plane,'TO-POST')
    
    fit_range = [8 80];
    ex_s = ["RF12","CB01"];
    
    [TURB,MOM,levels] = load_post_all(datapath);

end


% Plot overview of the segments

plot_seg_overview(MOM,levels);
title(plane)
print(gcf,[plotpath,'_overview'],'-dpng','-r300')



%% Calculate dissipation

% Constants

B_L = 2.0; B_T = 2.6;
C_L = 0.5; C_T = 0.66;


% Settings

sfc_method = "logmean";
sfc_fit_points = 6;

psd_method = "logmean";
psd_fit_points = 6;
psd_win_length = 1000; % m
psd_win_overlap = 500; % m

vars = {'W','UX','VY'};
B = [B_T B_L B_T];
C = [C_T C_L C_T];


% Compute

disp('Compute dissipation rate ...')

Nvar = numel(vars);
Nseg = size(MOM,1);
E = struct([]);

for i_v = 1:Nvar
    var = vars{i_v}; fprintf('%2s',var)
    
    for i_s = 1:Nseg
        fprintf(' %d',i_s)
        dr = MOM.dr(i_s);

        [MOM.(['edr_sfc_',var])(i_s),MOM.(['slp_sfc_',var])(i_s),es] = edr_sfc( detrend(TURB(i_s).(var)),...
            dr,fit_range,B(i_v),'Method',sfc_method,'FitPoints',sfc_fit_points );
        
        [MOM.(['edr_psd_',var])(i_s),MOM.(['slp_psd_',var])(i_s),ep] = edr_psd( detrend(TURB(i_s).(var)),...
            dr,fit_range,C(i_v),'Method',psd_method,'FitPoints',psd_fit_points,...
            'WindowLength',floor(psd_win_length/dr),'WindowOverlap',floor(psd_win_overlap/dr) );
        
        E(1).(['sfc_',var])(i_s) = es;
        E(1).(['psd_',var])(i_s) = es;
    end
    
    E.(['sfc_',var]) = struct2table(E.(['sfc_',var]));
    E.(['psd_',var]) = struct2table(E.(['psd_',var]));
    
    MOM.(['slp_psd_',var]) = - MOM.(['slp_psd_',var]);
    
    fprintf('\n')
end



% Dependent parameters

% sfc and psd prefactors
for i_v = 1:Nvar
    var = vars{i_v};
    
    MOM.(['off_sfc_',var]) = B(i_v)*MOM.(['edr_sfc_',var]).^(2/3);
    MOM.(['off_psd_',var]) = C(i_v)*MOM.(['edr_psd_',var]).^(2/3);
    
    E.(['sfc_',var]).off = MOM.(['off_sfc_',var]) .* E.(['sfc_',var]).offsetFixed;
end

% Anisotropy
MOM.ar_sfc_WU = MOM.off_sfc_W./MOM.off_sfc_UX;
MOM.ar_sfc_WV = MOM.off_sfc_W./MOM.off_sfc_VY;
MOM.ar_psd_WU = MOM.off_psd_W./MOM.off_psd_UX;
MOM.ar_psd_WV = MOM.off_psd_W./MOM.off_psd_VY;
MOM.ar_sfc_VU = MOM.off_sfc_VY./MOM.off_sfc_UX;
MOM.ar_psd_VU = MOM.off_psd_VY./MOM.off_psd_UX;

% Dissipation rates after reversal of longi/trans
MOM.edr_sfc_UY = MOM.edr_sfc_UX * (B_L/B_T).^(3/2);
MOM.edr_sfc_VX = MOM.edr_sfc_VY * (B_T/B_L).^(3/2);
MOM.edr_psd_UY = MOM.edr_psd_UX * (C_L/C_T).^(3/2);
MOM.edr_psd_VX = MOM.edr_psd_VY * (C_T/C_L).^(3/2);



%% Integral length scale

disp('Compute integral length scale ...')

var = 'W';
for i_s = 1:Nseg
    fprintf(' %d',i_s)
    MOM.(['ls_',var])(i_s) = int_ls_short(detrend(TURB(i_s).(var)),'Method','e-decay')*MOM.dr(i_s);
end
fprintf('\n')



%% PLOTS

dirs = {'along','cross'};


%% Exclude dirty segments

% ind_exc_psd = find( MOM.slp_psd_W > 7/3 | MOM.slp_psd_UX > 7/3 | MOM.slp_psd_VY > 7/3 | ...
%                     MOM.slp_psd_W < 2/3 | MOM.slp_psd_UX < 2/3 | MOM.slp_psd_VY < 2/3 );
% ind_exc_sfc = find( MOM.slp_sfc_W > 5/3 | MOM.slp_sfc_UX > 5/3 | MOM.slp_sfc_VY > 5/3 | ...
%                     MOM.slp_sfc_W < 0   | MOM.slp_sfc_UX < 0   | MOM.slp_sfc_VY < 0   );
%                   
% MOM ([ind_exc_psd;ind_exc_sfc],:) = [];
% TURB([ind_exc_psd;ind_exc_sfc],:) = [];



%% Examples of dissipation rate derivation

i_s = find( MOM.flight==ex_s(1) & MOM.name==ex_s(2) );
dr = MOM.dr(i_s);

for i_v = 1:Nvar
    var = vars{i_v};
    
    edr_sfc( detrend(TURB(i_s).(var)), dr,fit_range,B(i_v),'Method',sfc_method,...
        'FitPoints',sfc_fit_points,'Plot',true,'PlotXLim',[dr 400],'PlotYLim',[1e-2 0.3]);
    
    title(join([plane,MOM.flight(i_s),MOM.name(i_s),round(MOM.alt(i_s)),'m',var]))
    print(gcf,join([plotpath,'ex','sfc',var,string(i_s)],'_'),'-dpng','-r300')
end

%
for i_v = 1:Nvar
    var = vars{i_v};
    
    edr_psd( detrend(TURB(i_s).(var)), dr,fit_range,C(i_v),'Method',psd_method,...
        'FitPoints',psd_fit_points,'Plot',true,'PlotXLim',[2*dr 400],'PlotYLim',[1e-4 1],...
        'WindowLength',floor(psd_win_length/dr),'WindowOverlap',floor(psd_win_overlap/dr) );
    
    title(join([plane,MOM.flight(i_s),MOM.name(i_s),round(MOM.alt(i_s)),'m',var]))
    print(gcf,join([plotpath,'ex','psd',var,string(i_s)],'_'),'-dpng','-r300')
end


%% (u,v) offset

[fig,ax] = plot_xy(MOM,{'off_sfc_UX'},{'off_sfc_VY'},levels,'levels',1,1,{'cross3/4','cross1','cross4/3'});
legend(cat(2,levels,dirs),'Location','northwest')
xlabel('$B_L\epsilon_u^{2/3}\,\textrm{sfc}$','Interpreter','latex')
ylabel('$B_T\epsilon_v^{2/3}\,\textrm{sfc}$','Interpreter','latex')
title(plane)
print(fig,[plotpath,'_uv_sfc'],'-dpng','-r300')
% ax.XLim = ax.XLim(1) + [0 0.5*diff(ax.XLim)];
% ax.YLim = ax.YLim(1) + [0 0.5*diff(ax.YLim)];
% print(fig,[plotpath,'_uv23_sfc_zoom'],'-dpng','-r300')


[fig,ax] = plot_xy(MOM,{'off_psd_UX'},{'off_psd_VY'},levels,'levels',1,1,{'cross3/4','cross1','cross4/3'});
legend(cat(2,levels,dirs),'Location','northwest')
xlabel('$C_L\epsilon_u^{2/3}\,\textrm{psd}$','Interpreter','latex')
ylabel('$C_T\epsilon_v^{2/3}\,\textrm{psd}$','Interpreter','latex')
title(plane)
print(fig,[plotpath,'_uv_psd'],'-dpng','-r300')
% ax.XLim = ax.XLim(1) + [0 0.5*diff(ax.XLim)];
% ax.YLim = ax.YLim(1) + [0 0.5*diff(ax.YLim)];
% print(fig,[plotpath,'_uv23_psd_zoom'],'-dpng','-r300')


%% (u,w) offset

[fig,ax] = plot_xy(MOM,{'off_sfc_UX'},{'off_sfc_W'},levels,'levels',1,1,{'cross3/4','cross1','cross4/3'});
legend(cat(2,levels,dirs),'Location','northwest')
xlabel('$B_L\epsilon_u^{2/3}\,\textrm{sfc}$','Interpreter','latex')
ylabel('$B_T\epsilon_w^{2/3}\,\textrm{sfc}$','Interpreter','latex')
title(plane)
print(fig,[plotpath,'_uw_sfc'],'-dpng','-r300')
% ax.XLim = ax.XLim(1) + [0 0.5*diff(ax.XLim)];
% ax.YLim = ax.YLim(1) + [0 0.5*diff(ax.YLim)];
% print(fig,[plotpath,'_uw23_sfc_zoom'],'-dpng','-r300')

[fig,ax] = plot_xy(MOM,{'off_psd_UX'},{'off_psd_W'},levels,'levels',1,1,{'cross3/4','cross1','cross4/3'});
legend(cat(2,levels,dirs),'Location','northwest')
xlabel('$C_L\epsilon_u^{2/3}\,\textrm{psd}$','Interpreter','latex')
ylabel('$C_T\epsilon_w^{2/3}\,\textrm{psd}$','Interpreter','latex')
title(plane)
print(fig,[plotpath,'_uw_psd'],'-dpng','-r300')
% ax.XLim = ax.XLim(1) + [0 0.5*diff(ax.XLim)];
% ax.YLim = ax.YLim(1) + [0 0.5*diff(ax.YLim)];
% print(fig,[plotpath,'_uw23_psd_zoom'],'-dpng','-r300')


%% (w/u,w/v) offset

[fig,ax] = plot_xy(MOM,{'ar_sfc_WU'},{'ar_sfc_WV'},levels,'levels',1,1,...
    {'cross3/4','cross4/3','ver3/4','ver4/3','hor3/4','hor4/3'});
legend(cat(2,levels,dirs),'Location','southeast')
xlabel('$B_T\epsilon_w^{2/3}/(B_L\epsilon_u^{2/3})\,\textrm{sfc}$','Interpreter','latex')
ylabel('$B_T\epsilon_w^{2/3}/(B_T\epsilon_v^{2/3})\,\textrm{sfc}$','Interpreter','latex')
title(plane)
print(fig,[plotpath,'_ar_sfc'],'-dpng','-r300')
% ax.XLim = [0.5 2]; ax.YLim = [0.5 2];
% print(fig,[plotpath,'_ar23_sfc_zoom'],'-dpng','-r300')

[fig,ax] = plot_xy(MOM,{'ar_psd_WU'},{'ar_psd_WV'},levels,'levels',1,1,...
    {'cross3/4','cross4/3','ver3/4','ver4/3','hor3/4','hor4/3'});
legend(cat(2,levels,dirs),'Location','southeast')
xlabel('$C_T\epsilon_w^{2/3}/(C_L\epsilon_u^{2/3})\,\textrm{psd}$','Interpreter','latex')
ylabel('$C_T\epsilon_w^{2/3}/(C_T\epsilon_v^{2/3})\,\textrm{psd}$','Interpreter','latex')
title(plane)
print(fig,[plotpath,'_ar_psd'],'-dpng','-r300')
% ax.XLim = [0.4 1.6]; ax.YLim = [0.4 1.6];
% print(fig,[plotpath,'_ar23_psd_zoom'],'-dpng','-r300')


%% (psd,sfc) v/u

[fig,ax] = plot_xy(MOM,{'ar_psd_VU'},{'ar_sfc_VU'},levels,'levels',1,1,...
    {'ver3/4','ver4/3','hor3/4','hor4/3'});
legend(cat(2,levels,dirs),'Location','best')
xlabel('$P_v/P_u$','Interpreter','latex')
ylabel('$D_v/D_u$','Interpreter','latex')
title(plane)
print(fig,[plotpath,'_ar_uv'],'-dpng','-r300')
% ax.XLim = [0.4 1.6]; ax.YLim = [0.4 1.6];
% print(fig,[plotpath,'_ax23_zoom'],'-dpng','-r300')


%% (psd,sfc) w/u

[fig,ax] = plot_xy(MOM,{'ar_psd_WU'},{'ar_sfc_WU'},levels,'levels',1,1,...
    {'ver3/4','ver4/3','hor3/4','hor4/3'});
legend(cat(2,levels,dirs),'Location','best')
xlabel('$P_w/P_u$','Interpreter','latex')
ylabel('$D_w/D_u$','Interpreter','latex')
title(plane)
print(fig,[plotpath,'_ar_uw'],'-dpng','-r300')


%% (s,p)

fig = plot_xy(MOM,{'slp_sfc_W','slp_sfc_UX','slp_sfc_VY'},...
    {'slp_psd_W','slp_psd_UX','slp_psd_VY'},levels,'vars',1,0,{'ver2/3','hor5/3'},...
    'XLim',[0 5/3],'YLim',[2/3 7/3]);
% fig.PaperSize = [20 12]; fig.PaperPosition = [0 0 20 12];
legend(cat(2,vars,levels),'Location','best')
xlabel('$s\,\textrm{sfc}$','Interpreter','latex')
ylabel('$s\,\textrm{psd}$','Interpreter','latex')
title(plane)
print(fig,[plotpath,'_slp'],'-dpng','-r300')


%% edr stats

% along and cross wind

% h = plot_whisker(MOM,{'edr_sfc_W','edr_sfc_UX','edr_sfc_VY'},levels,1,...
%     'PrimaryLabels',{'W a','W c','UX a','UX c','VY a','VY c'});
% ylabel('$\epsilon\,[\mathrm{m^2s^{-3}}]\,\textrm{sfc}$','Interpreter','latex')
% title(plane)
% print(h.figure,[plotpath,filesep,'edr_wsk_',plane,'_sfc_dir'],'-dpng','-r300')
% 
% h = plot_whisker(MOM,{'edr_psd_W','edr_psd_UX','edr_psd_VY'},levels,1,...
%     'PrimaryLabels',{'W a','W c','UX a','UX c','VY a','VY c'});
% ylabel('$\epsilon\,[\mathrm{m^2s^{-3}}]\,\textrm{psd}$','Interpreter','latex')
% title(plane)
% print(h.figure,[plotpath,filesep,'edr_wsk_',plane,'_psd_dir'],'-dpng','-r300')

% before and after reversal

h = plot_whisker(MOM,{'edr_sfc_W','edr_sfc_UX','edr_sfc_VY','edr_sfc_UY','edr_sfc_VX'},...
    levels,0,'PrimaryLabels',{'W','UX','VY','UY','VX'});
ylabel('$\epsilon\,[\mathrm{m^2s^{-3}}]\,\textrm{sfc}$','Interpreter','latex')
title(plane)
print(h.figure,[plotpath,'_edr_wsk_sfc_rev'],'-dpng','-r300')

h = plot_whisker(MOM,{'edr_psd_W','edr_psd_UX','edr_psd_VY','edr_psd_UY','edr_psd_VX'},...
    levels,0,'PrimaryLabels',{'W','UX','VY','UY','VX'});
ylabel('$\epsilon\,[\mathrm{m^2s^{-3}}]\,\textrm{psd}$','Interpreter','latex')
title(plane)
print(h.figure,[plotpath,'_edr_wsk_psd_rev'],'-dpng','-r300')


%% slp stats

h = plot_whisker(MOM,{'slp_sfc_W','slp_sfc_UX','slp_sfc_VY'},...
    levels,0,'PrimaryLabels',{'W','UX','VY'});
hold on
xlim = h.axis.XLim;
plot(h.axis,xlim,2/3*[1 1],'Color','black','LineWidth',1)
h.axis.XLim = xlim;
ylabel('$s\,\textrm{sfc}$','Interpreter','latex')
title(plane)
print(h.figure,[plotpath,'_slp_wsk_sfc'],'-dpng','-r300')

h = plot_whisker(MOM,{'slp_psd_W','slp_psd_UX','slp_psd_VY'},...
    levels,0,'PrimaryLabels',{'W','UX','VY'});
hold on
xlim = h.axis.XLim;
plot(h.axis,xlim,5/3*[1 1],'Color','black','LineWidth',1)
h.axis.XLim = xlim;
ylabel('$s\,\textrm{psd}$','Interpreter','latex')
title(plane)
print(h.figure,[plotpath,'_slp_wsk_psd'],'-dpng','-r300')



%% integral length scale stats

h = plot_whisker(MOM,{'ls_W'},levels,0,'PrimaryLabels',{'W'});
ylabel('$L\,[\mathrm{m}]$','Interpreter','latex')
title(plane)
print(h.figure,[plotpath,'_ls_wsk_W'],'-dpng','-r300')


%% Length

MOM.length_km = MOM.length/1000;

h = plot_whisker(MOM,{'length_km'},levels,0,'PrimaryLabels',{''});
ylabel('Segment length $[\mathrm{km}]$','Interpreter','latex')
title(plane)
print(h.figure,[plotpath,'_length_wsk'],'-dpng','-r300')


%% Altitude

MOM.vstop = MOM.top - MOM.alt;

h = plot_whisker(MOM,{'alt','vstop'},levels,0,'PrimaryLabels',{'alt','vstop'});
ylabel('Altitude $[\mathrm{m}]$','Interpreter','latex')
hold on
xlim = h.axis.XLim;
plot(h.axis,xlim,prctile(MOM.top,50)*[1 1],'Color','red','LineWidth',1)
plot(h.axis,xlim,prctile(MOM.top,25)*[1 1],'Color','red','LineWidth',1,'LineStyle','--')
plot(h.axis,xlim,prctile(MOM.top,75)*[1 1],'Color','red','LineWidth',1,'LineStyle','--')
% plot(h.axis,xlim,min(MOM.top)*[1 1],'Color','red','LineWidth',1,'LineStyle',':')
% plot(h.axis,xlim,max(MOM.top)*[1 1],'Color','red','LineWidth',1,'LineStyle',':')
h.axis.XLim = xlim;
h.axis.View = [0 90];
title(plane)
print(h.figure,[plotpath,'_alt_wsk'],'-dpng','-r300')



%% Summmary of segments

sortrows(groupsummary(MOM,"level",["mean","std","min","max"],...
    ["alt","vstop","ls_W","length_km"]),"mean_alt",'descend')

groupsummary(MOM,"level",["mean","std"],...
    ["slp_sfc_W","slp_sfc_UX","slp_sfc_VY","slp_psd_W","slp_psd_UX","slp_psd_VY"])

groupsummary(MOM,"level",["mean","std"],...
    ["ar_sfc_VU","ar_psd_VU","ar_sfc_WU","ar_psd_WU","ar_sfc_WV","ar_psd_WV"])


% annotation(f,'textbox',[0 .9 .1 .1],'String',ts{i},...
%         'EdgeColor','none','FontSize',fontsize+2,'FontWeight','bold')
%     
% #!/bin/bash
% pdfjam --nup 3x1 --scale 0.95 --papersize '{36cm,10cm}' pfig08*.pdf --outfile fig08.pdf
% pdfjam --nup 3x1 --scale 0.95 --papersize '{36cm,10cm}' pfig_ex_08_x.pdf pfig_ex_08_y.pdf pfig08_z.pdf --outfile fig08ex.pdf
