
%% Settings

% Variables
vars_vel = {'UX','VY','W'};
vars_ar = {'ar_VU','ar_WU'};

% Example segments
examples = {["RF12","R2B"],["RF06","SC01"],["RF09","C6"],["RF12","CB01"]};


% Axes limits
lim.ar    = [0 1.6];
lim.s     = [0 1.3];
lim.p     = [0.4 2.5];
lim.ar_sbs= [0 2];
lim.sfc   = [5e-3 0.4];
lim.psd   = [1e-3 4];
lim.csfc  = [1e-3 3e-2];
lim.cpsd  = [4e-4 3e-2];
lim.cbth  = [4e-4 3e-2];
lim.e_ar  = [0 0.3];
lim.e_s   = [0 0.1];
lim.e_p   = [0 0.3];

% Axes labels
scale_lab_m.sfc = '$r\,[\mathrm{m}]$';
scale_lab_m.psd = '$\lambda\,[\mathrm{m}]$';
scale_lab_m.bth = '$r,\,\lambda\,[\mathrm{m}]$';
scale_lab_n.sfc = '$r/L$';
scale_lab_n.psd = '$\lambda/L$';
scale_lab_n.bth = '$r/L,\,\lambda/L$';
stat_lab.sfc  = '$D(r)\,[\mathrm{m^2\,s^{-2}}]$';
stat_lab.psd  = '$P(k)\,[\mathrm{m^3\,s^{-2}}]$';
stat_lab.csfc = '$D(r)r^{-2/3}$';
stat_lab.cpsd = '$P(k)k^{5/3}$';
stat_lab.cbth = '$D(r)r^{-2/3},\,P(k)k^{5/3}$';

leg_vel = lower(cellfun(@(x) x(1),vars_vel,'UniformOutput',false));

% Direction of the x-axis
xdir = 'reverse';

% Plot style
mks = 8;
lw = 2;
font = 16; % affects only legend
alpha = 0.25;

% Plot format
pfrm = '-dpng';

% Plot resolution
pres = '-r300';


% Output path
addpath(genpath(myprojectpath))
plotpath = [myprojectpath,filesep,'figures'];
if ~isfolder(plotpath), mkdir(plotpath), end

% Load results
load([myprojectpath,filesep,'results.mat'])
Npl = numel(planes);



%% Examples: sfc + psd

ifnorm = false;
ifraw = true;

plotpath_ex = [plotpath,filesep,'ex'];
if ~isfolder(plotpath_ex), mkdir(plotpath_ex), end


for i_p = 1:Npl    
    plane = planes{i_p};
    
    MOM = MOM_vec{i_p};

    if ifraw
        SFC = rawSFC_vec{i_p}; PSD = rawPSD_vec{i_p};
        pstyle = '.'; psuff = 'raw';
    else 
        SFC = SFC_vec{i_p}; PSD = PSD_vec{i_p};
        pstyle = 'o-'; psuff = 'avg';
    end
    
    ind_s = find( MOM.flight==examples{i_p}(1) & MOM.name==examples{i_p}(2) );
    plotpath_plane = [plotpath_ex,filesep,plane];
    
    
    for ii_s = 1:numel(ind_s)
        i_s = ind_s(ii_s);
        dr = MOM.dr(i_s);
        
        if ifnorm
            L = MOM.int_scale(i_s);
            xlab = scale_lab_n;
        else
            L = 1;
            xlab = scale_lab_m;
        end
        xlim_sfc = [dr r_max]/L;
        xlim_psd = [2*dr r_max]/L;
        tit = join([plane,MOM.flight(i_s),MOM.name(i_s),round(MOM.alt(i_s)),'m']);
        
        
        % D(r)
        
        [fig,~,co] = fig16x12('loglog',[1 1],'on','XDir',xdir,...
            'XLim',xlim_sfc,'YLim',lim.sfc);
        
        for i_v = 1:numel(vars_vel)
            var = vars_vel{i_v}; c = co(i_v,:);
            plot( SFC(i_s).r/L, SFC(i_s).(var), pstyle, ...
                'MarkerSize',mks,'Color',c,'MarkerFaceColor',c)
        end
        
        r1 = sfc_bot_factors(i_p)*dr/L;
        r2 = sfc_upp_factor*MOM.int_scale(i_s)/L;
        ind = SFC(i_s).r/L > r1 & SFC(i_s).r/L < r2;
        C = max(cellfun(@(v) mean(SFC(i_s).(v)(ind)./(SFC(i_s).r(ind)/L).^(2/3)),vars_vel));
        r12 = [r1*1.3 r2/1.3];
        plot(r12,1.5*C*r12.^(2/3),'Color',co(4,:),'LineWidth',lw)
        plot(r1*[1 1],lim.sfc,':','Color','black','LineWidth',lw)
        plot(r2*[1 1],lim.sfc,':','Color','black','LineWidth',lw)
        
        xlabel(xlab.sfc,'Interpreter','latex')
        ylabel(stat_lab.sfc,'Interpreter','latex')
        legend(leg_vel,'Location','northeast','Interpreter','latex','FontSize',font)
        title(tit)
        print(fig,join([plotpath_plane,psuff,'sfc',MOM.level(i_s),string(i_s)],'_'),pfrm,pres)
        
        
        % P(k)
        
        [fig,~,co] = fig16x12('loglog',[1 1],'on','XDir',xdir,...
            'XLim',xlim_psd,'YLim',lim.psd);
        
        for i_v = 1:numel(vars_vel)
            var = vars_vel{i_v}; c = co(i_v,:);
            plot( PSD(i_s).r/L, PSD(i_s).(var), pstyle, ...
                'MarkerSize',mks,'Color',c,'MarkerFaceColor',c)
        end
        
        r1 = psd_bot_factors(i_p)*dr/L;
        r2 = psd_upp_factor*MOM.int_scale(i_s)/L;
        ind = PSD(i_s).r/L > r1 & PSD(i_s).r/L < r2;
        C = max(cellfun(@(v) mean(PSD(i_s).(v)(ind)./(PSD(i_s).r(ind)/L).^(5/3)),vars_vel));
        r12 = [r1*1.3 r2/1.3];
        plot(r12,2*C*r12.^(5/3),'Color',co(4,:),'LineWidth',lw)
        plot(r1*[1 1],lim.psd,'--','Color','black','LineWidth',lw)
        plot(r2*[1 1],lim.psd,'--','Color','black','LineWidth',lw)
        
        xlabel(xlab.psd,'Interpreter','latex')
        ylabel(stat_lab.psd,'Interpreter','latex')
        legend(leg_vel,'Location','northeast','Interpreter','latex','FontSize',font)
        title(tit)
        print(fig,join([plotpath_plane,psuff,'psd',MOM.level(i_s),string(i_s)],'_'),pfrm,pres)
        
    end
end



%% Examples: compensated sfc + psd

ifnorm = false;
ifraw = false;

plotpath_ex = [plotpath,filesep,'ex'];
if ~isfolder(plotpath_ex), mkdir(plotpath_ex), end


for i_p = 1:Npl    
    plane = planes{i_p};
    
    MOM = MOM_vec{i_p};

    if ifraw
        SFC = rawSFC_vec{i_p}; PSD = rawPSD_vec{i_p};
        pstyle = '.'; psuff = 'raw';
    else 
        SFC = SFC_vec{i_p}; PSD = PSD_vec{i_p};
        pstyle = 'o-'; psuff = 'avg';
    end
    
    ind_s = find( MOM.flight==examples{i_p}(1) & MOM.name==examples{i_p}(2) );
    plotpath_plane = [plotpath_ex,filesep,plane];
    
    
    for ii_s = 1:numel(ind_s)
        i_s = ind_s(ii_s);
        dr = MOM.dr(i_s);
        
        if ifnorm
            L = MOM.int_scale(i_s);
            xlab = scale_lab_n;
        else
            L = 1;
            xlab = scale_lab_m;
        end
        xlim_sfc = [dr r_max]/L;
        xlim_psd = [2*dr r_max]/L;
        tit = join([plane,MOM.flight(i_s),MOM.name(i_s),round(MOM.alt(i_s)),'m']);
    
        
        % D(r)*r^-2/3 & P(k)*k^5/3
        
        [fig,~,co] = fig16x12('loglog',[1 1],'on','XDir',xdir,...
            'XLim',xlim_sfc,'YLim',lim.cbth);

        for i_v = 1:numel(vars_vel)
            var = vars_vel{i_v}; c = co(i_v,:);
            plot( PSD(i_s).r/L, PSD(i_s).(var).*PSD(i_s).k.^(5/3), pstyle, ...
                'MarkerSize',mks,'Color',c,'MarkerFaceColor',c)
        end
        for i_v = 1:numel(vars_vel)
            var = vars_vel{i_v}; c = co(i_v,:);
            plot( SFC(i_s).r/L, SFC(i_s).(var)./SFC(i_s).r.^(2/3), pstyle, ...
                'MarkerSize',mks,'Color',c)
        end

        plot(psd_upp_factor*MOM.int_scale(i_s)/L*[1 1], lim.cbth, '--','Color','black','LineWidth',lw)
        plot(psd_bot_factors(i_p)*dr/L*[1 1],           lim.cbth, '--','Color','black','LineWidth',lw)
        plot(sfc_upp_factor*MOM.int_scale(i_s)/L*[1 1], lim.cbth, ':','Color','black','LineWidth',lw)
        plot(sfc_bot_factors(i_p)*dr/L*[1 1],           lim.cbth, ':','Color','black','LineWidth',lw)

        xlabel(xlab.bth,'Interpreter','latex')
        ylabel(stat_lab.cbth,'Interpreter','latex')
        legend(leg_vel,'Location','best','Interpreter','latex','FontSize',font)
        title(tit)
        print(fig,join([plotpath_plane,psuff,'c_bth',MOM.level(i_s),string(i_s)],'_'),pfrm,pres)
        
        
        % D(r)*r^-2/3
        
%         [fig,~,co] = fig16x12('loglog',[1 1],'on','XDir',xdir,...
%             'XLim',xlim_sfc,'YLim',lim.csfc);
%         
%         for i_v = 1:numel(vars_vel)
%             var = vars_vel{i_v}; c = co(i_v,:);
%             plot( SFC(i_s).r/L, SFC(i_s).(var)./SFC(i_s).r.^(2/3), pstyle, ...
%                 'MarkerSize',mks,'Color',c,'MarkerFaceColor',c)
%         end
%         
%         plot(sfc_upp_factor*MOM.int_scale(i_s)/L*[1 1], lim.csfc,':','Color','black','LineWidth',lw)
%         plot(sfc_bot_factors(i_p)*dr/L*[1 1],           lim.csfc,':','Color','black','LineWidth',lw)
%         
%         xlabel(xlab.sfc,'Interpreter','latex')
%         ylabel(stat_lab.csfc,'Interpreter','latex')
%         legend(leg_vel,'Location','best','Interpreter','latex','FontSize',font)
%         title(tit)
%         print(fig,join([plotpath_plane,psuff,'c_sfc',MOM.level(i_s),string(i_s)],'_'),pfrm,pres)
        
        
        % P(k)*k^5/3
        
%         [fig,~,co] = fig16x12('loglog',[1 1],'on','XDir',xdir,...
%             'XLim',xlim_sfc,'YLim',lim.cpsd);
%         
%         for i_v = 1:numel(vars_vel)
%             var = vars_vel{i_v}; c = co(i_v,:);
%             plot( PSD(i_s).r/L, PSD(i_s).(var).*PSD(i_s).k.^(5/3), pstyle, ...
%                 'MarkerSize',mks,'Color',c,'MarkerFaceColor',c)
%         end
% 
%         plot(psd_upp_factor*MOM.int_scale(i_s)/L*[1 1], lim.cpsd,'--','Color','black','LineWidth',lw)
%         plot(psd_bot_factors(i_p)*dr/L*[1 1],           lim.cpsd,'--','Color','black','LineWidth',lw)
%         
%         xlabel(xlab.psd,'Interpreter','latex')
%         ylabel(stat_lab.cpsd,'Interpreter','latex')
%         legend(leg_vel,'Location','best','Interpreter','latex','FontSize',font)
%         title(tit)
%         print(fig,join([plotpath_plane,psuff,'c_psd',MOM.level(i_s),string(i_s)],'_'),pfrm,pres)
        
    end
end



%% Examples: scale-by-scale ratios

ifnorm = false;

plotpath_ex = [plotpath,filesep,'ex'];
if ~isfolder(plotpath_ex), mkdir(plotpath_ex), end


for i_p = 1:Npl    
    plane = planes{i_p};
    
    MOM = MOM_vec{i_p};
    SFC = SFC_vec{i_p};
    PSD = PSD_vec{i_p};
    
    ind_s = find( MOM.flight==examples{i_p}(1) & MOM.name==examples{i_p}(2) );
    plotpath_plane = [plotpath_ex,filesep,plane];
    
    
    for ii_s = 1:numel(ind_s)
        i_s = ind_s(ii_s);
        dr = MOM.dr(i_s);
        
        if ifnorm
            L = MOM.int_scale(i_s);
            xlab = scale_lab_n;
        else
            L = 1;
            xlab = scale_lab_m;
        end
        xlim_sfc = [dr r_max]/L;
        xlim_psd = [2*dr r_max]/L;
        tit = join([plane,MOM.flight(i_s),MOM.name(i_s),round(MOM.alt(i_s)),'m']);

        
        % Pv/Pu, Pw/Pu, Dv/Pu, Dw/Pu
        
        [fig,~,co] = fig16x12('loglin',[1 1],'on','XDir',xdir,...
            'XLim',xlim_sfc,'YLim',lim.ar_sbs+[0 1]);

        plot(PSD(i_s).r/L,PSD(i_s).ar_VU,  'o-','MarkerSize',mks,'Color',co(2,:),'MarkerFaceColor',co(2,:))
        plot(PSD(i_s).r/L,PSD(i_s).ar_WU,  'o-','MarkerSize',mks,'Color',co(3,:),'MarkerFaceColor',co(3,:))
        plot(SFC(i_s).r/L,SFC(i_s).ar_VU+1,'o-','MarkerSize',mks,'Color',co(2,:))
        plot(SFC(i_s).r/L,SFC(i_s).ar_WU+1,'o-','MarkerSize',mks,'Color',co(3,:))

        plot(xlim_sfc,[1 1]*4/3, '--','Color','black','LineWidth',lw)
        plot(xlim_sfc,[1 1]*4/3+1,':','Color','black','LineWidth',lw)

        xlabel(xlab.bth,'Interpreter','latex')
        legend({'$P_v/D_u$','$P_w/D_u$','$D_v/D_u+1$','$D_w/D_u+1$'},...
            'Location','northeast','Interpreter','latex','FontSize',font)
        title(tit)
        print(fig,join([plotpath_plane,'ar_bth',MOM.level(i_s),string(i_s)],'_'),pfrm,pres)
        
        
        % Dv/Du, Dw/Du
        
%         [fig,~,co] = fig16x12('loglin',[1 1],'on','XDir',xdir,...
%             'XLim',xlim_sfc,'YLim',lim.ar_sbs);
% 
%         plot(SFC(i_s).r/L,SFC(i_s).ar_VU,'o-','MarkerSize',mks,'Color',co(2,:),'MarkerFaceColor',co(2,:))
%         plot(SFC(i_s).r/L,SFC(i_s).ar_WU,'o-','MarkerSize',mks,'Color',co(3,:),'MarkerFaceColor',co(3,:))
% 
%         plot(xlim_sfc,[1 1]*4/3,'--','Color','black','LineWidth',lw)
%         plot(sfc_upp_factor*MOM.int_scale(i_s)/L*[1 1], lim.ar_sbs, '--','Color','black','LineWidth',lw)
%         plot(sfc_bot_factors(i_p)*dr/L*[1 1],           lim.ar_sbs, '--','Color','black','LineWidth',lw)
% 
%         xlabel(xlab.sfc,'Interpreter','latex')
%         legend({'$D_v/D_u$','$D_w/D_u$'},'Location','best','Interpreter','latex','FontSize',font)
%         title(tit)
%         print(fig,join([plotpath_plane,'ar_sfc',MOM.level(i_s),string(i_s)],'_'),pfrm,pres)

        
        % Pv/Pu, Pw/Pu
            
%         [fig,~,co] = fig16x12('loglin',[1 1],'on','XDir',xdir,...
%             'XLim',xlim_psd,'YLim',lim.ar_sbs);
% 
%         plot(PSD(i_s).r/L,PSD(i_s).ar_VU,'o-','MarkerSize',mks,'Color',co(2,:),'MarkerFaceColor',co(2,:))
%         plot(PSD(i_s).r/L,PSD(i_s).ar_WU,'o-','MarkerSize',mks,'Color',co(3,:),'MarkerFaceColor',co(3,:))
% 
%         plot(xlim_psd,[1 1]*4/3,'--','Color','black','LineWidth',lw)
%         plot(psd_upp_factor*MOM.int_scale(i_s)/L*[1 1], lim.ar_sbs, '--','Color','black','LineWidth',lw)
%         plot(psd_bot_factors(i_p)*dr/L*[1 1],           lim.ar_sbs, '--','Color','black','LineWidth',lw)
% 
%         xlabel(xlab.psd,'Interpreter','latex')
%         legend({'$P_v/P_u$','$P_w/P_u$'},'Location','best','Interpreter','latex','FontSize',font)
%         title(tit)
%         print(fig,join([plotpath_plane,'ar_psd',MOM.level(i_s),string(i_s)],'_'),pfrm,pres)
        
    end
end



%% Bulk ratios

plotpath_bulk = [plotpath,filesep,'bulk'];
if ~isfolder(plotpath_bulk), mkdir(plotpath_bulk), end


for i_p = 1:Npl    
    plane = planes{i_p};
    
    MOM = MOM_vec{i_p};
    
    levels = sortrows(groupsummary(MOM,{'level'},{'mean'},{'alt'}),'mean_alt','descend').level';
    for i_l = 1:numel(levels)
        MOM.level_id(MOM.level==levels(i_l)) = i_l;
    end
    
    if ismember(plane,{'ATR-EUREC4A','TO-POST'})
        dirvar = 'dir2';
    else
        dirvar = '';
    end

    
    % (Pv/Pu,Dv/Du)
    
    [fig,ax] = plot_xy_uni(MOM,{'ar_psd_VU'},{'ar_sfc_VU'},'level_id','',dirvar,true,...
        {'ver4/3','hor4/3','cross1'},[],'XLim',lim.ar,'YLim',lim.ar);
    plot(ax,4/3,4/3,'d','Color',"#77AC30",'MarkerFaceColor',"#77AC30",'MarkerSize',12)
    legend(horzcat(levels,{'HIT'}),'Location','northwest','Interpreter','latex')
    xlabel('$P_v/P_u$','Interpreter','latex')
    ylabel('$D_v/D_u$','Interpreter','latex')
    title(plane)
    print(fig,[plotpath_bulk,filesep,plane,'_ar_vu'],pfrm,pres)
    
    
    % (Pw/Pu,Dw/Du)
    
    [fig,ax] = plot_xy_uni(MOM,{'ar_psd_WU'},{'ar_sfc_WU'},'level_id','',dirvar,true,...
        {'ver4/3','hor4/3','cross1'},[],'XLim',lim.ar,'YLim',lim.ar);
    plot(ax,4/3,4/3,'d','Color',"#77AC30",'MarkerFaceColor',"#77AC30",'MarkerSize',12)
    legend(horzcat(levels,{'HIT'}),'Location','southeast','Interpreter','latex')
    xlabel('$P_w/P_u$','Interpreter','latex')
    ylabel('$D_w/D_u$','Interpreter','latex')
    title(plane)
    print(fig,[plotpath_bulk,filesep,plane,'_ar_wu'],pfrm,pres)
    
    
    % (p,s)

    [fig,ax] = plot_xy_uni(MOM,{'slp_psd_UX','slp_psd_VY','slp_psd_W'},...
        {'slp_sfc_UX','slp_sfc_VY','slp_sfc_W'},'yN','level_id',dirvar,false,...
        {'hor2/3','ver5/3'},[],'XLim',lim.p,'YLim',lim.s);
    plot(ax,5/3,2/3,'d','Color',"#77AC30",'MarkerFaceColor',"#77AC30",'MarkerSize',12)
    legend(horzcat({'u','v','w'},levels,{'K41'}),'Location','northwest','Interpreter','latex')
    xlabel('$p$','Interpreter','latex')
    ylabel('$s$','Interpreter','latex')
    title(plane)
    print(fig,[plotpath_bulk,filesep,plane,'_slp'],pfrm,pres)
    
end



%% Scale-by-scale ratios

ifnorm = true;
lvl_std = {'cloud-base','cloud-base','sub-cloud','near-surface'};

plotpath_sbs = [plotpath,filesep,'sbs'];
if ~isfolder(plotpath_sbs), mkdir(plotpath_sbs), end 


for i_p = 1:Npl  
    plane = planes{i_p};
    
    MOM   = MOM_vec{i_p};
    avSFC = avSFC_vec{i_p};
    avPSD = avPSD_vec{i_p};
    
    levels = [avSFC(:).level]';
    Nlvl = numel(levels);
    for i_l = 1:Nlvl
        MOM.level_id(MOM.level==levels(i_l)) = i_l;
    end
    
    if ifnorm
        xvar = 'r_iL';
        vars_ari = {'ar_VU_iL','ar_WU_iL'};
        xlim_sfc = r_L_range;
        xlim_psd = r_L_range;
        xlab = scale_lab_n;
    else
        xvar = 'r_i';
        vars_ari = {'ar_VU_i','ar_WU_i'};
        xlim_sfc = [min(MOM.dr) r_max];
        xlim_psd = [2*min(MOM.dr) r_max];
        xlab = scale_lab_m;
    end
    
    
    % Dv/Du + Pv/Du and Dw/Du + Pw/Du
    
    for i_v = 1:numel(vars_ari)
        var = vars_ari{i_v};
        
        [fig,~,co] = fig16x12('loglin',[1 1],'on','XDir',xdir,...
            'XLim',xlim_psd,'YLim',lim.ar_sbs+[0 1]);
        
        i_l = find(strcmp(lvl_std{i_p},levels));
        if ~isempty(i_l)
%         for i_l = 1:Nlvl
            c = co(i_l,:);
%             plot( avPSD(i_l).(xvar), avPSD(i_l).(var) + avPSD(i_l).([var,'_std']), '--', 'Color',c,'HandleVisibility','off')
%             plot( avPSD(i_l).(xvar), avPSD(i_l).(var) - avPSD(i_l).([var,'_std']), '--', 'Color',c,'HandleVisibility','off')
%             plot( avSFC(i_l).(xvar), avSFC(i_l).(var) + avSFC(i_l).([var,'_std']) +1, '--', 'Color',c,'HandleVisibility','off')
%             plot( avSFC(i_l).(xvar), avSFC(i_l).(var) - avSFC(i_l).([var,'_std']) +1, '--', 'Color',c,'HandleVisibility','off')
            x = avPSD(i_l).(xvar); y = avPSD(i_l).(var); dy = avPSD(i_l).([var,'_std']);
            x = x(~isnan(y)); dy = dy(~isnan(y)); y = y(~isnan(y));
            patch( [x; flip(x)], [y-dy; flip(y+dy)], c,'FaceAlpha',alpha,'LineStyle','none','HandleVisibility','off')
            x = avSFC(i_l).(xvar); y = avSFC(i_l).(var); dy = avSFC(i_l).([var,'_std']);
            x = x(~isnan(y)); dy = dy(~isnan(y)); y = y(~isnan(y));
            patch( [x; flip(x)], [y-dy; flip(y+dy)]+1, c,'FaceAlpha',alpha,'LineStyle','none','HandleVisibility','off')
        end
        
        for i_l = 1:Nlvl
            c = co(i_l,:);
            plot( avPSD(i_l).(xvar), avPSD(i_l).(var), 'o-', ...
                'MarkerSize',mks,'Color',c,'MarkerFaceColor',c)
        end
        for i_l = 1:Nlvl
            c = co(i_l,:);
            plot( avSFC(i_l).(xvar), avSFC(i_l).(var)+1, 'o-', ...
                'MarkerSize',mks,'Color',c)
        end
                    
        plot(xlim_sfc,4/3*[1 1], '--','Color','black','LineWidth',lw)
        plot(xlim_sfc,4/3*[1 1]+1,':','Color','black','LineWidth',lw)

        xlabel(xlab.bth,'Interpreter','latex')
        ylabel(['$P_',lower(var(4)),'/P_u,\,D_',lower(var(4)),'/D_u+1$'],'Interpreter','latex')
        legend(levels,'Location','northeast','Interpreter','latex') %,'FontSize',font)
        title(plane)
        print(fig,[plotpath_sbs,filesep,plane,'_bth_',lower(var(1:5))],pfrm,pres)
    end   

end



%% Uncertainties

plotpath_bulk = [plotpath,filesep,'bulk'];
if ~isfolder(plotpath_bulk), mkdir(plotpath_bulk), end


for i_p = 1:Npl    
    plane = planes{i_p};
    
    MOM = MOM_vec{i_p};

    levels = sortrows(groupsummary(MOM,{'level'},{'mean'},{'alt'}),'mean_alt','descend').level';
    for i_l = 1:numel(levels)
        MOM.level_id(MOM.level==levels(i_l)) = i_l;
    end

    
    % Ratios
    
    h = plot_whisker(MOM,{'e_ar_sfc_VU','e_ar_psd_VU','e_ar_sfc_WU','e_ar_psd_WU'},...
        levels,0,'PrimaryLabels',{'$D\,v/u$','$P\,v/u$','$D\,w/u$','$P\,w/u$'},...{'$D_v/D_u$','$P_v/P_u$','$D_w/D_u$','$P_w/P_u$'},...
        'DataLim',lim.e_ar);
    hold on
    h.axis.YLim = lim.e_ar;
    ylabel('Uncertainty','Interpreter','latex')
    title(plane,'FontSize',12)
    print(h.figure,[plotpath_bulk,filesep,plane,'_e_wsk_ar'],pfrm,pres)
    
    
    % s
    h = plot_whisker(MOM,{'e_slp_sfc_UX','e_slp_sfc_VY','e_slp_sfc_W'},...
        levels,0,'PrimaryLabels',{'$s\,u$','$s\,v$','$s\,w$'},...{'$s_u$','$s_v$','$s_w$'},...
        'DataLim',lim.e_s);
    hold on
    h.axis.YLim = lim.e_s;
    ylabel('Uncertainty','Interpreter','latex')
    title(plane,'FontSize',12)
    print(h.figure,[plotpath_bulk,filesep,plane,'_e_wsk_slp_sfc'],pfrm,pres)
    
    
    % p
    
    h = plot_whisker(MOM,{'e_slp_psd_UX','e_slp_psd_VY','e_slp_psd_W'},...
        levels,0,'PrimaryLabels',{'$p\,u$','$p\,v$','$p\,w$'},...{'$p_u$','$p_v$','$p_w$'},
        'DataLim',lim.e_p);
    hold on
    h.axis.YLim = lim.e_p;
    ylabel('Uncertainty','Interpreter','latex')
    title(plane,'FontSize',12)
    print(h.figure,[plotpath_bulk,filesep,plane,'_e_wsk_slp_psd'],pfrm,pres)

end



%% Scale-by-scale bundles

ifnorm = true;

plotpath_bundle = [plotpath,filesep,'bundles'];
if ~isfolder(plotpath_bundle), mkdir(plotpath_bundle), end


for i_p = 1:Npl    
    plane = planes{i_p};
    
    MOM = MOM_vec{i_p};
    SFC = SFC_vec{i_p};
    PSD = PSD_vec{i_p};
    Nseg = size(MOM,1);
    
    levels = sortrows(groupsummary(MOM,{'level'},{'mean'},{'alt'}),...
        'mean_alt','descend').level';
    Nlvl = numel(levels);
    for i_l = 1:Nlvl
        MOM.level_id(MOM.level==levels(i_l)) = i_l;
    end
    
    if ifnorm
        xlim_sfc = r_L_range;
        xlim_psd = r_L_range;
    else
        xlim_sfc = [min(MOM.dr) r_max];
        xlim_psd = [2*min(MOM.dr) r_max];
    end
    
    
    % Dv/Du and Dw/Du
    
    for i_v = 1:numel(vars_ar)
        var = vars_ar{i_v};
        
        [fig,~,co] = fig16x12('loglin',[1 1],'on','XDir',xdir,...
            'XLim',xlim_sfc,'YLim',lim.ar_sbs);
        
        for i_l = 1:Nlvl
            plot(nan,nan,'Color',co(i_l,:))
        end

        for i_s = 1:Nseg
            if ifnorm
                L = MOM.int_scale(i_s);
            else
                L = 1;
            end
            plot(SFC(i_s).r/L,SFC(i_s).(var),'-','Color',co(MOM.level_id(i_s),:))
        end

        plot(xlim_sfc,4/3*[1 1],'--','Color','black','LineWidth',lw)

        xlabel(xlab.sfc,'Interpreter','latex')
        ylabel(['$D_',lower(var(end-1)),'/D_u$'],'Interpreter','latex')
        legend(levels,'Location','best','Interpreter','latex')
        title(plane)
        print(fig,[plotpath_bundle,filesep,plane,'_sfc_',lower(var)],pfrm,pres)
    end
    
    
    % Pv/Pu and Pw/Pu
    
    for i_v = 1:numel(vars_ar)
        var = vars_ar{i_v};
        
        [fig,~,co] = fig16x12('loglin',[1 1],'on','XDir',xdir,...
            'XLim',xlim_psd,'YLim',lim.ar_sbs);
       
        for i_l = 1:Nlvl
            plot(nan,nan,'Color',co(i_l,:))
        end

        for i_s = 1:Nseg
            if ifnorm
                L = MOM.int_scale(i_s);
            else
                L = 1;
            end
            plot(PSD(i_s).r/L,PSD(i_s).(var),'-','Color',co(MOM.level_id(i_s),:))
        end

        plot(xlim_psd,4/3*[1 1],'--','Color','black','LineWidth',lw)

        xlabel(xlab.psd,'Interpreter','latex')
        ylabel(['$P_',lower(var(end-1)),'/P_u$'],'Interpreter','latex')
        legend(levels,'Location','best','Interpreter','latex')
        title(plane)
        print(fig,[plotpath_bundle,filesep,plane,'_psd_',lower(var)],pfrm,pres)
    end
    
end

