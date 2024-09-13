
scale_by_scale_2_avg


%% Settings

% Variables
vars_vel = {'UX','VY','W'};
vars_ar = {'ar_VU','ar_WU'};

% Fit range: bottom limit (= factor*dr)
sfc_bot_factors = [2 2 2 3]; % list of factors for each plane/experiment
psd_bot_factors = [4 4 4 6];

% Fit range: upper limit (= factor*L)
sfc_upp_factor = 1.0; % one number for all planes/experiments
psd_upp_factor = 2.0; 

% Example segments to plot
examples = {["RF12","R2B"],["RF06","SC01"],["RF09","C6"],["RF12","CB01"]};


% Direction of the x-axis
xdir = 'reverse';

% Axes lims
lim.sfc  = [5e-3 0.4];
lim.psd  = [1e-3 4];
lim.csfc = [1e-3 3e-2];
lim.cpsd = [4e-4 3e-2];
lim.cbth = [4e-4 3e-2];
lim.ratio= [0 2];
% liminf = cell2struct(repmat({[-inf inf]},6,1),fieldnames(lim));

% Axes labels
xlabm.sfc = '$r\,[\mathrm{m}]$';
xlabm.psd = '$\lambda\,[\mathrm{m}]$';
xlabm.bth = '$r,\,\lambda\,[\mathrm{m}]$';
xlabn.sfc = '$r/L$';
xlabn.psd = '$\lambda/L$';
xlabn.bth = '$r/L,\,\lambda/L$';
ylab.sfc  = '$D(r)\,[\mathrm{m^2\,s^{-2}}]$';
ylab.psd  = '$P(k)\,[\mathrm{m^3\,s^{-2}}]$';
ylab.csfc = '$D(r)r^{-2/3}$';
ylab.cpsd = '$P(k)k^{5/3}$';
ylab.cbth = '$D(r)r^{-2/3},\,P(k)k^{5/3}$';


% Plot style
mks = 8;
lw = 2;

% Plot format
pfrm = '-dpng';
% Plot resolutions
pres = '-r300';


% Plot output path
addpath(genpath(myprojectpath))
plotpath = [myprojectpath,filesep,'figures',filesep,'scale_by_scale'];
if ~isfolder(plotpath), mkdir(plotpath), end



%% Example segments: sfc + psd

ifnorm = false;
ifraw = true;


leg = lower(cellfun(@(x) x(1),vars_vel,'UniformOutput',false));

for i_p = 1:Npl    
    plane = planes{i_p};
    
    MOM = MOM_vec{i_p};
    Nseg = size(MOM,1);

    if ifraw
        SFC = rawSFC_vec{i_p};
        PSD = rawPSD_vec{i_p};
        pstyle = '.';
        plotpath_ex = [plotpath,filesep,'ex_raw'];
    else 
        SFC = SFC_vec{i_p};
        PSD = PSD_vec{i_p};
        pstyle = 'o-';
        plotpath_ex = [plotpath,filesep,'ex_avg'];
    end
    
    if ifnorm
        plotpath_ex = [plotpath_ex,'_norm'];
        xlab = xlabn;
    else
        xlab = xlabm;
    end
    if ~isfolder(plotpath_ex), mkdir(plotpath_ex), end
    
    ind_s = find( MOM.flight==examples{i_p}(1) & MOM.name==examples{i_p}(2) );
    plotpath_plane = [plotpath_ex,filesep,plane,'_'];
    
    
    for ii_s = 1:numel(ind_s)
        i_s = ind_s(ii_s);
        dr = MOM.dr(i_s);
        
        if ifnorm
            L = MOM.int_scale(i_s);
        else
            L = 1;
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
        ylabel(ylab.sfc,'Interpreter','latex')
        legend(leg,'Location','northeast','Interpreter','latex')
        title(tit)
        print(fig,join([[plotpath_plane,'sfc'],MOM.level(i_s),string(i_s)],'_'),pfrm,pres)
        
        
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
        ylabel(ylab.psd,'Interpreter','latex')
        legend(leg,'Location','northeast','Interpreter','latex')
        title(tit)
        print(fig,join([[plotpath_plane,'psd'],MOM.level(i_s),string(i_s)],'_'),pfrm,pres)
        
        
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
%         ylabel(ylab.csfc,'Interpreter','latex')
%         legend(leg,'Location','best','Interpreter','latex')
%         title(tit)
%         print(fig,join([[plotpath_plane,'c_sfc'],MOM.level(i_s),string(i_s)],'_'),pfrm,pres)
        
        
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
%         ylabel(ylab.cpsd,'Interpreter','latex')
%         legend(leg,'Location','best','Interpreter','latex')
%         title(tit)
%         print(fig,join([[plotpath_plane,'c_psd'],MOM.level(i_s),string(i_s)],'_'),pfrm,pres)
        
    end
end



%% Example segments: ratios

ifnorm = false;


for i_p = 1:Npl    
    plane = planes{i_p};
    
    MOM = MOM_vec{i_p};
    Nseg = size(MOM,1);
  
    SFC = SFC_vec{i_p};
    PSD = PSD_vec{i_p};
    pstyle = 'o-';
    plotpath_ex = [plotpath,filesep,'ex_avg'];
 
    if ifnorm
        plotpath_ex = [plotpath_ex,'_norm'];
        xlab = xlabn;
    else
        xlab = xlabm;
    end
    if ~isfolder(plotpath_ex), mkdir(plotpath_ex), end
    
    ind_s = find( MOM.flight==examples{i_p}(1) & MOM.name==examples{i_p}(2) );
    plotpath_plane = [plotpath_ex,filesep,plane,'_'];
    
    
    for ii_s = 1:numel(ind_s)
        i_s = ind_s(ii_s);
        dr = MOM.dr(i_s);
        
        if ifnorm
            L = MOM.int_scale(i_s);
        else
            L = 1;
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
        ylabel(ylab.cbth,'Interpreter','latex')
        leg = lower(cellfun(@(x) x(1),vars_vel,'UniformOutput',false));
        legend(leg,'Location','best','Interpreter','latex')
        title(tit)
        print(fig,join([[plotpath_plane,'c_bth'],MOM.level(i_s),string(i_s)],'_'),pfrm,pres)
        
        
        % Dv/Du, Dw/Du
        
%         [fig,~,co] = fig16x12('loglin',[1 1],'on','XDir',xdir,...
%             'XLim',xlim_sfc,'YLim',lim.ratio);
% 
%         plot(SFC(i_s).r/L,SFC(i_s).ar_VU,'o-','MarkerSize',mks,'Color',co(2,:),'MarkerFaceColor',co(2,:))
%         plot(SFC(i_s).r/L,SFC(i_s).ar_WU,'o-','MarkerSize',mks,'Color',co(3,:),'MarkerFaceColor',co(3,:))
% 
%         plot(xlim_sfc,[1 1]*4/3,'--','Color','black','LineWidth',lw)
%         plot(sfc_upp_factor*MOM.int_scale(i_s)/L*[1 1], lim.ratio, '--','Color','black','LineWidth',lw)
%         plot(sfc_bot_factors(i_p)*dr/L*[1 1],           lim.ratio, '--','Color','black','LineWidth',lw)
% 
%         xlabel(xlab.sfc,'Interpreter','latex')
%         legend({'$D_v/D_u$','$D_w/D_u$'},'Location','best','Interpreter','latex')
%         title(tit)
%         print(fig,join([[plotpath_plane,'ar_sfc'],MOM.level(i_s),string(i_s)],'_'),pfrm,pres)

        
        % Pv/Pu, Pw/Pu
            
%         [fig,~,co] = fig16x12('loglin',[1 1],'on','XDir',xdir,...
%             'XLim',xlim_psd,'YLim',lim.ratio);
% 
%         plot(PSD(i_s).r/L,PSD(i_s).ar_VU,'o-','MarkerSize',mks,'Color',co(2,:),'MarkerFaceColor',co(2,:))
%         plot(PSD(i_s).r/L,PSD(i_s).ar_WU,'o-','MarkerSize',mks,'Color',co(3,:),'MarkerFaceColor',co(3,:))
% 
%         plot(xlim_psd,[1 1]*4/3,'--','Color','black','LineWidth',lw)
%         plot(psd_upp_factor*MOM.int_scale(i_s)/L*[1 1], lim.ratio, '--','Color','black','LineWidth',lw)
%         plot(psd_bot_factors(i_p)*dr/L*[1 1],           lim.ratio, '--','Color','black','LineWidth',lw)
% 
%         xlabel(xlab.psd,'Interpreter','latex')
%         legend({'$P_v/P_u$','$P_w/P_u$'},'Location','best','Interpreter','latex')
%         title(tit)
%         print(fig,join([[plotpath_plane,'ar_psd'],MOM.level(i_s),string(i_s)],'_'),pfrm,pres)
            

        % Pv/Pu, Pw/Pu, Dv/Pu, Dw/Pu
        
        [fig,~,co] = fig16x12('loglin',[1 1],'on','XDir',xdir,...
            'XLim',xlim_sfc,'YLim',lim.ratio+[0 1]);

        plot(PSD(i_s).r/L,PSD(i_s).ar_VU,  'o-','MarkerSize',mks,'Color',co(2,:),'MarkerFaceColor',co(2,:))
        plot(PSD(i_s).r/L,PSD(i_s).ar_WU,  'o-','MarkerSize',mks,'Color',co(3,:),'MarkerFaceColor',co(3,:))
        plot(SFC(i_s).r/L,SFC(i_s).ar_VU+1,'o-','MarkerSize',mks,'Color',co(2,:))
        plot(SFC(i_s).r/L,SFC(i_s).ar_WU+1,'o-','MarkerSize',mks,'Color',co(3,:))

        plot(xlim_sfc,[1 1]*4/3, '--','Color','black','LineWidth',lw)
        plot(xlim_sfc,[1 1]*4/3+1,':','Color','black','LineWidth',lw)

        xlabel(xlab.bth,'Interpreter','latex')
        legend({'$P_v/D_u$','$P_w/D_u$','$D_v/D_u+1$','$D_w/D_u+1$'},...
            'Location','northeast','Interpreter','latex')
        title(tit)
        print(fig,join([[plotpath_plane,'ar_bth'],MOM.level(i_s),string(i_s)],'_'),pfrm,pres)
       
    end
    
end



%% Level averages

ifnorm = true;


if ifnorm
    plotpath_sbs = [plotpath,filesep,'lvl_norm'];
    xlab = xlabn;
else
    plotpath_sbs = [plotpath,filesep,'lvl'];
    xlab = xlabm;
end
if ~isfolder(plotpath_sbs), mkdir(plotpath_sbs), end 


for i_p = 1:Npl    
    plane = planes{i_p};
    
    MOM = MOM_vec{i_p};
    
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
    else
        xvar = 'r_i';
        vars_ari = {'ar_VU_i','ar_WU_i'};
        xlim_sfc = [min(MOM.dr) r_max];
        xlim_psd = [2*min(MOM.dr) r_max];
    end
    
    
    % Dv/Du and Dw/Du
    
%     for i_v = 1:numel(vars_ari)
%         var = vars_ari{i_v};
%         
%         [fig,~,co] = fig16x12('loglin',[1 1],'on','XDir',xdir,...
%             'XLim',xlim_sfc,'YLim',lim.ratio);
% 
%         for i_l = 1:Nlvl
%             c = co(i_l,:);
%             plot( avSFC(i_l).(xvar), avSFC(i_l).(var), 'o-', ...
%                 'MarkerSize',mks,'Color',c,'MarkerFaceColor',c)
%         end
% 
%         plot(xlim_sfc,4/3*[1 1],'--','Color','black','LineWidth',lw)
% 
%         xlabel(xlab.sfc,'Interpreter','latex')
%         ylabel(['$D_',lower(var(4)),'/D_u$'],'Interpreter','latex')
%         legend(levels,'Location','best','Interpreter','latex')
%         title(plane)
%         print(fig,[plotpath_sbs,filesep,plane,'_sfc_',lower(var(1:5))],pfrm,pres)
%     end
    
    
    % Pv/Du and Pw/Du
    
%     for i_v = 1:numel(vars_ari)
%         var = vars_ari{i_v};
%         
%         [fig,~,co] = fig16x12('loglin',[1 1],'on','XDir',xdir,...
%             'XLim',xlim_psd,'YLim',lim.ratio);
% 
%         for i_l = 1:Nlvl
%             c = co(i_l,:);
%             plot( avPSD(i_l).(xvar), avPSD(i_l).(var), 'o-', ...
%                 'MarkerSize',mks,'Color',c,'MarkerFaceColor',c)
%         end
% 
%         plot(xlim_psd,4/3*[1 1],'--','Color','black','LineWidth',lw)
% 
%         xlabel(xlab.psd,'Interpreter','latex')
%         ylabel(['$P_',lower(var(4)),'/P_u$'],'Interpreter','latex')
%         legend(levels,'Location','best','Interpreter','latex')
%         title(plane)
%         print(fig,[plotpath_sbs,filesep,plane,'_psd_',lower(var(1:5))],pfrm,pres)
%     end
    
    
    % Dv/Du + Pv/Du and Dw/Du + Pw/Du
    
    for i_v = 1:numel(vars_ari)
        var = vars_ari{i_v};
        
        [fig,~,co] = fig16x12('loglin',[1 1],'on','XDir',xdir,...
            'XLim',xlim_psd,'YLim',lim.ratio+[0 1]);
        
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
        legend(levels,'Location','northeast','Interpreter','latex')
        title(plane)
        print(fig,[plotpath_sbs,filesep,plane,'_bth_',lower(var(1:5))],pfrm,pres)
    end
    
end
 


%% Bundles

ifnorm = true;


if ifnorm
    plotpath_sbs = [plotpath,filesep,'bundle_norm'];
    xlab = xlabn;
else
    plotpath_sbs = [plotpath,filesep,'bundle'];
    xlab = xlabm;
end
if ~isfolder(plotpath_sbs), mkdir(plotpath_sbs), end 


for i_p = 1:Npl    
    plane = planes{i_p};
    
    MOM = MOM_vec{i_p};
    Nseg = size(MOM,1);
    
    SFC = SFC_vec{i_p};
    PSD = PSD_vec{i_p};
    
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
            'XLim',xlim_sfc,'YLim',lim.ratio);
        
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
        print(fig,[plotpath_sbs,filesep,plane,'_sfc_',lower(var)],pfrm,pres)
    end
    
    
    % Pv/Pu and Pw/Pu
    
    for i_v = 1:numel(vars_ar)
        var = vars_ar{i_v};
        
        [fig,~,co] = fig16x12('loglin',[1 1],'on','XDir',xdir,...
            'XLim',xlim_psd,'YLim',lim.ratio);
       
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
        print(fig,[plotpath_sbs,filesep,plane,'_psd_',lower(var)],pfrm,pres)
    end
    
end