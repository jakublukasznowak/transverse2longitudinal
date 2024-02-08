
% Detect horizontal flight segments based on the small derivatives of
% altitude and true heading with respect ot distance.
% See Appendix B of the manuscript for description of the criteria.


function SEG = calc_seg(DATA,options)

arguments
    DATA (:,1) struct
    options.Plots (1,1) logical = false
    options.AltAvScale (1,1) {mustBePositive, mustBeFinite, mustBeNonempty} = 4e3 % m
    options.AltDrvLimit (1,1) {mustBePositive, mustBeFinite, mustBeNonempty} = 0.01 % m/m
    options.HdgAvScale (1,1) {mustBePositive, mustBeFinite, mustBeNonempty} = 20e3 % m
    options.HdgDrvLimit (1,1) {mustBePositive, mustBeFinite, mustBeNonempty} = 0.003 % deg/m
    options.TASLimit (1,1) {mustBePositive, mustBeFinite, mustBeNonempty} = 0.8 % 
    options.MinLength (1,1) {mustBePositive, mustBeFinite, mustBeNonempty} = 20e3 % m
    options.MaxTrend (1,1) {mustBePositive, mustBeFinite, mustBeNonempty} = 0.002 % m/m
    options.Title (1,1) string = "EXPERIMENT"
end


Nf = size(DATA,1);

SEG = cell(Nf,1);

for i_f = 1:Nf
    
    fsamp = DATA(i_f).fsamp;
    tas  = median(DATA(i_f).TAS,'omitnan');
    
    L = length(DATA(i_f).time);
    x = (1:L)'*tas/fsamp/1e3; % km


    % Condition (1): small derivative of altitude
    dzdx = diff(DATA(i_f).ALT)*fsamp/tas;
    dzdx_m = movmean(dzdx,floor(options.AltAvScale/tas*fsamp));
    mask_alt = [false; abs(dzdx_m) < options.AltDrvLimit ];
    
    % Condition (2): small derivative of heading
    dhdx = diff(unwrap(DATA(i_f).THDG))*fsamp/tas;
    dhdx_m = movmean(dhdx,floor(options.HdgAvScale/tas*fsamp));
    mask_hdg = [false; abs(dhdx_m) < options.HdgDrvLimit];

    % Condition (3): minimum TAS
    tas_m = movmean(DATA(i_f).TAS,floor(options.AltAvScale/tas*fsamp));
    mask_tas = tas_m > options.TASLimit*tas;
    
    mask_2 = mask_alt & mask_hdg;
    mask_3 = mask_alt & mask_hdg & mask_tas;
    
    % Condition (4): minimum segment length
    ind_3 = mask2ind(mask_3);
    ind_len = ind_3( diff(ind_3,1,2)*tas/fsamp >= options.MinLength, :);
    mask_len = ind2mask(ind_len,L);
    
    % Condition (5): small altitude trend within the segment
    slopes = zeros(size(ind_len,1),1);
    for i = 1:size(ind_len,1)
        p = polyfit(x(ind_len(i,1):ind_len(i,2)),DATA(i_f).ALT(ind_len(i,1):ind_len(i,2)),1);
        slopes(i) = p(1);
    end
    ind_trend = ind_len( abs(slopes) < options.MaxTrend*1e3, :);
    mask_trend = ind2mask(ind_trend,L);

    ind_final = ind_trend;

    
    if options.Plots
        
        set(groot,'defaultLineMarkerSize',10);
        font = 14;
        
        f = figure('Color','white','PaperUnits','centimeters',...
            'PaperSize',[36 20],'PaperPosition',[0 0 36 20]);

        t = tiledlayout(3,1,'Parent',f);
        
        ax1 = nexttile([2 1]); hold on, grid on
        plot(x,DATA(i_f).ALT,'.')
        plot(x(mask_alt),DATA(i_f).ALT(mask_alt),'.')
        plot(x(mask_2),DATA(i_f).ALT(mask_2)+100,'.')
        plot(x(mask_3),DATA(i_f).ALT(mask_3)+200,'.')
        plot(x(mask_len),DATA(i_f).ALT(mask_len)+300,'.')
        plot(x(mask_trend),DATA(i_f).ALT(mask_trend)+400,'.')
        legend({'$z$','$dz/dx$','+$d\psi/dx$','+TAS','+length','+trend'},'Interpreter','latex') %,'Location','best')
        ylabel('Altitude [m]','Interpreter','latex')
        
        ax2 = nexttile; hold on, grid on
        plot(x,DATA(i_f).THDG,'.')
        plot(x(mask_hdg),DATA(i_f).THDG(mask_hdg),'.')
        legend({'$\psi$','$d\psi/dx$'},'Interpreter','latex') %,'Location','best')
        ylabel('Heading [deg]','Interpreter','latex')
        
        t.TileSpacing = 'compact';
        t.Padding = 'compact';
        
        tit = title(t,join([options.Title,DATA(i_f).flight]));
        xl = xlabel(t,'Distance [km]','Interpreter','latex');
        
        ax1.Box = 'on'; ax2.Box = 'on';
        ax1.FontSize = font; ax2.FontSize = font;
        tit.FontSize = font; xl.FontSize = font;
%         ax1.TickLabelInterpreter = 'latex'; ax2.TickLabelInterpreter = 'latex';
        ax2.YLim = [0 360];
        
        print(f,strcat('seg_',options.Title,'_',DATA(i_f).flight),'-dpng','-r300')

%         figure, hold on, grid on
% %         plot(x(2:end),   1e3*dzdx)
%         plot(x(2:end),   1e3*dzdx_m)
%         plot(x([2,end]), 1e3*[1 1]*options.AltDrvLimit,'Color','black')
%         plot(x([2,end]),-1e3*[1 1]*options.AltDrvLimit,'Color','black')
%         title(DATA(i_f).flight)
%         xlabel('x [km]')
%         ylabel('dz/dx [m/km]')
%         
%         figure, hold on, grid on
% %         plot(x(2:end),   1e3*dhdx)
%         plot(x(2:end),   1e3*dhdx_m)
%         plot(x([2,end]), 1e3*[1 1]*options.HdgDrvLimit,'Color','black')
%         plot(x([2,end]),-1e3*[1 1]*options.HdgDrvLimit,'Color','black')
%         title(DATA(i_f).flight)
%         xlabel('x [km]')
%         ylabel('dh/dx [deg/km]')
        
%         figure, hold on, grid on
%         plot(x,tas_m)
%         plot(x([1,end]),[1 1]*options.TASLimit*tas,'Color','black')
%         title(DATA(i_f).flight)
%         xlabel('x [km]')
%         ylabel('TAS [m/s]')
    end
        
    
    
    seg = table;
    
    Nseg = size(ind_final,1);
    seg.flight = repmat(DATA(i_f).flight,Nseg,1);
    seg.name = string(num2str((1:Nseg)','%03d'));
    
    seg.start_idx = ind_final(:,1);
    seg.end_idx   = ind_final(:,2);
    seg.start = DATA(i_f).time(seg.start_idx);
    seg.end   = DATA(i_f).time(seg.end_idx);
    
    SEG{i_f} = seg;
    
end

SEG = cat(1,SEG{:});


frontfields = {'flight','name','start','end'};
SEG = movevars(SEG,frontfields,'Before',1);


end


    