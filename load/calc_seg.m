
function SEG = calc_seg(DATA,showfigures)

if nargin<2 || isempty(showfigures)
    showfigures = false;
end


Nf = size(DATA,1);

SEG = cell(Nf,1);

for i_f = 1:Nf
    
    samp = DATA(i_f).fsamp;
    tas  = mean(DATA(i_f).TAS,'omitnan');
    
    L = length(DATA(i_f).time);
    x = (1:L)'*tas/samp/1e3;

    
    % Condition (1): small time derivative of averaged altitude
    dz = movmean(diff(DATA(i_f).ALT),floor(1e3/tas*samp));
    mask_alt = [false; (abs(dz)<0.025)];
    
    % Condition (2): small time derivative of averaged heading
    dhdg = movmean(diff(unwrap(DATA(i_f).THDG)),60*samp);
    mask_hdg = [false; abs(dhdg)<0.005];

    % Condition (3): minimum TAS 50 m/s
    mask_tas = (movmean(DATA(i_f).TAS,floor(1e3/tas*samp))>50);
    
    mask_and = mask_alt & mask_hdg & mask_tas;
    
    % Condition (4): minimum length of 10 km
    ind_and = mask2ind(mask_and);
    ind_10 = ind_and(diff(ind_and,1,2)*tas/samp>=10e3,:);
    mask_10 = ind2mask(ind_10,L);
    
    % Condition (5): small inclination within the segment
    slopes = zeros(size(ind_10,1),1);
    for i = 1:size(ind_10,1)
        p = polyfit(x(ind_10(i,1):ind_10(i,2)),DATA(i_f).ALT(ind_10(i,1):ind_10(i,2)),1);
        slopes(i) = p(1);
    end
    ind_slopes = ind_10(abs(slopes)<3,:);
    mask_slopes = ind2mask(ind_slopes,L);

    ind_final = ind_slopes;

    if showfigures
        figure, hold on, grid on
        plot(x,DATA(i_f).ALT)
        plot(x,DATA(i_f).THDG,'.')
        plot(x(mask_alt),DATA(i_f).ALT(mask_alt),'.')
        plot(x(mask_hdg),DATA(i_f).ALT(mask_hdg)+20,'.')
        plot(x(mask_10),DATA(i_f).ALT(mask_10)+40,'.')
        plot(x(mask_slopes),DATA(i_f).ALT(mask_slopes)+60,'.')
%         plot(get(gca,'XLim'),[1 1]*DATA(i_f).cloud_base,'--','Color','b')
%         plot(get(gca,'XLim'),[1 1]*DATA(i_f).cloud_top,'--','Color','b')
        title(DATA(i_f).flight)
        legend({'alt','hdg','small dz','small dhdg','>10 km','small slope'})

%         figure, hold on, grid on
%         plot(x(2:end),dz)
%         plot(x(mask_alt(2:end)),dz(mask_alt(2:end)),'.')
% 
%         figure, hold on, grid on
%         plot(x(2:end),dhdg)
%         plot(x(mask_hdg(2:end)),dhdg(mask_hdg(2:end)),'.')
    end
    
    
    seg = table;
    
    Nseg = size(ind_final,1);
    seg.flight = repmat(DATA(i_f).flight,Nseg,1);
    
    seg.start_idx = ind_final(:,1);
    seg.end_idx   = ind_final(:,2);
    seg.start = DATA(i_f).time(seg.start_idx);
    seg.end   = DATA(i_f).time(seg.end_idx);
    
    SEG{i_f} = seg;
    
end

SEG = cat(1,SEG{:});


frontfields = {'flight','start','end'};
SEG = movevars(SEG,frontfields,'Before',1);


end


    