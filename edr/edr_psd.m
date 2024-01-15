
function [edrFixed,slpFree,e,fig] = edr_psd (x,dr,fit_range,C,options)

arguments
    x (:,1) {mustBeReal, mustBeFinite, mustBeNonempty}
    dr (1,1) {mustBePositive, mustBeFinite, mustBeNonempty} % = TAS/samp
    fit_range (1,2) {mustBePositive, mustBeFinite, mustBeNonempty, mustBeValidRange(fit_range,x,dr)}
    C (1,1) {mustBeReal, mustBeFinite, mustBeNonempty} = 0.5
    options.Slope (1,1) {mustBeReal, mustBeFinite, mustBeNonempty} = -5/3
    options.Method (1,1) string {mustBeMember(options.Method,{'direct','logmean'})} = 'logmean'
    options.FitPoints (1,1) {mustBeInteger, mustBePositive, mustBeFinite, mustBeNonempty} = 20
    options.Plot (1,1) logical = false
    options.PlotXLim (1,2) {mustBePositive, mustBeFinite, mustBeNonempty, mustBeValidRange(options.PlotXLim,x,dr)} = fit_range
    options.PlotYLim (1,2) {mustBeReal, mustBeNonempty} = [-inf inf];
    options.WindowLength (1,1) {mustBeInteger, mustBePositive, mustBeFinite, mustBeNonempty} = floor(length(x)/2)
    options.WindowOverlap (1,1) {mustBeInteger, mustBeNonnegative, mustBeFinite, mustBeNonempty} = ceil(length(x)/4)
end


% Fitting range bounds

r1 = fit_range(1); 
r2 = fit_range(2);
if r2>dr*options.WindowLength
    r2 = dr*options.WindowLength;
    fprintf('Warning in EDR_PSD: Fitting range modified to [%.2f %.2f] to comply with the pwelch window length.\n',r1,r2)
end
w1 = 2*pi*dr/r2;
w2 = 2*pi*dr/r1; 


% Calculate power spectrum

Nfft = 2^nextpow2(options.WindowLength);

[psd,wv] = pwelch(x,options.WindowLength,options.WindowOverlap,Nfft);
psd = psd(wv>0); wv = wv(wv>0);


% Select fitting range

ind1 = find(wv>=w1,1,'first');
ind2 = find(wv<=w2,1,'last');
wv_fit = wv(ind1:ind2);
psd_fit = psd(ind1:ind2);


% Average PSD in log-equal bins (for LOGMEAN method)

if strcmp(options.Method,'logmean')
    [wv_fit,psd_fit] = logmean(wv_fit,psd_fit,options.FitPoints);
end


% Fit (1): fixed slope

slpFixed = options.Slope;

offsetFixed = mean(log(psd_fit)-slpFixed*log(wv_fit));
edrFixed = (exp(offsetFixed)/C)^(3/2)/dr;

e.offsetFixed = std(log(psd_fit)-slpFixed*log(wv_fit))/sqrt(length(wv_fit)); % standard error from LS fit
e.edrFixed = 3/2*edrFixed*e.offsetFixed; % error propagation


% Fit (2): free slope

[p,S] = polyfit(log(wv_fit),log(psd_fit),1);

slpFree = p(1);
offsetFree = p(2);
edrFree = (exp(offsetFree)/C)^(3/2)/dr;

covarM = (inv(S.R)*inv(S.R)')*S.normr^2/S.df; % covariance matix
e.slopeFree  = sqrt(covarM(1,1));
e.offsetFree = sqrt(covarM(2,2));
e.edrFree = 3/2*edrFree*e.offsetFree; % error propagation


% Linear correlation

corrM = corrcoef(log(wv_fit),log(psd_fit));
e.R2 = corrM(1,2);


% Diagnostic plot

if options.Plot
    
    rv = 2*pi*dr./wv;
    rv_fit = 2*pi*dr./wv_fit;
  
    [fig,~,co] = fig16x12('loglog',[1 1],'XLim',options.PlotXLim,'YLim',options.PlotYLim);
    
    plot(rv,psd,'.','Color',co(1,:))
    plot(rv_fit,psd_fit,'^','MarkerFaceColor',co(2,:))
    
    plot(rv_fit,C*(dr*edrFree).^(2/3).*wv_fit.^slpFree,'-','Color',co(4,:),'LineWidth',1)
    plot(rv_fit,C*(dr*edrFixed).^(2/3).*wv_fit.^slpFixed,'-','Color',co(5,:),'LineWidth',1)

    xlabel('$r\,[\mathrm{m}]$','Interpreter','latex')
    ylabel('$P\,[\mathrm{m^2\,s^{-2}\,rad^{-1}}]$','Interpreter','latex')
    
    legend({'power spectrum','fit points',...
        ['fit slope ',num2str(slpFree,'%.2f')],...
        ['slope ','-5/3']},...
        'Location','northwest')
%     text(0.05,0.10,['$\epsilon = ',sprintf('%.2f',edrFixed/10^floor(log10(edrFixed))),'\cdot10^',...
%         sprintf('{%d}',floor(log10(edrFixed))),'\,\mathrm{m^2\,s^{-3}}$',...
%         newline,'$R = ',sprintf('%.3f',R2),'$'],...
%         'FontSize',12,'Units','Normalized',...
%         'HorizontalAlignment','left','Interpreter','latex')
    text(0.66,0.10,['$\epsilon = ',sprintf('%.2f',edrFixed/10^floor(log10(edrFixed))),'\cdot10^',...
        sprintf('{%d}',floor(log10(edrFixed))),'\,\mathrm{m^2\,s^{-3}}$',...
        newline,'$R = ',sprintf('%.3f',abs(e.R2)),'$'],...
        'FontSize',12,'Units','Normalized',...
        'HorizontalAlignment','left','Interpreter','latex')
    
else
    fig = [];
end


end


function mustBeValidRange(a,x,dr)
    if ~ge(a(1),dr*2)
        eid = 'Range:firstTooLow';
        msg = sprintf('Fitting range must be within [dr*2 dr*length(x)] = [%.2f %.2f].',dr*2,dr*length(x));
        throwAsCaller(MException(eid,msg))
    end
    if ~le(a(2),length(x)*dr)
        eid = 'Range:lastTooHigh';
        msg = sprintf('Fitting range must be within [dr*2 dr*length(x)] = [%.2f %.2f].',dr*2,dr*length(x));
        throwAsCaller(MException(eid,msg))
    end
    if ge(a(1),a(2))
        eid = 'Range:notIncreasing';
        msg = 'Fitting range must be of nonzero length.';
        throwAsCaller(MException(eid,msg))
    end
end