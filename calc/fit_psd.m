
% Calculate power spectral density for the input vector X representing wind
% velocity fluctuations with spatial distance between sample points DR [in
% meters] and in the range os scales defined by the 2-element vector
% FIT_RANGE [in meters] fit the two formulas:
%       (1) P(w) = O*w^{-5/3}    to get O
%       (2) P(w) = Ostar*w{-slp} to get slp
% where w is normalized angular frequency w = 2*pi*f/fs = 2*pi*dr/r where f
% is frequency, fs is sampling frequency, r is separation distance
% Plus provide error info E and plot optional diagnostic figure FIG.
%
% Optional Name-Value arguments:
% 'Method': 
%       'direct'    - evaluates power spectrum within the maximum
%                     available range of normalized frequencies W, then
%                     selects the points which belong to the FIT_RANGE and
%                     performs a linear fit in log-log coordinates.
%       'logmean'   - (default) evaluates power spectrum within the maximum available
%                     range of normalized frequencies W, then averages the
%                     values in several log-equally-distributed bins covering
%                     the FITTING_RANGE and performs a linear fit in log-log
%                     coordinates using the averaged points.
% 'FitPoints': number of fit points for the logmean method
% 'WindowLength': window length used in the calculation of power spectrum,
% it is passed directly to PWELCH
% 'WindowOverlap': overlap of the windows used in the calculation of power
% spectrum, it is passed directly to PWELCH
% 'Slope': use another fixed exponent than 2/3
% 'Plot': selects whether to show a diagnostic plot (true/false)
% 'PlotXLim': plot axes XLim [in meters]
% 'PlotYLim': plot axes YLim [in meters]


function [O,slp,e,fig] = fit_psd (x,dr,fit_range,options)

arguments
    x (:,1) {mustBeReal, mustBeFinite, mustBeNonempty}
    dr (1,1) {mustBePositive, mustBeFinite, mustBeNonempty} % = TAS/samp
    fit_range (1,2) {mustBePositive, mustBeFinite, mustBeNonempty} 
    options.Method (1,1) string {mustBeMember(options.Method,{'direct','logmean'})} = 'logmean'
    options.FitPoints (1,1) {mustBeInteger, mustBeFinite, mustBeNonempty, mustBeGreaterThan(options.FitPoints,1)} = 20
    options.WindowLength (1,1) {mustBeInteger, mustBePositive, mustBeFinite, mustBeNonempty} = floor(length(x)/2)
    options.WindowOverlap (1,1) {mustBeInteger, mustBeNonnegative, mustBeFinite, mustBeNonempty} = ceil(length(x)/4)
    options.Slope (1,1) {mustBeReal, mustBeFinite, mustBeNonempty} = -5/3
    options.Plot (1,1) logical = false
    options.PlotXLim (1,2) {mustBePositive, mustBeFinite, mustBeNonempty} = fit_range
    options.PlotYLim (1,2) {mustBeReal, mustBeNonempty} = [-inf inf];
    
end


% Check if the fit range is valid

if fit_range(1)<2*dr || fit_range(2)>dr*options.WindowLength
    if fit_range(1)<2*dr
        fit_range(1) = 2*dr;
    else
        fit_range(2) = dr*options.WindowLength;
    end
    warning('FIT_PSD:InvalidFitRange','Invalid fit range was changed to [%.2f %.2f].',fit_range(1),fit_range(2))
end


% Fit range bounds

r1 = fit_range(1); 
r2 = fit_range(2);
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
Li = length(wv_fit);

if Li<2
    throw(MException('FIT_PSD:TooFewFitPoints','Number of fit points must be at least 2.'))
end


% Average PSD in log-equal bins (for LOGMEAN method)

if strcmp(options.Method,'logmean')
    [wv_fit,psd_fit] = logmean(wv_fit,psd_fit,options.FitPoints);
    e.N = length(psd_fit);
end


% Fit (1): fixed slope

slpFixed = options.Slope;

logO = mean(log(psd_fit)-slpFixed*log(wv_fit));
e.logO = std(log(psd_fit)-slpFixed*log(wv_fit))/sqrt(length(wv_fit)); % standard error

O = exp(logO);
e.O = O*e.logO; % error propagation


% Fit (2): free slope

[p,S] = polyfit(log(wv_fit),log(psd_fit),1);

slp = p(1);
logOstar = p(2);

covarM = (inv(S.R)*inv(S.R)')*S.normr^2/S.df; % covariance matix
e.slp  = sqrt(covarM(1,1));
% e.logOstar = sqrt(covarM(2,2));

Ostar = exp(logOstar);
% e.Ostar = Ostar*e_logOstar;


% Linear correlation

corrM = corrcoef(log(wv_fit),log(psd_fit));
e.R2 = corrM(1,2);


% Diagnostic plot

if options.Plot
    
    rv = 2*pi*dr./wv;
    rv_fit = 2*pi*dr./wv_fit;
  
    [fig,~,co] = fig16x12('loglog',[1 1],'on','XLim',options.PlotXLim,'YLim',options.PlotYLim);
    
    plot(rv,psd,'.','Color',co(1,:),'MarkerSize',8)
    plot(rv_fit,psd_fit,'^','MarkerFaceColor',co(2,:),'MarkerSize',8)
    
    plot(rv_fit,Ostar*wv_fit.^slp,'-','Color',co(4,:),'LineWidth',2)
    plot(rv_fit,O*wv_fit.^slpFixed,'-','Color',co(5,:),'LineWidth',2)

    xlabel('$\lambda\,[\mathrm{m}]$','Interpreter','latex')
    ylabel('$P\,[\mathrm{m^2\,s^{-2}\,rad^{-1}}]$','Interpreter','latex')
    
    legend({'$P$','fit points',...
        ['$p=$ ',num2str(abs(slp),'%.2f')],...
        '$p=$ 5/3'},...
        'Location','northwest','Interpreter','latex')
    
else
    fig = [];
end


end