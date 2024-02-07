
% Calculate 2nd order structure function for the input vector X representing
% wind velocity fluctuations with spatial distance between sample points DR 
% [in meters] in the range of scales defined by the 2-element vector FIT_RANGE
% [in meters] and fit the two formulas in log-log coordinates
%       (1) D(r) = O*r^{2/3}    to get O
%       (2) D(r) = Ostar*r^slp  to get slp
% plus provide error info E and plot optional diagnostic figure FIG.
%
% Optional Name-Value arguments:
% 'Method': 
%       'direct'    - evaluates structure function at all
%                     possible displacements which belong to the
%                     FIT_RANGE and use those values in the fit.
%       'logmean'   - (default) evaluates structure function at all possible
%                     displacements which belong to the FIT_RANGE but
%                     averages those values in several log-equally-distributed 
%                     bins covering the FIT_RANGE before the fit is performed
% 'FitPoints': number of fit points for the logmean method
% 'Slope': use another fixed exponent than 2/3
% 'Plot': selects whether to show a diagnostic plot (true/false)
% 'PlotXLim': plot axes XLim [in meters]
% 'PlotYLim': plot axes YLim [in meters]


function [O,slp,e,fig] = fit_sfc (x,dr,fit_range,options)

arguments
    x (:,1) {mustBeReal, mustBeFinite, mustBeNonempty}
    dr (1,1) {mustBePositive, mustBeFinite, mustBeNonempty} % = TAS/samp
    fit_range (1,2) {mustBePositive, mustBeFinite, mustBeNonempty, mustBeValidRange(fit_range,x,dr)}
    options.Method (1,1) string {mustBeMember(options.Method,{'direct','logmean'})} = 'logmean'
    options.FitPoints (1,1) {mustBeInteger, mustBePositive, mustBeFinite, mustBeNonempty} = 10
    options.Slope (1,1) {mustBeReal, mustBeFinite, mustBeNonempty} = 2/3
    options.Plot (1,1) logical = false
    options.PlotXLim (1,2) {mustBePositive, mustBeFinite, mustBeNonempty, mustBeValidRange(options.PlotXLim,x,dr)} = fit_range
    options.PlotYLim (1,2) {mustBeReal, mustBeNonempty} = [-inf inf];
end


% Prepare the list of displacements

iv = ( ceil(fit_range(1)/dr) : fit_range(2)/dr )';
rv = iv*dr;
Li = length(iv);
Lx = length(x);


% Calculate structure function values for the displacements from the list

sfc = nan(Li,1);
for i = 1:Li
    sfc(i) = mean( ( x(iv(i)+1:Lx) - x(1:Lx-iv(i)) ).^2 );
end


% Average SFC in log-equal bins (for LOGMEAN method)

if strcmp(options.Method,'logmean')
    [rv_fit,sfc_fit] = logmean(rv,sfc,options.FitPoints);
else
    rv_fit = rv;
    sfc_fit = sfc;
end
Li_fit = length(rv_fit);

if Li_fit<options.FitPoints
    if ~strcmp(options.Method,'direct')
        fprintf('Warning in EDR_SFC: Number of fitting points was reduced to %d.\n',Li_fit)
    end
end


% Fit (1): fixed slope

slpFixed = options.Slope;

logO = mean(log(sfc_fit)-slpFixed*log(rv_fit));
e.logO = std(log(sfc_fit)-slpFixed*log(rv_fit))/sqrt(length(rv_fit)); % standard error

O = exp(logO);
e.O = O*e.logO; % error propagation


% Fit (2): free slope

[p,S] = polyfit(log(rv_fit),log(sfc_fit),1);

slp = p(1);
logOstar = p(2);

covarM = (inv(S.R)*inv(S.R)')*S.normr^2/S.df; % covariance matix
e.slp  = sqrt(covarM(1,1));
% e.logOstar = sqrt(covarM(2,2));

Ostar = exp(logOstar);
% e.Ostar = Ostar*e_logOstar;


% Linear correlation

corrM = corrcoef(log(rv_fit),log(sfc_fit)); % correlation matix
e.R2 = corrM(1,2);


% Diagnostic plot

if options.Plot
    
    iv = ( ceil(options.PlotXLim(1)/dr) : options.PlotXLim(2)/dr )';
    rv = iv*dr;
    Li = length(iv);
    Lx = length(x);

    sfc = nan(Li,1);
    for i = 1:Li
        sfc(i) = mean( ( x(iv(i)+1:Lx) - x(1:Lx-iv(i)) ).^2 );
    end

    
    [fig,~,co] = fig16x12('loglog',[1 1],'on','XLim',options.PlotXLim,'YLim',options.PlotYLim);
    
    plot(rv,sfc,'.','Color',co(1,:),'MarkerSize',8)
    plot(rv_fit,sfc_fit,'^','MarkerFaceColor',co(2,:),'MarkerSize',8)
    
    plot(rv_fit,Ostar*rv_fit.^slp,'-','Color',co(4,:),'LineWidth',2)
    plot(rv_fit,O*rv_fit.^slpFixed,'-','Color',co(5,:),'LineWidth',2)

    xlabel('$r\,[\mathrm{m}]$','Interpreter','latex')
    ylabel('$D\,[\mathrm{m^2\,s^{-2}}]$','Interpreter','latex')
    
    legend({'$D$','fit points',...
        ['$s=$ ',num2str(slp,'%.2f')],...
        '$s=$ 2/3'},...
        'Location','northwest','Interpreter','latex')
    
else
    fig = [];
end


end


function mustBeValidRange(a,x,dr)
    if ~ge(a(1),dr)
        eid = 'Range:firstTooLow';
        msg = sprintf('Fitting range must be within [dr dr*length(x)] = [%.2f %.2f].',dr,dr*length(x));
        throwAsCaller(MException(eid,msg))
    end
    if ~le(a(2),length(x)*dr)
        eid = 'Range:lastTooHigh';
        msg = sprintf('Fitting range must be within [dr dr*length(x)] = [%.2f %.2f].',dr,dr*length(x));
        throwAsCaller(MException(eid,msg))
    end
    if ge(a(1),a(2))
        eid = 'Range:notIncreasing';
        msg = 'Fitting range must be of nonzero length.';
        throwAsCaller(MException(eid,msg))
    end
end