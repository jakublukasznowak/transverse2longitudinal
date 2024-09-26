% Fit the two formulas in log-log coordinates
%       (1) y = O*x^s        to get O
%       (2) y = Ostar*x^slp  to get slp
% plus provide error info E

function [O,slp,e] = fit_uni (x,y,s)

arguments
    x (:,1) {mustBeReal, mustBeFinite, mustBeNonempty}
    y (:,1) {mustBeReal, mustBeFinite, mustBeNonempty}
    s (1,1) {mustBeReal, mustBeFinite, mustBeNonempty} = 2/3
end

e.N = length(y);


% Fit (1): fixed slope

logO = mean(log(y)-s*log(x));
e.logO = std(log(y)-s*log(x))/sqrt(length(x)); % standard error

O = exp(logO);
e.O = O*e.logO; % error propagation


% Fit (2): free slope

[p,S] = polyfit(log(x),log(y),1);

slp = p(1);
% logOstar = p(2);

covarM = (inv(S.R)*inv(S.R)')*S.normr^2/S.df; % covariance matix
e.slp  = sqrt(covarM(1,1));
% e.logOstar = sqrt(covarM(2,2));

% Ostar = exp(logOstar);
% e.Ostar = Ostar*e.logOstar;


% Linear correlation

corrM = corrcoef(log(x),log(y)); % correlation matix
e.R2 = corrM(1,2);



end
