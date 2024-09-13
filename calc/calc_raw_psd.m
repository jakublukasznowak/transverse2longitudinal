
% Calculate power spectral density for the input vector X representing wind
% velocity fluctuations with spatial distance between sample points DR [in
% meters]
%
% Optional Name-Value arguments:
% 'WindowLength': window length used in the calculation of power spectrum,
% it is passed directly to PWELCH
% 'WindowOverlap': overlap of the windows used in the calculation of power
% spectrum, it is passed directly to PWELCH


function [psd,kv] = calc_raw_psd (x,dr,options)

arguments
    x (:,1) {mustBeReal, mustBeFinite, mustBeNonempty}
    dr (1,1) {mustBePositive, mustBeFinite, mustBeNonempty} % = TAS/samp
    options.WindowLength (1,1) {mustBeInteger, mustBePositive, mustBeFinite, mustBeNonempty} = floor(length(x)/2)
    options.WindowOverlap (1,1) {mustBeInteger, mustBeNonnegative, mustBeFinite, mustBeNonempty} = ceil(length(x)/4)
end


ks = 2*pi/dr;


% Calculate power spectrum

Nfft = 2^nextpow2(options.WindowLength);

[psd,kv] = pwelch(x,options.WindowLength,options.WindowOverlap,Nfft,ks);
psd = psd(kv>0); kv = kv(kv>0);


end