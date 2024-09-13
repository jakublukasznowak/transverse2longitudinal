
function [xm,ym,xc] = logmean (x,y,Nbin,range)

% LOGMEAN averages the values from two input vectors X and Y in NBIN bins
% which are log-equally distributed across the range of X values.

if nargin<4 || isempty(range)
    range = [min(x) max(x)];
end

bin_edge = exp( linspace(log(range(1)),log(range(2)),Nbin+1)' );
bin_edge(1) = range(1);
bin_edge(end) = range(2);
bin_ind = discretize(x,bin_edge);

xc = mean([bin_edge(1:end-1) bin_edge(2:end)],2);

nnind = ~isnan(bin_ind);
xm = accumarray(bin_ind(nnind),x(nnind),[],@mean);
ym = accumarray(bin_ind(nnind),y(nnind),[],@mean);

emptyBins = ~ismember((1:Nbin)',bin_ind);
Neb = sum(emptyBins);
if Neb > 0
    warning('LOGMEAN:EmptyBins','%d bins are empty.',Neb)
    xc = xm(~emptyBins);
    xm = xm(~emptyBins);
    ym = ym(~emptyBins);
end

end