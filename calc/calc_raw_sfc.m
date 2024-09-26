
% Calculate 2nd order structure function for the input vector X representing
% wind velocity fluctuations with spatial distance between sample points DR 
% [in meters] in the range of scales defined by the 2-element vector RANGE
% [in meters]

function [sfc,rv] = calc_raw_sfc (x,dr,range)

arguments
    x (:,1) {mustBeReal, mustBeFinite, mustBeNonempty}
    dr (1,1) {mustBePositive, mustBeFinite, mustBeNonempty} % = TAS/samp
    range (1,2) {mustBePositive, mustBeFinite, mustBeNonempty}
end


Lx = length(x);

% Check if the range is valid

if range(1)<dr || range(2)>Lx*dr
    if range(1)<dr
        range(1) = dr;
    else
        range(2) = Lx*dr;
    end
    warning('CALC_RAW_SFC:InvalidRange','Invalid range was changed to [%.2f %.2f].',range(1),range(2))
end


% Prepare the list of displacements

iv = ( floor(range(1)/dr) : ceil(range(2)/dr) )';
% iv = ( ceil(range(1)/dr) : range(2)/dr )';
rv = iv*dr;
Li = length(iv);


% Calculate structure function values for the displacements from the list

sfc = nan(Li,1);
for i = 1:Li
    sfc(i) = mean( ( x(iv(i)+1:Lx) - x(1:Lx-iv(i)) ).^2 );
end


end
