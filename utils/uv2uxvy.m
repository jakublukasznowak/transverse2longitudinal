
% Rotate horizontal wind velocity component from geographical
% eastward/northward (u/v) to longitudinal/lateral (ux/vy) with respect to
% the aircraft using true heading angle


function TURB = uv2uxvy (TURB)

% avscale = 1e3; % m

Nseg = size(TURB,1);

for i_s = 1:Nseg
    
    hdg = TURB(i_s).THDG;
    
%     hdg = movmean(unwrap(hdg),avscale/mean(TURB(i_s).TAS)*TURB(i_s).fsamp);
    
    ang = 90-mod( hdg+360, 360);
    
    TURB(i_s).UX = ( cosd(ang) .* TURB(i_s).U + sind(ang) .* TURB(i_s).V );
    TURB(i_s).VY = (-sind(ang) .* TURB(i_s).U + cosd(ang) .* TURB(i_s).V );
    
end

end