
function MOM = calc_mom(TURB)

ff = fieldnames(TURB)';
vars = setdiff(ff,union({'flight','name','level'},...
    ff(startsWith(ff,'fsamp') | startsWith(ff,'time')),'stable'));


Nseg = size(TURB,1);
Nvar = numel(vars);

MOM = struct2table(TURB);
MOM = MOM(:,intersect({'flight','name','level'},ff,'stable'));

for i_s = 1:Nseg
    
    MOM.start(i_s) = TURB(i_s).time(1);
    MOM.end(i_s)   = TURB(i_s).time(end);
    
    for i_v = 1:Nvar
        var = vars{i_v};
        MOM.(['MEAN_',var])(i_s) = mean(TURB(i_s).(var));
    end
    
    if isfield(TURB,'THDG')
        MOM.MEAN_THDG(i_s) = mod( mean(unwrap(TURB(i_s).THDG))+360, 360);
    end
    
end


if ismember('MEAN_TAS',MOM.Properties.VariableNames)
    MOM.length = seconds(MOM.end-MOM.start).*MOM.MEAN_TAS;
end


if all(ismember({'MEAN_U','MEAN_V'},MOM.Properties.VariableNames))
    MOM.MEAN_UU = sqrt( MOM.MEAN_U.^2 + MOM.MEAN_V.^2 );
    MOM.MEAN_WDIR = mod(90-atan2d(-MOM.MEAN_V,-MOM.MEAN_U)+360,360);

    ang = 90-MOM.MEAN_WDIR;
    MOM.UL = -( cosd(ang) .* MOM.MEAN_U + sind(ang) .* MOM.MEAN_V );
    MOM.VT = -(-sind(ang) .* MOM.MEAN_U + cosd(ang) .* MOM.MEAN_V );
end


if ismember('MEAN_ALT',MOM.Properties.VariableNames)
    MOM.alt = MOM.MEAN_ALT;
end


if all(ismember({'MEAN_THDG','MEAN_WDIR'},MOM.Properties.VariableNames))
    [MOM.dir2,MOM.dir4] = dev_angle(MOM.MEAN_THDG,MOM.MEAN_WDIR);
end


end

