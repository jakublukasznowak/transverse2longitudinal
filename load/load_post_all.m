
function [TURB,MOM,levels,turb_info] = load_post_all(datapath)

% List of levels
levels  = {'cloud-top','cloud-base','sub-cloud','near-surface'};
tags = {'CT','CB','SC','S'};

% List of variables from turbulence dataset
turb_vars = {'time','Time';
             'ALT','RADALT';
             'TAS','TAS';
             'THDG','GTRK'; 
             'U','WX';
             'V','WY';
             'W','WZ'};


         
% Read data files

% High-rate signals
disp('Load turbulence data ...')
[DATA,turb_info] = read_turb([datapath,filesep,'TURBULENCE'],turb_vars(:,2),turb_vars(:,1));
for i_s = 1:size(DATA,1)
    DATA(i_s).ALT = movmean(DATA(i_s).ALT,DATA(i_s).fsamp*2);
end

% Table with cloud top and base heights
disp('Load table of average BL levels ...')
CT = readtable([datapath,filesep,'cloud_tops.txt']);



% Process

% Algorithmic detection of horizontal segments
disp('Detect horizontal segments ...')
SEG = calc_seg(DATA,'MinLength',20e3,'Plots',false,...
    'AltAvScale',2e3,'AltDrvLimit',0.012,...
    'HdgAvScale',2e3,'HdgDrvLimit',0.005,...
    'TASLimit',0.9); 

% Cut signals to the segments
disp('Cut turbulence data to the segments ...')
TURB = calc_turb(SEG,DATA);

% Wind rotation from U,V to UX,VY
TURB = uv2uxvy(TURB);

% Calculate mean segment values
disp('Calculate mean values ...')
MOM = calc_mom(TURB);

% Append cloud heights to MOM
for i_s = 1:size(MOM,1)
    ind_f = find(strcmp(CT.flight,MOM.flight(i_s)));
    MOM.cloud_base(i_s) = CT.cloud_base(ind_f);
    MOM.cloud_top(i_s)  = CT.cloud_top(ind_f);
    MOM.cloud_top_std(i_s) = CT.cloud_top_std(ind_f);
end
MOM.cloud_mid = mean([MOM.cloud_base MOM.cloud_top],2);

% Classify segments according to height
disp('Classify segments ...')
MOM.level = repmat("",size(MOM,1),1);
MOM.level(MOM.alt<60) = "near-surface";
MOM.level(MOM.alt>=60 & MOM.alt<MOM.cloud_base) = "sub-cloud";
MOM.level(MOM.alt>=MOM.cloud_base & MOM.alt<MOM.cloud_mid) = "cloud-base";
MOM.level(MOM.alt>=MOM.cloud_mid & MOM.alt<MOM.cloud_top+MOM.cloud_top_std) = "cloud-top";
MOM.level(MOM.alt>=MOM.cloud_top+MOM.cloud_top_std) = "free-troposphere";
MOM = movevars(MOM,{'flight','level','alt','length',...
    'cloud_base','cloud_top','cloud_top_std'},'Before',1);
MOM.top = MOM.cloud_top;

% Select segments
ind_s = ismember(MOM.level,levels) & MOM.length>=20e3;
TURB = TURB(ind_s,:);
MOM  = MOM(ind_s,:);

% Label segments
flights = unique(MOM.flight);
for i_f = 1:numel(flights)
    for i_l = 1:numel(levels)
        ind = MOM.flight==flights(i_f) & MOM.level==levels{i_l};
        MOM.name(ind) = num2str((1:sum(ind))',[tags{i_l},'%02d']);
    end
end

% Exclude segments
exc_seg = ["RF05","CT03";  % periodic signal oscillations
           "RF07","CT03";  % noise
           "RF04","CT02"]; % noise
for i_e = 1:size(exc_seg,1)
    ind_s = find(MOM.flight==exc_seg(i_e,1) & MOM.name==exc_seg(i_e,2));
    MOM(ind_s,:) = [];
    TURB(ind_s) = [];
end


MOM.dr = MOM.MEAN_TAS./[TURB.fsamp]';


end