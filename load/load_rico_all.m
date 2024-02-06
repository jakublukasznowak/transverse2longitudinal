
function [TURB,MOM,levels,turb_info] = load_rico_all(datapath)


% List of levels
levels  = {'cloud-layer','cloud-base','sub-cloud','near-surface'};
tags = {'CL','CB','SC','S'};

% List of variables from turbulence dataset
turb_vars = {'time','Time';
             'P','PSXC';
             'ALT','GGALTC';
             'TAS','TASX';
             'THDG','THDG';
             'U','UIC';
             'V','VIC';
             'W','WIC';
             'UX','UXC';
             'VY','VYC'};

         

% Read data files

disp('Load turbulence data ...')
[DATA,turb_info] = read_turb([datapath,filesep,'TURBULENCE'],turb_vars(:,2),turb_vars(:,1));



% Process

% Algorithmic detection of horizontal segments
disp('Detect horizontal segments ...')
SEG  = calc_seg(DATA,'MinLength',30e3,'Plots',false,'Title','C130-RICO',...
    'AltAvScale',4e3,'AltDrvLimit',0.01,...
    'HdgAvScale',20e3,'HdgDrvLimit',0.003);

% Cut signals to the segments
disp('Cut turbulence data to the segments ...')
TURB = calc_turb(SEG,DATA);   

% Calculate mean segment values
disp('Calculate mean values ...')
MOM  = calc_mom(TURB);                       

% Classify segments according to height 
disp('Classify segments ...')
MOM.level = repmat("",size(MOM,1),1);
MOM.level(MOM.MEAN_P>=990) = "near-surface";
MOM.level(MOM.MEAN_P<990 & MOM.MEAN_P>=950) = "sub-cloud";
MOM.level(MOM.MEAN_P<950 & MOM.MEAN_P>=900) = "cloud-base";
MOM.level(MOM.MEAN_P<900 & MOM.MEAN_P>=800) = "cloud-layer";
MOM = movevars(MOM,{'flight','level','alt','length'},'Before',1);
MOM.top = repmat(1000,size(MOM,1),1);

% Select segments
ind_s = ismember(MOM.level,levels);
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
exc_seg = ["RF11","S02";    % spikes
           "RF16","SC09";   % strange maneuvers
           "RF03","SC03"];  % strange maneuvers
for i_e = 1:size(exc_seg,1)
    ind_s = find(MOM.flight==exc_seg(i_e,1) & MOM.name==exc_seg(i_e,2));
    MOM(ind_s,:) = [];
    TURB(ind_s) = [];
end


MOM.dr = MOM.MEAN_TAS./[TURB.fsamp]';


end