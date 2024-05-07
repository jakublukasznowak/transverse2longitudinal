
function [TURB,MOM,levels,turb_info,seg_info] = load_vocals_all(datapath)


% List of levels
levels  = {'in-cloud','cloud-base','sub-cloud'};

% List of variables from turbulence dataset
turb_vars = {'time','Time';
             'ALT','ALTX';
             'TAS','TASX';
             'THDG','THDG';
%              'U','UIC';
%              'V','VIC';
             'W','WIC';
             'UX','UXC';
             'VY','VYC'};

         

% Read data files

% 25 Hz signals
disp('Load turbulence data ...')
[DATA,turb_info] = read_turb([datapath,filesep,'TURBULENCE'],turb_vars(:,2),turb_vars(:,1));

% Segmentation from auxiliary dataset
disp('Load flight segmentation ...')
[SEG,seg_info] = load_vocals_seg(datapath);
SEG = SEG(ismember(SEG.level,levels),:);


% Process

% Cut signals to the segments
disp('Cut turbulence data to the segments ...')
TURB = calc_turb(SEG,DATA); 

% Calculate mean segment values
disp('Calculate mean values ...')
MOM = calc_mom(TURB);                    
MOM.top = SEG.cloud_top;
MOM.cloud_base = SEG.cloud_base;
MOM.LCL = SEG.LCL;

% Exclude segments
exc_seg = ["RF12","C5"]; % drifting TAS
for i_e = 1:size(exc_seg,1)
    ind_s = find(MOM.flight==exc_seg(i_e,1) & MOM.name==exc_seg(i_e,2));
    MOM(ind_s,:) = [];
    TURB(ind_s) = [];
end


MOM.dr = MOM.MEAN_TAS./[TURB.fsamp]';


end