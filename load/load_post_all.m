
function [TURB,MOM,levels,turb_info,pvm_info] = load_post_all(datapath)

lwc_threshold = 0.02;


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
%              'DP','DT';
%              'DPlicor','tdl'};


         
% Read data files

% High-rate signals
disp('Load turbulence data ...')
[DATA,turb_info] = read_turb([datapath,filesep,'TURBULENCE'],turb_vars(:,2),turb_vars(:,1));
for i_s = 1:size(DATA,1)
    DATA(i_s).ALT = movmean(DATA(i_s).ALT,DATA(i_s).fsamp*2);
end
[PVM_DATA,pvm_info] = read_turb([datapath,filesep,'PVM_100Hz'],{'Time','LWC3'},{'time','LWC'});

% Table with cloud top and base heights
disp('Load table of average BL levels ...')

carman = readtable([datapath,filesep,'Carman2012_table.txt']);
carman.cloud_depth = carman.cloud_top-carman.cloud_base;
carman.cloud_middle = 0.5*(carman.cloud_base+carman.cloud_top);

gerber = readtable([datapath,filesep,'Gerber2013_table.txt']);
gerber.cloud_base = gerber.cloud_top - gerber.cloud_depth;
gerber.cloud_middle = 0.5*(gerber.cloud_base+gerber.cloud_top);

CT = outerjoin(carman,gerber,'Keys','flight','MergeKeys',true,...
    'LeftVariables',{'flight','cloud_base','cloud_middle','cloud_top','cloud_top_std'},...
    'RightVariables',{'cloud_base','cloud_middle','cloud_top'});
CT{:,{'cloud_base','cloud_middle','cloud_top'}} = CT{:,{'cloud_base_carman','cloud_middle_carman','cloud_top_carman'}};
CT{isnan(CT.cloud_base),{'cloud_base','cloud_middle','cloud_top'}} = CT{isnan(CT.cloud_base),{'cloud_base_gerber','cloud_middle_gerber','cloud_top_gerber'}};
CT.flight = string(CT.flight);


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
PVM  = calc_turb(SEG,PVM_DATA);

% Wind rotation from U,V to UX,VY
TURB = uv2uxvy(TURB);

% Calculate mean segment values
disp('Calculate mean values ...')
MOM = calc_mom(TURB);
Nseg = size(MOM,1);
MOM.MEAN_LWC = nan(Nseg,1);
MOM.MEDIAN_LWC = nan(Nseg,1);
MOM.CF_LWC = nan(Nseg,1);
for i_s = 1:Nseg
    MOM.MEAN_LWC(i_s)   = mean  (PVM(i_s).LWC,'omitnan');
    MOM.MEDIAN_LWC(i_s) = median(PVM(i_s).LWC,'omitnan');
    MOM.CF_LWC(i_s)     = sum(PVM(i_s).LWC>lwc_threshold)/length(PVM(i_s).LWC);
end


% Append cloud heights to MOM
for i_s = 1:size(MOM,1)
    ind_f = find(strcmp(CT.flight,MOM.flight(i_s)));
    MOM.cloud_base(i_s) = CT.cloud_base(ind_f);
    MOM.cloud_top(i_s)  = CT.cloud_top(ind_f);
    MOM.cloud_middle(i_s)  = CT.cloud_middle(ind_f);
    MOM.cloud_top_std(i_s) = CT.cloud_top_std(ind_f);
end


% Classify segments according to height
disp('Classify segments ...')
MOM.level = repmat("",size(MOM,1),1);

% MOM.level(MOM.alt<60) = "near-surface";
% MOM.level(MOM.alt>=60 & MOM.alt<MOM.cloud_base) = "sub-cloud";
% MOM.level(MOM.alt>=MOM.cloud_base & MOM.alt<MOM.cloud_middle) = "cloud-base";
% MOM.level(MOM.alt>=MOM.cloud_middle & MOM.alt<MOM.cloud_top+MOM.cloud_top_std) = "cloud-top";
% MOM.level(MOM.alt>=MOM.cloud_top+MOM.cloud_top_std) = "free-troposphere";

MOM.level(MOM.alt<60) = "near-surface";
MOM.level(MOM.alt>=60 & MOM.alt<MOM.cloud_middle & MOM.CF_LWC<0.5) = "sub-cloud";
MOM.level(MOM.alt>=60 & MOM.alt<MOM.cloud_middle & MOM.CF_LWC>=0.5) = "cloud-base";
MOM.level(MOM.alt>=MOM.cloud_middle & MOM.CF_LWC>=0.5) = "cloud-top";
MOM.level(MOM.alt>=MOM.cloud_middle & MOM.CF_LWC<0.5) = "cloud-top-uncertain";
MOM.level(MOM.alt>=MOM.cloud_top & MOM.CF_LWC<0.05) = "free-troposphere";

MOM = movevars(MOM,{'flight','level','alt','length',...
    'cloud_base','cloud_top','cloud_top_std'},'Before',1);
MOM.top = MOM.cloud_top;


% Select segments
ind_s = ismember(MOM.level,levels) & MOM.length>=20e3;
TURB = TURB(ind_s,:);
MOM  = MOM(ind_s,:);

% Label segments
MOM.name_old = MOM.name;
flights = unique(MOM.flight);
for i_f = 1:numel(flights)
    for i_l = 1:numel(levels)
        ind = MOM.flight==flights(i_f) & MOM.level==levels{i_l};
        MOM.name(ind) = num2str((1:sum(ind))',[tags{i_l},'%02d']);
    end
end

% Exclude segments
exc_seg = ["RF05","006","CT03";  % periodic oscillations in UX, drift in TAS and ALT
           "RF07","006","CT03";  % noise in velocities
           "RF04","004","CT02"]; % noise in velocities
for i_e = 1:size(exc_seg,1)
    ind_s = find(MOM.flight==exc_seg(i_e,1) & MOM.name_old==exc_seg(i_e,2));
    MOM(ind_s,:) = [];
    TURB(ind_s) = [];
end


MOM.dr = MOM.MEAN_TAS./[TURB.fsamp]';



% Inspection plot
%
% [~,~,co,~,mk] = fig16x12; close gcf
% lw = 2;
% figure, hold on, grid on
% lv = unique(MOM.level);
% for i_l = 1:numel(lv)
%     ind_l = find(MOM.level==lv(i_l));
%     scatter(MOM.start(ind_l),MOM.alt(ind_l),50,MOM.CF_LWC(ind_l),mk{i_l},'filled')
% end
% colormap copper
% cb = colorbar;
% cb.Label.String = 'cloud fraction';
% for i_f = 1:size(CT,1)
%     ind_f = find(MOM.flight==CT.flight(i_f));
%     xv = [min(MOM.start(ind_f)) max(MOM.start(ind_f))];
%     plot(xv,[1 1]*CT.cloud_base_carman(i_f),'LineWidth',lw,'Color',co(1,:))
%     plot(xv,[1 1]*CT.cloud_base_gerber(i_f),'LineWidth',lw,'Color',co(2,:))
%     plot(xv,[1 1]*CT.cloud_middle_carman(i_f),'--','LineWidth',lw,'Color',co(1,:))
%     plot(xv,[1 1]*CT.cloud_middle_gerber(i_f),'--','LineWidth',lw,'Color',co(2,:))
%     plot(xv,[1 1]*CT.cloud_top_carman(i_f),'LineWidth',lw,'Color',co(1,:))
%     plot(xv,[1 1]*CT.cloud_top_gerber(i_f),'LineWidth',lw,'Color',co(2,:))
% end
% legend([lv',"carman","gerber"])


end