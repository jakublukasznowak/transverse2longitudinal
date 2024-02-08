
function [TURB,MOM,levels,turb_info,mom_info] = load_eureca_all(datapath)


% List of levels
levels  = {'cloud-base','top-subcloud','mid-subcloud','near-surface'};

% List of flights
flight_ids = num2cell(num2str((9:19)','RF%02d'),2); % RF09 - RF19

% List of variables from turbulent moments dataset
mom_vars = {'alt','MEAN_WDIR','MEAN_TAS','MEAN_THDG'};

% List of variables from turbulent fluctuations dataset
turb_vars = {'W','W_DET';
             'UX','UX_DET';
             'VY','VY_DET'};


         
% Read data files

% Flight segmentation
disp('Load flight segmentation ...')
SEG = load_atr_seg(datapath,'v1.9','longlegs');                          
SEG = SEG(ismember(SEG.flight,flight_ids) & ismember(SEG.level,levels),:);

% Mean values and moments
disp('Load mean values and moments ...')
[MOM,mom_info] = load_atr_mom(datapath,'L3','v1.9','longlegs',mom_vars); 
MOM = join(SEG,MOM,'Keys',{'start','end'});

% Stratification table
% disp('Load table of average BL levels ...')
% CT = readtable([datapath,filesep,'stratification.xlsx'], ...
%     detectImportOptions([datapath,filesep,'stratification.xlsx'],...
%     'TextType','string','VariableNamesRange','A1:L1','DataRange','A3:L20'));
% for i_s = 1:size(MOM,1)
%     ind_f = find(strcmp(CT.Flight,MOM.flight(i_s)));
%     MOM.LCL(i_s) = CT.zLCL(ind_f);
%     MOM.inversion(i_s)  = CT.zINV(ind_f);
%     MOM.top(i_s) = CT.zSC(ind_f);
% end

% Turbulent fluctuations
disp('Load turbulence data ...')
[TURB,turb_info] = load_atr_turb(MOM,datapath,'L3','v1.9',turb_vars(:,2),turb_vars(:,1));


MOM.dr = MOM.MEAN_TAS./[TURB.fsamp]';


end