
function [SEG,seg_info,DATA] = load_c130_seg(datapath)

% Level labels
level_dict = [
    "S",  "porpoise";
    "P",  "profile";
    "AC", "above-cloud";
    "SC", "sub-cloud";
    "C",  "in-cloud";
    "CB", "cloud-base";
    "BC", "cloud-base";
    "CT", "cloud-top";
    "FT", "above-cloud"];

% List of variables to read from nc files
vars = {'time','Time';
        'alt','ALTX';
        'LegType','LegType';
        'leg','leg';
        'cloud_top','cloudtop';
        'cloud_base','cloudbase';
        'cloud_base_source','cloudbaseSource';
        'LCL','LCL'};

% List of variables to average in segments
avg_vars = {'LegType','alt','cloud_top','cloud_base','LCL'};



% Read lidar files

d = dir([datapath,filesep,'LIDAR',filesep,'*.nc']);
Nf = numel(d);

DATA = cell(Nf,1);
seg_info = cell(Nf,1);

for i_f = 1:Nf
    
    fprintf('Load %s\n',d(i_f).name)
    datafile = [d(i_f).folder,filesep,d(i_f).name];

    [DATA{i_f},seg_info{i_f}] = load_nc(datafile,vars(:,2),vars(:,1));
    
    flight = textscan(d(i_f).name,'VOCALS_%4s_',1);
    DATA{i_f}.flight = string(flight{1}{1});
    
    epoch = textscan(ncreadatt(datafile,'Time','units'),'seconds since %q %q',1);
    epoch = datetime([epoch{1}{1},' ',epoch{2}{1}]);
    DATA{i_f}.time = datetime(DATA{i_f}.time,'ConvertFrom','epochtime',...
        'Epoch',epoch,'Format','yyyy-MM-dd HH:mm:ss.SS','TimeZone','UTC');  

end

DATA = cat(1,DATA{:});
seg_info = cat(1,seg_info{:});



% Read segment info

SEG = cell(Nf,1);

for i_f = 1:Nf
    
    seg = table;
    
    % Attributes of the variable LEG
    leginfo = seg_info(i_f).Variables(strcmp({seg_info(i_f).Variables.Name},'leg'));
    
    
    % Segment names from the atrribute LEGNAME of variable LEG
    leglist = textscan(leginfo.Attributes(2).Value,'%s','Delimiter',',');
    seg.name = upper(string(leglist{1})); 
    Nseg = numel(seg.name);
    
    % Segment names from the further attributes of variable LEG
    leglist2 = {leginfo.Attributes(4:2:end).Name}';
    if ~isempty(leglist2)
        seg.name2 = upper(string(leglist2));
    else
        seg.name2 = repmat("",Nseg,1);
    end
    
    % Segment start/end from variable LEG
    seg.start_idx = nan(Nseg,1);
    seg.end_idx = nan(Nseg,1);
    for i_s = 1:Nseg
        seg.start_idx(i_s) = find(DATA(i_f).leg==i_s-1,1,'first');
        seg.end_idx(i_s)   = find(DATA(i_f).leg==i_s-1,1,'last');
    end
    seg.start = DATA(i_f).time(seg.start_idx);
    seg.end   = DATA(i_f).time(seg.end_idx);
       
    % Segment start/end from the atrributes of variable LEG
    idxlist2 = cat(1,leginfo.Attributes(5:2:end).Value);
    if ~isempty(idxlist2)
        seg.start_idx2 = double(idxlist2(:,1))+1;
        seg.end_idx2   = double(idxlist2(:,2))+1;
    else
        seg.start_idx2 = nan(Nseg,1);
        seg.end_idx2   = nan(Nseg,1);
    end

    
    % Average auxiliary variables in segments
    for i_v = 1:numel(avg_vars)
        seg.(avg_vars{i_v}) = nan(Nseg,1);
        for i_s = 1:Nseg
            seg.(avg_vars{i_v})(i_s) = ...
                mean(DATA(i_f).(avg_vars{i_v})(seg.start_idx(i_s):seg.end_idx(i_s)),'omitnan');
        end
    end
    
    
    % Segment levels from translating names
    seg.level = repmat("",Nseg,1);
    for i_d = 1:size(level_dict,1)
        seg.level(startsWith(seg.name,level_dict(i_d,1))) = level_dict(i_d,2);
    end
    
    % Segment levels from variable LEGTYPE
    seg.level2 = repmat("",Nseg,1);
    seg.level2(round(seg.LegType)==0) = "sub-cloud";
    seg.level2(round(seg.LegType)==1) = "in-cloud";
    seg.level2(round(seg.LegType)==2) = "above-cloud";
    
    
    % Flight number
    seg.flight = repmat(DATA(i_f).flight,Nseg,1);
    
    
    SEG{i_f} = seg;
       
end
    
SEG = cat(1,SEG{:});


frontfields = {'flight','name','level','start','end','alt','name2','level2'};
SEG = movevars(SEG,frontfields,'Before',1);

end