
function [DATA,info] = read_turb(datapath,vars,varsnew)

if nargin<3
    varsnew = {};
end
if nargin<2
    vars = {};
end
if isempty(varsnew)
    varsnew = vars;
end
if ~isempty(vars) && ~ismember('Time',vars)
    vars    = union({'Time'},vars,'stable');
    varsnew = union({'time'},varsnew,'stable');
end


d = dir([datapath,filesep,'*.nc']);

Nf = size(d,1);
Nvar = numel(varsnew);

DATA = cell(Nf,1);
info = cell(Nf,1);


for i_f = 1:Nf
    
%     fprintf('Load %s\n',d(i_f).name)
    datafile = [d(i_f).folder,filesep,d(i_f).name];

    [DATA{i_f},info{i_f}] = load_nc(datafile,vars,varsnew);
    
    
    % Flight number
    
%     flight = textscan(d(i_f).name,'%4s',1);
%     DATA{i_f}.flight = string(flight{1}{1});
    irf = strfind(d(i_f).name,'RF');
    flight = d(i_f).name(irf:irf+3);
    DATA{i_f}.flight = string(flight);
    
    
    % Time
    
    epoch = textscan(ncreadatt(datafile,'Time','units'),'seconds since %q %q',1);
    epoch = datetime([epoch{1}{1},' ',epoch{2}{1}]);
    
    dims = {info{i_f}.Dimensions.Name};
    samps = erase(dims(startsWith(dims,'sps')),'sps');
    DATA{i_f}.fsampvec = [info{i_f}.Dimensions(startsWith(dims,'sps')).Length];
    
    sectime = double(DATA{i_f}.time');
    
    Nsamp = numel(samps);
    for i_s = 1:Nsamp
        [X,Y] = meshgrid(sectime,(0:1/DATA{i_f}.fsampvec(i_s):1-1e-5)');
        DATA{i_f}.(['time',samps{i_s}]) = datetime(X(:)+Y(:),'ConvertFrom','epochtime',...
            'Epoch',epoch,'Format','yyyy-MM-dd HH:mm:ss.SS','TimeZone','UTC');
    end
    
    [~,i_s_max] = max(DATA{i_f}.fsampvec);
    DATA{i_f}.fsamp = DATA{i_f}.fsampvec(i_s_max);
    DATA{i_f}.time = DATA{i_f}.(['time',samps{i_s_max}]);
    
    
    % Reshape variables

    for i_v = 1:Nvar
        var = varsnew{i_v};
        if isfield(DATA{i_f},var)
            DATA{i_f}.(var) = DATA{i_f}.(var)(:);
        end
    end

end


% Concatenate structures

fnames = cellfun(@fieldnames,DATA,'UniformOutput',false);
flist  = unique(cat(1,fnames{:}));
im = cellfun(@(x) ismember(flist,x)',fnames,'UniformOutput',false);
flist = flist(all(cat(1,im{:}),1));
DATA = cellfun(@(x) rmfield(x,setdiff(fieldnames(x),flist,'stable')),...
    DATA,'UniformOutput',false);

DATA = cat(1,DATA{:});
info = cat(1,info{:});


end