
function TURB = calc_turb(SEG,DATA)

vars = setdiff(fieldnames(DATA),{'flight','fsamp','fsampvec'},'stable');


Nseg = size(SEG,1);
Nvar = numel(vars);

TURB = cell(Nseg,1);

for i_s = 1:Nseg
    
    ind_f = find([DATA(:).flight]==SEG.flight(i_s));
    
    TURB{i_s}.flight = SEG.flight(i_s);
    if ismember('name',SEG.Properties.VariableNames)
        TURB{i_s}.name   = SEG.name(i_s);
    end
    if ismember('level',SEG.Properties.VariableNames)
        TURB{i_s}.level  = SEG.level(i_s);
    end
    TURB{i_s}.fsamp  = DATA(ind_f).fsamp;
    TURB{i_s}.fsampvec = DATA(ind_f).fsampvec;
    
    
    samps = cellfun(@num2str,num2cell(DATA(ind_f).fsampvec),'UniformOutput',false);
    tvec = intersect(fieldnames(DATA)',cellfun(@(x) ['time',x],samps,'UniformOutput',false));
    Lvec = cellfun(@(x) length(DATA(ind_f).(x)),tvec);
    
    for i_v = 1:Nvar
        var = vars{i_v};
        
        ind_s = find(length(DATA(ind_f).(var))==Lvec);
        
        ind1 = find(DATA(ind_f).(['time',samps{ind_s}])>=SEG.start(i_s),1,'first');
        ind2 = find(DATA(ind_f).(['time',samps{ind_s}])<=SEG.end(i_s),1,'last');
        
        TURB{i_s}.(var) = DATA(ind_f).(var)(ind1:ind2);
    end
    
end

TURB = cat(1,TURB{:});

frontfields = intersect({'flight','name','level','fsamp','time'},fieldnames(TURB)','stable');
TURB = orderfields(TURB,cat(2,frontfields,setdiff(fieldnames(TURB)',frontfields,'stable')));


end