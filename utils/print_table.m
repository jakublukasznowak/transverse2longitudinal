
% Print latex-formatted summary table of statistics for levels


function print_table(TAB,groupvars,vars,ifstd,ifcount,prec)

stats = ["mean","std"];

if nargin<6 || isempty(prec)
    prec = 2;
end
if nargin<5 || isempty(ifcount)
    ifcount = false;
end
if nargin<4 || isempty(ifstd)
    ifstd = true;
end

if iscell(groupvars)
    groupvars = string(groupvars);
end
if iscell(vars)
    vars = string(vars);
end


if ismember('length',TAB.Properties.VariableNames)
    TAB.length_km = TAB.length/1000;
end
% MOM.vstop = MOM.top - MOM.alt;


if ismember('alt',TAB.Properties.VariableNames)
    G_level = groupsummary(TAB,groupvars,stats,union(vars,{'alt'}));
    G_level = sortrows(G_level,"mean_alt",'descend');
    G_total = groupsummary(TAB,[],stats,union(vars,{'alt'}));
else
    G_level = groupsummary(TAB,groupvars,stats,vars);
    G_total = groupsummary(TAB,[],stats);
end

for i_gv = 1:numel(groupvars)
    G_total{:,groupvars(i_gv)} = "all";
end
G = [G_level;G_total];


Ns = 15;

fprintf(' %*s',Ns,groupvars(1))
for i_gv = 2:numel(groupvars)
    fprintf(' & %*s',Ns,groupvars(i_gv))
end
if ifcount
    fprintf(' & %5s','N')
end
for i_v = 1:numel(vars)
    fprintf(' & %*s',Ns,vars(i_v))
end
fprintf(' \\\\ \n')

for i_g = 1:size(G,1)
    fprintf(' %*s',Ns,G{i_g,groupvars(1)})
    for i_gv = 2:numel(groupvars)
        fprintf(' & %*s',Ns,G{i_g,groupvars(i_gv)})
    end
    if ifcount
        fprintf(' & %5d',G.GroupCount(i_g))
    end
    for i_v = 1:numel(vars)
        if ifstd
            fprintf(' & %*.*f (%.*f)',Ns-prec-5,prec,G{i_g,"mean_"+vars(i_v)},...
                prec,G{i_g,"std_"+vars(i_v)})
        else
            fprintf(' & %*.*f',Ns,prec,G{i_g,"mean_"+vars(i_v)})
        end
    end
    fprintf(' \\\\ \n')
end


end