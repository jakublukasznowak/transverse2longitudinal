
function h = plot_whisker(MOM,xvars,levels,ifdirs,varargin)

if nargin<4 || isempty(ifdirs)
    ifdirs = true;
end
if nargin<3 || isempty(levels)
    levels = unique(MOM.level);
end


dirs   = {'along','cross'};

Nvar = numel(xvars);
Nlvl = numel(levels);
Nseg = size(MOM,1);

[~,~,co] = fig16x12;
close gcf


if ifdirs
    bginput = repmat({nan(Nseg,Nlvl)},1,Nvar*2);
    co = reshape(repmat(co,1,2)',3,[])';
else
    bginput = repmat({nan(Nseg,Nlvl)},1,Nvar);
end

for i_v = 1:Nvar
    for i_l = 1:Nlvl
        if ifdirs
            ind = MOM.level==levels{i_l} & MOM.dir2==dirs{1};
            bginput{2*i_v-1}(ind,i_l) = MOM.(xvars{i_v})(ind);
            ind = MOM.level==levels{i_l} & MOM.dir2==dirs{2};
            bginput{2*i_v}(ind,i_l) = MOM.(xvars{i_v})(ind);
        else
            ind = MOM.level==levels{i_l};
            bginput{i_v}(ind,i_l) = MOM.(xvars{i_v})(ind);
        end
    end
end


h = boxplotGroup(bginput,'Whisker',Inf,...
    'SecondaryLabels',cellfun(@(x) [x,'   '],levels,'UniformOutput',false),...
    'InterGroupSpace',2,'groupLabelType','vertical',...
    'Colors',co,'GroupType','betweenGroups','Widths',0.8,...
    varargin{:});
if ifdirs
    for i_g = 1:2:numel(h.boxplotGroup)
        ind_box = find(strcmp({h.boxplotGroup(i_g).Children(:).Tag},'Box'));
        for i_b = 1:numel(ind_box)
            h.boxplotGroup(i_g).Children(ind_box(i_b)).LineStyle = '--';
        end
    end
    for i_g = 2:2:numel(h.boxplotGroup)
        ind_wsk = find(ismember({h.boxplotGroup(i_g).Children(:).Tag},{'Upper Whisker','Lower Whisker'}));          
        for i_w = 1:numel(ind_wsk)
            h.boxplotGroup(i_g).Children(ind_wsk(i_w)).LineStyle = '-';
        end
    end
end
h.axis.TickLabelInterpreter = 'latex';
h.axis.View = [90 90];
h.axis.Clipping = 'off';

end