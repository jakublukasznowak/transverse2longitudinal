
function [fig,ax] = plot_xy_uni(T,xfields,yfields,cfield,mfield,bifield,ifequal,blines,varargin)

mks = 6;
lw = 2;

if nargin<8
    blines = {};
end
if nargin<7 || isempty(ifequal)
    ifequal = true;
end
if nargin<6 || isempty(bifield)
    bifield = "bNull";
elseif ischar(bifield)
    bifield = string(bifield);
end
if nargin<5 || isempty(mfield)
    mfield = "mNull";
elseif ischar(mfield)
    mfield = string(mfield);
end
if nargin<4 || isempty(cfield)
    cfield = "yN";
elseif ischar(cfield)
    cfield = string(cfield);
end

    

Ny = numel(yfields);

TC = cell(1,Ny);
for i_y = 1:Ny
    TC{i_y} = T(:,intersect(horzcat(cfield,mfield,bifield),T.Properties.VariableNames));
    TC{i_y}.X = T{:,xfields(i_y)};
    TC{i_y}.Y = T{:,yfields(i_y)};
    TC{i_y}.yN(:) = i_y;
    TC{i_y}.bNull(:) = 0;
    TC{i_y}.mNull(:) = 0;
    TC{i_y}.cNull(:) = 0;
end
TC = vertcat(TC{:});

cvals = unique(TC.(cfield))';
mvals = unique(TC.(mfield))';
bvals = unique(TC.(bifield))';
if numel(bvals)>2
    disp('Warning! BIFIELD can take only 2 values.')
    bvals = bvals(1:2);
end

Nc = numel(cvals);
Nm = numel(mvals);
Nb = numel(bvals);



[fig,ax,co] = fig16x12('linlin',[1 1],'on',varargin{:});
mk = repmat({'o','^','d','s','v'},1,10);


if Nc>1
    for i_c = 1:Nc
        c = co(i_c,:);
        plot(nan,nan,'Color',c,'MarkerFaceColor',c,...
            'Marker',mk{1},'MarkerSize',mks,'LineStyle','none')
    end
end
if Nm>1
    for i_m = 1:Nm
        plot(nan,nan,'Color','black','Marker',mk{i_m},'MarkerSize',mks,'LineStyle','none')
    end
end
% if Nb>1
%     plot(nan,nan,'Color','black','MarkerFaceColor','black',...
%         'Marker',mk{1},'MarkerSize',mks,'LineStyle','none')
%     plot(nan,nan,'Color','black','MarkerFaceColor','none',...
%         'Marker',mk{1},'MarkerSize',mks,'LineStyle','none')
% end


for i_c = 1:Nc
    c = co(i_c,:);
    
    for i_m = 1:Nm
        m = mk{i_m};
        
        ind = TC.(cfield)==cvals(i_c) & TC.(mfield)==mvals(i_m) & TC.(bifield)==bvals(1);
            
        plot( TC.X(ind), TC.Y(ind),...
            'Color',c,'MarkerFaceColor',c,...
            'Marker',m,'MarkerSize',mks,...
            'LineStyle','none','HandleVisibility','off')  
        
        if Nb>1
            ind = TC.(cfield)==cvals(i_c) & TC.(mfield)==mvals(i_m) & TC.(bifield)==bvals(2);
            
            plot( TC.X(ind), TC.Y(ind),...
                'Color',c,'MarkerFaceColor','none',...
                'Marker',m,'MarkerSize',mks,...
                'LineStyle','none','HandleVisibility','off')  
        end
        
    end
    
end


% axis tight
if ifequal
    ax.XLim = [min([ax.XLim(1) ax.YLim(1)]) max([ax.XLim(2) ax.YLim(2)])];
    ax.YLim = ax.XLim;
end

xlim=ax.XLim; ylim=ax.YLim;
for i_b = 1:numel(blines)
    if startsWith(blines{i_b},'hor')
        plot(xlim,eval(erase(blines{i_b},'hor'))*[1 1],'LineWidth',lw,'LineStyle','--',...
            'Color','black','HandleVisibility','off')
    end
    if startsWith(blines{i_b},'ver')
        plot(eval(erase(blines{i_b},'ver'))*[1 1],ylim,'LineWidth',lw,'LineStyle','--',...
            'Color','black','HandleVisibility','off')
    end
    if startsWith(blines{i_b},'cross')
        plot(xlim,eval(erase(blines{i_b},'cross'))*xlim,'LineWidth',lw,'LineStyle','--',...
            'Color','black','HandleVisibility','off')
    end
end
ax.XLim=xlim; ax.YLim=ylim;


end

