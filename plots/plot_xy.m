
function [fig,ax] = plot_xy(MOM,xvars,yvars,levels,coloring,ifdirs,ifequal,blines,varargin)

if nargin<8
    blines = {};
end
if nargin<7 || isempty(ifequal)
    ifequal = true;
end
if nargin<6 || isempty(ifdirs)
    ifdirs = true;
end
if nargin<5 || isempty(coloring)
    coloring = 'vars';
end
if nargin<4 || isempty(levels)
    levels = unique(MOM.level);
end


mks = 6;
lw = 1.5;

dirs   = {'along','cross'};

Nvar = numel(xvars);
Nlvl = numel(levels);


[fig,ax,co,~,mk] = fig16x12('linlin',[1 1],varargin{:});


if strcmp(coloring,'vars')
%     co = co(5:end,:);
    if Nvar>1
        for i_v = 1:Nvar
            c = co(i_v,:);
            plot(nan,nan,'Color',c,'MarkerFaceColor',c,...
                'Marker',mk{1},'MarkerSize',mks,'LineStyle','none')
        end
    end
    if Nlvl>1
        for i_l = 1:Nlvl
            plot(nan,nan,'Color','black','Marker',mk{i_l},'MarkerSize',mks,'LineStyle','none')
        end
    end
%     if ifdirs
%         plot(nan,nan,'Color','black','MarkerFaceColor','black',...
%             'Marker',mk{1},'MarkerSize',mks,'LineStyle','none')
%         plot(nan,nan,'Color','black','MarkerFaceColor','none',...
%             'Marker',mk{1},'MarkerSize',mks,'LineStyle','none')
%     end
elseif strcmp(coloring,'levels')
    if Nlvl>1
        for i_l = 1:Nlvl
            c = co(i_l,:);
            plot(nan,nan,'Color',c,'MarkerFaceColor',c,...
                'Marker',mk{1},'MarkerSize',mks,'LineStyle','none')
        end
    end
    if Nvar>1
        for i_v = 1:Nvar
            plot(nan,nan,'Color','black','Marker',mk{i_v},'MarkerSize',mks,'LineStyle','none')
        end
    end
%     if ifdirs
%         plot(nan,nan,'Color','black','MarkerFaceColor','black',...
%             'Marker',mk{1},'MarkerSize',mks,'LineStyle','none')
%         plot(nan,nan,'Color','black','MarkerFaceColor','none',...
%             'Marker',mk{1},'MarkerSize',mks,'LineStyle','none')
%     end
end


for i_l = 1:Nlvl
    for i_v = 1:Nvar
        
        if strcmp(coloring,'vars')
            c = co(i_v,:);
            m = mk{i_l};
        elseif strcmp(coloring,'levels')
            c = co(i_l,:);
            m = mk{i_v};
        end
        
        if ifdirs
            ind = MOM.level==levels{i_l} & MOM.dir2==dirs{1};
            plot( MOM.(xvars{i_v})(ind),MOM.(yvars{i_v})(ind),...
                'Color',c,'MarkerFaceColor',c,...
                'Marker',m,'MarkerSize',mks,'LineStyle','none','HandleVisibility','off')     

            ind = MOM.level==levels{i_l} & MOM.dir2==dirs{2};
            plot( MOM.(xvars{i_v})(ind),MOM.(yvars{i_v})(ind),...
                'Color',c,'MarkerFaceColor','none',...
                'Marker',m,'MarkerSize',mks,'LineStyle','none','HandleVisibility','off')
        else
            ind = (MOM.level==levels{i_l});
            plot( MOM.(xvars{i_v})(ind),MOM.(yvars{i_v})(ind),...
                'Color',c,'MarkerFaceColor',c,...
                'Marker',m,'MarkerSize',mks,'LineStyle','none','HandleVisibility','off')
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

