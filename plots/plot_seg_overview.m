
function fig = plot_seg_overview(MOM,levels)

if nargin<2
    levels = unique(MOM.level);
end


[fig,ax,co] = fig16x12; % = figure('Units','normalized','Position',[0 0 0.6 0.3]);
% hold on, grid on
% co = get(gca,'ColorOrder');
fig.PaperSize = [36 10]; fig.PaperPosition = [0 0 36 10];
ax.Position = [0.07 0.14 1-1.5*0.07 1-1.5*0.14];

Nlvl = numel(levels);
for i_l = 1:Nlvl
    ind_l = (MOM.level==levels{i_l});
    plot(MOM.start(ind_l), MOM.alt(ind_l),...
        'Marker','o','MarkerSize',10,'LineStyle','none',...
        'Color',co(i_l,:),'MarkerFaceColor',co(i_l,:))
    if ismember('name',MOM.Properties.VariableNames)
        text(MOM.start(ind_l), MOM.alt(ind_l), MOM.name(ind_l))
    end
end

legend(levels,'Location','best')
ylabel('Altitude [m]')

end