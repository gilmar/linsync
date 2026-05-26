function applyCohortAxesStyle(ax, P)
%APPLYCOHORTAXESSTYLE  White background, x-grid only, box off, subtle bottom spine.
if nargin < 2 || isempty(P)
    P = cohortFigStyle();
end
set(ax, 'Box', 'off', 'TickLength', [0 0], ...
    'XColor', P.mute, 'YColor', P.mute, ...
    'GridColor', P.grid, 'GridAlpha', 1, ...
    'XGrid', 'on', 'YGrid', 'off');
ax.XAxis.LineWidth = 0.5;
ax.YAxis.LineWidth = 0.5;
ax.YAxis.Color = 'none';
end
