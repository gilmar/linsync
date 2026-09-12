function applyCohortAxesStyle(ax, P)
%APPLYCOHORTAXESSTYLE  Minimal axes styling for cohort summary figures.
%
%   applyCohortAxesStyle(ax)
%   applyCohortAxesStyle(ax, cohortFigStyle())
%
%   Strips the box, the tick marks and the y-axis line, and leaves only a
%   light vertical grid. For a horizontal bar chart whose categories are
%   already named in the row labels, axis furniture adds nothing and
%   competes with the data for attention.
%
%   See also COHORTFIGSTYLE, PLOTCOHORTOUTLIERSUMMARY.

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
