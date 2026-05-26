function addCohortRowLabelAnnotations(fig, ax, yPos, ytl, P)
%ADDCOHORTROWLABELANNOTATIONS  Row labels in figure space, left of the data axes.
%
%   Keeps region / mouse labels out of the bar plot and survives PNG export.

if nargin < 5 || isempty(P)
    P = cohortFigStyle();
end

axPos = ax.Position;
labelW = 0.28;
gap = 0.012;
labelRight = axPos(1) - gap;
labelLeft = max(labelRight - labelW, 0.02);
boxH = min(0.09, 0.85 * axPos(4) / max(numel(yPos), 1));
yLim = ax.YLim;

for k = 1:numel(yPos)
    frac = (yPos(k) - yLim(1)) / diff(yLim);
    if isequal(ax.YDir, 'reverse')
        frac = 1 - frac;
    end
    yNorm = axPos(2) + frac * axPos(4);
    annotation(fig, 'textbox', ...
        [labelLeft, yNorm - boxH / 2, labelRight - labelLeft, boxH], ...
        'String', ytl{k}, ...
        'EdgeColor', 'none', ...
        'Color', P.ink, ...
        'FontSize', 9, ...
        'FontName', 'Arial', ...
        'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'middle', ...
        'FitBoxToText', 'off', ...
        'Margin', 0, ...
        'Interpreter', 'none');
end
end
