function fig = plotRegionSuscOutliers(andersonSummary, varargin)
%PLOTREGIONSUSCOUTLIERS  Horizontal bar chart of Bonferroni-significant D(->i) outliers.
%
%   fig = plotRegionSuscOutliers(andersonSummary)
%   fig = plotRegionSuscOutliers(andersonSummary, 'OutputDir', 'comparison_figures', 'Save', true)

p = inputParser;
addParameter(p, 'OutputDir', '', @(s) ischar(s) || isstring(s));
addParameter(p, 'Save', false, @islogical);
parse(p, varargin{:});
opts = p.Results;

P = cohortFigStyle();
rows = collectSuscOutlierRows(andersonSummary);
if isempty(rows)
    warning('plotRegionSuscOutliers:Empty', 'No Bonferroni-significant susceptibility outliers.');
    fig = [];
    return;
end

baseNames = {rows.baseName};
allBases = unique(baseNames, 'stable');
robustBases = {};
for b = 1:numel(allBases)
    bn = allBases{b};
    miceHere = unique({rows(strcmp(baseNames, bn)).mouseId});
    if numel(miceHere) >= 2
        robustBases{end+1} = bn; %#ok<AGROW>
    end
end
isRobust = ismember(baseNames, robustBases);
ordRobust = find(isRobust);
ordOther  = find(~isRobust);
[~, subR] = sort([rows(ordRobust).eff], 'descend');
[~, subO] = sort([rows(ordOther).eff],  'descend');
ord = [ordRobust(subR), ordOther(subO)];
rows = rows(ord);
isRobust = ismember({rows.baseName}, robustBases);
nRobust = sum(isRobust);
nRows = numel(rows);

fig = figure('Units', 'inches', 'Position', [1 1 7.5 4.94], ...
    'Color', 'w', 'Visible', 'off', 'MenuBar', 'none', 'ToolBar', 'none');
ax = axes(fig, 'Position', [0.22 0.14 0.72 0.62]);
hold(ax, 'on');
applyCohortAxesStyle(ax, P);

yPos = 1:nRows;
xMax = max([rows.eff], [], 'omitnan');
xMin = min(0, min([rows.eff], [], 'omitnan'));
xPad = 0.08 * max(xMax - xMin, 1);
xLim = [xMin - 0.02 * (xMax - xMin), xMax + xPad + 0.25 * (xMax - xMin)];

if nRobust > 0
    yRobust = yPos(isRobust);
    patch(ax, [xLim(1) xLim(2) xLim(2) xLim(1)], ...
        [min(yRobust)-0.55, min(yRobust)-0.55, max(yRobust)+0.55, max(yRobust)+0.55], ...
        P.cream, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    text(ax, xLim(1) + 0.02 * diff(xLim), max(yRobust) + 0.65, ...
        'flagged in BOTH mice (cohort-robust)', ...
        'FontWeight', 'bold', 'Color', P.red, 'FontSize', 9, ...
        'VerticalAlignment', 'bottom', 'Interpreter', 'none');
end
if nRobust < nRows
    yOther = yPos(~isRobust);
    text(ax, xLim(1) + 0.02 * diff(xLim), max(yOther) + 0.65, ...
        'one mouse only', ...
        'Color', P.mute, 'FontSize', 9, ...
        'VerticalAlignment', 'bottom', 'Interpreter', 'none');
end

hA1 = []; hA2 = [];
ytl = cell(nRows, 1);
for k = 1:nRows
    r = rows(k);
    if strcmp(r.mouseId, 'Anderson_1')
        col = P.a1;
        if isempty(hA1)
            hA1 = barh(ax, yPos(k), r.eff, 0.62, 'FaceColor', col, 'EdgeColor', 'none');
        else
            barh(ax, yPos(k), r.eff, 0.62, 'FaceColor', col, 'EdgeColor', 'none', ...
                'HandleVisibility', 'off');
        end
    else
        col = P.red;
        if isempty(hA2)
            hA2 = barh(ax, yPos(k), r.eff, 0.62, 'FaceColor', col, 'EdgeColor', 'none');
        else
            barh(ax, yPos(k), r.eff, 0.62, 'FaceColor', col, 'EdgeColor', 'none', ...
                'HandleVisibility', 'off');
        end
    end
    if ~isempty(r.hemi)
        regDisp = sprintf('%s (%s)', prettyRegionName(r.baseName), r.hemi);
    else
        regDisp = prettyRegionName(r.baseName);
    end
    ytl{k} = sprintf('%s\n%s', regDisp, r.mouseId);
    text(ax, r.eff + 0.01 * diff(xLim), yPos(k), ...
        sprintf('+%.1f%%   (z %.1f)', r.eff, r.z), ...
        'Color', P.ink, 'FontSize', 9, 'VerticalAlignment', 'middle', ...
        'Interpreter', 'none');
end

set(ax, 'YTick', yPos, 'YTickLabel', ytl, 'YDir', 'reverse', 'FontSize', 9);
xlim(ax, xLim);
ylim(ax, [0.4, nRows + 0.9]);

xlabel(ax, 'excess over control mean,  D(\rightarrow i) = (D_{And} - D_{Arn})/D_{Arn}', ...
    'Interpreter', 'tex', 'Color', P.ink);
title(ax, { ...
    'Where Anderson exceeds the control band — susceptibility D(\rightarrow i)', ...
    'column scheme, n = 2 Anderson vs 4 Arnold; Bonferroni-significant'}, ...
    'FontWeight', 'bold', 'Color', P.ink, 'Interpreter', 'tex');
ax.Subtitle.Color = P.mute;
ax.Subtitle.FontWeight = 'normal';

legH = [];
legL = {};
if ~isempty(hA1)
    legH(end+1) = hA1; legL{end+1} = 'Anderson_1'; %#ok<AGROW>
end
if ~isempty(hA2)
    legH(end+1) = hA2; legL{end+1} = 'Anderson_2'; %#ok<AGROW>
end
if ~isempty(legH)
    legend(ax, legH, legL, 'Location', 'southeast', 'Box', 'off');
end

if opts.Save
    outDir = resolveFigureOutDir(opts.OutputDir, andersonSummary);
    saveCohortFigure(fig, outDir, 'region_susc_outliers');
end
end

%% ------------------------------------------------------------------
function s = prettyRegionName(baseName)
s = strrep(char(baseName), '_', ' ');
s = lower(s);
end
