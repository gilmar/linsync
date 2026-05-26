function fig = plotRegionSuscOutliers(caseSummary, varargin)
%PLOTREGIONSUSCOUTLIERS  Horizontal bar chart of Bonferroni-significant D(->i) outliers.
%
%   fig = plotRegionSuscOutliers(caseSummary)
%   fig = plotRegionSuscOutliers(caseSummary, 'OutputDir', 'comparison_figures', 'Save', true)

p = inputParser;
addParameter(p, 'OutputDir', '', @(s) ischar(s) || isstring(s));
addParameter(p, 'Save', false, @islogical);
parse(p, varargin{:});
opts = p.Results;

info = cohortGroupInfo(caseSummary);
P = cohortFigStyle();
colors = caseMouseColors(info.caseIds, P);

rows = collectSuscOutlierRows(caseSummary);
if isempty(rows)
    warning('plotRegionSuscOutliers:Empty', 'No Bonferroni-significant susceptibility outliers.');
    fig = [];
    return;
end

robustMin = 2;
if info.nCase < robustMin
    robustMin = inf;
end

baseNames = {rows.baseName};
allBases = unique(baseNames, 'stable');
robustBases = {};
for b = 1:numel(allBases)
    bn = allBases{b};
    miceHere = unique({rows(strcmp(baseNames, bn)).mouseId});
    if numel(miceHere) >= robustMin
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
ax = axes(fig, 'Position', [0.34 0.14 0.56 0.62]);
hold(ax, 'on');
applyCohortAxesStyle(ax, P);

yPos = 1:nRows;
xMax = max([rows.eff], [], 'omitnan');
xMin = min(0, min([rows.eff], [], 'omitnan'));
xPad = 0.08 * max(xMax - xMin, 1);
barOrigin = max(xMin, 0);
xLim = [barOrigin, xMax + xPad + 0.25 * (xMax - xMin)];

if nRobust > 0 && isfinite(robustMin)
    yRobust = yPos(isRobust);
    patch(ax, [barOrigin xLim(2) xLim(2) barOrigin], ...
        [min(yRobust)-0.55, min(yRobust)-0.55, max(yRobust)+0.55, max(yRobust)+0.55], ...
        P.cream, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    nRobustMice = numel(unique({rows(isRobust).mouseId}));
    if nRobustMice >= info.nCase
        robustLbl = sprintf('flagged in all %d case mice (cohort-robust)', info.nCase);
    else
        robustLbl = sprintf('flagged in %d case mice (cohort-robust)', nRobustMice);
    end
    text(ax, barOrigin + 0.02 * diff(xLim), min(yRobust) - 0.42, robustLbl, ...
        'FontWeight', 'bold', 'Color', P.red, 'FontSize', 9, ...
        'VerticalAlignment', 'bottom', 'Interpreter', 'none');
end
if nRobust < nRows && info.nCase > 1
    yOther = yPos(~isRobust);
    text(ax, barOrigin + 0.02 * diff(xLim), min(yOther) - 0.42, ...
        'one case mouse only', ...
        'Color', P.mute, 'FontSize', 9, ...
        'VerticalAlignment', 'bottom', 'Interpreter', 'none');
end

legH = [];
legL = {};
legSeen = containers.Map('KeyType', 'char', 'ValueType', 'logical');
ytl = cell(nRows, 1);
for k = 1:nRows
    r = rows(k);
    col = colors(r.mouseId);
    if ~legSeen.isKey(r.mouseId)
        bh = barh(ax, yPos(k), r.eff, 0.62, 'FaceColor', col, 'EdgeColor', 'none');
        legH(end+1) = bh; %#ok<AGROW>
        legL{end+1} = r.mouseId; %#ok<AGROW>
        legSeen(r.mouseId) = true;
    else
        barh(ax, yPos(k), r.eff, 0.62, 'FaceColor', col, 'EdgeColor', 'none', ...
            'HandleVisibility', 'off');
    end
    if ~isempty(r.hemi)
        regDisp = sprintf('%s (%s)', prettyRegionName(r.baseName), r.hemi);
    else
        regDisp = prettyRegionName(r.baseName);
    end
    ytl{k} = sprintf('%s\n%s', regDisp, r.mouseId);
    textX = max(r.eff, 0) + 0.02 * diff(xLim);
    text(ax, textX, yPos(k), ...
        sprintf('%s   (z %.1f)', formatPercentEffect(r.eff), r.z), ...
        'Color', P.ink, 'FontSize', 9, 'VerticalAlignment', 'middle', ...
        'HorizontalAlignment', 'left', 'Interpreter', 'none');
end

yLim = [0.15, nRows + 0.9];
set(ax, 'YTick', yPos, 'YTickLabel', [], 'YDir', 'reverse', 'FontSize', 9);
xlim(ax, xLim);
ylim(ax, yLim);
ax.YAxis.Visible = 'off';
addCohortRowLabelAnnotations(fig, ax, yPos, ytl, P);

xlabel(ax, sprintf(['excess over %s mean,  D(\\rightarrow i) = ' ...
    '(D_{%s} - D_{%s})/D_{%s}'], info.controlLabel, ...
    info.caseTex, info.controlTex, info.controlTex), ...
    'Interpreter', 'tex', 'Color', P.ink);
title(ax, { ...
    sprintf('Where %s exceeds the %s band — susceptibility D(\\rightarrow i)', ...
        info.caseLabel, info.controlLabel), ...
    cohortCountSubtitle(info)}, ...
    'FontWeight', 'bold', 'Color', P.ink, 'Interpreter', 'tex');
ax.Subtitle.Color = P.mute;
ax.Subtitle.FontWeight = 'normal';

if ~isempty(legH)
    legend(ax, legH, legL, 'Location', 'southeast', 'Box', 'off', 'Interpreter', 'none');
end

if opts.Save
    outDir = resolveFigureOutDir(opts.OutputDir, caseSummary);
    saveCohortFigure(fig, outDir, 'region_susc_outliers');
end
end

%% ------------------------------------------------------------------
function s = prettyRegionName(baseName)
s = strrep(char(baseName), '_', ' ');
s = lower(s);
end

%% ------------------------------------------------------------------
function line = cohortCountSubtitle(info)
if ~isempty(info.scheme)
    line = sprintf('%s scheme, n = %d %s vs %d %s; Bonferroni-significant', ...
        info.scheme, info.nCase, info.caseLabel, info.nControl, info.controlLabel);
else
    line = sprintf('n = %d %s vs %d %s; Bonferroni-significant', ...
        info.nCase, info.caseLabel, info.nControl, info.controlLabel);
end
end
