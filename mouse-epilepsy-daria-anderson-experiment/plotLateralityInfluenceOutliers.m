function fig = plotLateralityInfluenceOutliers(lrSummary, varargin)
%PLOTLATERALITYINFLUENCEOUTLIERS  Significant L-R influence D(k->) pairs (signed LI).

p = inputParser;
addParameter(p, 'OutputDir', '', @(s) ischar(s) || isstring(s));
addParameter(p, 'Save', false, @islogical);
addParameter(p, 'Alpha', [], @(x) isempty(x) || isscalar(x));
parse(p, varargin{:});
opts = p.Results;

if isempty(opts.Alpha)
    if isfield(lrSummary, 'alpha')
        alpha = lrSummary.alpha;
    else
        alpha = 0.05;
    end
else
    alpha = opts.Alpha;
end

P = cohortFigStyle();
rows = collectLateralityDinflRows(lrSummary, alpha);
if isempty(rows)
    warning('plotLateralityInfluenceOutliers:Empty', ...
        'No Bonferroni-significant D(k->) laterality pairs.');
    fig = [];
    return;
end

miceWithRows = unique({rows.mouseId});
[~, ord] = sort([rows.raw], 'descend');
rows = rows(ord);

fig = figure('Units', 'inches', 'Position', [1 1 7.5 4.94], ...
    'Color', 'w', 'Visible', 'off', 'MenuBar', 'none', 'ToolBar', 'none');
ax = axes(fig, 'Position', [0.22 0.18 0.72 0.58]);
hold(ax, 'on');
applyCohortAxesStyle(ax, P);

nRows = numel(rows);
yPos = 1:nRows;
xVals = [rows.raw];
muVals = [rows.arn_mean];
sdVals = [rows.arn_std];
xMax = max([xVals, muVals + sdVals], [], 'omitnan');
xMin = min([0, xVals, muVals - sdVals], [], 'omitnan');
xPad = 0.12 * max(xMax - xMin, 0.01);
xLim = [xMin - xPad, xMax + xPad + 0.2 * max(xMax, 0.01)];

yline(ax, 0, '--', 'Color', P.mute, 'LineWidth', 1, 'HandleVisibility', 'off');

hBar = [];
hDiamond = [];
ytl = cell(nRows, 1);
for k = 1:nRows
    r = rows(k);
    if isempty(hBar)
        hBar = barh(ax, yPos(k), r.raw, 0.55, 'FaceColor', P.red, 'EdgeColor', 'none');
    else
        barh(ax, yPos(k), r.raw, 0.55, 'FaceColor', P.red, 'EdgeColor', 'none', ...
            'HandleVisibility', 'off');
    end
    if isempty(hDiamond)
        hDiamond = errorbar(ax, r.arn_mean, yPos(k), sdVals(k), sdVals(k), ...
            'horizontal', 'Color', P.grey, 'LineStyle', 'none', ...
            'Marker', 'd', 'MarkerSize', 7, 'MarkerFaceColor', P.grey, ...
            'CapSize', 0);
    else
        errorbar(ax, r.arn_mean, yPos(k), sdVals(k), sdVals(k), ...
            'horizontal', 'Color', P.grey, 'LineStyle', 'none', ...
            'Marker', 'd', 'MarkerSize', 7, 'MarkerFaceColor', P.grey, ...
            'CapSize', 0, 'HandleVisibility', 'off');
    end
    ytl{k} = sprintf('%s\n%s', prettyRegionName(r.region), r.mouseId);
    text(ax, r.raw + 0.02 * diff(xLim), yPos(k), ...
        sprintf('+%.2f   (z = %.1f)', r.raw, r.z), ...
        'Color', P.ink, 'FontSize', 9, 'VerticalAlignment', 'middle', ...
        'Interpreter', 'none');
end

set(ax, 'YTick', yPos, 'YTickLabel', ytl, 'YDir', 'reverse', 'FontSize', 9);
xlim(ax, xLim);
ylim(ax, [0.4, nRows + 0.6]);

xlabel(ax, 'D(k\rightarrow)_L - D(k\rightarrow)_R   (raw difference, column)', ...
    'Interpreter', 'tex', 'Color', P.ink);
subText = deriveLateralitySubtitle(lrSummary, alpha, miceWithRows);
title(ax, { ...
    'L–R asymmetry per region pair: influence D(k\rightarrow)', ...
    subText}, ...
    'FontWeight', 'bold', 'Color', P.ink, 'Interpreter', 'tex');
ax.Subtitle.Color = P.mute;
ax.Subtitle.FontWeight = 'normal';
ax.Subtitle.Interpreter = 'none';

text(ax, xLim(1), 0.02, ...
    'A single-mouse signal at this sample size — suggestive, not yet replicated.', ...
    'Units', 'normalized', 'FontAngle', 'italic', 'Color', P.mute, ...
    'FontSize', 7.5, 'VerticalAlignment', 'bottom', 'Interpreter', 'none');

legEntries = {};
legHandles = [];
if ~isempty(hBar)
    legHandles(end+1) = hBar; legEntries{end+1} = miceWithRows{1}; %#ok<AGROW>
end
if ~isempty(hDiamond)
    legHandles(end+1) = hDiamond; legEntries{end+1} = 'Arnold mean \pm SD'; %#ok<AGROW>
end
if ~isempty(legHandles)
    legend(ax, legHandles, legEntries, 'Location', 'southeast', 'Box', 'off', ...
        'Interpreter', 'none');
end

if opts.Save
    outDir = resolveFigureOutDir(opts.OutputDir, lrSummary);
    saveCohortFigure(fig, outDir, 'laterality_influence_outliers');
end
end

%% ------------------------------------------------------------------
function subText = deriveLateralitySubtitle(lrSummary, alpha, miceWithRows)
andersonNames = lrSummary.mice(lrSummary.isAnderson);
rows = collectLateralityDinflRows(lrSummary, alpha);
miceWithSig = unique({rows.mouseId});
noSurvivors = setdiff(andersonNames(:)', miceWithSig, 'stable');
if numel(miceWithRows) == 1 && numel(noSurvivors) == 1
    subText = sprintf('%s only — no homotopic pair survives in %s', ...
        miceWithRows{1}, noSurvivors{1});
elseif isempty(miceWithRows)
    subText = 'column scheme, n = 2 Anderson vs 4 Arnold';
else
    subText = sprintf('%s — column scheme', strjoin(miceWithRows, ', '));
end
end

%% ------------------------------------------------------------------
function s = prettyRegionName(baseName)
s = strrep(char(baseName), '_', ' ');
s = lower(s);
end
