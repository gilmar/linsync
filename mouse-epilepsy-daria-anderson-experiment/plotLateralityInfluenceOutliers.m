function fig = plotLateralityInfluenceOutliers(lrSummary, varargin)
%PLOTLATERALITYINFLUENCEOUTLIERS  Significant L-R influence D(k->) pairs (signed LI).

p = inputParser;
addParameter(p, 'OutputDir', '', @(s) ischar(s) || isstring(s));
addParameter(p, 'Save', false, @islogical);
addParameter(p, 'Alpha', [], @(x) isempty(x) || isscalar(x));
parse(p, varargin{:});
opts = p.Results;

info = cohortGroupInfo(lrSummary);
if isempty(opts.Alpha)
    alpha = info.alpha;
else
    alpha = opts.Alpha;
end

P = cohortFigStyle();
colors = caseMouseColors(info.caseIds, P);

rows = collectLateralityDinflRows(lrSummary, alpha);
if isempty(rows)
    warning('plotLateralityInfluenceOutliers:Empty', ...
        'No Bonferroni-significant D(k->) laterality pairs.');
    fig = [];
    return;
end

miceWithRows = unique({rows.mouseId}, 'stable');
[~, ord] = sort([rows.raw], 'descend');
rows = rows(ord);

fig = figure('Units', 'inches', 'Position', [1 1 7.5 4.94], ...
    'Color', 'w', 'Visible', 'off', 'MenuBar', 'none', 'ToolBar', 'none');
ax = axes(fig, 'Position', [0.34 0.24 0.56 0.52]);
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

legH = [];
legL = {};
legSeen = containers.Map('KeyType', 'char', 'ValueType', 'logical');
hDiamond = [];
ytl = cell(nRows, 1);
for k = 1:nRows
    r = rows(k);
    col = colors(r.mouseId);
    if ~legSeen.isKey(r.mouseId)
        bh = barh(ax, yPos(k), r.raw, 0.55, 'FaceColor', col, 'EdgeColor', 'none');
        legH(end+1) = bh; %#ok<AGROW>
        legL{end+1} = r.mouseId; %#ok<AGROW>
        legSeen(r.mouseId) = true;
    else
        barh(ax, yPos(k), r.raw, 0.55, 'FaceColor', col, 'EdgeColor', 'none', ...
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
    annot = sprintf('%+.2f   (z = %.1f)', r.raw, r.z);
    textX = r.raw + 0.03 * diff(xLim);
    if r.raw < 0
        textX = r.raw - 0.03 * diff(xLim);
        ha = 'right';
    else
        ha = 'left';
    end
    text(ax, textX, yPos(k), annot, ...
        'Color', P.ink, 'FontSize', 9, 'VerticalAlignment', 'middle', ...
        'HorizontalAlignment', ha, 'Interpreter', 'none');
end

yLim = [0.4, nRows + 0.6];
set(ax, 'YTick', yPos, 'YTickLabel', [], 'YDir', 'reverse', 'FontSize', 9);
xlim(ax, xLim);
ylim(ax, yLim);
ax.YAxis.Visible = 'off';
addCohortRowLabelAnnotations(fig, ax, yPos, ytl, P);

schemeTag = info.scheme;
if isempty(schemeTag)
    xlabSuffix = '(raw signed difference)';
else
    xlabSuffix = sprintf('(raw difference, %s)', schemeTag);
end
xlabel(ax, ['D(k\rightarrow)_L - D(k\rightarrow)_R   ' xlabSuffix], ...
    'Interpreter', 'tex', 'Color', P.ink);

title(ax, { ...
    'L–R asymmetry per region pair: influence D(k\rightarrow)', ...
    deriveLateralitySubtitle(lrSummary, info, alpha, miceWithRows)}, ...
    'FontWeight', 'bold', 'Color', P.ink, 'Interpreter', 'tex');
ax.Subtitle.Color = P.mute;
ax.Subtitle.FontWeight = 'normal';
ax.Subtitle.Interpreter = 'none';

footnote = deriveLateralityFootnote(info, miceWithRows);
if ~isempty(footnote)
    annotation(fig, 'textbox', [0.34 0.06 0.56 0.10], ...
        'String', footnote, 'EdgeColor', 'none', 'Color', P.mute, ...
        'FontAngle', 'italic', 'FontSize', 7.5, 'FontName', 'Arial', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
        'FitBoxToText', 'off', 'Interpreter', 'none');
end

if ~isempty(hDiamond)
    legH(end+1) = hDiamond; %#ok<AGROW>
    legL{end+1} = [info.controlLabel ' mean ' char(177) ' SD']; %#ok<AGROW>
end
if ~isempty(legH)
    legend(ax, legH, legL, 'Location', 'southeast', 'Box', 'off', ...
        'Interpreter', 'none');
end

if opts.Save
    outDir = resolveFigureOutDir(opts.OutputDir, lrSummary);
    saveCohortFigure(fig, outDir, 'laterality_influence_outliers');
end
end

%% ------------------------------------------------------------------
function subText = deriveLateralitySubtitle(lrSummary, info, alpha, miceWithRows)
allRows = collectLateralityDinflRows(lrSummary, alpha);
miceWithSig = unique({allRows.mouseId}, 'stable');
caseIds = info.caseIds;
noSig = setdiff(caseIds, miceWithSig, 'stable');

if numel(miceWithSig) == 1 && numel(noSig) == 1 && info.nCase > 1
    subText = sprintf('%s only — no homotopic pair survives in %s', ...
        miceWithSig{1}, noSig{1});
elseif numel(miceWithSig) == 1 && info.nCase > 1
    subText = sprintf('%s only — %d of %d case mice', ...
        miceWithSig{1}, 1, info.nCase);
elseif isempty(miceWithSig)
    subText = cohortCountSubtitle(info);
else
    subText = sprintf('%s — %s', strjoin(miceWithSig, ', '), cohortCountSubtitle(info));
end
end

%% ------------------------------------------------------------------
function footnote = deriveLateralityFootnote(info, miceWithRows)
nSig = numel(miceWithRows);
if info.nCase <= 1
    footnote = '';
elseif nSig == 1
    footnote = 'A single case-mouse signal at this sample size — suggestive, not yet replicated.';
elseif nSig < info.nCase
    footnote = sprintf('Significant in %d of %d case mice — interpret with replication in mind.', ...
        nSig, info.nCase);
else
    footnote = '';
end
end

%% ------------------------------------------------------------------
function line = cohortCountSubtitle(info)
if ~isempty(info.scheme)
    line = sprintf('%s scheme, n = %d %s vs %d %s', ...
        info.scheme, info.nCase, info.caseLabel, info.nControl, info.controlLabel);
else
    line = sprintf('n = %d %s vs %d %s', ...
        info.nCase, info.caseLabel, info.nControl, info.controlLabel);
end
end

%% ------------------------------------------------------------------
function s = prettyRegionName(baseName)
s = strrep(char(baseName), '_', ' ');
s = lower(s);
end
