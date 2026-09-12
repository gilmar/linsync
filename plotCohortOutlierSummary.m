function fig = plotCohortOutlierSummary(rows, varargin)
%PLOTCOHORTOUTLIERSUMMARY  Publication-style bar chart of significant cohort findings.
%
%   fig = plotCohortOutlierSummary(rows, 'Title', 'Regions more susceptible ...')
%   fig = plotCohortOutlierSummary(rows, 'ValueField', 'value', 'XLabel', 'L - R')
%
%   Takes the flattened finding rows produced by COLLECTCOHORTOUTLIERROWS
%   or COLLECTLATERALITYOUTLIERROWS and draws one horizontal bar per
%   finding, coloured by case subject.
%
%   Findings whose region was flagged in at least 'RobustMin' different
%   case subjects are grouped at the top on a tinted band. Separating them
%   matters: with a handful of case subjects, a region flagged once is a
%   lead, while a region flagged in every case subject is a result, and a
%   single ranked list would hide that distinction.
%
%   Required row fields: subjectId, region, baseName, hemisphere, and the
%   field named by 'ValueField' (default 'effectPct').
%
%   Name-value options
%     'Title'      : figure title
%     'Subtitle'   : line under the title (e.g. cohort sizes)
%     'Footnote'   : small print under the axes
%     'XLabel'     : x-axis caption
%     'ValueField' : row field to plot (default 'effectPct')
%     'ValueFormat': 'percent' (default) | 'numeric'
%     'RobustMin'  : subjects a region must appear in to count as robust.
%                    Default 2, disabled automatically with one case subject.
%     'Style'      : struct from cohortFigStyle()
%     'Visible'    : false (default) -- build off-screen for batch export
%
%   Returns [] when `rows` is empty, so a caller can skip saving.
%
%   See also COLLECTCOHORTOUTLIERROWS, COLLECTLATERALITYOUTLIERROWS,
%   RENDERCOHORTCOMPARISONFIGURES.

p = inputParser;
addParameter(p, 'Title',       '', @(s) ischar(s) || isstring(s));
addParameter(p, 'Subtitle',    '', @(s) ischar(s) || isstring(s));
addParameter(p, 'Footnote',    '', @(s) ischar(s) || isstring(s));
addParameter(p, 'XLabel',      'effect', @(s) ischar(s) || isstring(s));
addParameter(p, 'ValueField',  'effectPct', @(s) ischar(s) || isstring(s));
addParameter(p, 'ValueFormat', 'percent', @(s) ischar(s) || isstring(s));
addParameter(p, 'RobustMin',   2, @isscalar);
addParameter(p, 'Style',       [], @(x) isempty(x) || isstruct(x));
addParameter(p, 'Visible',     false, @islogical);
parse(p, varargin{:});
opts = p.Results;

if isempty(rows)
    fig = [];
    return;
end

P = opts.Style;
if isempty(P)
    P = cohortFigStyle();
end

valueField = char(opts.ValueField);
values = [rows.(valueField)];
subjectIds = unique({rows.subjectId}, 'stable');
colors = cohortSubjectColors(sort(subjectIds), P);

%% ---- Order: robust findings first, each block by effect size ---------
robustMin = opts.RobustMin;
if numel(subjectIds) < 2
    robustMin = inf;   % robustness is meaningless with one case subject
end

baseNames = {rows.baseName};
uniqueBases = unique(baseNames, 'stable');
robustBases = {};
for b = 1:numel(uniqueBases)
    subjectsHere = unique({rows(strcmp(baseNames, uniqueBases{b})).subjectId});
    if numel(subjectsHere) >= robustMin
        robustBases{end+1} = uniqueBases{b}; %#ok<AGROW>
    end
end

isRobust = ismember(baseNames, robustBases);
idxRobust = find(isRobust);
idxOther = find(~isRobust);
[~, subR] = sort(values(idxRobust), 'descend');
[~, subO] = sort(values(idxOther), 'descend');
rows = rows([idxRobust(subR), idxOther(subO)]);

values = [rows.(valueField)];
isRobust = ismember({rows.baseName}, robustBases);
nRows = numel(rows);

%% ---- Canvas ----------------------------------------------------------
if opts.Visible
    visible = 'on';
else
    visible = 'off';
end
fig = figure('Units', 'inches', 'Position', [1 1 7.5 4.94], ...
    'Color', 'w', 'Visible', visible, 'MenuBar', 'none', 'ToolBar', 'none');
ax = axes(fig, 'Position', [0.34 0.14 0.56 0.62]);
hold(ax, 'on');
applyCohortAxesStyle(ax, P);

% Leave a blank row between the robust block and the rest, so the caption
% for the lower block cannot collide with the last bar of the upper one.
blockGap = 1.2;
yPos = 1:nRows;
if any(isRobust) && any(~isRobust) && isfinite(robustMin)
    yPos(~isRobust) = yPos(~isRobust) + blockGap;
end

% Bars grow from zero in both directions: a finding can be a deficit as
% well as an excess, and clipping the axis at zero would hide the negative
% ones entirely.
vMax = max(values, [], 'omitnan');
vMin = min(values, [], 'omitnan');
if ~isfinite(vMax); vMax = 1; end
if ~isfinite(vMin); vMin = 0; end
lo = min(0, vMin);
hi = max(0, vMax);
span = max(hi - lo, eps);
barOrigin = 0;
% Value labels sit outside the end of their bar, so each side that carries
% bars needs room for a label.
if lo < 0
    leftPad = 0.30 * span;
else
    leftPad = 0.06 * span;
end
xLim = [lo - leftPad, hi + 0.30 * span];

%% ---- Robustness band and captions ------------------------------------
if any(isRobust) && isfinite(robustMin)
    yR = yPos(isRobust);
    patch(ax, [xLim(1) xLim(2) xLim(2) xLim(1)], ...
        [min(yR) - 0.55, min(yR) - 0.55, max(yR) + 0.55, max(yR) + 0.55], ...
        P.cream, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    nRobustSubjects = numel(unique({rows(isRobust).subjectId}));
    if nRobustSubjects >= numel(subjectIds)
        bandLabel = sprintf('flagged in all %d case subjects (cohort-robust)', ...
            numel(subjectIds));
    else
        bandLabel = sprintf('flagged in %d case subjects (cohort-robust)', nRobustSubjects);
    end
    text(ax, xLim(1) + 0.02 * diff(xLim), min(yR) - 0.42, bandLabel, ...
        'FontWeight', 'bold', 'Color', P.red, 'FontSize', 9, ...
        'VerticalAlignment', 'bottom', 'Interpreter', 'none');
end
if any(~isRobust) && numel(subjectIds) > 1 && isfinite(robustMin)
    yO = yPos(~isRobust);
    text(ax, xLim(1) + 0.02 * diff(xLim), min(yO) - 0.42, ...
        'one case subject only', 'Color', P.mute, 'FontSize', 9, ...
        'VerticalAlignment', 'bottom', 'Interpreter', 'none');
end

%% ---- Bars -------------------------------------------------------------
legendHandles = [];
legendLabels = {};
seen = containers.Map('KeyType', 'char', 'ValueType', 'logical');
rowLabels = cell(nRows, 1);

for k = 1:nRows
    sid = rows(k).subjectId;
    if isKey(colors, sid)
        c = colors(sid);
    else
        c = P.grey;
    end
    v = values(k);
    if ~isfinite(v)
        v = 0;
    end

    h = barh(ax, yPos(k), v - barOrigin, 0.62, 'BaseValue', barOrigin, ...
        'FaceColor', c, 'EdgeColor', 'none');
    if ~isKey(seen, sid)
        seen(sid) = true;
        legendHandles(end+1) = h;              %#ok<AGROW>
        legendLabels{end+1} = strrep(sid, '_', ' '); %#ok<AGROW>
    else
        set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
    end

    if v >= barOrigin
        labelX = v + 0.015 * diff(xLim);
        labelAlign = 'left';
    else
        labelX = v - 0.015 * diff(xLim);
        labelAlign = 'right';
    end
    text(ax, labelX, yPos(k), formatValue(values(k), opts.ValueFormat), ...
        'FontSize', 8.5, 'Color', P.ink, 'VerticalAlignment', 'middle', ...
        'HorizontalAlignment', labelAlign, 'Interpreter', 'none');

    rowLabels{k} = rowLabel(rows(k));
end

if lo < 0
    xline(ax, 0, '-', 'Color', P.mute, 'LineWidth', 0.75, 'HandleVisibility', 'off');
end

set(ax, 'YTick', yPos, 'YTickLabel', repmat({''}, nRows, 1), ...
    'YDir', 'reverse');
ylim(ax, [0.3, max(yPos) + 0.7]);
xlim(ax, xLim);
xlabel(ax, char(opts.XLabel), 'Color', P.ink, 'FontSize', 9.5);

if ~isempty(legendHandles)
    legend(ax, legendHandles, legendLabels, 'Location', 'southeast', ...
        'Box', 'off', 'FontSize', 8.5, 'Interpreter', 'none');
end

addRowLabels(fig, ax, yPos, rowLabels, P);
addTitles(fig, opts, P);
end

%% ------------------------------------------------------------------
function label = rowLabel(row)
%ROWLABEL  "REGION NAME (L)" -- readable region plus hemisphere when known.
name = strrep(row.baseName, '_', ' ');
if isempty(row.hemisphere)
    label = name;
else
    label = sprintf('%s (%s)', name, row.hemisphere);
end
end

%% ------------------------------------------------------------------
function s = formatValue(v, valueFormat)
if ~isfinite(v)
    s = 'n/a';
    return;
end
switch lower(char(valueFormat))
    case 'percent'
        if v >= 0
            s = sprintf('+%.1f%%', v);
        else
            s = sprintf('%.1f%%', v);
        end
    otherwise
        s = sprintf('%+.3g', v);
end
end

%% ------------------------------------------------------------------
function addRowLabels(fig, ax, yPos, labels, P)
%ADDROWLABELS  Row names in figure coordinates, left of the axes.
% Drawn as annotations rather than tick labels so long region names can
% overflow the axes box without being clipped on export.
axPos = ax.Position;
labelWidth = 0.28;
gap = 0.012;
labelRight = axPos(1) - gap;
labelLeft = max(labelRight - labelWidth, 0.02);
boxHeight = min(0.09, 0.85 * axPos(4) / max(numel(yPos), 1));
yLim = ax.YLim;

for k = 1:numel(yPos)
    frac = (yPos(k) - yLim(1)) / diff(yLim);
    if strcmp(ax.YDir, 'reverse')
        frac = 1 - frac;
    end
    yNorm = axPos(2) + frac * axPos(4);
    annotation(fig, 'textbox', ...
        [labelLeft, yNorm - boxHeight / 2, labelRight - labelLeft, boxHeight], ...
        'String', labels{k}, ...
        'EdgeColor', 'none', ...
        'Color', P.ink, ...
        'FontSize', 9, ...
        'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'middle', ...
        'FitBoxToText', 'off', ...
        'Margin', 0, ...
        'Interpreter', 'none');
end
end

%% ------------------------------------------------------------------
function addTitles(fig, opts, P)
if ~isempty(char(opts.Title))
    annotation(fig, 'textbox', [0.06 0.885 0.88 0.07], ...
        'String', char(opts.Title), ...
        'EdgeColor', 'none', 'Color', P.ink, ...
        'FontSize', 12.5, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
        'Interpreter', 'none');
end
if ~isempty(char(opts.Subtitle))
    annotation(fig, 'textbox', [0.06 0.825 0.88 0.055], ...
        'String', char(opts.Subtitle), ...
        'EdgeColor', 'none', 'Color', P.mute, ...
        'FontSize', 9.5, ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
        'Interpreter', 'none');
end
if ~isempty(char(opts.Footnote))
    annotation(fig, 'textbox', [0.06 0.015 0.88 0.06], ...
        'String', char(opts.Footnote), ...
        'EdgeColor', 'none', 'Color', P.mute, ...
        'FontSize', 8, ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', ...
        'Interpreter', 'none');
end
end
