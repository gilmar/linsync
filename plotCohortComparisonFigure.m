function fig = plotCohortComparisonFigure(spec)
%PLOTCOHORTCOMPARISONFIGURE  One case subject's per-node metric against a control band.
%
%   fig = plotCohortComparisonFigure(spec)
%
%   Draws the standard per-node comparison used throughout the connectome
%   pipeline: the control cohort as a mean +- SD band with its individual
%   subjects as dots, the case subject overlaid as markers, and the nodes
%   that survived multiple-comparison correction ringed, z-labelled and
%   highlighted in the axis tick labels. Showing the individual control
%   values alongside the band matters -- with a handful of controls, a band
%   alone hides how thin the evidence for an "outlier" really is.
%
%   Nodes flagged trivial are dropped from the axis entirely: they are
%   excluded from the statistics already, and keeping them only crowds the
%   region labels.
%
%   Required spec fields
%     caseValues    : N x 1 values for the case subject
%     controlValues : N x nControl matrix of control values
%     controlMean   : N x 1 control mean
%     controlStd    : N x 1 control SD
%     z             : N x 1 z-scores
%     isOutlier     : N x 1 logical, already corrected
%     labels        : N x 1 cellstr of node labels
%
%   Optional spec fields
%     isTrivial   : N x 1 logical, nodes to omit (default: none)
%     yLabel      : TeX label for the y axis (default: 'value')
%     panelTitle  : axes title
%     figureTitle : figure-level title
%     figureName  : window name
%     caseLabel   : legend entry for the case subject (default 'case')
%     controlLabel: legend entry for the control cohort (default 'control')
%     xLabel      : caption under the tick labels (default 'Region')
%
%   See also COMPARECOHORTGROUPS, SAVEFIGUREBOTH.

spec = fillDefaults(spec);

keep = ~spec.isTrivial(:);
if ~any(keep)
    keep = true(numel(spec.caseValues), 1);
end
caseVals    = spec.caseValues(keep);
controlVals = spec.controlValues(keep, :);
ctrlMean    = spec.controlMean(keep);
ctrlStd     = spec.controlStd(keep);
zScore      = spec.z(keep);
outlierMask = spec.isOutlier(keep);
labels      = spec.labels(keep);

N = numel(caseVals);
x = 1:N;

% Legend entries go through the TeX interpreter, so an underscore in a
% subject ID would silently become a subscript.
caseLegend = texEscape(spec.caseLabel);
controlLegend = texEscape(spec.controlLabel);

figW = max(1600, min(2200, 1100 + 14 * N));
figH = max(640,  min(920,   500 +  5 * N));
fig = figure('Name', spec.figureName, 'Position', [60 60 figW figH], 'Color', 'w');
ax = axes(fig); hold(ax, 'on');

yl = [min([caseVals; ctrlMean - ctrlStd; controlVals(:)], [], 'omitnan'), ...
      max([caseVals; ctrlMean + ctrlStd; controlVals(:)], [], 'omitnan')];
if ~all(isfinite(yl)) || diff(yl) <= 0
    yl = [0 1];
end
ySpan = max(diff(yl), eps);
yl = yl + [-1 1] * 0.06 * ySpan;

% Control mean +- SD band.
valid = isfinite(ctrlMean) & isfinite(ctrlStd);
if any(valid)
    xv = x(valid);
    yLo = ctrlMean - ctrlStd;
    yHi = ctrlMean + ctrlStd;
    fill(ax, [xv fliplr(xv)], [yLo(valid).' fliplr(yHi(valid).')], ...
        [0.80 0.84 0.95], 'EdgeColor', 'none', 'FaceAlpha', 0.65, ...
        'DisplayName', sprintf('%s mean \\pm SD', controlLegend));
end

plot(ax, x, ctrlMean, '-', 'Color', [0.20 0.30 0.70], 'LineWidth', 1.2, ...
    'DisplayName', sprintf('%s mean', controlLegend));

for j = 1:size(controlVals, 2)
    plot(ax, x, controlVals(:, j), '.', 'Color', [0.55 0.60 0.85], ...
        'MarkerSize', 6, 'HandleVisibility', 'off');
end

plot(ax, x, caseVals, 'o', 'Color', [0.65 0.10 0.10], ...
    'MarkerFaceColor', [0.85 0.33 0.10], 'MarkerSize', 5, ...
    'DisplayName', caseLegend);

outlierIdx = find(outlierMask & isfinite(zScore));
if ~isempty(outlierIdx)
    peakY = max([caseVals(outlierIdx); ...
                 ctrlMean(outlierIdx) + ctrlStd(outlierIdx)], [], 'omitnan');
    if isfinite(peakY)
        yl(2) = max(yl(2), peakY + 0.12 * ySpan);
    end
end
xlim(ax, [0.5 N + 0.5]);
ylim(ax, yl);

if ~isempty(outlierIdx)
    plot(ax, x(outlierIdx), caseVals(outlierIdx), 'o', 'Color', 'k', ...
        'MarkerSize', 10, 'LineWidth', 1.5, 'HandleVisibility', 'off');
    for k = 1:numel(outlierIdx)
        i = outlierIdx(k);
        if ~isfinite(caseVals(i))
            continue;
        end
        text(ax, i, caseVals(i) + 0.025 * ySpan, sprintf('z=%+.1f', zScore(i)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
            'FontSize', 8, 'Color', [0.15 0.15 0.15], ...
            'Interpreter', 'none', 'Clipping', 'on');
    end
end

tickStyle = drawRegionTicks(ax, labels, outlierIdx);
grid(ax, 'on'); box(ax, 'on');
ylabel(ax, spec.yLabel, 'Interpreter', 'tex', 'FontSize', 10);
title(ax, spec.panelTitle, 'Interpreter', 'tex', 'FontSize', 10);
lg = legend(ax, 'Location', 'northeastoutside', 'FontSize', 8, 'Box', 'off');
sg = sgtitle(fig, spec.figureTitle, 'Interpreter', 'tex', 'FontSize', 11);
layoutAxes(ax, N, lg, sg, tickStyle);
placeXLabel(fig, ax, spec.xLabel);
end

%% ------------------------------------------------------------------
function spec = fillDefaults(spec)
required = {'caseValues', 'controlValues', 'controlMean', 'controlStd', ...
    'z', 'isOutlier', 'labels'};
for k = 1:numel(required)
    if ~isfield(spec, required{k})
        error('plotCohortComparisonFigure:MissingField', ...
            'spec.%s is required.', required{k});
    end
end
spec.caseValues = spec.caseValues(:);
spec.controlMean = spec.controlMean(:);
spec.controlStd = spec.controlStd(:);
spec.z = spec.z(:);
spec.isOutlier = logical(spec.isOutlier(:));
spec.labels = cellstr(spec.labels);
spec.labels = spec.labels(:);

defaults = { ...
    'isTrivial',    false(numel(spec.caseValues), 1); ...
    'yLabel',       'value'; ...
    'panelTitle',   ''; ...
    'figureTitle',  ''; ...
    'figureName',   'Cohort comparison'; ...
    'caseLabel',    'case'; ...
    'controlLabel', 'control'; ...
    'xLabel',       'Region'};
for k = 1:size(defaults, 1)
    if ~isfield(spec, defaults{k, 1}) || isempty(spec.(defaults{k, 1}))
        spec.(defaults{k, 1}) = defaults{k, 2};
    end
end
spec.isTrivial = logical(spec.isTrivial(:));
end

%% ------------------------------------------------------------------
function tickStyle = drawRegionTicks(ax, labels, outlierIdx)
%DRAWREGIONTICKS  Hand-placed tick labels so outliers can be bold and coloured.
% Built-in tick labels cannot be styled per tick on every MATLAB release,
% so the labels are drawn as text objects instead.
N = numel(labels);
isOutlier = false(N, 1);
isOutlier(outlierIdx) = true;

if N > 40
    tickAngle = 90;
    tickFont = 6;
else
    tickAngle = 45;
    tickFont = 7;
end

xticks(ax, 1:N);
xticklabels(ax, repmat({''}, N, 1));
ax.TickLabelInterpreter = 'none';
if isprop(ax, 'XAxis')
    ax.XAxis.TickLabelGapMultiplier = 1.8;
end

yBase = ax.YLim(1);
for i = 1:N
    if isOutlier(i)
        fw = 'bold';
        fc = [0.78 0.10 0.05];
    else
        fw = 'normal';
        fc = [0.32 0.32 0.32];
    end
    text(ax, i, yBase, labels{i}, ...
        'Rotation', tickAngle, ...
        'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'top', ...
        'FontSize', tickFont, ...
        'FontWeight', fw, ...
        'Color', fc, ...
        'Interpreter', 'none', ...
        'Clipping', 'off', ...
        'Tag', 'cohortComparisonTick');
end

tickStyle = struct('tickAngle', tickAngle, 'tickFont', tickFont, 'nRegions', N);
end

%% ------------------------------------------------------------------
function layoutAxes(ax, nRegions, legendObj, sgTitleObj, tickStyle)
%LAYOUTAXES  Leave room for rotated tick labels, titles and the legend.
drawnow;
ax.Units = 'normalized';
ti = ax.TightInset;

left = max(ti(1), 0.08);
bottom = max(ti(2), min(0.42, 0.14 + 0.004 * nRegions));
if tickStyle.tickAngle == 90
    tickPad = min(0.16, 0.04 + 0.0020 * tickStyle.nRegions);
else
    tickPad = min(0.12, 0.03 + 0.0014 * tickStyle.nRegions);
end
bottom = max(bottom, ti(2) + tickPad + 0.038);

top = max(ti(4), 0.10);
if ~isempty(sgTitleObj) && isgraphics(sgTitleObj)
    top = top + 0.05;
end

right = max(ti(3), 0.06);
if ~isempty(legendObj) && isgraphics(legendObj)
    legendObj.Location = 'northeastoutside';
    drawnow;
    right = max(right, ax.TightInset(3) + 0.12);
end

ax.Position = [left, bottom, max(0.45, 1 - left - right), max(0.30, 1 - bottom - top)];
drawnow;
end

%% ------------------------------------------------------------------
function placeXLabel(fig, ax, labelText)
%PLACEXLABEL  Caption below the hand-drawn tick labels, in figure coordinates.
% An axes xlabel would collide with the text-object ticks, so the caption
% is an annotation positioned under the lowest tick label.
drawnow;
fig.Units = 'normalized';
ax.Units = 'normalized';
if isprop(ax, 'XLabel')
    ax.XLabel.String = '';
end
delete(findobj(fig, 'Tag', 'cohortComparisonXLabel'));

gap = 0.016;
labelH = 0.026;
minYFloor = 0.012;

yPos = lowestTickY(fig, ax) - gap - labelH;
if yPos < minYFloor
    shift = minYFloor - yPos;
    pos = ax.Position;
    ax.Position = [pos(1), pos(2) + shift, pos(3), max(pos(4) - 0.5 * shift, 0.28)];
    drawnow;
    yPos = max(minYFloor, lowestTickY(fig, ax) - gap - labelH);
end

ap = ax.Position;
annotation(fig, 'textbox', [ap(1), yPos, ap(3), labelH], ...
    'String', labelText, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'FontSize', 10, ...
    'Interpreter', 'none', ...
    'Tag', 'cohortComparisonXLabel');
end

%% ------------------------------------------------------------------
function yMin = lowestTickY(fig, ax)
yMin = ax.Position(2);
ht = findobj(ax, 'Tag', 'cohortComparisonTick', 'Type', 'text');
if isempty(ht)
    return;
end
fig.Units = 'normalized';
ax.Units = 'normalized';
ap = ax.Position;
yLim = ax.YLim;
ySpan = diff(yLim);
if ySpan <= 0
    return;
end
for k = 1:numel(ht)
    ext = ht(k).Extent;
    for yd = [ext(2), ext(2) + ext(4)]
        yMin = min(yMin, ap(2) + ((yd - yLim(1)) / ySpan) * ap(4));
    end
end
end

%% ------------------------------------------------------------------
function s = texEscape(str)
s = strrep(char(str), '_', '\_');
end
