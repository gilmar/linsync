function summary = compareCentralityMeasures(varargin)
%COMPARECENTRALITYMEASURES  Stability centralities against classical centralities.
%
%   summary = compareCentralityMeasures('ResultsDir', resultsDir, ...
%                'Normalisation', 'column', 'CorrType', 'Spearman')
%
%   Asks whether the stability centralities D(->i) and D(k->) say anything
%   that the classical graph centralities do not. For each subject it
%   correlates every pair of node-level measures across nodes, then
%   aggregates across the cohort. A measure that correlates near 1 with
%   D(->i) is redundant with it; one that does not is carrying independent
%   structure and is worth reporting alongside.
%
%   Rank (Spearman) correlation is the default because these measures live
%   on wildly different scales and are heavy-tailed -- Pearson would be
%   dominated by a few hub nodes.
%
%   Reads the per-subject files written by RUNCOHORTSTABILITYCENTRALITIES;
%   each must carry the `centralities` struct from COMPUTENETWORKCENTRALITIES.
%
%   Name-value options
%     'ResultsDir'    : ''
%     'Normalisation' : 'column'
%     'CorrType'      : 'Spearman' | 'Pearson'
%     'Metrics'       : {} -- default: the full catalog that is present
%     'PerSubject'    : true  -- also emit one heatmap per subject
%     'Reorder'       : 'cluster' | 'fixed' -- row/column ordering of the
%                       heatmaps; 'cluster' groups measures with similar
%                       correlation profiles (needs the Statistics Toolbox,
%                       falls back to 'fixed' without it)
%     'SaveResults'   : true
%     'Plot'          : true
%     'RunParameters' : []
%
%   Outputs (under <ResultsDir>/centrality_correlation/)
%     centrality_corr_<subject>_<scheme>.{fig,png}
%     centrality_corr_mean_<scheme>.{fig,png}
%     centrality_corr_bars_<scheme>.{fig,png}
%     centrality_corr_<scheme>.mat
%
%   See also COMPUTENETWORKCENTRALITIES, CONNECTOMEMETRICCATALOG,
%   PAIRWISECORRELATION.

p = inputParser;
addParameter(p, 'ResultsDir',    '',         @(s) ischar(s) || isstring(s));
addParameter(p, 'Normalisation', 'column',   @(s) ischar(s) || isstring(s));
addParameter(p, 'CorrType',      'Spearman', @(s) ischar(s) || isstring(s));
addParameter(p, 'Metrics',       {},         @(c) iscell(c) || isstring(c));
addParameter(p, 'PerSubject',    true,  @islogical);
addParameter(p, 'Reorder',       'cluster', @(s) ischar(s) || isstring(s));
addParameter(p, 'SaveResults',   true,  @islogical);
addParameter(p, 'Plot',          true,  @islogical);
addParameter(p, 'CasePrefix',    'Case',    @(s) ischar(s) || isstring(s));
addParameter(p, 'ControlPrefix', 'Control', @(s) ischar(s) || isstring(s));
addParameter(p, 'RunParameters', [],    @(x) isempty(x) || isstruct(x));
parse(p, varargin{:});
opts = p.Results;

reorder = lower(char(opts.Reorder));
if ~ismember(reorder, {'cluster', 'fixed'})
    error('compareCentralityMeasures:BadReorder', ...
        'Reorder must be ''cluster'' or ''fixed'' (got ''%s'').', reorder);
end

scheme = normalisationSchemeInfo(opts.Normalisation).scheme;
resultsDir = resolveExperimentResultsDir(opts.ResultsDir);
outDir = experimentResultsDir(resultsDir, 'centrality_correlation');

cohort = loadCohortResults(resultsDir, scheme, ...
    'CasePrefix', opts.CasePrefix, 'ControlPrefix', opts.ControlPrefix);

%% ---- Choose the measures actually present ---------------------------
if isempty(opts.Metrics)
    catalog = connectomeMetricCatalog();
    % x0_crit is an outcome variable, not a centrality: exclude it here.
    catalog = catalog(~strcmp({catalog.field}, 'x0_crit'));
else
    catalog = connectomeMetricCatalog(cellstr(opts.Metrics));
end

nSubjects = numel(cohort.subjects);
values = cell(numel(catalog), 1);
present = false(numel(catalog), 1);
for k = 1:numel(catalog)
    values{k} = packCohortMetric(cohort, catalog(k).field);
    present(k) = any(isfinite(values{k}(:)));
end
catalog = catalog(present);
values = values(present);
nMeasures = numel(catalog);
if nMeasures < 2
    error('compareCentralityMeasures:TooFewMeasures', ...
        ['Fewer than two node measures are available in the results under %s. ' ...
         'Was the Brain Connectivity Toolbox on the path during the run?'], resultsDir);
end

measureLabels = {catalog.label};
fprintf('compareCentralityMeasures: %d subject(s), %d measure(s), %s correlation.\n', ...
    nSubjects, nMeasures, char(opts.CorrType));

%% ---- Per-subject correlation matrices -------------------------------
corrPerSubject = NaN(nMeasures, nMeasures, nSubjects);
for s = 1:nSubjects
    X = NaN(cohort.N, nMeasures);
    for k = 1:nMeasures
        X(:, k) = values{k}(:, s);
    end
    % Trivial nodes would inject a block of identical zeros and inflate
    % every correlation, so they are dropped before correlating.
    X = X(~cohort.isTrivialAny, :);
    corrPerSubject(:, :, s) = pairwiseCorrelation(X, opts.CorrType);
end

meanCorr = mean(corrPerSubject, 3, 'omitnan');
sdCorr = std(corrPerSubject, 0, 3, 'omitnan');

%% ---- Ordering --------------------------------------------------------
if strcmp(reorder, 'cluster')
    [order, orderNote] = clusterOrder(meanCorr);
else
    order = 1:nMeasures;
    orderNote = 'fixed catalog order';
end

%% ---- Figures ---------------------------------------------------------
if opts.Plot
    if opts.PerSubject
        for s = 1:nSubjects
            fig = plotCorrelationHeatmap(corrPerSubject(order, order, s), ...
                measureLabels(order), ...
                sprintf('%s  --  %s %s correlation across nodes', ...
                    strrep(cohort.subjects{s}, '_', '\_'), scheme, char(opts.CorrType)), ...
                nMeasures <= 16);
            if opts.SaveResults
                saveFigureBoth(fig, fullfile(outDir, ...
                    sprintf('centrality_corr_%s_%s', cohort.subjects{s}, scheme)));
            end
            close(fig);
        end
    end

    fig = plotCorrelationHeatmap(meanCorr(order, order), measureLabels(order), ...
        sprintf('Mean across %d subject(s)  --  %s %s correlation  [%s]', ...
            nSubjects, scheme, char(opts.CorrType), orderNote), ...
        nMeasures <= 16);
    if opts.SaveResults
        saveFigureBoth(fig, fullfile(outDir, sprintf('centrality_corr_mean_%s', scheme)));
    end
    close(fig);

    fig = plotStabilityVsClassicalBars(catalog, corrPerSubject, scheme, char(opts.CorrType));
    if ~isempty(fig)
        if opts.SaveResults
            saveFigureBoth(fig, fullfile(outDir, sprintf('centrality_corr_bars_%s', scheme)));
        end
        close(fig);
    end
end

%% ---- Summary ---------------------------------------------------------
summary = struct();
summary.scheme         = scheme;
summary.resultsDir     = resultsDir;
summary.subjects       = cohort.subjects;
summary.isCase         = cohort.isCase;
summary.isControl      = cohort.isControl;
summary.measures       = catalog;
summary.measureLabels  = measureLabels;
summary.corrType       = char(opts.CorrType);
summary.corrPerSubject = corrPerSubject;
summary.meanCorr       = meanCorr;
summary.sdCorr         = sdCorr;
summary.order          = order;
summary.orderNote      = orderNote;

patch = experimentRunParameters('fromOptions', opts, ...
    'Analysis', 'compareCentralityMeasures', 'ResultsDir', resultsDir, ...
    'Scheme', scheme, 'Subjects', cohort.subjects);
summary.runParameters = experimentRunParameters('merge', opts.RunParameters, patch);

if opts.SaveResults
    save(fullfile(outDir, sprintf('centrality_corr_%s.mat', scheme)), 'summary');
end
end

%% ------------------------------------------------------------------
function fig = plotCorrelationHeatmap(C, labels, titleText, withText)
%PLOTCORRELATIONHEATMAP  Diverging heatmap of a correlation matrix, fixed to [-1, 1].
% The colour limits are fixed so heatmaps from different subjects and
% schemes stay comparable by eye.
fig = figure('Position', [60 60 840 740], 'Color', 'w');
ax = axes(fig);
imagesc(ax, C, [-1 1]);
axis(ax, 'square');
colormap(ax, divergingColormap(256));
cb = colorbar(ax);
cb.Label.String = 'correlation';

M = size(C, 1);
set(ax, 'XTick', 1:M, 'XTickLabel', labels, ...
        'YTick', 1:M, 'YTickLabel', labels, ...
        'TickLabelInterpreter', 'tex', 'FontSize', 9);
xtickangle(ax, 45);
title(ax, titleText, 'Interpreter', 'tex', 'FontSize', 11);

if withText
    for i = 1:M
        for j = 1:M
            if ~isfinite(C(i, j))
                continue;
            end
            if abs(C(i, j)) > 0.6
                textColor = [1 1 1];
            else
                textColor = [0.1 0.1 0.1];
            end
            text(ax, j, i, sprintf('%.2f', C(i, j)), ...
                'HorizontalAlignment', 'center', 'FontSize', 7, 'Color', textColor);
        end
    end
end
end

%% ------------------------------------------------------------------
function fig = plotStabilityVsClassicalBars(catalog, corrPerSubject, scheme, corrType)
%PLOTSTABILITYVSCLASSICALBARS  Mean +- SD correlation of each stability
% measure with every classical measure, across subjects.
fields = {catalog.field};
stabilityIdx = find(ismember(fields, {'D_susceptibility', 'D_influence'}));
classicalIdx = find(~ismember(fields, {'D_susceptibility', 'D_influence'}));
if isempty(stabilityIdx) || isempty(classicalIdx)
    fig = [];
    return;
end

nStab = numel(stabilityIdx);
nClass = numel(classicalIdx);
means = NaN(nClass, nStab);
sds = NaN(nClass, nStab);
for a = 1:nStab
    for b = 1:nClass
        v = squeeze(corrPerSubject(stabilityIdx(a), classicalIdx(b), :));
        means(b, a) = mean(v, 'omitnan');
        sds(b, a) = std(v, 0, 'omitnan');
    end
end

fig = figure('Position', [80 80 max(760, 90 * nClass) 500], 'Color', 'w');
ax = axes(fig);
hold(ax, 'on');
hBar = bar(ax, means, 'grouped');
for a = 1:nStab
    hBar(a).DisplayName = catalog(stabilityIdx(a)).label;
end

% Error bars centred on each grouped bar.
for a = 1:nStab
    if isprop(hBar(a), 'XEndPoints')
        x = hBar(a).XEndPoints;
    else
        x = (1:nClass) + hBar(a).XOffset;
    end
    errorbar(ax, x, means(:, a), sds(:, a), 'k', 'linestyle', 'none', ...
        'CapSize', 6, 'HandleVisibility', 'off');
end

set(ax, 'XTick', 1:nClass, 'XTickLabel', {catalog(classicalIdx).label}, ...
    'TickLabelInterpreter', 'tex');
xtickangle(ax, 30);
ylim(ax, [-1.05 1.05]);
yline(ax, 0, '-', 'Color', [0.6 0.6 0.6], 'HandleVisibility', 'off');
ylabel(ax, sprintf('%s correlation (mean \\pm SD across subjects)', corrType), ...
    'Interpreter', 'tex');
title(ax, sprintf('Stability centralities vs classical centralities  --  %s', scheme), ...
    'Interpreter', 'tex');
legend(ax, 'Location', 'best', 'Interpreter', 'tex');
grid(ax, 'on');
box(ax, 'on');
end

%% ------------------------------------------------------------------
function cmap = divergingColormap(n)
%DIVERGINGCOLORMAP  Blue - white - red, so sign is readable at a glance.
if nargin < 1 || isempty(n)
    n = 256;
end
half = floor(n / 2);
lower = [linspace(0.10, 1, half).', linspace(0.25, 1, half).', linspace(0.60, 1, half).'];
upper = [linspace(1, 0.75, n - half).', linspace(1, 0.10, n - half).', linspace(1, 0.10, n - half).'];
cmap = [lower; upper];
end

%% ------------------------------------------------------------------
function [order, note] = clusterOrder(C)
%CLUSTERORDER  Order measures so similar correlation profiles sit together.
% Distance is 1 - correlation (signed), so anti-correlated measures end up
% far apart. Average linkage (UPGMA), with optimal leaf ordering when it is
% available. Falls back to the catalog order without the Statistics Toolbox.
M = size(C, 1);
if M < 3 || exist('linkage', 'file') ~= 2 || exist('squareform', 'file') ~= 2
    order = 1:M;
    note = 'fixed catalog order (Statistics Toolbox unavailable)';
    return;
end

Cf = C;
Cf(~isfinite(Cf)) = 0;
Cf(1:M+1:end) = 1;
D = 1 - Cf;
D(1:M+1:end) = 0;
D = (D + D.') / 2;

try
    d = squareform(D, 'tovector');
    Z = linkage(d, 'average');
    if exist('optimalleaforder', 'file') == 2
        order = optimalleaforder(Z, d);
        note = 'average linkage, optimal leaf ordering';
    else
        h = figure('Visible', 'off');
        [~, ~, order] = dendrogram(Z, 0);
        close(h);
        note = 'average linkage, dendrogram ordering';
    end
catch ME
    warning('compareCentralityMeasures:Cluster', ...
        'Hierarchical clustering failed (%s); using the catalog order.', ME.message);
    order = 1:M;
    note = 'fixed catalog order (clustering failed)';
end
order = order(:).';
end
