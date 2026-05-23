function summary = compareCentralityMeasures(varargin)
%COMPARECENTRALITYMEASURES  Compare the stability-based centralities
% D(->i) and D(k->) with the classical graph-theoretic centralities
% computed by computeNetworkCentralities, in the style of the
% correlation analyses in Liao & Lizier (2026), Section IIIC and
% Appendix H.
%
% For each mouse, computes the cross-node Spearman (rank) correlation
% between every pair of centralities. Produces three figure types:
%
%   1) Per-mouse correlation heatmap   (Fig 3 / Fig 10 of the paper)
%   2) Mean correlation heatmap across mice  (aggregate of (1))
%   3) Bar chart -- mean +- SD correlation of D(->i) and D(k->) with
%      each classical centrality across mice  (Fig 4 style)
%
% Reads stabilityCentralities_<mouseId>_<scheme>_results.mat under results/. Each
% file must contain the `results.centralities` field added by
% runMouseStabilityCentralities (re-run runAllMiceStabilityCentralities if missing).
%
% Outputs (under results/):
%   centrality_corr_<mouseId>_<scheme>.{fig,png}   -- per-mouse heatmap
%   centrality_corr_mean_<scheme>.{fig,png}        -- mean heatmap
%   centrality_corr_bars_<scheme>.{fig,png}        -- bar chart
%   centrality_corr_<scheme>.mat                   -- packed summary
%
% Usage:
%   compareCentralityMeasures;
%   compareCentralityMeasures('Normalisation', 'tvb');
%   compareCentralityMeasures('CorrType', 'Pearson', 'PerMouse', false);
%
% Name-value options:
%   Normalisation : 'parkes' (default) -- which scheme to load
%   CorrType      : 'Spearman' (default) | 'Pearson' | 'Kendall'
%   PerMouse      : true (default)  -- emit per-mouse heatmap figures
%   SaveResults   : true (default)  -- write figures + summary mat
%   Reorder       : 'cluster' (default) | 'fixed'
%                    'cluster' applies hierarchical clustering on the
%                    mean correlation matrix to reorder rows/columns so
%                    measures with similar correlation profiles sit next
%                    to each other (matches Fig 12/13 of the paper).
%                    'fixed' keeps the canonical ordering from
%                    defineMeasures.

setupMousePaths();

p = inputParser;
addParameter(p, 'Normalisation', 'column',   @(s) ischar(s) || isstring(s));
addParameter(p, 'CorrType',      'Spearman', @(s) ischar(s) || isstring(s));
addParameter(p, 'PerMouse',      true,       @islogical);
addParameter(p, 'SaveResults',   true,       @islogical);
addParameter(p, 'Reorder',       'cluster',  @(s) ischar(s) || isstring(s));
addParameter(p, 'ResultsDir',    '',         @(s) ischar(s) || isstring(s));
addParameter(p, 'RunParameters', [],         @(x) isempty(x) || isstruct(x));
parse(p, varargin{:});
opts     = p.Results;
scheme   = char(opts.Normalisation);
corrType = char(opts.CorrType);
reorder  = lower(char(opts.Reorder));
if ~ismember(reorder, {'cluster', 'fixed'})
    error('compareCentralityMeasures:BadReorder', ...
        'Reorder must be ''cluster'' or ''fixed'' (got ''%s'').', reorder);
end

resultsDir = resolveMouseResultsDir(opts.ResultsDir);

%% Discover and load per-mouse result files
pattern = sprintf('stabilityCentralities_*_%s_results.mat', scheme);
files = dir(fullfile(resultsDir, pattern));
if isempty(files)
    error('compareCentralityMeasures:NoResults', ...
        'No files matching %s found under %s. Run runAllMiceStabilityCentralities first.', ...
        pattern, resultsDir);
end

[measures, displayLabels] = defineMeasures();
M = numel(measures);

mice   = cell(0, 1);
loaded = cell(0, 1);
for k = 1:numel(files)
    s = load(fullfile(resultsDir, files(k).name));
    if ~isfield(s, 'results'); continue; end
    r = s.results;
    if ~isfield(r, 'centralities') || isempty(fieldnames(r.centralities))
        warning('compareCentralityMeasures:NoCentralities', ...
            'Mouse %s has no centralities field -- skipped. Re-run runMouseStabilityCentralities.', ...
            r.mouseId);
        continue;
    end
    mice{end+1, 1}   = char(r.mouseId); %#ok<AGROW>
    loaded{end+1, 1} = r;                %#ok<AGROW>
end
if isempty(loaded)
    error('compareCentralityMeasures:Empty', ...
        'No mice with the centralities field. Re-run runAllMiceStabilityCentralities.');
end
nMice = numel(loaded);
N     = loaded{1}.N;

fprintf('compareCentralityMeasures: %d mouse(/mice), %d centrality measures, %d nodes.\n', ...
    nMice, M, N);

%% Per-mouse correlation matrices
corrAll = NaN(M, M, nMice);
for m = 1:nMice
    X = NaN(N, M);
    for i = 1:M
        X(:, i) = extractMeasure(loaded{m}, measures{i}, N);
    end
    if isfield(loaded{m}, 'isTrivial')
        X(loaded{m}.isTrivial, :) = NaN;  % drop trivial nodes
    end
    corrAll(:, :, m) = safeCorr(X, corrType);
end

meanCorr = mean(corrAll, 3, 'omitnan');
sdCorr   = std(corrAll, 0, 3, 'omitnan');

%% Hierarchical-clustering reorder
% Compute the ordering once on the mean correlation matrix so per-mouse
% heatmaps and the mean heatmap share the same axes (matches the Fig 12/13
% convention in the paper). Falls back to fixed ordering if the Statistics
% Toolbox is unavailable.
if strcmp(reorder, 'cluster')
    [order, reorderNote] = clusterOrder(meanCorr);
else
    order = 1:M;
    reorderNote = 'fixed canonical order';
end
displayOrdered = displayLabels(order);
fprintf('Heatmap reorder = %s (%s).\n', reorder, reorderNote);

%% Per-mouse heatmaps
if opts.PerMouse
    for m = 1:nMice
        nUsed = N;
        if isfield(loaded{m}, 'isTrivial')
            nUsed = sum(~loaded{m}.isTrivial);
        end
        ttl = sprintf('Centrality correlations: %s   (%s, n_{nodes}=%d, scheme=%s, %s)', ...
            strrep(mice{m}, '_', '\_'), corrType, nUsed, scheme, reorder);
        fig = plotCorrHeatmap(corrAll(order, order, m), displayOrdered, ttl, true);
        if opts.SaveResults
            baseName = sprintf('centrality_corr_%s_%s', mice{m}, scheme);
            saveFigBoth(fig, fullfile(resultsDir, baseName));
        end
    end
end

%% Mean-across-mice heatmap (analogous to Fig 3)
ttl = sprintf('Mean centrality correlations across %d mice   (%s, scheme=%s, %s)', ...
    nMice, corrType, scheme, reorder);
fig = plotCorrHeatmap(meanCorr(order, order), displayOrdered, ttl, true);
if opts.SaveResults
    baseName = sprintf('centrality_corr_mean_%s', scheme);
    saveFigBoth(fig, fullfile(resultsDir, baseName));
end

%% Bar chart -- D(->i) and D(k->) vs each classical centrality
%   (analogous to Fig 4 of the paper, but aggregated across mice
%   instead of swept across a randomization parameter)
classicalIdx    = 3:M;            % skip the two stability measures
classicalLabels = displayLabels(classicalIdx);
nClassical      = numel(classicalIdx);

corrSusc = squeeze(corrAll(1, classicalIdx, :)).';   % nMice x nClassical
corrInfl = squeeze(corrAll(2, classicalIdx, :)).';
if nMice == 1
    corrSusc = corrSusc(:).';   % keep row orientation when a singleton
    corrInfl = corrInfl(:).';
end

meanSusc = mean(corrSusc, 1, 'omitnan');
sdSusc   = std(corrSusc,  0, 1, 'omitnan');
meanInfl = mean(corrInfl, 1, 'omitnan');
sdInfl   = std(corrInfl,  0, 1, 'omitnan');

fig = figure('Name', 'Centrality correlations -- D(->i), D(k->) vs classical', ...
             'Position', [80 80 1100 520]);
ax = axes(fig); hold(ax, 'on');
b = bar(ax, [meanSusc; meanInfl].');
b(1).FaceColor = [0.20 0.40 0.80];
b(2).FaceColor = [0.85 0.33 0.10];
b(1).EdgeColor = 'none';
b(2).EdgeColor = 'none';
% Error bars at the actual XEndPoints of each bar group
errorbar(ax, b(1).XEndPoints, meanSusc, sdSusc, 'k.', 'LineWidth', 0.8, ...
    'CapSize', 6, 'HandleVisibility', 'off');
errorbar(ax, b(2).XEndPoints, meanInfl, sdInfl, 'k.', 'LineWidth', 0.8, ...
    'CapSize', 6, 'HandleVisibility', 'off');
yline(ax, 0, '--', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
set(ax, 'XTick', 1:nClassical, 'XTickLabel', classicalLabels, ...
    'TickLabelInterpreter', 'tex');
ylim(ax, [-1 1]);
ylabel(ax, sprintf('%s correlation', corrType));
legend(ax, {'D(\rightarrow i)', 'D(k \rightarrow)'}, 'Location', 'best');
title(ax, sprintf(['Mean \\pm SD correlation of D(\\rightarrow i) and ' ...
    'D(k \\rightarrow) with classical centralities (%d mice, scheme=%s)'], ...
    nMice, scheme), 'Interpreter', 'tex');
grid(ax, 'on'); box(ax, 'on');

if opts.SaveResults
    baseName = sprintf('centrality_corr_bars_%s', scheme);
    saveFigBoth(fig, fullfile(resultsDir, baseName));
end

%% Pack summary
summary = struct();
summary.scheme        = scheme;
summary.corrType      = corrType;
summary.mice          = mice;
summary.measureNames  = measures;
summary.displayLabels = {displayLabels{:}}.';  %#ok<CCAT1>
summary.corrAll       = corrAll;     % M x M x nMice
summary.meanCorr      = meanCorr;    % M x M
summary.sdCorr        = sdCorr;      % M x M
summary.bar.classical = classicalLabels;
summary.bar.meanSusc  = meanSusc;
summary.bar.sdSusc    = sdSusc;
summary.bar.meanInfl  = meanInfl;
summary.bar.sdInfl    = sdInfl;
summary.reorder       = reorder;
summary.clusterOrder  = order;
summary.orderedLabels = {displayOrdered{:}}.';  %#ok<CCAT1>

if isempty(opts.RunParameters)
    summary.runParameters = mouseExperimentRunParameters('buildFromCompareCentrality', ...
        opts, scheme, resultsDir);
else
    summary.runParameters = mouseExperimentRunParameters('merge', opts.RunParameters, ...
        mouseExperimentRunParameters('buildFromCompareCentrality', opts, scheme, resultsDir));
end

if opts.SaveResults
    save(fullfile(resultsDir, sprintf('centrality_corr_%s.mat', scheme)), ...
        '-struct', 'summary');
end
end

%% ------------------------------------------------------------------
function [measures, labels] = defineMeasures()
%DEFINEMEASURES  Ordered list of centrality measures to include.
% First two entries MUST be the stability centralities D(->i) and D(k->),
% so the bar chart can pick them up by fixed index.
measures = { ...
    'D_susceptibility', ...
    'D_influence', ...
    'centralities.in_strength', ...
    'centralities.out_strength', ...
    'centralities.eigenvector_left', ...
    'centralities.eigenvector_right', ...
    'centralities.pagerank', ...
    'centralities.katz', ...
    'centralities.self_comm', ...
    'centralities.betweenness', ...
    'centralities.closeness_in', ...
    'centralities.closeness_out'};
labels = { ...
    'D(\rightarrow i)', ...
    'D(k \rightarrow)', ...
    'in-DC', ...
    'out-DC', ...
    'EC^L', ...
    'EC^R', ...
    'PR', ...
    'KZ', ...
    'SelfC', ...
    'BC', ...
    'CC^{in}', ...
    'CC^{out}'};
end

%% ------------------------------------------------------------------
function v = extractMeasure(r, dottedName, N)
%EXTRACTMEASURE  Resolve a dotted field path like 'centralities.pagerank'
% from a results struct, returning an N-by-1 vector. Returns NaN(N,1) if
% any field is missing.
parts = strsplit(dottedName, '.');
cur = r;
for k = 1:numel(parts)
    if isstruct(cur) && isfield(cur, parts{k})
        cur = cur.(parts{k});
    else
        v = NaN(N, 1);
        return;
    end
end
if isnumeric(cur)
    v = cur(:);
    if numel(v) ~= N
        v = NaN(N, 1);
    end
else
    v = NaN(N, 1);
end
end

%% ------------------------------------------------------------------
function C = safeCorr(X, corrType)
%SAFECORR  Pairwise correlation that tolerates all-NaN columns.
M = size(X, 2);
C = NaN(M, M);
finiteCol = false(1, M);
for j = 1:M
    finiteCol(j) = any(isfinite(X(:, j)) & ~isnan(X(:, j)));
end
if ~any(finiteCol); return; end
sub = X(:, finiteCol);
try
    Csub = corr(sub, 'Type', corrType, 'Rows', 'pairwise');
catch ME
    warning('compareCentralityMeasures:Corr', ...
        'corr() failed (%s) -- correlation matrix is NaN.', ME.message);
    return;
end
C(finiteCol, finiteCol) = Csub;
end

%% ------------------------------------------------------------------
function fig = plotCorrHeatmap(C, labels, ttl, withText)
%PLOTCORRHEATMAP  Render a correlation matrix as a diverging heatmap.
fig = figure('Position', [60 60 820 720]);
ax = axes(fig);
imagesc(ax, C, [-1, 1]);
axis(ax, 'square');
colormap(ax, divergingBlueWhiteRed(256));
cb = colorbar(ax);
cb.Label.String = 'correlation';
M = size(C, 1);
set(ax, 'XTick', 1:M, 'XTickLabel', labels, 'YTick', 1:M, ...
    'YTickLabel', labels, 'TickLabelInterpreter', 'tex');
xtickangle(ax, 45);
if withText
    for i = 1:M
        for j = 1:M
            v = C(i, j);
            if isnan(v)
                txt = '\bullet';
                color = [0.4 0.4 0.4];
            else
                txt = sprintf('%.2f', v);
                if abs(v) > 0.55
                    color = [1 1 1];
                else
                    color = [0 0 0];
                end
            end
            text(ax, j, i, txt, 'HorizontalAlignment', 'center', ...
                'FontSize', 7, 'Color', color, 'Interpreter', 'tex');
        end
    end
end
title(ax, ttl, 'Interpreter', 'tex');
end

%% ------------------------------------------------------------------
function cmap = divergingBlueWhiteRed(n)
%DIVERGINGBLUEWHITERED  Smooth blue -> white -> red colormap on [-1, 1].
t = linspace(-1, 1, n).';
r = clamp(1 + 0.7 * t, 0, 1);
g = clamp(1 - 0.9 * abs(t), 0, 1);
b = clamp(1 - 0.7 * t, 0, 1);
cmap = [r, g, b];
end

function y = clamp(x, lo, hi)
y = min(max(x, lo), hi);
end

%% ------------------------------------------------------------------
function saveFigBoth(fig, basePath)
%SAVEFIGBOTH  Save a figure as both .fig and .png.
savefig(fig, [basePath '.fig']);
try
    exportgraphics(fig, [basePath '.png'], 'Resolution', 200);
catch
    saveas(fig, [basePath '.png']);
end
end

%% ------------------------------------------------------------------
function [order, note] = clusterOrder(C)
%CLUSTERORDER  Hierarchical-clustering permutation that places measures
% with similar correlation profiles next to each other in the heatmap
% (matches Liao & Lizier Fig 12/13). Distance = 1 - C (signed), so
% strongly anti-correlated measures end up far apart, mirroring the
% paper's clustering. Linkage uses the average-linkage UPGMA rule, with
% optimal leaf ordering applied when the Statistics Toolbox supports it.
% Falls back to the canonical order if the toolbox is unavailable.
M = size(C, 1);
if M < 2 || exist('linkage', 'file') ~= 2
    order = 1:M;
    note = 'fallback (Statistics Toolbox unavailable)';
    return;
end
Cf = C; Cf(isnan(Cf)) = 0;
Cf(1:M+1:end) = 1;                  % ensure self-corr is 1
D  = 1 - Cf;                        % distance in [0, 2]
D(1:M+1:end) = 0;
D = (D + D.') / 2;                  % numerical symmetry
try
    d = squareform(D, 'tovector');
    Z = linkage(d, 'average');
    if exist('optimalleaforder', 'file') == 2
        order = optimalleaforder(Z, d);
        note  = 'average linkage, optimal leaf ordering';
    else
        h = figure('Visible', 'off');
        [~, ~, order] = dendrogram(Z, 0);
        close(h);
        note = 'average linkage, dendrogram default ordering';
    end
catch ME
    warning('compareCentralityMeasures:Cluster', ...
        'Hierarchical clustering failed (%s) -- using fixed order.', ME.message);
    order = 1:M;
    note  = sprintf('fallback (linkage failed: %s)', ME.message);
end
order = order(:).';
end
