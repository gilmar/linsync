function fig = plotConnectomeHeatmaps(varargin)
%PLOTCONNECTOMEHEATMAPS  Quality-control grid of every subject's raw connectome.
%
%   fig = plotConnectomeHeatmaps('DataRoot', dataRoot)
%   fig = plotConnectomeHeatmaps('DataRoot', dataRoot, 'Layout', 'hemisphere')
%
%   Renders each subject's weighted adjacency matrix as a heatmap in one
%   figure. This is the first thing to look at in a new study: a subject
%   whose matrix is transposed, mis-parcellated, empty or on a wildly
%   different scale is obvious here and invisible three analysis steps
%   later.
%
%   Colour scaling runs from 0 to the 99th percentile of the off-diagonal
%   entries, so a few heavy tracts (and any self-loops) cannot saturate the
%   image and hide the structure.
%
%   Name-value options
%     'DataRoot'       : '' -- subject folders (default: pwd/data)
%     'ConnectomeFile' : 'connectome.csv'
%     'Subjects'       : {} -- default: discover under DataRoot
%     'Layout'         : 'file' (default) keeps the file's node order;
%                        'hemisphere' groups by hemisphere prefix and
%                        sorts alphabetically within each, which makes
%                        homotopic structure and inter-hemispheric blocks
%                        visible as blocks
%     'LeftPrefix'     : 'L-'
%     'RightPrefix'    : 'R-'
%     'TrivialLabelSuffixes' : {}
%     'ResultsDir'     : '' -- when saving, the run folder
%     'SaveResults'    : true
%
%   Output (under <ResultsDir>/qc/)
%     connectome_heatmaps_<layout>.{fig,png}
%
%   See also LOADCONNECTOME, LISTCONNECTOMESUBJECTS.

p = inputParser;
addParameter(p, 'DataRoot',       '', @(s) ischar(s) || isstring(s));
addParameter(p, 'ConnectomeFile', 'connectome.csv', @(s) ischar(s) || isstring(s));
addParameter(p, 'Subjects',       {}, @(c) iscell(c) || isstring(c));
addParameter(p, 'Layout',         'file', @(s) ischar(s) || isstring(s));
addParameter(p, 'LeftPrefix',     'L-', @(s) ischar(s) || isstring(s));
addParameter(p, 'RightPrefix',    'R-', @(s) ischar(s) || isstring(s));
addParameter(p, 'TrivialLabelSuffixes', {}, @(c) iscell(c) || isstring(c) || ischar(c));
addParameter(p, 'ResultsDir',     '', @(s) ischar(s) || isstring(s));
addParameter(p, 'SaveResults',    true, @islogical);
parse(p, varargin{:});
opts = p.Results;

layout = validatestring(char(opts.Layout), {'file', 'hemisphere'});

dataRoot = char(opts.DataRoot);
if isempty(dataRoot)
    dataRoot = fullfile(pwd, 'data');
end

subjects = opts.Subjects;
if isempty(subjects)
    subjects = listConnectomeSubjects(dataRoot, opts.ConnectomeFile);
else
    subjects = cellstr(subjects);
end
if isempty(subjects)
    error('plotConnectomeHeatmaps:NoSubjects', ...
        'No subjects with a "%s" file found under %s.', opts.ConnectomeFile, dataRoot);
end

nSubjects = numel(subjects);
nCols = min(3, nSubjects);
nRows = ceil(nSubjects / nCols);

fig = figure('Name', sprintf('Connectome QC (%s layout)', layout), ...
             'Position', [50 50 560 * nCols, 480 * nRows], 'Color', 'w');

for s = 1:nSubjects
    [K, labels, info] = loadConnectome(subjects{s}, ...
        'DataRoot', dataRoot, 'FileName', opts.ConnectomeFile, ...
        'TrivialLabelSuffixes', opts.TrivialLabelSuffixes, ...
        'WarnAsymmetry', false);

    fprintf('%-16s N=%3d  symmetry error=%.3g  trivial nodes=%d\n', ...
        subjects{s}, info.N, info.symmetryError, sum(info.isTrivial));

    if strcmp(layout, 'hemisphere')
        perm = hemispherePermutation(labels, opts.LeftPrefix, opts.RightPrefix);
        K = K(perm, perm);
        labels = labels(perm);
    end

    ax = subplot(nRows, nCols, s);
    imagesc(ax, K, [0, offDiagonalCeiling(K)]);
    axis(ax, 'square');
    colormap(ax, parula);
    colorbar(ax);
    title(ax, sprintf('%s  (N = %d)', subjects{s}, info.N), 'Interpreter', 'none');

    tickStep = max(1, round(numel(labels) / 20));
    ticks = 1:tickStep:numel(labels);
    set(ax, 'XTick', ticks, 'XTickLabel', labels(ticks), ...
            'YTick', ticks, 'YTickLabel', labels(ticks), ...
            'TickLabelInterpreter', 'none', 'FontSize', 6);
    xtickangle(ax, 90);
end

sgtitle(sprintf('Raw connectomes  --  %s node order', layout), 'Interpreter', 'none');

if opts.SaveResults
    resultsDir = resolveExperimentResultsDir(opts.ResultsDir);
    outDir = experimentResultsDir(resultsDir, 'qc');
    saveFigureBoth(fig, fullfile(outDir, sprintf('connectome_heatmaps_%s', layout)));
    fprintf('Wrote connectome_heatmaps_%s.{fig,png} to %s\n', layout, outDir);
end
end

%% ------------------------------------------------------------------
function ceiling = offDiagonalCeiling(K)
%OFFDIAGONALCEILING  99th percentile of the off-diagonal weights.
v = K(~eye(size(K)));
v = v(isfinite(v));
if isempty(v)
    ceiling = 1;
    return;
end
v = sort(v);
ceiling = v(max(1, min(numel(v), round(0.99 * numel(v)))));
if ~(ceiling > 0)
    ceiling = max(max(v), 1);
end
end

%% ------------------------------------------------------------------
function perm = hemispherePermutation(labels, leftPrefix, rightPrefix)
%HEMISPHEREPERMUTATION  Group nodes by hemisphere, alphabetical within each.
labelStr = string(labels(:));
isLeft  = startsWith(labelStr, string(leftPrefix));
isRight = startsWith(labelStr, string(rightPrefix));
other = ~isLeft & ~isRight;

perm = [sortedIdx(labelStr, isLeft); ...
        sortedIdx(labelStr, isRight); ...
        sortedIdx(labelStr, other)];
end

function idx = sortedIdx(labelStr, mask)
idx = find(mask);
if isempty(idx)
    idx = zeros(0, 1);
    return;
end
[~, order] = sort(labelStr(idx));
idx = idx(order);
idx = idx(:);
end
