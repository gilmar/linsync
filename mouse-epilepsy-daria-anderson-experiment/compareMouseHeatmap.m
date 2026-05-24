%% compareMouseHeatmap.m
% Sanity check: load each mouse's coarse connectome (66x66) and render
% it with imagesc; lay it out next to the reference PNG that ships in
% data/<mouseId>/coarse_connectome_<mouseid>.png so we can confirm
% loadMouseConnectome agrees with the source.
%
% Usage:
%   compareMouseHeatmap                     % default = 'reference' style
%   compareMouseHeatmap('reference')        % match reference PNG layout
%   compareMouseHeatmap('csv')              % keep CSV row/column order
%   compareMouseHeatmap('ResultsDir', path) % write figures under path
%
% Layout styles:
%   'csv'        - imagesc(K) verbatim. Row 1 (L-CORTEX_VISUAL) at the top,
%                  row 66 (R-BACKGROUND) at the bottom. Matches the order
%                  in fine_family_labelled_coarse.csv.
%   'reference'  - permutes rows/columns to match the supplied PNGs:
%                    * R-* hemisphere first, then L-*
%                    * alphabetical sort within each hemisphere
%                    * y-axis flipped (set(gca,'YDir','normal')) so the
%                      first plotted row ends up at the BOTTOM, exactly as
%                      in data/<mouseId>/coarse_connectome_*.png.
%
% Color scaling is linear from 0 to the 99th percentile of the
% non-diagonal entries to suppress saturation by self-loops.

function compareMouseHeatmap(varargin)
setupMousePaths();
style = 'reference';
nvArgs = varargin;
if ~isempty(varargin) && (ischar(varargin{1}) || isstring(varargin{1})) ...
        && ismember(lower(char(varargin{1})), {'csv', 'reference'})
    style = char(varargin{1});
    nvArgs = varargin(2:end);
end
p = inputParser;
addParameter(p, 'ResultsDir', '', @(s) ischar(s) || isstring(s));
parse(p, nvArgs{:});
style = validatestring(style, {'csv', 'reference'});
resultsDir = resolveMouseResultsDir(p.Results.ResultsDir);

mice = listAvailableMice();
nMice = numel(mice);
if nMice == 0
    error('compareMouseHeatmap: no mouse coarse CSVs found under data/.');
end
nCols = 3;
nRows = ceil(nMice / nCols);

figName = sprintf('Mouse coarse connectomes (%s layout)', style);
fig = figure('Name', figName, ...
             'Position', [50 50 600 * nCols, 500 * nRows]);

for m = 1:nMice
    mouseId = mice{m};
    [K, labels, info] = loadMouseConnectome(mouseId);

    fprintf('%s : %dx%d, symmetryError = %.3g, trivial nodes = %d\n', ...
        mouseId, info.N, info.N, info.symmetryError, sum(info.isTrivial));
    if any(info.isTrivial)
        idx = find(info.isTrivial);
        for k = 1:numel(idx)
            fprintf('   trivial node %2d: %s\n', idx(k), labels{idx(k)});
        end
    end

    Knd = K - diag(diag(K));
    posVals = Knd(Knd > 0);
    if isempty(posVals)
        topVal = 1;
    else
        topVal = prctile(posVals, 99);
    end

    switch style
        case 'csv'
            perm        = 1:numel(labels);
            flipYAxis   = false;
        case 'reference'
            perm        = referencePermutation(labels);
            flipYAxis   = true;
    end
    Kp     = K(perm, perm);
    labelsP = labels(perm);

    ax = subplot(nRows, nCols, m);
    imagesc(Kp, [0 topVal]);
    axis(ax, 'square');
    colormap(ax, parula);
    colorbar(ax);
    title(ax, sprintf('%s (clim 0..%.0f)', mouseId, topVal), ...
          'Interpreter', 'none');
    xlabel(ax, 'to (column)');
    ylabel(ax, 'from (row)');

    set(ax, 'XTick', 1:numel(labelsP), 'XTickLabel', labelsP, ...
            'YTick', 1:numel(labelsP), 'YTickLabel', labelsP, ...
            'TickLabelInterpreter', 'none', 'FontSize', 5);
    xtickangle(ax, 90);

    if flipYAxis
        set(ax, 'YDir', 'normal');   % row 1 at bottom, like the PNG
    end
end
sgtitle(sprintf(['Mouse coarse connectomes (%s layout) -- compare with ', ...
                 'data/<mouseId>/coarse\\_connectome\\_*.png'], style));

baseName = sprintf('mouse_heatmaps_overview_%s', style);
savefig(fig, fullfile(resultsDir, [baseName '.fig']));
try
    exportgraphics(fig, fullfile(resultsDir, [baseName '.png']), 'Resolution', 200);
catch
    saveas(fig, fullfile(resultsDir, [baseName '.png']));
end
fprintf('Saved %s\n', fullfile(resultsDir, [baseName '.png']));
close(fig);
end


function perm = referencePermutation(labels)
%REFERENCEPERMUTATION  Permutation that reorders the CSV labels to match
% the layout used by the supplied coarse_connectome_*.png files:
%   - R-* hemisphere first, then L-*
%   - alphabetical sort by region name inside each hemisphere
% Combined with set(gca,'YDir','normal') in the caller, this places R-*
% labels in the upper half of the plot exactly like the reference image.

labels = string(labels(:));
isR    = startsWith(labels, "R-");
isL    = startsWith(labels, "L-");

idxR = find(isR);
idxL = find(isL);
idxOther = find(~(isR | isL));   % defensive; expected to be empty

[~, ordR] = sort(labels(idxR));
[~, ordL] = sort(labels(idxL));

perm = [idxR(ordR); idxL(ordL); idxOther];
end
