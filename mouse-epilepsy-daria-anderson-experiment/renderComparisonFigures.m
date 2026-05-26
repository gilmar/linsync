function paths = renderComparisonFigures(andersonSummary, lrSummary, varargin)
%RENDERCOMPARISONFIGURES  Export cohort comparison summary figures (column, in-memory).
%
%   paths = renderComparisonFigures(andersonSummary, lrSummary)
%   paths = renderComparisonFigures(..., 'OutputDir', fullfile(resultsDir,'comparison_figures'))

p = inputParser;
addParameter(p, 'OutputDir', '', @(s) ischar(s) || isstring(s));
addParameter(p, 'OverrideCfg', struct(), @isstruct);
addParameter(p, 'CloseFigures', true, @islogical);
parse(p, varargin{:});
opts = p.Results;

if nargin < 2
    lrSummary = struct();
end
outDir = resolveFigureOutDir(opts.OutputDir, andersonSummary);
paths = struct();

f1 = plotRegionSuscOutliers(andersonSummary, 'OutputDir', outDir, 'Save', true);
paths.regionSusc = fullfile(outDir, 'region_susc_outliers.png');

f2 = plotRobustVsSuggestiveSummary(andersonSummary, lrSummary, ...
    'OutputDir', outDir, 'Save', true, 'OverrideCfg', opts.OverrideCfg);
paths.robustSuggestive = fullfile(outDir, 'robust_vs_suggestive_summary.png');

f3 = plotLateralityInfluenceOutliers(lrSummary, 'OutputDir', outDir, 'Save', true);
paths.lateralityInfluence = fullfile(outDir, 'laterality_influence_outliers.png');

if opts.CloseFigures
    closeIfValid(f1);
    closeIfValid(f2);
    closeIfValid(f3);
end

fprintf('Comparison figures -> %s\n', outDir);
end

%% ------------------------------------------------------------------
function closeIfValid(fig)
if ~isempty(fig) && isgraphics(fig, 'figure')
    close(fig);
end
end
