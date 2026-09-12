function outDir = experimentResultsDir(resultsDir, category, createDir)
%EXPERIMENTRESULTSDIR  Resolve one results category to a path, creating it.
%
%   outDir = experimentResultsDir(resultsDir, 'stability')
%   outDir = experimentResultsDir(resultsDir, 'cohort_comparison', false)
%
%   Category names are matched case-insensitively with '-' and ' ' treated
%   as '_', and each accepts a few aliases so calling code can read
%   naturally. Pass createDir = false to resolve a path without touching
%   the filesystem.
%
%   Categories (aliases in brackets)
%     provenance
%     qc                    [heatmap]
%     stability             [stability_centralities, stabilitycentralities]
%     centrality_correlation[centrality_corr, centralitycorrelation]
%     cohort_comparison     [cohortcomparison, case_vs_control, comparison_groups]
%     hemispheric_asymmetry [lr_asymmetry, lr, laterality]
%     dst_cohort            [dst, dstcohort]
%     comparison_figures    [comparison, comparisonfigures]
%
%   See also EXPERIMENTRESULTSLAYOUT, LOCATEEXPERIMENTRESULTSFILE.

if nargin < 3 || isempty(createDir)
    createDir = true;
end

layout = experimentResultsLayout(resultsDir);
key = lower(strrep(strrep(char(category), '-', '_'), ' ', '_'));

switch key
    case 'provenance'
        outDir = layout.provenance;
    case {'qc', 'heatmap'}
        outDir = layout.qc;
    case {'stability', 'stability_centralities', 'stabilitycentralities'}
        outDir = layout.stability;
    case {'centrality_correlation', 'centralitycorrelation', 'centrality_corr'}
        outDir = layout.centralityCorrelation;
    case {'cohort_comparison', 'cohortcomparison', 'case_vs_control', 'comparison_groups'}
        outDir = layout.cohortComparison;
    case {'hemispheric_asymmetry', 'hemisphericasymmetry', 'lr_asymmetry', 'lr', 'laterality'}
        outDir = layout.hemisphericAsymmetry;
    case {'dst_cohort', 'dstcohort', 'dst'}
        outDir = layout.dstCohort;
    case {'comparison_figures', 'comparisonfigures', 'comparison'}
        outDir = layout.comparisonFigures;
    otherwise
        error('experimentResultsDir:UnknownCategory', ...
            'Unknown results category "%s".', category);
end

if createDir && ~exist(outDir, 'dir')
    mkdir(outDir);
end
end
