function layout = experimentResultsLayout(resultsDir, createDirs)
%EXPERIMENTRESULTSLAYOUT  Standard subfolder paths for one experiment run.
%
%   layout = experimentResultsLayout(resultsDir)
%   layout = experimentResultsLayout(resultsDir, true)   % also create them
%
%   One run of the connectome pipeline emits several families of artefacts.
%   Keeping each family in its own subfolder means a run directory stays
%   readable, and a reader can find the provenance record without wading
%   through a few hundred figures.
%
%   Fields (each an absolute path)
%     root                  : the run directory itself
%     provenance            : config snapshot, parameters, manifest, logs
%     qc                    : input connectome quality-control figures
%     stability             : per-subject and cohort stability centralities
%     centralityCorrelation : stability vs classical centrality correlations
%     cohortComparison      : case-vs-control per-node comparisons
%     hemisphericAsymmetry  : left-right laterality analyses
%     dstCohort             : network-level D_st across the cohort
%     comparisonFigures     : publication-style summary figures
%
%   See also EXPERIMENTRESULTSDIR, LOCATEEXPERIMENTRESULTSFILE.

if nargin < 2 || isempty(createDirs)
    createDirs = false;
end

resultsDir = char(resultsDir);
layout = struct();
layout.root                  = resultsDir;
layout.provenance            = fullfile(resultsDir, 'provenance');
layout.qc                    = fullfile(resultsDir, 'qc');
layout.stability             = fullfile(resultsDir, 'stability_centralities');
layout.centralityCorrelation = fullfile(resultsDir, 'centrality_correlation');
layout.cohortComparison      = fullfile(resultsDir, 'cohort_comparison');
layout.hemisphericAsymmetry  = fullfile(resultsDir, 'hemispheric_asymmetry');
layout.dstCohort             = fullfile(resultsDir, 'dst_cohort');
layout.comparisonFigures     = fullfile(resultsDir, 'comparison_figures');

if createDirs
    names = fieldnames(layout);
    for k = 1:numel(names)
        d = layout.(names{k});
        if ~exist(d, 'dir')
            mkdir(d);
        end
    end
end
end
