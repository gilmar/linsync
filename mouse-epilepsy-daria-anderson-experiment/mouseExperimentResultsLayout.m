function layout = mouseExperimentResultsLayout(resultsDir)
%MOUSEEXPERIMENTRESULTSLAYOUT  Subfolder paths for one experiment run under results/.
%
%   layout = mouseExperimentResultsLayout(resultsDir)
%
%   Fields:
%     root, provenance, qc, stability, centralityCorrelation,
%     andersonVsArnold, lrAsymmetry, dstCohort, comparisonFigures

resultsDir = char(resultsDir);
layout = struct();
layout.root                  = resultsDir;
layout.provenance            = fullfile(resultsDir, 'provenance');
layout.qc                    = fullfile(resultsDir, 'qc');
layout.stability             = fullfile(resultsDir, 'stability_centralities');
layout.centralityCorrelation = fullfile(resultsDir, 'centrality_correlation');
layout.andersonVsArnold      = fullfile(resultsDir, 'anderson_vs_arnold');
layout.lrAsymmetry           = fullfile(resultsDir, 'lr_asymmetry');
layout.dstCohort             = fullfile(resultsDir, 'dst_cohort');
layout.comparisonFigures     = fullfile(resultsDir, 'comparison_figures');
end
