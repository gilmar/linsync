function outDir = mouseResultsDir(resultsDir, category)
%MOUSERESULTSDIR  Path for one results category; creates the subfolder if needed.
layout = mouseExperimentResultsLayout(resultsDir);
outDir = resolveResultsCategory(layout, category);
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
end

%% ------------------------------------------------------------------
function d = resolveResultsCategory(layout, category)
category = lower(strrep(strrep(category, '-', '_'), ' ', '_'));
switch category
    case 'provenance'
        d = layout.provenance;
    case {'qc', 'heatmap'}
        d = layout.qc;
    case {'stability', 'stability_centralities', 'stabilitycentralities'}
        d = layout.stability;
    case {'centrality_correlation', 'centralitycorrelation', 'centrality_corr'}
        d = layout.centralityCorrelation;
    case {'anderson_vs_arnold', 'andersonvsarnold', 'anderson'}
        d = layout.andersonVsArnold;
    case {'lr_asymmetry', 'lrasymmetry', 'lr'}
        d = layout.lrAsymmetry;
    case {'dst_cohort', 'dstcohort', 'dst'}
        d = layout.dstCohort;
    case {'comparison_figures', 'comparisonfigures', 'comparison'}
        d = layout.comparisonFigures;
    otherwise
        error('mouseResultsDir:UnknownCategory', 'Unknown category "%s".', category);
end
end
