function outDir = resolveFigureOutDir(outputDir, summary)
%RESOLVEFIGUREOUTDIR  Default comparison_figures/ under results dir when available.
if nargin >= 1 && ~isempty(outputDir)
    outDir = char(outputDir);
    return;
end
if nargin >= 2 && isstruct(summary) && isfield(summary, 'runParameters') ...
        && isstruct(summary.runParameters) && isfield(summary.runParameters, 'resultsDir')
    outDir = fullfile(summary.runParameters.resultsDir, 'comparison_figures');
else
    outDir = 'comparison_figures';
end
if ~isfolder(outDir)
    mkdir(outDir);
end
end
