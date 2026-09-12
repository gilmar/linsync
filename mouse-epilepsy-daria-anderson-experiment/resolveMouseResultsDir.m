function resultsDir = resolveMouseResultsDir(resultsDirArg)
%RESOLVEMOUSERESULTSDIR  Resolve output directory for mouse experiment scripts.
%
%   If resultsDirArg is non-empty, use it (create if missing). Otherwise
%   return setupMousePaths() default (flat results/).

if nargin < 1
    resultsDirArg = '';
end
if ~isempty(resultsDirArg)
    resultsDir = char(resultsDirArg);
    if ~exist(resultsDir, 'dir')
        mkdir(resultsDir);
    end
    ensureMouseExperimentResultsDirs(resultsDir);
else
    resultsDir = setupMousePaths();
end
end
