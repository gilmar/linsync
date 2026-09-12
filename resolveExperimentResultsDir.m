function resultsDir = resolveExperimentResultsDir(resultsDirArg, projectRoot)
%RESOLVEEXPERIMENTRESULTSDIR  Normalise the ResultsDir option used by every writer.
%
%   resultsDir = resolveExperimentResultsDir(opts.ResultsDir)
%   resultsDir = resolveExperimentResultsDir(opts.ResultsDir, projectRoot)
%
%   When a directory is supplied it is created (with its standard category
%   subfolders) and returned. When it is empty the default is
%   <projectRoot>/results, with projectRoot defaulting to the current
%   folder -- so an interactive call with no options still writes somewhere
%   sensible instead of erroring.
%
%   See also EXPERIMENTRESULTSLAYOUT, SETUPCONNECTOMEPATHS.

if nargin < 1
    resultsDirArg = '';
end
if nargin < 2 || isempty(projectRoot)
    projectRoot = pwd;
end

resultsDir = char(resultsDirArg);
if isempty(resultsDir)
    resultsDir = fullfile(char(projectRoot), 'results');
end
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end
experimentResultsLayout(resultsDir, true);
end
