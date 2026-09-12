function fpath = locateExperimentResultsFile(resultsDir, category, filename)
%LOCATEEXPERIMENTRESULTSFILE  Find a results file in its category, or at the run root.
%
%   fpath = locateExperimentResultsFile(resultsDir, 'stability', 'summary.mat')
%
%   Searches the category subfolder first, then the run directory itself so
%   that flat result folders written by hand, or by older versions of the
%   pipeline, are still readable. Returns '' when nothing matches.
%
%   See also EXPERIMENTRESULTSDIR.

searchDirs = {experimentResultsDir(resultsDir, category, false), char(resultsDir)};
fpath = '';
for k = 1:numel(searchDirs)
    candidate = fullfile(searchDirs{k}, filename);
    if isfile(candidate)
        fpath = candidate;
        return;
    end
end
end
