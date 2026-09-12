function files = listSubjectResultFiles(resultsDir, scheme)
%LISTSUBJECTRESULTFILES  Discover the per-subject result .mat files of a run.
%
%   files = listSubjectResultFiles(resultsDir, 'column')
%
%   Returns a dir() struct array with one entry per subject, for the files
%   runStabilityCentralities wrote as
%     <prefix>_<subjectId>_<scheme>_results.mat
%
%   Downstream comparisons discover the cohort this way rather than being
%   handed a subject list, so a comparison always describes exactly the
%   results present in the folder.
%
%   Both the stability_centralities/ subfolder and the run root are
%   searched; when a subject appears in both, the subfolder copy wins.
%
%   See also EXPERIMENTRESULTPREFIX, RUNCOHORTSTABILITYCENTRALITIES.

if nargin < 2 || isempty(scheme)
    error('listSubjectResultFiles:Args', 'resultsDir and scheme are required.');
end

scheme = char(scheme);
prefix = experimentResultPrefix();
primaryDir = experimentResultsDir(resultsDir, 'stability', false);
searchDirs = {primaryDir, char(resultsDir)};

allFiles = [];
for d = 1:numel(searchDirs)
    found = dir(fullfile(searchDirs{d}, sprintf('%s_*_%s_results.mat', prefix, scheme)));
    allFiles = [allFiles; found]; %#ok<AGROW>
end
if isempty(allFiles)
    files = allFiles;
    return;
end

bySubject = containers.Map('KeyType', 'char', 'ValueType', 'any');
for k = 1:numel(allFiles)
    subjectId = subjectIdFromResultFile(allFiles(k).name, prefix, scheme);
    if isempty(subjectId)
        continue;
    end
    if ~isKey(bySubject, subjectId) || ...
            (strcmp(allFiles(k).folder, primaryDir) && ...
             ~strcmp(bySubject(subjectId).folder, primaryDir))
        bySubject(subjectId) = allFiles(k);
    end
end

vals = bySubject.values;
if isempty(vals)
    files = allFiles(1:0);
else
    files = vertcat(vals{:});
    files = files(:);
end
end

%% ------------------------------------------------------------------
function subjectId = subjectIdFromResultFile(fname, prefix, scheme)
subjectId = '';
[~, base, ext] = fileparts(fname);
if ~strcmp(ext, '.mat')
    return;
end
head = [prefix '_'];
tail = ['_' scheme '_results'];
if ~startsWith(base, head) || ~endsWith(base, tail)
    return;
end
subjectId = base(numel(head) + 1:end - numel(tail));
end
