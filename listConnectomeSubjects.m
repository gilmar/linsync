function subjects = listConnectomeSubjects(dataRoot, fileName, varargin)
%LISTCONNECTOMESUBJECTS  Discover the subjects present under a study data folder.
%
%   subjects = listConnectomeSubjects(dataRoot)
%   subjects = listConnectomeSubjects(dataRoot, 'connectome.csv')
%   subjects = listConnectomeSubjects(dataRoot, fileName, 'Warn', false)
%
%   Returns the sorted IDs of every immediate subfolder of `dataRoot` that
%   actually contains the named connectome file. Subjects are discovered
%   from the filesystem rather than declared in a list, so a cohort that
%   grows or loses a subject needs no code or config change -- and a folder
%   whose data has not arrived yet is skipped instead of failing a run
%   halfway through.
%
%   Inputs
%     dataRoot : folder holding one subfolder per subject
%     fileName : connectome filename inside each subject folder.
%                Default 'connectome.csv'.
%
%   Name-value options
%     'Warn' : true -- warn once per skipped folder, naming the missing file
%
%   Output
%     subjects : Nx1 cellstr of subject IDs (folder names), sorted
%
%   See also LOADCONNECTOME, COHORTGROUPMASK.

if nargin < 1 || isempty(dataRoot)
    error('listConnectomeSubjects:NoDataRoot', 'dataRoot is required.');
end
if nargin < 2 || isempty(fileName)
    fileName = 'connectome.csv';
end

p = inputParser;
addParameter(p, 'Warn', true, @islogical);
parse(p, varargin{:});

dataRoot = char(dataRoot);
if ~isfolder(dataRoot)
    error('listConnectomeSubjects:NoSuchFolder', ...
        'Data folder does not exist: %s', dataRoot);
end

d = dir(dataRoot);
candidates = {d([d.isdir]).name};
candidates = candidates(~ismember(candidates, {'.', '..'}));
candidates = sort(candidates);

subjects = cell(0, 1);
for k = 1:numel(candidates)
    csvFile = fullfile(dataRoot, candidates{k}, char(fileName));
    if isfile(csvFile)
        subjects{end+1, 1} = candidates{k}; %#ok<AGROW>
    elseif p.Results.Warn
        warning('listConnectomeSubjects:NoConnectome', ...
            'Skipping "%s": no %s in that folder.', candidates{k}, fileName);
    end
end
end
