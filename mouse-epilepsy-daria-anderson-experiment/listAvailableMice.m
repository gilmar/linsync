function mice = listAvailableMice()
%LISTAVAILABLEMICE  Return the list of mouse IDs under data/ that have a
% coarse-grained connectome CSV (fine_family_labelled_coarse.csv).
% Mice with no such file (e.g. Arnold_2 in this snapshot) are skipped.

experimentRoot = fileparts(mfilename('fullpath'));
dataRoot = fullfile(experimentRoot, 'data');
d = dir(dataRoot);
candidates = {d([d.isdir]).name};
candidates = candidates(~ismember(candidates, {'.', '..'}));
candidates = sort(candidates);

mice = {};
for k = 1:numel(candidates)
    csvFile = fullfile(dataRoot, candidates{k}, 'fine_family_labelled_coarse.csv');
    if exist(csvFile, 'file')
        mice{end+1, 1} = candidates{k}; %#ok<AGROW>
    else
        warning('listAvailableMice:NoCSV', ...
            'Skipping %s -- no fine_family_labelled_coarse.csv', candidates{k});
    end
end
end
