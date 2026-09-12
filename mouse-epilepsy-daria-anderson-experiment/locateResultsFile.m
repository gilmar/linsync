function fpath = locateResultsFile(resultsDir, category, filename)
%LOCATERESULTSFILE  First match under category subfolder or legacy flat root.
searchDirs = mouseResultsSearchDirs(resultsDir, category);
fpath = '';
for k = 1:numel(searchDirs)
    candidate = fullfile(searchDirs{k}, filename);
    if isfile(candidate)
        fpath = candidate;
        return;
    end
end
end
