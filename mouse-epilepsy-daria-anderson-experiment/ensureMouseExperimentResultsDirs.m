function ensureMouseExperimentResultsDirs(resultsDir)
%ENSUREMOUSEEXPERIMENTRESULTSDIRS  Create all standard results subfolders.
layout = mouseExperimentResultsLayout(resultsDir);
cats = fieldnames(layout);
for k = 1:numel(cats)
    if strcmp(cats{k}, 'root')
        continue;
    end
    d = layout.(cats{k});
    if ~exist(d, 'dir')
        mkdir(d);
    end
end
end
