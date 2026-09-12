function searchDirs = mouseResultsSearchDirs(resultsDir, category)
%MOUSERESULTSSEARCHDIRS  Search category subfolder first, then legacy flat root.
layout = mouseExperimentResultsLayout(resultsDir);
searchDirs = {mouseResultsDir(resultsDir, category), layout.root};
end
