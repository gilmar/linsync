function combineSyncResults(parameters)
%
% Combine the sync results for randomly generated networks as specified in the 
% parameters object (or file)
%
% Inputs:
% - parameters - an object containing the expected properties, or a string
% describing the filename to run load this object in.
% Importantly, the folder specified in parameters.combineResultsFrom
%  will point to where the sub-folders all contain results that should be
%  combined.
%
%% Linear Sync Toolkit (linsync)
% Copyright (C) 2023 Joseph T. Lizier
% Distributed under GNU General Public License v3

if ischar(parameters)
    % Assume that this string contains a filename which when run will load
    % a properties object for this run
    eval(['run ', parameters]);
end
% Postcondition: parameters are in the parameters object

addpath(parameters.syncToolkitPath);

fprintf('Beginning combining sync results from dir %s into %s\n', parameters.combineResultsFrom, parameters.folder);

% Generate string for the boolean arguments ready for file names
if (parameters.undirected)
    undirString = 'un';
else
    undirString = 'dir';
end
if (parameters.discretized)
    discString = 'disc';
else
    discString = 'cont';
end

% Pull out the array of parameters that we will sweep through here:
paramsToRunThrough = parameters.(parameters.tosweep);
% But replace the swept parameter with its final value (for generating the
% results file names)
parameters.(parameters.tosweep) = parameters.(parameters.tosweep)(end);

[networkType, netTypeSuffix] = generateNetworkTypeStrings(parameters);

for s = parameters.SRangeToPlot
    tic
    
    % Get filename prefix sorted:
    % This will generate filename prefix with the *final* value of the swept
    %  parameter, but parameter being swept will be indicated
    fileNameSuffix = sprintf('N%d-%s%s-b%.2f-c%.2f-p%.4f-sweep_%s-%s-k%d-%s-S%d-repeats*', ...
            parameters.N, networkType, netTypeSuffix, parameters.b, parameters.c, parameters.p, ...
            parameters.tosweep, undirString, parameters.maxMotifLength, discString, s);
    
    totalRepeats = 0;
    allSyncWidths = [];
    allSyncWidthApproxes = [];
    allSyncWidthEmpirical = [];
    allDominantEigenvalues = [];
    allSecondEigenvalues = [];
    
    filesInResultsDir = dir(parameters.combineResultsFrom);
    for index = 3:length(filesInResultsDir) % Skip . and ..
        file = filesInResultsDir(index);
        if file.isdir
            resultsFilenameTemplate = [parameters.combineResultsFrom, '/', file.name, '/', fileNameSuffix, '.mat'];
            try
                resultsFilename = strtrim(ls(resultsFilenameTemplate)); % match the correct number of repeats
            catch ME
                % ls will return an error if there are no matches
                fprintf('No results file found in folder %s\n', file.name);
                continue;
            end
            if ~exist(resultsFilename, 'file')
                fprintf('No results file found in folder %s (should not happen here)\n', file.name);
                continue;
            end
            % Make sure we don't overwrite the parameters variable, just
            % load in the others to be re-saved later
            load(resultsFilename, 'N', 'd', 'b', 'c', 'p', 'undirected', 'maxMotifLength', 'discretized', ...
                'networkType', 'paramsToRunThrough', 'S', 'repeats', 'syncWidths', 'syncWidthApproxes', ...
                'syncWidthEmpirical', 'dominantEigenvalues', 'secondEigenvalues');
            totalRepeats = totalRepeats + repeats;
            fprintf('Loaded %d repeats from folder %s (running total %d)\n', repeats, file.name, totalRepeats);
            allSyncWidths = [allSyncWidths, syncWidths];
            allSyncWidthApproxes = cat(3, allSyncWidthApproxes, syncWidthApproxes);
            allSyncWidthEmpirical = [allSyncWidthEmpirical, syncWidthEmpirical];
            allDominantEigenvalues = [allDominantEigenvalues, dominantEigenvalues];
            allSecondEigenvalues = [allSecondEigenvalues, secondEigenvalues];
        end
    end
    if (totalRepeats == 0)
        error('No results found - did you set parameters.combineResultsFrom correctly (%s)?\nFilename suffix sought is %s', ...
            parameters.combineResultsFrom, fileNameSuffix);
    end
    % All results have been sorted. So now rename them and save into the
    % new file. Also change the final filename for the correct number of
    % repeats
    fileNameSuffix = strrep(fileNameSuffix, sprintf('-repeats*', repeats), sprintf('-repeats%d', totalRepeats));
    repeats = totalRepeats;
    syncWidths = allSyncWidths;
    syncWidthApproxes = allSyncWidthApproxes;
    syncWidthEmpirical = allSyncWidthEmpirical;
    dominantEigenvalues = allDominantEigenvalues;
    secondEigenvalues = allSecondEigenvalues;
    
    % Ready to save - let's make sure the folder exists first:
    if (~exist(parameters.folder, 'dir'))
        fprintf('Creating folder %s as it did not exist\n', parameters.folder);
        mkdir(parameters.folder);
    end
    combinedFilename = [parameters.folder, '/', fileNameSuffix, '.mat'];
    save(combinedFilename, '-mat', 'N', 'd', 'b', 'c', 'p', 'undirected', 'maxMotifLength', 'discretized', ...
         'networkType', 'paramsToRunThrough', 'S', 'repeats', ...
         'syncWidths', 'syncWidthApproxes', 'syncWidthEmpirical', ...
         'dominantEigenvalues', 'secondEigenvalues', 'parameters');
    fprintf('Saving results from %d total repeats to %s\n', repeats, combinedFilename);

    toc
end


end
