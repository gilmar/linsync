function resultsDir = setupMousePaths()
%SETUPMOUSEPATHS  Add linsync toolkit, the epileptor experiment helpers,
%and this experiment folder to the MATLAB path.
%
%   The mouse-epilepsy experiment reuses the 1-D Epileptor dynamics from
%   epileptor-experiment/ (oneDepileptor, CouplingMatrix, normal,
%   covariancesGaussianNet from linsync/), so we add both folders to the
%   path here.
%
%   Returns resultsDir -- full path to the results/ subfolder
%   (created if missing).

experimentRoot = fileparts(mfilename('fullpath'));
toolkitRoot    = fileparts(experimentRoot);
epileptorRoot  = fullfile(toolkitRoot, 'epileptor-experiment');

addpath(toolkitRoot);
if exist(epileptorRoot, 'dir')
    addpath(epileptorRoot);
else
    warning('setupMousePaths:NoEpileptorExperiment', ...
        'epileptor-experiment folder not found at %s. oneDepileptor / CouplingMatrix / normal will not resolve.', ...
        epileptorRoot);
end
addpath(experimentRoot);

resultsDir = fullfile(experimentRoot, 'results');
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end
end
