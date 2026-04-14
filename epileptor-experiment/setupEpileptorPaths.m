function resultsDir = setupEpileptorPaths()
%SETUPEPILEPTORPATHS  Add linsync toolkit root and this experiment folder to the MATLAB path.
%
%   Call once at the start of each experiment script so covariancesGaussianNet,
%   con2cov, etc. resolve from the parent linsync directory, and local
%   CouplingMatrix, loadPatientWeights, etc. resolve from epileptor-experiment/.
%
%   Returns resultsDir — full path to the results/ subfolder (created if missing).

experimentRoot = fileparts(mfilename('fullpath'));
toolkitRoot = fileparts(experimentRoot);
addpath(toolkitRoot);
addpath(experimentRoot);

resultsDir = fullfile(experimentRoot, 'results');
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end
end
