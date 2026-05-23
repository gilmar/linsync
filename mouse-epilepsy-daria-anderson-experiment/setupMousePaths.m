function resultsDir = setupMousePaths()
%SETUPMOUSEPATHS  Add linsync toolkit, BCT, and this experiment folder to the MATLAB path.
%
%   Shared Epileptor dynamics and §4.5 helpers (oneDepileptor, CouplingMatrix,
%   normal, findCriticalX0, computeStabilityCentralities, etc.) live in the
%   linsync root and are added via toolkitRoot below.
%
%   Returns resultsDir -- full path to the results/ subfolder
%   (created if missing).

experimentRoot = fileparts(mfilename('fullpath'));
toolkitRoot    = fileparts(experimentRoot);
bctRoot        = fullfile(toolkitRoot, '2019_03_03_BCT');

addpath(toolkitRoot);
if exist(bctRoot, 'dir')
    addpath(bctRoot);
else
    warning('setupMousePaths:NoBCT', ...
        ['Brain Connectivity Toolbox not found at %s. betweenness_wei, ' ...
         'distance_wei, and other BCT functions will be unavailable. ' ...
         'Download BCT from https://sites.google.com/site/bctnet/ and ' ...
         'place it at that path.'], bctRoot);
end
addpath(experimentRoot);

resultsDir = fullfile(experimentRoot, 'results');
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end
end
