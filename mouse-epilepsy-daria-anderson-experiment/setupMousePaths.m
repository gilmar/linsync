function resultsDir = setupMousePaths(varargin)
%SETUPMOUSEPATHS  Add linsync toolkit, BCT, and this experiment folder to the MATLAB path.
%
%   Shared Epileptor dynamics and stability-centrality helpers (oneDepileptor, CouplingMatrix,
%   normal, findCriticalX0, computeStabilityCentralities, etc.) live in the
%   linsync root and are added via toolkitRoot below.
%
%   resultsDir = setupMousePaths()
%       -> .../mouse-epilepsy-daria-anderson-experiment/results/
%
%   resultsDir = setupMousePaths('ExperimentName', 'initial_column')
%       -> .../results/initial_column/
%
%   Name-value options:
%     ExperimentName  -- subfolder under results/ for this experiment run

experimentRoot = fileparts(mfilename('fullpath'));
toolkitRoot    = fileparts(experimentRoot);
bctRoot        = fullfile(toolkitRoot, '2019_03_03_BCT');

addpath(toolkitRoot);
if exist(bctRoot, 'dir')
    addpath(bctRoot);
else
    warning('setupMousePaths:NoBCT', ...
        ['Brain Connectivity Toolbox not found at %s. betweenness_wei and ' ...
         'distance_wei will be unavailable (stability centralities still run). ' ...
         'Install instructions: linsync/docs/BCT.md'], bctRoot);
end
addpath(experimentRoot);

p = inputParser;
addParameter(p, 'ExperimentName', '', @(s) ischar(s) || isstring(s));
parse(p, varargin{:});
expName = char(p.Results.ExperimentName);

baseResults = fullfile(experimentRoot, 'results');
if ~exist(baseResults, 'dir')
    mkdir(baseResults);
end

if isempty(expName)
    resultsDir = baseResults;
else
    resultsDir = fullfile(baseResults, expName);
    if ~exist(resultsDir, 'dir')
        mkdir(resultsDir);
    end
end
end
