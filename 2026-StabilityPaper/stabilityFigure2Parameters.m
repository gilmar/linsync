%% Script to set up the parameters object for stability Figure 2 (Liao & Lizier 2026):
%   D_st vs. rewiring probability p on a Watts-Strogatz ring network.
% You can not only assign values here, but have differential processing
% (e.g. to do different things on your desktop or cluster).
% Required members are described as they appear below.
%
%% Linear Sync Toolkit (linsync)
% Copyright (C) 2023 Joseph T. Lizier
% Distributed under GNU General Public License v3

clear parameters; % In case parameters previously held the name of the parameters file

% You can use this boolean to flag different folders for local and cluster runs.
if exist('/home/joseph/', 'dir')
    isCluster = false;
elseif ~isempty(getenv('PBS_JOBID'))
    isCluster = true;
else
    isCluster = false;
end

% Set the location of the code library - only required for
% runComputeSyncResults which is designed for cluster use
if isCluster
    parameters.syncToolkitPath = 'syncToolkit';
else
    parameters.syncToolkitPath = '..'; % Toolkit is in the parent folder
end

% - N - Network size
parameters.N = 100;

% - b - Self / row-sum coupling weight (see weightNetworkStandard).
parameters.b = 0.95;

% - c - Total cross-coupling strength per node (see weightNetworkStandard).
parameters.c = 0.7;

% - d - Degree for randRing: d/2 links on either side of each node.
parameters.d = 4;

% - undirected - whether the network is undirected.
parameters.undirected = false;

% - discretized - discrete-time AR (true) or continuous-time OU (false).
parameters.discretized = false;

% - computationMode - 'stability' for D_st via covariancesGaussianNet (full Omega).
parameters.computationMode = 'stability';

% - p - rewiring probability for randRing (can be an array to sweep).
parameters.p = [0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0];

% - tosweep - parameter field name to sweep in computeSyncResults
parameters.tosweep = 'p';
% - tosweep_label - axis label for plots
parameters.tosweep_label = 'p';

% - generateNetworkFunction / weightTheNetworkFunction
parameters.generateNetworkFunction = 'generateNewRandomRingMatrix';
parameters.weightTheNetworkFunction = 'weightNetworkStandard';

% - repeats - networks per parameter value (paper uses 2000; smaller for tests).
parameters.repeats = 20;

% - S - empirical sample count; 0 skips empirical covariance.
parameters.S = 0;
parameters.SRangeToPlot = 0;

% - maxMotifLength / motifLengthsToCheck (used for filenames and parseParameters;
%   motif approximations are skipped when computationMode is 'stability').
parameters.maxMotifLength = 50;
parameters.motifLengthsToCheck = [2, 10, 50];

parameters.MaxK = 100000;
parameters.randSeed = 1;
parameters.checkDiagonalizable = false;

parameters.dt = 1;

% - folder - where computeSyncResults writes the .mat file
if isCluster
    parameters.folder = '2026-StabilityPaper/results/N100-randRing-d4-b0.95-c0.70-dir-k50-cont/[@P1]';
else
    parameters.folder = '2026-StabilityPaper/results';
end

% combineResultsFrom - for cluster combineSyncResults workflows
parameters.combineResultsFrom = './results/clusterResults/';
