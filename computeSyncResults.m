function computeSyncResults(varargin)
%% function computeSyncResults(properties)
% function computeSyncResults(N, b, c, undirected, discretized, repeats, networkType, p, d, S, maxMotifLength, folder, MaxK, dt)
%
% Generate networks and compute the synchronizability for a number of networks assuming gaussian dynamics
%  (discrete time AR or continuous time Ornstein-Uhlenbeck as specified in input).
% Save the results for each network to a file.
%
% Inputs:
% See the specs of inputs in the parametersTemplate.m script.
% Parameters can be supplied in one of two ways:
% - Option 1 -- (a filename string or object)
%   - parameters - an object containing the expected properties (as outlined below), or a string
%     describing the filename to run load this object in
% - Option 2 -- DEPRECTATED -- (all supplied individually)
%   - definitions of individual variables are as specified in the parametersTemplate.m script.
%     We require the following to be defined in order:
%     N, b, c, undirected, discretized, repeats, networkType, p, d, S, maxMotifLength, folder
%     The following are optional:
%     MaxK, dt, randSeed
%
%% Linear Sync Toolkit (linsync)
% Copyright (C) 2023 Joseph T. Lizier
% Distributed under GNU General Public License v3

%% Preliminaries:

% Parse the vargin to have each parameter in its own variable (See list
% above), and to have parameters defined.
parseParameters;
% Note: be careful not to refer to parameters until referring to a member
% of it, else on the cluster if we have parameters.m defined, Matlab can
% think that we're referring to that (weird error)

tic;

I = eye(N);
G = ones(N)./N;
U = I - G;

% Pull out the array of parameters that we will sweep through here:
paramsToRunThrough = parameters.(parameters.tosweep);
indices = length(paramsToRunThrough);

% Create storage for the means
syncWidths = zeros(indices, repeats);
syncWidthApproxes = zeros(indices, length(parameters.motifLengthsToCheck), repeats);
syncWidthEmpirical = zeros(indices, repeats);
dominantEigenvalues = zeros(indices, repeats);
secondEigenvalues = zeros(indices, repeats);
diagonalizable = zeros(indices, repeats);

% Seed the random number generator
rng(randSeed);

%% Main loop over parameter we are sweeping whilst running experiments:
for paramIndex = 1 : indices
    fprintf('Updating parameter %s(%d)=%.4f:\n===============\n\n', ...
        parameters.tosweep, paramIndex, paramsToRunThrough(paramIndex));

    % Assign this particular parameter in the parameters object:
    parameters.(parameters.tosweep) = paramsToRunThrough(paramIndex);

    % Now run experiments for each sample network for this parameter
    for r = 1 : repeats
        % Generate a new sample network.
        err = 2; % So the loop runs until a connected network is sampled
        while (err == 2)
            % Generate a new unweighted adjacency matrix A representing the network:
            fprintf('Generating matrix %d\n', r);

            % Generate an unweighted NxN network structure, without self-connections,
            %  which we ensure is not *undirectionally* connected:
            % (hard coding these instead of allowing user to set them)
            parameters.allowSelf = false;
            parameters.ensureConnected = true;
            parameters.randRing.includeSelf = false;
            A = feval(parameters.generateNetworkFunction, parameters);

            % Generate the weighted update matrix C from A (in row-vector form following Barnett).
            % For the no-delay case, C is NxN, with delays it will be NxNx(tau+1)
            C = feval(parameters.weightTheNetworkFunction, A, parameters);

            % Generate the projected covariance matrix (UcovarianceU == U /Omega U) for it
            [sortedLambdasBU, UcovarianceXU, B, err] = covarianceUGaussianNet(C, discretized, MaxK, false, 1);
            % if err==2 then the UcovarianceU matrix failed to converge and we'll loop again
        end

        % compute <\sigma^2> (syncWidth) for this sample network
        syncWidth = synchronizability(UcovarianceXU) % Letting this print to std out for logging
        syncWidths(paramIndex, r) = syncWidth;

        % Now check what the sync width approximation by low order motifs looks like here:
        % Motif length means:
        %  -- for discrete, the length u of both walks (so a walk length of 2 means size 4 motif);
        %  -- for continuous, the actual motif size m (composed of walk
        %       length u and m-u)
        for kIndex = 1 : length(parameters.motifLengthsToCheck)
            k = parameters.motifLengthsToCheck(kIndex);
            % verbose = -1 to suppress warnings on lack of convergence, since we're only asking for a low order approximation
            [~, UcovarianceUApprox, ~, ~] = covarianceUGaussianNet(C, discretized, k, true, -1, true);
            syncWidthApproxes(paramIndex, kIndex, r) = synchronizability(UcovarianceUApprox);
        end
        
        % Compute sync width empirically (<\sigma^2>_E) from simulations of the dynamics
        %  on this network structure
        if (S ~= 0)
            % The time samples are taken dt time units apart (for
            % continuous time)
            empiricalCovarianceProjected = empiricalCovariancesProjected(B, N, discretized, S, dt);
        else
            % Don't run any empirical calculations:
            empiricalCovarianceProjected = 0;
        end
        empiricalSyncWidth = trace(empiricalCovarianceProjected) ./ N;
        syncWidthEmpirical(paramIndex, r) = empiricalSyncWidth;
        
        % Grab the dominant eigenvalue component here, depending on the criteria (magnitude or real component):
        if (discretized)
            % For discrete time we want that with the largest magnitude. (should be < 1 for stability)
            % Not sure how the sort function treated complex values before, so re-sort:
            sortedEigsByCriteria = sort(abs(sortedLambdasBU));
        else
            % For continuous time we want that with the largest real component. (should be < 1 for stability)
            % Not sure how the sort function treated complex values before, so re-sort:
            sortedEigsByCriteria = sort(real(sortedLambdasBU));
        end
        dominantEigenvalue = sortedEigsByCriteria(length(sortedEigsByCriteria));
        dominantEigenvalues(paramIndex, r) = dominantEigenvalue;
        % Retaining second eigenvalue also in case more is seen here:
        secondEigenvalue = sortedEigsByCriteria(length(sortedEigsByCriteria) - 1);
        secondEigenvalues(paramIndex, r) = secondEigenvalue;
        
        if (parameters.checkDiagonalizable)
            % And finally check whether the weighted adjacency matrix was
            % diagonlizable:
            diagonalizable(paramIndex, r) = isdiagonalizable(B);
        else
            diagonalizable(paramIndex, r) = nan;
        end

        fprintf('Repeat %d for parameter(%d)=%.4f: sync_width=%.1f, empirical=%.1f, approx(1)=%.1f, approx(2)=%.1f, approx(3)=%.1f, lambda_1=%.4f, lambda_2=%.4f, isdiag=%d\n\n', ...
            r, paramIndex, paramsToRunThrough(paramIndex), syncWidth, empiricalSyncWidth, syncWidthApproxes(paramIndex, 1, r), syncWidthApproxes(paramIndex, 2, r), ...
            syncWidthApproxes(paramIndex, 3, r), dominantEigenvalue, secondEigenvalue, diagonalizable(paramIndex, r));
    end
    
    fprintf('====\n**All repeats for parameter %.4f completed: <sync_width>=%.3f, <empirical>=%.3f, <approx(1)>=%.3f, <approx(2)>=%.3f, <approx(3)>=%.3f, <lambda_1>=%.4f, <lambda_2>=%.4f, <isdiag>=%.3f\n\n', ...
        paramsToRunThrough(paramIndex), mean(syncWidths(paramIndex,:)), ...
        mean(syncWidthEmpirical(paramIndex,:)), mean(syncWidthApproxes(paramIndex, 1, :)), ...
        mean(syncWidthApproxes(paramIndex, 2, :)), mean(syncWidthApproxes(paramIndex, 3, :)), ...
        mean(dominantEigenvalues(paramIndex,:)), mean(secondEigenvalues(paramIndex,:)), ...
        mean(diagonalizable(paramIndex, :)));
    toc
end                

%% Save the processed results here:
[networkType, netTypeSuffix] = generateNetworkTypeStrings(parameters);
% This will generate filename prefix with the *final* value of the swept
%  parameter, but parameter being swept will be indicated
fileNamePrefix = sprintf('%s/N%d-%s%s-b%.2f-c%.2f-p%.4f-sweep_%s-%s-k%d-%s-S%d-repeats%d', ...
                folder, parameters.N, networkType, netTypeSuffix, parameters.b, parameters.c, parameters.p, ...
                parameters.tosweep, undirString, parameters.maxMotifLength, discString, S, repeats);

% Convert booleans to integers so Matlab can save them
if (discretized)
    discretized = 1;
else
    discretized = 0;
end
if (undirected)
    undirected = 1;
else
    undirected = 0;
end

parameters = originalParameters; % Restore the original parameters

% Ready to save - let's make sure the folder exists first:
if (~exist(folder, 'dir'))
    fprintf('Creating folder %s as it did not exist\n', folder);
    mkdir(folder);
end
save([fileNamePrefix, '.mat'], '-mat', 'N', 'd', 'b', 'c', 'p', 'undirected', ...
    'motifLengthsToCheck', 'maxMotifLength', 'discretized', ...
    'networkType', 'paramsToRunThrough', 'S', 'repeats', ...
    'syncWidths', 'syncWidthApproxes', 'syncWidthEmpirical', ...
    'dominantEigenvalues', 'secondEigenvalues', 'diagonalizable', ...
    'parameters'); % This is the important one, the other inputs are legacy
toc

end
