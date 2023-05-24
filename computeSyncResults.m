function computeSyncResults(varargin)
%% function computeSyncResults(N, b, c, undirected, discretized, repeats, networkType, p, d, S, maxMotifLength, folder, MaxK, dt)
% function computeSyncResults(properties)
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

varyingP = false;
if (length(p) > 1)
    % we're varying p
    varyingP = true;
    indices = length(p);
    paramsToRunThrough = p;
else
    % assume we're varying c
    indices = length(c);
    paramsToRunThrough = c;
end
if (size(paramsToRunThrough, 2) > size(paramsToRunThrough,1)) % More columns than rows
    % Make it a column vector
    paramsToRunThrough = paramsToRunThrough';
end

% Create storage for the means
syncWidths = zeros(indices, repeats);
syncWidthApproxes = zeros(indices, maxMotifLength, repeats);
syncWidthEmpirical = zeros(indices, repeats);
dominantEigenvalues = zeros(indices, repeats);
secondEigenvalues = zeros(indices, repeats);
diagonalizable = zeros(indices, repeats);

% Seed the random number generator
rng(randSeed);

%% Main loop over parameters to run experiments for (p or c array):
for paramIndex = 1 : indices
    fprintf('Param index %d:\n===============\n\n', paramIndex);

    if (varyingP)
        parameters.p = paramsToRunThrough(paramIndex);
    else
        c = paramsToRunThrough(paramIndex);
        parameters.c = paramsToRunThrough(paramIndex);
    end

    % Now run experiments for each sample network for this parameter
    for r = 1 : repeats
        % Generate a new sample network.
        err = 2; % So the loop runs until a connected network is sampled
        while (err == 2)
            % Generate a new unweighted adjacency matrix A representing the network:
            fprintf('Generating matrix %d\n', r);

            % Generate an unweighted network structure, without self-connections,
            %  which we ensure is not *undirectionally* connected:
            % (hard coding these instead of allowing user to set them)
            parameters.allowSelf = false;
            parameters.ensureConnected = true;
            parameters.randRing.includeSelf = false;
            A = feval(parameters.generateNetworkFunction, parameters);

            % Generate the weighted update matrix C from A (in row-vector form following Barnett).
            C = feval(parameters.weightTheNetworkFunction, A, parameters);

            % Generate the projected covariance matrix (UcovarianceU == U /Omega U) for it
            [sortedLambdasCU, UcovarianceU, err] = covarianceUGaussianNet(C, discretized, MaxK, false, 1);
            % if err==2 then the UcovarianceU matrix failed to converge and we'll loop again
        end

        % compute <\sigma^2> (syncWidth) for this sample network
        syncWidth = synchronizability(UcovarianceU) % Letting this print to std out for logging
        syncWidths(paramIndex, r) = syncWidth;

        % Now check what the sync width approximation by low order motifs looks like here:
        % Max motif length means:
        %  -- for discrete, the length u of both walks (so a walk length of 2 means size 4 motif);
        %  -- for continuous, the actual motif size m (composed of walk
        %       length u and m-u)
        for k = 1 : maxMotifLength
            % verbose = -1 to suppress warnings on lack of convergence, since we're only asking for a low order approximation
            [~, UcovarianceUApprox, ~] = covarianceUGaussianNet(C, discretized, k, true, -1);
            syncWidthApproxes(paramIndex, k, r) = synchronizability(UcovarianceUApprox);
        end
        
        % Compute sync width empirically (<\sigma^2>_E) from simulations of the dynamics
        %  on this network structure
        if (S ~= 0)
            % The time samples are taken dt time units apart (for
            % continuous time)
            empiricalCovarianceProjected = empiricalCovariancesProjected(C, discretized, S, dt);
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
            sortedEigsByCriteria = sort(abs(sortedLambdasCU));
        else
            % For continuous time we want that with the largest real component. (should be < 1 for stability)
            % Not sure how the sort function treated complex values before, so re-sort:
            sortedEigsByCriteria = sort(real(sortedLambdasCU));
        end
        dominantEigenvalue = sortedEigsByCriteria(length(sortedEigsByCriteria));
        dominantEigenvalues(paramIndex, r) = dominantEigenvalue;
        % Retaining second eigenvalue also in case more is seen here:
        secondEigenvalue = sortedEigsByCriteria(length(sortedEigsByCriteria) - 1);
        secondEigenvalues(paramIndex, r) = secondEigenvalue;
        
        % And finally check whether the weighted adjacency matrix was
        % diagonlizable:
        diagonalizable(paramIndex, r) = isdiagonalizable(C);

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
end                

parameters = originalParameters; % Restore the original parameters

%% Save the processed results here:
[networkType, netTypeSuffix] = generateNetworkTypeStrings(parameters);
if (varyingP)
    % Put full range of p back into the variable p
    p = paramsToRunThrough;
    fileNamePrefix = sprintf('%s/N%d-%s%s-b%.2f-c%.2f-%s-k%d-%s-S%d-repeats%d', ...
                folder, N, networkType, netTypeSuffix, b, c, undirString, maxMotifLength, discString, S, repeats);
else
    % Put full range of c back into the variable c
    c = paramsToRunThrough;
    fileNamePrefix = sprintf('%s/N%d-%s%s-b%.2f-p%.4f-%s-k%d-%s-S%d-repeats%d', ...
                folder, N, networkType, netTypeSuffix, b, p, undirString, maxMotifLength, discString, S, repeats);
end
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

% Ready to save - let's make sure the folder exists first:
if (~exist(folder, 'dir'))
    fprintf('Creating folder %s as it did not exist\n', folder);
    mkdir(folder);
end
save([fileNamePrefix, '.mat'], '-mat', 'N', 'd', 'b', 'c', 'p', 'undirected', 'maxMotifLength', 'discretized', ...
    'networkType', 'paramsToRunThrough', 'S', 'repeats', ...
    'syncWidths', 'syncWidthApproxes', 'syncWidthEmpirical', ...
    'dominantEigenvalues', 'secondEigenvalues', 'diagonalizable', ...
    'parameters'); % This is the important one, the other inputs are legacy
toc

end
