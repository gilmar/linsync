function [sortedLambdas, covariance, B, err] = covarianceUGaussianNet(C, discreteTime, maximumIterations, forcePowerSeriesForSymmetric, verbose, skipEigsCalculateAndCheck)
%
% Computes the eigenvalues of C, and the covariance matrix (\Omega) for the given network.
%
% Inputs
% - C - Either:
%     - an NxN connectivity matrix, including self-connection weights (in the Cij
%       = i -> j format, i.e. for row-vectors), for standard connections at the
%       standard delay (being regular delay 1 for discrete time or no delay
%       for continuous time).
%     - or an NxNx(tau+1) matrix representing coupling matrices of C(n)
%       from discrete time steps 1+n in the past (for n=0:tau).
%       Each C(i,j,n) matrix represents connections Cij = i->j as before but
%       over lag 1+n. Can have non-zero weights at multiple n for each
%       i->j.
%       Only valid for discreteTime=true at present.
% - discreteTime - whether the time-series dynamics from a discrete-time AR process (true)
%    or continuous time Ornstein-Uhlenbeck process.
% - maximumIterations - the number of components to add into the power series for UcovarianceU. Default is 1000
% - forcePowerSeriesForSymmetric - force the use of power series for covariance even if the
%     connectivity matrix is symmetric
% - verbose - level of verbosity for discreteCon2Cov and con2cov
%     (set to -1 if only using a small number of iterations, to avoid
%     warnings)
% - skipEigsCalculateAndCheck - to skip computing eigenvaluces and checking conditions.
%      Default is false. Should only be used for debugging, or calls where
%      this connectivity matrix has already been checked (e.g. for low
%      order approximations in computeSyncResults)
% 
% Outputs
% - sortedLambdas - sorted eigenvalues of C (if C is NxN) or 
%    B (if including delays), (from smallest to largest magnitude if
%    complex).
% - covariance - covariance matrix \Omega_X, of system X
%     (whether we are using standard delay only, or multiple delays -- 
%      i.e. this is not the project covariance for embedded system Z
%      if we are dealing with delays)
%     \Omega is the covariance matrix of X (not Z), assuming the Ornstein-Uhlenbeck or VAR process.
% - B - the matrix is just the C matrix is the case of no-delays (C is
%     NxN), else if C is NxNx(tau+1), then each NxN matrix is embedded into
%     B which is the connectivity matrix for the embedded system Z.
% - err - error status from the call to work out the covariance matrix here
%
%% Linear Sync Toolkit (linsync)
% Copyright (C) 2023 Joseph T. Lizier
% Distributed under GNU General Public License v3

N = size(C,1);
hasDelays = (length(size(C)) == 3);
if (hasDelays)
    tauPlus1 = size(C,3); % Number of delays we consider. If > 1 we're looking beyond standard case
    if (~discreteTime)
        error('We do not handle delays for continuous time yet');
    end
else
    tauPlus1 = 1; % Only the standard delay here
end
% For system X:
I = eye(N);

if (nargin < 3)
    maximumIterations = 1000; % This will be on the low side, should be setting something higher for proper experiments.
end
tol = 100; % Multiples of machine epsilon within which we want to consider a value to be zero.

if (nargin < 4)
    forcePowerSeriesForSymmetric = false;
end

if (nargin < 5)
    verbose = 0;
end

if (nargin < 6)
    skipEigsCalculateAndCheck = false;
end

% Construct the embedded matrix B (if no delays B = C)
B = embedNetworkWithDelays(C);
if hasDelays
    zLength = size(B,1);
    IZ = eye(zLength);
    % IX = zeros(zLength); % not required later
    % IX(1:N,1:N) = I; % not required later
end

if (~skipEigsCalculateAndCheck)
    % Compute eigenvalues of C (no delays) or B (delayed case):
    %  (since B is square, the eigenvalues are the same as B^T - i.e. it doesn't matter that they correspond to row/column vectors)
    if (hasDelays)
        % Could optimise using eigs to select only top 2, but this function
        % seems unstable; asking for all instead (this isn't the bottleneck
        % anymore)
        lambdasB = eig(B);
    else
        lambdasB = eig(C); % Would be equivalent to eig(B); but stepping it out explicitly
    end
    
    sortedLambdas = sort(lambdasB); % Sorts the elements by magnitude
    % Find the maximum eigenvalue
    lambdaMax = sortedLambdas(size(sortedLambdas,1));
    
    % Check for stationarity for B (generalised discrete) or C (continuous time)
    if (discreteTime)
        if (abs(lambdaMax) >= 1)
            % Non-stationary
            save('nonconvergentNetwork.mat', '-mat', 'C'); % save for later investigation
            error('Discrete system with |\\lambda_B_max| >= 1 (%.6f) i.e. non-stationary of C or B\n', abs(lambdaMax));
        end
    else
        sortedRealPartsOfLambdas = sort(real(lambdasB));
        realLambdaBMax = sortedRealPartsOfLambdas(size(sortedRealPartsOfLambdas,1));
        if (realLambdaBMax >= 1)
            % Non-stationary
            save('nonconvergentNetwork.mat', '-mat', 'C'); % save for later investigation
            error('Continuous system with Max(Re(\\lambda_B)) >= 1 (%.6f), i.e. non-stationary of B\n', realLambdaBMax);
        end
    end
    
    % Check for convergence of the projected covariance matrix: (same condition for both continuous and discrete):
    if (abs(lambdaMax) >= 1)
        save('nonconvergentNetwork.mat', 'C'); % save for later investigation
        error('|\\lambda_B_max| >= 1 (%.2f) implies that the \\Omega matrix will not converge.\n', abs(lambdaMax));
    end
else
    sortedLambdas = []; % Not required for this call
end

% Check whether matrix B is symmetric or not, to help speed up the
% projected covariance calculation:
symmetric = isempty(find(B - B' > tol*eps, 1));

% Now compute the UcovarianceU matrix
if (discreteTime)
    if (hasDelays)
        % Case with delays
        error('The discreteCon2CovWithDelays function is not yet defined')
        [covariance,err] = discreteCon2CovWithDelays(B,maximumIterations,tol,verbose,N);
    else
        % Standard case without delays:
        [covariance,err] = discreteCon2Cov(B,symmetric && ~forcePowerSeriesForSymmetric,maximumIterations,tol,verbose);
    end
else
    [covariance,err] = con2cov(C,symmetric && ~forcePowerSeriesForSymmetric,maximumIterations,tol,verbose);
end

% Can handle err == 2 here, or just let it go through to the caller

end

