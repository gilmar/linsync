function [sortedLambdasAfterProj, UcovarianceXU, B, err] = covarianceUGaussianNet(C, discreteTime, maximumIterations, forcePowerSeriesForSymmetric, verbose, skipPowerSeriesConvergenceCheck)
%
% Computes the eigenvalues of C, and the relevant projected covariance matrix (U^T \Omega U) for the given network.
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
% - forcePowerSeriesForSymmetric - force the use of power series for UcovarianceU even if the
%     connectivity matrix is symmetric
% - verbose - level of verbosity for discreteCon2CovProjected and contCon2CovProjected
%     (set to -1 if only using a small number of iterations, to avoid
%     warnings)
% - skipPowerSeriesConvergenceCheck - to skip checking that max eigenvalue
%     lies within unit circle (only has an impact for continuous time
%     dynamics). Default is false.
% 
% Outputs
% - sortedLambdasAfterProj - sorted eigenvalues of CU (if C is NxN) or 
%    BU (if including delays), (from smallest to largest magnitude if complex)
% - UcovarianceXU - projected covariance matrix U^T \Omega_X U, of system X
%     (whether we are using standard delay only, or multiple delays -- 
%      i.e. this is not the project covariance for embedded system Z
%      if we are dealing with delays)
%     for U = I - 1/N and
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
    if (~discreteTime)
        error('We do not handle delays for continuous time yet');
    end
end
% For system X:
I = eye(N);
G = ones(N) / N;
U = I - G;

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
    skipPowerSeriesConvergenceCheck = false;
end

% Construct the embedded matrix B (if no delays B = C)
B = embedNetworkWithDelays(C);
if hasDelays
    zLength = size(B,1);
    IZ = eye(zLength);
    % IX = zeros(zLength); % not required later
    % IX(1:N,1:N) = I; % not required later
    GX = [repmat(G, 1, tauPlus1); zeros(N*(tauPlus1-1), N*tauPlus1)];
    UX = IZ - GX; % not required later
end

% Compute eigenvalues of C * U (no delays) or B * UX (delayed case):
%  (since B is square, the eigenvalues are the same as B^T - i.e. it doesn't matter that they correspond to row/column vectors)
if (hasDelays)
    lambdasBU = eig(B * UX);
else
    lambdasBU = eig(C * U); % Would be equivalent to eig(B * UX); but stepping it out explicitly
end

sortedLambdasAfterProj = sort(lambdasBU); % Sorts the elements by magnitude
% Find the minimum and maximum eigenvalues
lambdaAfterProjMax = sortedLambdasAfterProj(size(sortedLambdasAfterProj,1));

% Check for stationarity for B*U (generalised discrete) or C*U (continuous time)
if (discreteTime)
    if (abs(lambdaAfterProjMax) >= 1)
        % Non-stationary
        save('nonconvergentNetwork.mat', '-mat', 'C'); % save for later investigation
        error('Discrete system with |\\lambda_BU_max| >= 1 (%.6f) i.e. non-stationary of C*U or B*U\n', abs(lambdaAfterProjMax));
    end
else
    sortedRealPartsOfLambdasCU = sort(real(lambdasBU));
    realLambdaCUMax = sortedRealPartsOfLambdasCU(size(sortedRealPartsOfLambdasCU,1));
    if (realLambdaCUMax >= 1)
        % Non-stationary
        save('nonconvergentNetwork.mat', '-mat', 'C'); % save for later investigation
        error('Continuous system with Max(Re(\\lambda_CU)) >= 1 (%.6f), i.e. non-stationary of C*U\n', realLambdaCUMax);
    end
end

% Check for convergence of the projected covariance matrix: (same condition for both continuous and discrete):
if  (~skipPowerSeriesConvergenceCheck && (abs(lambdaAfterProjMax) >= 1))
    save('nonconvergentNetwork.mat', 'C'); % save for later investigation
    error('|\\lambda_CU_max| >= 1 (%.2f) implies that the U^T \\Omega U matrix will not converge.\n', abs(lambdaAfterProjMax));
end

% Check whether matrix B is symmetric or not, to help speed up the
% projected covariance calculation:
symmetric = isempty(find(B - B' > tol*eps, 1));

% Now compute the UcovarianceU matrix
if (discreteTime)
    if (hasDelays)
        % Case with delays
        [UcovarianceXU,err] = discreteCon2CovProjectedWithDelays(B,maximumIterations,tol,verbose,N);
    else
        % Standard case without delays:
        [UcovarianceXU,err] = discreteCon2CovProjected(B,symmetric && ~forcePowerSeriesForSymmetric,maximumIterations,tol,verbose);
    end
else
    [UcovarianceXU,err] = contCon2CovProjected(C,symmetric && ~forcePowerSeriesForSymmetric,maximumIterations,tol,verbose);
end

% Can handle err == 2 here, or just let it go through to the caller

end

