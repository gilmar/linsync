function [UcovarianceXU,err] = discreteCon2CovProjectedWithDelays(B,maxi,tol,verbose,N)
%
% This code computes the projected covariance matrix for the discrete-time
% case, with delays (without delays is handled in
% discreteCon2CovProjected).
% The code converts con2cov from ncomp from Barnett et al to the discrete case
%  and for the projection, and for delays.
%  See reference for the Barnett code:
%   - ncomp_tools: L. Barnett, C. L. Buckley and S. Bullock (2009)
%      (licensed to "use as you wish")
%   - L. Barnett, C. L. Buckley and S. Bullock (2009)
%       On Neural Complexity and Structural Connectivity, Physical Review E, 79, 051914
%
% Calculates the projected covariance matrix
%    \Omega_XU = U \Omega_X U
%              = [U_X \Omega_Z U_X]_[1,1]
% for the multivariate
% discrete AR process (where Z may have eigenvector \psi_0 but is otherwise stationary):
%
%     Z(t+1) = Z(t)*B + R_X(t)
%
% where X(t) is the first N components of Z(t) (rest are embedded values of
% X(t)), and
% where R_X(t) is uncorrelated mean-zero unit-variance Gaussian noise
%  for the first N terms and zero otherwise,
%  and B is the update matrix for the system Z.
% The projection is orthogonal to the fully synchronized state vector psi_Z0.
%
% B embeds the delayed update matrices C into an update for system Z.
%  Use of row vectors follows the format of Barnett et al..
%
% Inputs
% - B - connectivity matrix for the embedded system Z (including self-connection weights) for row-vectors
% - maxi - maximum number of iterations.
% - tol - tolerance for concluding convergence has been reached.
%  These parameters are used the same way as it the ncomp toolbox of Barnett et al:
% "
%     'maxi' gives the maximum iterations for the non-symmetric version.
%     Convergence criterion is to stop when the ratios of the added terms
%     approach the machine epsilon with given tolerance 'tol'. If maximum
%     iterations are exceeded a warning is issued. Failure to converge
%     is treated as an error. You can play off accuracy against no. of iterations.
%     Try maxi = 1000, tol = 10 for starters. If it converges you should get a max.
%     relative error of order 1e-15 or better.
% "
% - verbose - whether to report warnings and the number of iterations and convergence results
%   - -1 - minimal, output errors only, not warnings or number of iterations and convergence results
%        (useful for batch runs with low order approximations)
%   - 0 (or false) - default, output errors and warnings but not number of iterations and convergence results
%   - 1 - verbose, output errors, warnings, and number of iterations and convergence results
% - N - number of nodes of the original system X that is embedded in system Z
%
% Outputs
% - UcovarianceU - covariance matrix of the X process (sub-component of Z) projected into the orthonormal
%    space to the fully synchronized state vector psi_0
% - err - 0: no error, 1: warning that maximum iterations were reached, 2: covergence failed
%
%% Linear Sync Toolkit (linsync)
% Copyright (C) 2023 Joseph T. Lizier
% Distributed under GNU General Public License v3

err = 0;

% For system X:
I = eye(N);
% Averaging projection operator:
G = 1 ./ N * ones(N);
% Unaveraging operator:
U = I - G;

% For embedded system Z:
zLength = size(B,1);
tauPlus1 = zLength / N;
IZ = eye(zLength);
IX = zeros(zLength);
IX(1:N,1:N) = I;
GX = [repmat(G, 1, tauPlus1); zeros(N*(tauPlus1-1), N*tauPlus1)];
UX = IZ - GX;

% UcovarianceXU = (UX' * IX * UX)(1:N,1:N); % Initialised for u=0 in the sum
% So: UcovarianceXU = (U' * IX(1:N,1:N) * U);
UcovarianceXU = (U' * U);
% leftComponent = UX';
% rightComponent = UX;
rightComponentFirstBlockCol = UX(:,1:N); % First N-block column of right component
% leftMultiplier = UX' * B';
% rightMultiplier = B * UX;

% B should have a good sparse representation because of the many zeroes re the lags
Bsparse = sparse(B);
% BTsparse = sparse(B');

for i = 1:maxi
    % d_UcovarianceU holds the previous term added in (for i > 1)
    % UcovarianceXU holds the sum of previous projected terms
    %
    % Try this a faster way:
    %  leftComponent should be U^T (B^i)^T = U^T B^T^i = U^T (B^T)^(i-1) B^T
    %  So we can just right sparse multiply by B, and no need to apply U
    %  again.
    % leftComponent = leftComponent * BTsparse;
    % rightComponent = rightComponent * rightMultiplier; % It's just leftComponent'
    %
    % Faster again: we only need to maintain the first N-block column of rightComponent
    %  or first N-block row of leftComponent, because these are all that is
    %  ever used to update the top left NxN block (which is all we need for
    %  UcovarianceXU:
    rightComponentFirstBlockCol = Bsparse * rightComponentFirstBlockCol;

    % We want: d_UcovarianceXU = (leftComponent * IX * rightComponent)(1:N,1:N); so
    % So: d_UcovarianceXU = (leftComponent(1:N,:) * IX * rightComponent(:,1:N));
    % But only IX(1:N,1:N) is non-zero so: d_UcovarianceXU = (leftComponent(1:N,1:N) * rightComponent(1:N,1:N));
    % And leftComponent is just the transpose of rightComponent:
    d_UcovarianceXU = (rightComponentFirstBlockCol(1:N,:)' * rightComponentFirstBlockCol(1:N,:));
    if (mod(i, 100) == 0) && negligible(d_UcovarianceXU,UcovarianceXU,tol)
        % Save time on calls to negligible by only checking every 100 steps
        break
    end
    UcovarianceXU = UcovarianceXU + d_UcovarianceXU;
end

mre = maxrelerr(d_UcovarianceXU,UcovarianceXU);

if isnan(mre) || isinf(mre)
    fprintf(2,'ERROR (discreteCon2CovProjectedWithDelays): failed to converge\n');
    err = 2;
    return;
end

if (i == maxi) && (verbose >= 0)
    mac = mean(mean(abs(UcovarianceXU)));
    fprintf(1,'WARNING (discreteCon2CovProjectedWithDelays): exceeded maximum iterations (%d): mean abs. covariance = %.3f, max. relative error = %.3e (both in projected space)\n',maxi,mac,mre);
    err = 1;
    return;
end

if (verbose > 0)
    fprintf(1,'discreteCon2CovProjectedWithDelays: iterations = %i, max. relative error = %d\n',i,mre);
end
