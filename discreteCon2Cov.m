function [M,err] = discreteCon2Cov(C,useSymmetricForm,maxi,tol,verbose)
%
% This code computes the covariance matrix for the discrete-time case.
% converting con2cov from ncomp from Barnett et al to the discrete case.
%  See reference for the Barnett code:
%   - ncomp_tools: L. Barnett, C. L. Buckley and S. Bullock (2009)
%      (licensed to "use as you wish")
%   - L. Barnett, C. L. Buckley and S. Bullock (2009)
%       On Neural Complexity and Structural Connectivity, Physical Review E, 79, 051914
%
% Calculates the projected covariance matrix U \Omega U for the multivariate
% discrete AR process (where X may have eigenvector \psi_0 but is otherwise stationary):
%
%     X(t+1) = X(t)*C + R(t)
%
% where R(t) is uncorrelated mean-zero unit-variance Gaussian noise
%  and C is the update matrix.
%
% C(i,j) represents the update effect of i -> j ;
%  Note this follows the format of Barnett et al, where activity is in row vectors.
%
%
% Inputs
% - C - connectivity matrix, including self-connection weights (in the Cij = i -> j format, i.e. for row-vectors)
% - useSymmetricForm - whether C is symmetric or not AND whether we should use the
%       analytic shortcut for symmetric C's. Should set this to false if
%       you want to see approximations to a limited order.
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
%
% Outputs
% - M - covariance matrix of the process
% - err - 0: no error, 1: warning that maximum iterations were reached, 2: covergence failed
%
%% Linear Sync Toolkit (linsync)
% Copyright (C) 2023 Joseph T. Lizier
% Distributed under GNU General Public License v3

err = 0;
N = size(C,1);
I = eye(n);

if (useSymmetricForm)
    M = inv(I - C*C);
    return;
end

M = I;
dM = M;
for i = 1:maxi
    dM = C'*dM*C;
    if negligible(dM,M,tol)
       break
    end
    M = M + dM;
end

mre = maxrelerr(dM,M);

if isnan(mre) || isinf(mre)
    fprintf(2,'ERROR (discreteCon2Cov): failed to converge\n');
    err = 2;
    return;
end

if i == maxi
    mac = mean(mean(abs(M)));
    fprintf(1,'WARNING (discreteCon2Cov): exceeded maximum iterations: mean abs. covariance = %d, max. relative error = %d\n',mac,mre);
    err = 1;
    return;
end

if verbose
    fprintf(1,'discreteCon2Cov: iterations = %i, max. relative error = %d\n',i,mre);
end

