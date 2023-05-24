function C = weightNetworkStandardWithFixedDelay(A, parameters)
%
% Weight the given adjacency matrix with standard self and cross-coupling
% strengths, with the cross-couplings applied at a fixed delay
%
% - Inputs
%   - A - the unweighted adjacency matrix for an NxN matrix (assumed to
%   have entries 0 or 1, and none on diagonal).
%   - parameters - object containing the experimental parameters, including:
%       - parameters.b - total weight summed into each target
%       - parameters.c - total cross-connection weight
%       - parameters.delay_fixedCross - fixed cross coupling delay
%
% Outputs
% - C - NxNxtauPlus1 weighted connectivity matrix (can be directed; C(i,j) means a link exists from i->j)
%       where the cross-coupling is applied with additional delay parameters.delay.fixedCross
%       (compared to standard delay of "0")
%
%% Linear Sync Toolkit (linsync)
% Copyright (C) 2023 Joseph T. Lizier
% Distributed under GNU General Public License v3

I = eye(parameters.N);
tau = parameters.delay_fixedCross; % Additional delay beyond standard delay 1
tauPlus1 = tau + 1;

if (tau == 0)
    C = weightNetworkStandard(A, parameters);
    return;
end

% Else we have a non-trivial delay here

% Compute the in-degrees for each node (take a column sum)
D = diag(sum(A));

% Compute embedded connectivity matrix:
%  (notice how b-c is the self-weight, and c is the total
%   weight from d other inputs, which each have c/d (where d might be different for each node).
%   Equal weights c/d are not required by the maths, but used for simply experiments here.)
C = zeros(N,N,tauPlus1);
C(:,:,1) = (parameters.b - parameters.c) .* I; % Self-connections at standard lag 1
C(:,:,tauPlus1) = parameters.c .* A  / D;

end
