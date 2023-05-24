function C = weightNetworkStandard(A, parameters)
%
% Weight the given adjacency matrix with standard self and cross-coupling
% strengths
%
% - Inputs
%   - A - the unweighted adjacency matrix for an NxN matrix (assumed to
%   have entries 0 or 1, and none on diagonal).
%   - parameters - object containing the experimental parameters, including:
%       - parameters.b - total weight summed into each target
%       - parameters.c - total cross-connection weight
%
% Outputs
% - C - NxN weighted connectivity matrix (can be directed; C(i,j) means a link exists from i->j)
%
%% Linear Sync Toolkit (linsync)
% Copyright (C) 2023 Joseph T. Lizier
% Distributed under GNU General Public License v3

I = eye(parameters.N);

% Compute the in-degrees for each node (take a column sum)
D = diag(sum(A));
% Compute connectivity matrix:
%  (notice how b-c is the self-weight, and c is the total
%   weight from d other inputs, which each have c/d (where d might be different for each node).
%   Equal weights c/d are not required by the maths, but used for simply experiments here.)
C = (parameters.b - parameters.c) .* I + parameters.c .* A  / D;

end
