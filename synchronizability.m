function exp_sigma_sqr = synchronizability(UcovarianceXU)
% 
% Compute the expectation value of \sigma^2, meaning the expected average
% square deviation of nodes X from the network mean.
%
% Inputs
% - UcovarianceXU - covariance matrix, projected into space orthogonal to the
%    synchronized state psi_0, computed from the coupling matrix C
%    (or several C's over multiple delays).
%    If using this to compute exp_sigma_sqr for dynamics with delays, then
%    only pass in the NxN matrix UcovarianceU corresponding to
%    the nodes of X (not the matrix for the embedded system Z).
%
% Outputs
% - exp_sigma_sqr - expectation of the average square distance of nodes
%    from the network mean, a.k.a the width of the sync landscape
%
%% Linear Sync Toolkit (linsync)
% Copyright (C) 2023 Joseph T. Lizier
% Distributed under GNU General Public License v3

N = size(UcovarianceXU,1);

exp_sigma_sqr = trace(UcovarianceXU) ./ N;

end

