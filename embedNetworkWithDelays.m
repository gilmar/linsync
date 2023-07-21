function B = embedNetworkWithDelays(C)
%% function B = embedNetworkWithDelays(C)
%
% Converts a network connectivity matrix, which potentially has delayed
% connectivity, into an embedded network which has the standard delay only
% but with dummy nodes.
%
% Inputs:
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
%
% Outputs
% - B - the matrix is just the C matrix is the case of no-delays (C is
%     NxN), else if C is NxNx(tau+1), then each NxN matrix is embedded into
%     B which is the connectivity matrix for the embedded system Z.
%
%% Linear Sync Toolkit (linsync)
% Copyright (C) 2023 Joseph T. Lizier
% Distributed under GNU General Public License v3

N = size(C,1);
hasDelays = (length(size(C)) == 3);
if (hasDelays)
    tauPlus1 = size(C,3); % Number of delays we consider. If > 1 we're looking beyond standard case
else
    tauPlus1 = 1; % Only the standard delay here
end

% Construct the embedded matrix B (if no delays B = C)
B = zeros(N * tauPlus1);
if hasDelays
    for n = 1 : tauPlus1
        % Insert the delayed coupling matrix C_{n-1} into B
        B(1+(n-1)*N:n*N,1:N) = C(:,:,n);
    end
    % Insert I along the upper block diagonal
    B(1:(tauPlus1-1)*N, N+1:end) = eye((tauPlus1-1) * N);
    % zLength = size(B,1);
    % IZ = eye(zLength);
    % IX = zeros(zLength); % not required later
    % IX(1:N,1:N) = I; % not required later
    % GX = [repmat(G, 1, tauPlus1); zeros(N*(tauPlus1-1), N*tauPlus1)];
    % UX = IZ - GX; % not required later
else
    B = C; % standard case B collapses to C
end

