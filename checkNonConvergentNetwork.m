function checkNonConvergentNetwork(network)
%% function checkNonConvergentNetwork(network)
%
% Quick check of why this network was not convergent
%
% Inputs:
% - network (optional) - either:
%    - a matrix being the connectivity matrix to check, OR
%    - a filename to load the network from, OR
%    - (default) if no provided the network is loaded from the local file
%    nonconvergentNetwork.mat
%
%% Linear Sync Toolkit (linsync)
% Copyright (C) 2023 Joseph T. Lizier
% Distributed under GNU General Public License v3

% Load the network into matrix C
if (nargin > 0)
    if (~isstring(network))
        C = network;
    else
        % network is a filename
        load(network); % Assume this loads a network into matrix C
    end
else
    load('nonconvergentNetwork.mat');
end

% Postcondition: network is loaded into variable C
if (~exist('C', 'var'))
    error('Matrix C was not found in the given file (nonconvergentNetwork.mat if none provided)');
end

N = size(C,1);
% Construct the embedded matrix B (if no delays B = C)
hasDelays = (length(size(C)) == 3);
B = embedNetworkWithDelays(C);
% For system X:
% I = eye(N); % not required below
G = ones(N) / N;
% U = I - G; % not required below
if hasDelays
    tauPlus1 = size(C,3); % Number of delays we consider. If > 1 we're looking beyond standard case
    zLength = size(B,1);
    IZ = eye(zLength);
    % IX = zeros(zLength); % not required later
    % IX(1:N,1:N) = I; % not required later
    GX = [repmat(G, 1, tauPlus1); zeros(N*(tauPlus1-1), N*tauPlus1)];
    UX = IZ - GX;
%else
%    tauPlus1 = 1; % Only the standard delay here
end

fprintf('Network B is diagonalizable? %s (the following power iteration can only give direct sole eigenvector interpretation if so)\n', mat2str(isdiagonalizable(B)));

% Run power iteration method for B*U to see where it lands (if not zeros then we
% have an eigenvalue with |lambda| > 1
v = randn(1, N * tauPlus1);
maxIterations = 10000;
v = v * (B * UX)^maxIterations;

figure(); imagesc(B); colorbar; title('Embedded connectivity matrix B');
figure(); imagesc(B * UX); colorbar; title('Projected embedded connectivity matrix B');
figure(); imagesc((B * UX)^maxIterations); colorbar; title('Projected embedded connectivity matrix B, many powers');
figure(); plot(v, 'rx'); title('Elements of iterated vector v');
