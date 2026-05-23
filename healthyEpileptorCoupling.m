function [z_fixed, C_healthy, rho_healthy] = healthyEpileptorCoupling(K, x0Base, tau0, errorId)
%HEALTHYEPILEPTORCOUPLING  Healthy-state fixed point and coupling matrix C.
%
%   Solve the 1-D Epileptor fixed point with x0_i = x0Base for all nodes,
%   form the effective coupling matrix C, and error if rho(C) >= 1.
%
%   Optional errorId prefixes error identifiers (e.g. 'runEZ1_stabilityCentralities').

if nargin < 4 || isempty(errorId)
    errorId = 'healthyEpileptorCoupling';
end

N = size(K, 1);
x0_healthy = x0Base + zeros(N, 1);
Z0 = 3 + zeros(N, 1);
opt = optimset('TolFun', 1e-14, 'TolX', 1e-14, 'Display', 'off');

[z_fixed, ~, exitflag] = fsolve( ...
    @(z) oneDepileptor(z, x0_healthy, K, tau0), Z0, opt);
if exitflag <= 0
    error('%s:HealthyFP', ...
        'Failed to solve the healthy-state fixed point (exitflag=%d).', exitflag);
end

C_healthy = CouplingMatrix(z_fixed, K, tau0);
rho_healthy = max(abs(eig(C_healthy)));
if rho_healthy >= 1
    error('%s:Unstable', ...
        'Healthy configuration is already unstable (rho(C) = %.4f >= 1).', rho_healthy);
end
end
