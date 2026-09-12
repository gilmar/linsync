function [stable, rho_val, z_out] = evalEpileptorStability(x0_i, i, x0_vec_base, K, tau0, Z0, opt)
%EVALEPILEPTORSTABILITY  Test stability for a trial x0_i on one node.
%
%   Solve the 1-D Epileptor fixed point with x0(i) = x0_i (others held at
%   baseline) and test rho(C) < 1. Returns stability flag, spectral radius
%   of C, and the solved fixed point (NaN vector if fsolve failed).

x0_vec = x0_vec_base;
x0_vec(i) = x0_i;

try
    [z_out, ~, exitflag] = fsolve(@(z) oneDepileptor(z, x0_vec, K, tau0), Z0, opt);
catch
    stable = false;
    rho_val = Inf;
    z_out = NaN(size(Z0));
    return;
end

if exitflag <= 0 || any(~isfinite(z_out))
    stable = false;
    rho_val = Inf;
    z_out = NaN(size(Z0));
    return;
end

C_test = CouplingMatrix(z_out, K, tau0);
lambdas = eig(C_test);
if any(~isfinite(lambdas))
    stable = false;
    rho_val = Inf;
    return;
end
rho_val = max(abs(lambdas));
stable = rho_val < 1;
end
