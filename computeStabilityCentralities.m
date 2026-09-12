function out = computeStabilityCentralities(C_healthy, MaxK, warnId)
%COMPUTESTABILITYCENTRALITIES  D(->i), D(k->), and D_st from coupling C.
%
%   Computes the stationary covariance Omega via covariancesGaussianNet for
%   C and C', then returns per-node stability centralities and D_st =
%   trace(Omega)/N (see README Section 1.6).
%
%   Optional warnId prefixes warning identifiers (e.g. 'runEZ1_stabilityCentralities').

if nargin < 3 || isempty(warnId)
    warnId = 'computeStabilityCentralities';
end

N = size(C_healthy, 1);

out = struct();
out.D_susceptibility = NaN(N, 1);
out.D_influence      = NaN(N, 1);
out.Omega            = NaN(N);
out.OmegaTranspose   = NaN(N);
out.err_fwd          = NaN;
out.err_trans        = NaN;
out.D_st_healthy     = NaN;

try
    [~, out.Omega, ~, out.err_fwd] = covariancesGaussianNet(C_healthy, false, MaxK, false, 1);
    out.D_susceptibility = diag(out.Omega);
catch ME
    warning('%s:covariancesGaussianNet', warnId, '%s', ME.message);
end
if ~isnan(out.err_fwd) && out.err_fwd ~= 0
    warning('%s:covariancesGaussianNet', warnId, ...
        'covariancesGaussianNet(C) reported err = %d', out.err_fwd);
end

try
    [~, out.OmegaTranspose, ~, out.err_trans] = covariancesGaussianNet(C_healthy', false, MaxK, false, 1);
    out.D_influence = diag(out.OmegaTranspose);
catch ME
    warning('%s:covariancesGaussianNetTranspose', warnId, '%s', ME.message);
end
if ~isnan(out.err_trans) && out.err_trans ~= 0
    warning('%s:covariancesGaussianNetTranspose', warnId, ...
        'covariancesGaussianNet(C'') reported err = %d', out.err_trans);
end

out.D_st_healthy = trace(out.Omega) / N;
end
