%% runEZ1_section45.m
% Reproduce Chapter 4, Section 4.5 ("Measuring Susceptible Nodes in
% Epileptic Networks") of Jieru Liao's thesis. In particular, this
% script reproduces Figure 4.2:
%
%   For Patient 1's structural connectivity K, set a *fully healthy*
%   Epileptor configuration (x0_i = -2.3 for every node i), compute the
%   linearised effective coupling matrix C, and from it the per-node
%   stability-based centralities:
%
%       D(-> i) = Omega_ii         (stability susceptibility centrality)
%       D(k ->) = (Omega_{C'})_kk  (stability influence centrality)
%
%   where Omega is the stationary covariance of the linear OU process
%   with system matrix C (D_st = trace(Omega)/N, see README
%   Section 1.6). Omega is obtained by covariancesGaussianNet for C and
%   for C' respectively, matching the convention in JL_playing.m.
%
%   In parallel, quantify each node's "epileptogenicity" ground truth:
%   keep all other nodes at x0 = -2.3 and increase x0_i of the target
%   node toward the epileptogenic regime until the network becomes
%   unstable. The threshold value x0_i^c (critical excitability) at
%   which the spectral radius rho(C) first reaches 1 is recorded. This
%   is the stability criterion from Chapter 3 Section II.B / Section
%   4.3 (|lambda_v| < 1, equivalently rho(C) < 1).
%
%   The three outputs are saved and used to reproduce Figure 4.2 as a
%   two-panel scatter plot of D(-> i) and D(k ->) versus x0_i^c.

resultsDir = setupEpileptorPaths();

%% Basic parameters (align with runEZ1.m / Section 4.4 of the thesis)
MaxK      = 100000000;
tau0      = 6667;
patient   = 'P1';
x0_base   = -2.3;          % healthy baseline excitability for every node
x0_upper  = -0.5;          % upper bound of the per-node critical search
x0_step   = 0.005;         % coarse sweep step
bisectTol = 1e-4;          % bisection tolerance on x0

K = loadPatientWeights(patient);
K = normal(K);             % 95th-percentile-truncated, [0,1] rescaled
[N, ~] = size(K);

fprintf('Patient %s, N = %d nodes.\n', patient, N);

%% Healthy-state fixed point and effective coupling matrix C
x0_healthy = x0_base + zeros(N, 1);
Z0 = 3 + zeros(N, 1);
opt = optimset('TolFun', 1e-14, 'TolX', 1e-14, 'Display', 'off');

[z_fixed_healthy, ~, exitflag_h] = fsolve( ...
    @(z) oneDepileptor(z, x0_healthy, K, tau0), Z0, opt);
if exitflag_h <= 0
    error('runEZ1_section45: failed to solve the healthy-state fixed point (exitflag=%d).', exitflag_h);
end

C_healthy = CouplingMatrix(z_fixed_healthy, K, tau0);
rho_healthy = max(abs(eig(C_healthy)));
fprintf('Healthy rho(C) = %.6f\n', rho_healthy);
if rho_healthy >= 1
    error(['runEZ1_section45: healthy configuration is already unstable ' ...
           '(rho(C) = %.4f >= 1); cannot compute stability centralities.'], ...
           rho_healthy);
end

%% Stability-based centralities from healthy C
%   D(-> i) = Omega_ii for dynamics with coupling C
%   D(k ->) = Omega_kk for dynamics with coupling C' (transpose)
D_susceptibility = NaN(N, 1);
D_influence      = NaN(N, 1);
Omega  = NaN(N);
OmegaT = NaN(N);
err_fwd = NaN;
err_trans = NaN;

try
    [~, Omega, ~, err_fwd] = covariancesGaussianNet(C_healthy, false, MaxK, false, 1);
    D_susceptibility = diag(Omega);
catch ME
    warning('runEZ1_section45:covariancesGaussianNet', '%s', ME.message);
end
if (~isnan(err_fwd) && err_fwd ~= 0)
    warning('runEZ1_section45:covariancesGaussianNet', ...
        'covariancesGaussianNet(C) reported err = %d', err_fwd);
end

try
    [~, OmegaT, ~, err_trans] = covariancesGaussianNet(C_healthy', false, MaxK, false, 1);
    D_influence = diag(OmegaT);
catch ME
    warning('runEZ1_section45:covariancesGaussianNetTranspose', '%s', ME.message);
end
if (~isnan(err_trans) && err_trans ~= 0)
    warning('runEZ1_section45:covariancesGaussianNetTranspose', ...
        'covariancesGaussianNet(C'') reported err = %d', err_trans);
end

D_st_healthy = trace(Omega) / N;
fprintf('D_st (healthy) = %.6f (via trace(Omega)/N)\n', D_st_healthy);

%% Per-node critical excitability x0_i^c
%   For each node i: keep x0_j = -2.3 for j ~= i, and sweep x0_i upward
%   from -2.3 until the network is unstable (rho(C) >= 1). Refine with
%   bisection between the last stable and the first unstable sample.
x0_crit = NaN(N, 1);
rho_crit = NaN(N, 1);
tic
for i = 1:N
    [x0_crit(i), rho_crit(i)] = findCriticalX0(i, N, K, tau0, x0_base, ...
        x0_upper, x0_step, bisectTol, opt);
    fprintf('node %3d/%d: x0^c = %+.4f  rho(C)@crit = %.4f\n', ...
        i, N, x0_crit(i), rho_crit(i));
end
sweepTime = toc;
fprintf('Critical x0 sweep finished in %.1f s.\n', sweepTime);

%% Reproduce Figure 4.2 (two-panel scatter)
fig = figure('Name', sprintf('Figure 4.2 (%s, healthy)', patient), ...
             'Position', [100 100 1200 480]);

subplot(1, 2, 1);
scatter(x0_crit, D_susceptibility, 48, 'filled');
xlabel('critical excitability  x^c_{0,i}');
ylabel('D(\rightarrow i)');
title('Stability susceptibility centrality');
grid on;

subplot(1, 2, 2);
scatter(x0_crit, D_influence, 48, 'filled', 'MarkerFaceColor', [0.85 0.33 0.10]);
xlabel('critical excitability  x^c_{0,i}');
ylabel('D(k \rightarrow)');
title('Stability influence centrality');
grid on;

sgtitle(sprintf(['Figure 4.2 - Patient %s, fully healthy network ' ...
                 '(x_{0,j\\neqi} = %.2f)'], patient, x0_base));

%% Persist results
results.patient          = patient;
results.N                = N;
results.tau0             = tau0;
results.x0_base          = x0_base;
results.x0_upper         = x0_upper;
results.x0_step          = x0_step;
results.bisectTol        = bisectTol;
results.K                = K;
results.z_fixed_healthy  = z_fixed_healthy;
results.C_healthy        = C_healthy;
results.Omega            = Omega;
results.OmegaTranspose   = OmegaT;
results.D_susceptibility = D_susceptibility;   % D(-> i)
results.D_influence      = D_influence;        % D(k ->)
results.D_st_healthy     = D_st_healthy;
results.x0_crit          = x0_crit;
results.rho_at_crit      = rho_crit;
results.err_fwd          = err_fwd;
results.err_trans        = err_trans;
results.sweepTime        = sweepTime;

baseName = sprintf('section45_%s_healthy', patient);
save(fullfile(resultsDir, [baseName '_results.mat']), 'results');
savefig(fig, fullfile(resultsDir, [baseName '_figure42.fig']));
try
    exportgraphics(fig, fullfile(resultsDir, [baseName '_figure42.png']), 'Resolution', 200);
catch
    saveas(fig, fullfile(resultsDir, [baseName '_figure42.png']));
end
fprintf('Saved results to %s\n', fullfile(resultsDir, [baseName '_results.mat']));

%% ------------------------------------------------------------------
function [x0_c, rho_c] = findCriticalX0(i, N, K, tau0, x0_base, x0_upper, ...
                                        x0_step, bisectTol, opt)
%FINDCRITICALX0 Find the smallest x0_i for node i at which the linearised
% effective coupling matrix of the 1-D Epileptor network becomes unstable
% (rho(C) >= 1), with all other nodes held at x0_j = x0_base.
%
% Strategy: coarse forward sweep with step x0_step until either
% (a) rho(C) >= 1 is detected, or (b) fsolve fails (treated as seizure
% onset). Then bisect between the last stable x0_i and the first
% unstable x0_i to the specified tolerance bisectTol.

x0_vec_base = x0_base + zeros(N, 1);
Z0 = 3 + zeros(N, 1);

last_stable     = x0_base;
last_stable_rho = NaN;
last_stable_z   = Z0;
first_unstable  = NaN;
first_unstable_rho = NaN;

coarse_grid = x0_base:x0_step:x0_upper;
prev_z = Z0;
fsolve_failure_streak = 0;
for k = 1:length(coarse_grid)
    [stable, rho_val, z_out] = evalStability(coarse_grid(k), i, x0_vec_base, ...
        K, tau0, prev_z, opt);
    if stable
        last_stable     = coarse_grid(k);
        last_stable_rho = rho_val;
        last_stable_z   = z_out;
        prev_z          = z_out;
        fsolve_failure_streak = 0;
    elseif ~isfinite(rho_val) && fsolve_failure_streak < 1
        % Tolerate a single isolated fsolve failure: skip without declaring
        % instability so a transient numerical glitch does not truncate the sweep.
        fsolve_failure_streak = fsolve_failure_streak + 1;
        continue;
    else
        first_unstable     = coarse_grid(k);
        first_unstable_rho = rho_val;
        break;
    end
end

if isnan(first_unstable)
    x0_c = NaN;
    rho_c = last_stable_rho;
    return;
end

lo = last_stable;
hi = first_unstable;
lo_rho = last_stable_rho;
hi_rho = first_unstable_rho;
warm_z = last_stable_z;
while (hi - lo) > bisectTol
    mid = 0.5 * (lo + hi);
    [stable, rho_val, z_out] = evalStability(mid, i, x0_vec_base, K, tau0, warm_z, opt);
    if stable
        lo = mid;
        lo_rho = rho_val;
        warm_z = z_out;
    else
        hi = mid;
        hi_rho = rho_val;
    end
end

x0_c = hi;
rho_c = hi_rho;
end

function [stable, rho_val, z_out] = evalStability(x0_i, i, x0_vec_base, K, tau0, Z0, opt)
%EVALSTABILITY Solve fixed point with x0(i) = x0_i (others held at baseline)
% and test rho(C) < 1. Returns stability flag, spectral radius of C, and
% the solved fixed point (NaN vector if fsolve failed).

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
