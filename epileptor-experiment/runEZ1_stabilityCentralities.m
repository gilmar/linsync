%% runEZ1_stabilityCentralities.m
% Stability centralities and per-node critical excitability for Patient 1.
% Reproduces the analysis behind Liao (2026) thesis Fig. 4.2 (formerly runEZ1_section45):
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
%   for C' respectively.
%
%   In parallel, quantify each node's "epileptogenicity" ground truth:
%   keep all other nodes at x0 = -2.3 and increase x0_i of the target
%   node toward the epileptogenic regime until the network becomes
%   unstable. The threshold value x0_i^c (critical excitability) at
%   which the spectral radius rho(C) first reaches 1 is recorded.
%
%   The outputs are saved and plotted as a two-panel scatter of
%   D(-> i) and D(k ->) versus x0_i^c.

resultsDir = setupEpileptorPaths();

%% Basic parameters (align with runEZ1.m)
MaxK      = 100000000;
tau0      = 6667;
patient   = 'P1';
x0_base   = -2.3;
x0_upper  = -0.5;
x0_step   = 0.005;
bisectTol = 1e-4;

K = loadPatientWeights(patient);
K = normal(K);
[N, ~] = size(K);

fprintf('Patient %s, N = %d nodes.\n', patient, N);

%% Healthy-state fixed point and effective coupling matrix C
opt = optimset('TolFun', 1e-14, 'TolX', 1e-14, 'Display', 'off');
[z_fixed_healthy, C_healthy, rho_healthy] = healthyEpileptorCoupling(K, x0_base, tau0, 'runEZ1_stabilityCentralities');
fprintf('Healthy rho(C) = %.6f\n', rho_healthy);

%% Stability-based centralities from healthy C
stab = computeStabilityCentralities(C_healthy, MaxK, 'runEZ1_stabilityCentralities');
D_susceptibility = stab.D_susceptibility;
D_influence      = stab.D_influence;
Omega            = stab.Omega;
OmegaT           = stab.OmegaTranspose;
err_fwd          = stab.err_fwd;
err_trans        = stab.err_trans;
D_st_healthy     = stab.D_st_healthy;
fprintf('D_st (healthy) = %.6f (via trace(Omega)/N)\n', D_st_healthy);

%% Per-node critical excitability x0_i^c
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

%% Two-panel scatter (Fig. 4.2 style)
fig = plotStabilityScatterFigure(x0_crit, D_susceptibility, D_influence, ...
    sprintf('Figure 4.2 - Patient %s', patient), x0_base, ...
    sprintf('Stability centralities — Patient %s (healthy)', patient));

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
results.D_susceptibility = D_susceptibility;
results.D_influence      = D_influence;
results.D_st_healthy     = D_st_healthy;
results.x0_crit          = x0_crit;
results.rho_at_crit      = rho_crit;
results.err_fwd          = err_fwd;
results.err_trans        = err_trans;
results.sweepTime        = sweepTime;

baseName = sprintf('stabilityCentralities_%s_healthy', patient);
save(fullfile(resultsDir, [baseName '_results.mat']), 'results');
savefig(fig, fullfile(resultsDir, [baseName '_figure42.fig']));
try
    exportgraphics(fig, fullfile(resultsDir, [baseName '_figure42.png']), 'Resolution', 200);
catch
    saveas(fig, fullfile(resultsDir, [baseName '_figure42.png']));
end
fprintf('Saved results to %s\n', fullfile(resultsDir, [baseName '_results.mat']));
