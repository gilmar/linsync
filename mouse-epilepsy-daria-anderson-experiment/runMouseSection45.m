function results = runMouseSection45(mouseId, varargin)
%RUNMOUSESECTION45  Apply the §4.5 (Liao thesis) workflow to one mouse.
%
%   For the chosen mouse:
%     1. Load the 66x66 weighted coarse connectome K.
%     2. Normalise K (default: same TVB-style 95th percentile truncation
%        + rescale to [0,1] used for the human Epileptor experiment).
%     3. Set a fully healthy 1-D Epileptor (x0_i = -2.3 for every node),
%        solve the fixed point and form the effective coupling matrix C.
%     4. Compute per-node stability susceptibility centrality D(-> i)
%        and stability influence centrality D(k ->) by reading the
%        diagonal of Omega and Omega' (covariancesGaussianNet).
%     5. For each *non-trivial* node i, sweep x0_i upward from -2.3
%        with all others held at -2.3 and find the smallest x0_i
%        at which rho(C) >= 1 (critical excitability x^c_{0,i}),
%        coarse + bisection refinement.
%     6. Render the §4.5 two-panel scatter and persist results.
%
%   Usage
%     results = runMouseSection45('Anderson_1');
%     results = runMouseSection45('Arnold_2', ...
%                  'Normalisation', 'tvb', ...
%                  'x0Upper', -1.0, ...
%                  'Plot', true, 'SaveResults', true);
%
%   Name-value options (all optional):
%     Normalisation : 'tvb' (default) | 'none' | 'parkes'
%                       'parkes' applies the nctpy.utils.matrix_normalization
%                       scheme from Parkes et al. (Nat Protoc 2024,
%                       https://doi.org/10.1038/s41596-024-01023-w):
%                           A_norm = A / (|lambda(A)|_max + ParkesC).
%                       Linsync's con2cov implicitly adds the (-I) term
%                       used by the Parkes continuous-time dynamics
%                       (dx/dt = -X (I - A_norm) + dW), so we feed
%                       A_norm straight in as the coupling matrix C.
%                       The Epileptor fixed point and per-node x^c_0
%                       sweep are skipped for this scheme -- only D(->i)
%                       and D(k->) are produced.
%     ParkesC       : 1.0    c parameter in the Parkes normalisation.
%     x0Base        : -2.3   baseline (healthy) excitability (Epileptor schemes)
%     x0Upper       : -1.0   upper bound of per-node search    (Epileptor schemes)
%     x0Step        : 0.01   coarse sweep step                 (Epileptor schemes)
%     BisectTol     : 1e-4   bisection tolerance (in x0)        (Epileptor schemes)
%     MaxK          : 1e8    covariancesGaussianNet maxIterations
%     tau0          : 6667   1-D Epileptor slow timescale       (Epileptor schemes)
%     SaveResults   : true   write results .mat / figures to results/
%     Plot          : true   produce the figure
%     Verbose       : true   per-node fprintf

setupMousePaths();

p = inputParser;
addRequired(p,  'mouseId', @(s) ischar(s) || isstring(s));
addParameter(p, 'Normalisation', 'tvb', @(s) ischar(s) || isstring(s));
addParameter(p, 'ParkesC', 1.0, @isscalar);
addParameter(p, 'x0Base',  -2.3,  @isscalar);
addParameter(p, 'x0Upper', -1.0,  @isscalar);
addParameter(p, 'x0Step',   0.01, @isscalar);
addParameter(p, 'BisectTol', 1e-4, @isscalar);
addParameter(p, 'MaxK',  100000000, @isscalar);
addParameter(p, 'tau0',  6667,   @isscalar);
addParameter(p, 'SaveResults', true, @islogical);
addParameter(p, 'Plot',        true, @islogical);
addParameter(p, 'Verbose',     true, @islogical);
parse(p, mouseId, varargin{:});
opts = p.Results;
mouseId = char(opts.mouseId);

resultsDir = setupMousePaths();

%% Load + normalise
[K_raw, labels, info] = loadMouseConnectome(mouseId);
N = info.N;

if opts.Verbose
    fprintf('\n=== %s : N = %d, symmetryError = %.3g, trivial nodes = %d ===\n', ...
        mouseId, N, info.symmetryError, sum(info.isTrivial));
end

K = applyMouseNormalisation(K_raw, opts.Normalisation, opts.ParkesC);

isParkes = strcmpi(char(opts.Normalisation), 'parkes');

%% Healthy fixed point + coupling matrix C
fopt = optimset('TolFun', 1e-14, 'TolX', 1e-14, 'Display', 'off');

if isParkes
    % Parkes' A_norm = A / (rho(A) + c) is itself the coupling matrix for
    % the linsync continuous-time OU process (con2cov implicitly applies
    % the -I term used in the Parkes continuous-time formulation). No
    % Epileptor fixed point needed.
    z_fixed_healthy = NaN(N, 1);
    C_healthy       = K;
    rho_healthy     = max(abs(eig(C_healthy)));
    if opts.Verbose
        fprintf('Parkes normalisation -- skipping Epileptor solve. rho(C) = %.6f (target < 1)\n', ...
            rho_healthy);
    end
    if rho_healthy >= 1
        error('runMouseSection45:ParkesUnstable', ...
            ['Parkes-normalised %s has rho(C) = %.4f >= 1. ' ...
             'Try a larger ParkesC (default 1.0).'], mouseId, rho_healthy);
    end
else
    x0_healthy = opts.x0Base + zeros(N, 1);
    Z0 = 3 + zeros(N, 1);
    [z_fixed_healthy, ~, exitflag_h] = fsolve( ...
        @(z) oneDepileptor(z, x0_healthy, K, opts.tau0), Z0, fopt);
    if exitflag_h <= 0
        error('runMouseSection45:HealthyFP', ...
            'Failed to solve healthy-state fixed point for %s (exitflag=%d).', ...
            mouseId, exitflag_h);
    end

    C_healthy   = CouplingMatrix(z_fixed_healthy, K, opts.tau0);
    rho_healthy = max(abs(eig(C_healthy)));
    if opts.Verbose
        fprintf('Healthy rho(C) = %.6f\n', rho_healthy);
    end
    if rho_healthy >= 1
        error('runMouseSection45:Unstable', ...
            'Healthy %s already unstable (rho(C) = %.4f >= 1)', mouseId, rho_healthy);
    end
end

%% D(-> i), D(k ->) from healthy C
D_susceptibility = NaN(N, 1);
D_influence      = NaN(N, 1);
Omega  = NaN(N);
OmegaT = NaN(N);
err_fwd = NaN; err_trans = NaN;

try
    [~, Omega, ~, err_fwd] = covariancesGaussianNet(C_healthy, false, opts.MaxK, false, 1);
    D_susceptibility = diag(Omega);
catch ME
    warning('runMouseSection45:covariancesGaussianNet', '%s', ME.message);
end
if ~isnan(err_fwd) && err_fwd ~= 0
    warning('runMouseSection45:covariancesGaussianNet', ...
        'covariancesGaussianNet(C) reported err = %d', err_fwd);
end

try
    [~, OmegaT, ~, err_trans] = covariancesGaussianNet(C_healthy', false, opts.MaxK, false, 1);
    D_influence = diag(OmegaT);
catch ME
    warning('runMouseSection45:covariancesGaussianNetTranspose', '%s', ME.message);
end
if ~isnan(err_trans) && err_trans ~= 0
    warning('runMouseSection45:covariancesGaussianNetTranspose', ...
        'covariancesGaussianNet(C'') reported err = %d', err_trans);
end

D_st_healthy = trace(Omega) / N;
if opts.Verbose
    fprintf('D_st (healthy) = %.6f (via trace(Omega)/N)\n', D_st_healthy);
end

%% Classical graph-theoretic centralities (for comparison with D(->i), D(k->))
% Computed on the same C_healthy used above so the comparison is
% apples-to-apples with the stability centralities. See Liao & Lizier
% (2026), Appendix G, and Oldham et al. (2019).
try
    classicalCentralities = computeNetworkCentralities(C_healthy);
    if opts.Verbose
        fprintf('Classical centralities computed on C_healthy.\n');
    end
catch ME
    warning('runMouseSection45:Centralities', ...
        'computeNetworkCentralities failed: %s', ME.message);
    classicalCentralities = struct();
end

%% Per-node critical excitability x^c_{0,i}
x0_crit  = NaN(N, 1);
rho_crit = NaN(N, 1);
sweepTime = 0;
if isParkes
    if opts.Verbose
        fprintf('Parkes normalisation -- skipping per-node x0^c sweep (no Epileptor model).\n');
    end
else
    tic;
    for i = 1:N
        if info.isTrivial(i)
            % Decoupled node -> sweep cannot destabilise the network here
            if opts.Verbose
                fprintf('node %3d/%d %-30s : skipped (trivial / no edges)\n', ...
                    i, N, labels{i});
            end
            continue;
        end
        [x0_crit(i), rho_crit(i)] = findCriticalX0Mouse(i, N, K, opts.tau0, ...
            opts.x0Base, opts.x0Upper, opts.x0Step, opts.BisectTol, fopt);
        if opts.Verbose
            fprintf('node %3d/%d %-30s : x0^c = %+.4f  rho(C)@crit = %.4f\n', ...
                i, N, labels{i}, x0_crit(i), rho_crit(i));
        end
    end
    sweepTime = toc;
    if opts.Verbose
        fprintf('Critical x0 sweep finished in %.1f s.\n', sweepTime);
    end
end

%% Plot
fig = [];
if opts.Plot
    if isParkes
        % D(k->)=D(->i) for symmetric K under Parkes; show only D(->i).
        fig = figure('Name', sprintf('§4.5 mouse %s (Parkes)', mouseId), ...
                     'Position', [100 100 800 520]);

        bar(D_susceptibility, 'FaceColor', [0.20 0.40 0.80]);
        xlabel('node index');
        ylabel('D(\rightarrow i)');
        xlim([0.5 N + 0.5]); grid on;

        sgtitle({sprintf('Mouse %s  --  Parkes normalisation (c = %.2f), D_{st} = %.4f', ...
                         strrep(mouseId, '_', '\_'), opts.ParkesC, D_st_healthy), ...
                 'D(\rightarrow i)  (susceptibility);  D(k\rightarrow) \equiv D(\rightarrow i) for symmetric K'}, ...
                'Interpreter', 'tex');
    else
        fig = figure('Name', sprintf('§4.5 mouse %s', mouseId), ...
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

        sgtitle(sprintf(['Mouse %s, fully healthy (x_{0,j\\neqi} = %.2f), ' ...
                         'normalisation = %s'], strrep(mouseId, '_', '\_'), ...
                         opts.x0Base, opts.Normalisation), ...
                'Interpreter', 'tex');
    end
end

%% Pack results
results = struct();
results.mouseId          = mouseId;
results.labels           = labels;
results.isTrivial        = info.isTrivial;
results.symmetryError    = info.symmetryError;
results.N                = N;
results.tau0             = opts.tau0;
results.x0_base          = opts.x0Base;
results.x0_upper         = opts.x0Upper;
results.x0_step          = opts.x0Step;
results.bisectTol        = opts.BisectTol;
results.normalisation    = char(opts.Normalisation);
results.parkesC          = opts.ParkesC;
results.isParkes         = isParkes;
results.K_raw            = K_raw;
results.K                = K;
results.z_fixed_healthy  = z_fixed_healthy;
results.C_healthy        = C_healthy;
results.rho_healthy      = rho_healthy;
results.Omega            = Omega;
results.OmegaTranspose   = OmegaT;
results.D_susceptibility = D_susceptibility;
results.D_influence      = D_influence;
results.D_st_healthy     = D_st_healthy;
results.centralities     = classicalCentralities;
results.x0_crit          = x0_crit;
results.rho_at_crit      = rho_crit;
results.err_fwd          = err_fwd;
results.err_trans        = err_trans;
results.sweepTime        = sweepTime;

%% Persist
if opts.SaveResults
    baseName = sprintf('section45_%s_%s', mouseId, opts.Normalisation);
    save(fullfile(resultsDir, [baseName '_results.mat']), 'results');
    if opts.Plot && ~isempty(fig)
        savefig(fig, fullfile(resultsDir, [baseName '_figure.fig']));
        try
            exportgraphics(fig, fullfile(resultsDir, [baseName '_figure.png']), 'Resolution', 200);
        catch
            saveas(fig, fullfile(resultsDir, [baseName '_figure.png']));
        end
    end
    if opts.Verbose
        fprintf('Saved results to %s\n', ...
            fullfile(resultsDir, [baseName '_results.mat']));
    end
end
end

%% ------------------------------------------------------------------
function K = applyMouseNormalisation(K_raw, scheme, parkesC)
%APPLYMOUSENORMALISATION  Dispatch over normalisation strategies.
%   'tvb'    replicates normal() from epileptor-experiment (95th
%            percentile truncation, zero diagonal, rescale to [0,1]).
%   'none'   returns K_raw.
%   'parkes' applies the nctpy.utils.matrix_normalization scheme of
%            Parkes et al. (Nat Protoc 2024,
%            https://doi.org/10.1038/s41596-024-01023-w):
%                A_norm = A / (|lambda(A)|_max + parkesC)
%            (with parkesC = 1 by default). The continuous-time
%            -I subtraction described in the paper is applied
%            implicitly inside linsync's con2cov, so only the rescaling
%            is performed here. The diagonal is preserved as-is.

if nargin < 3 || isempty(parkesC)
    parkesC = 1.0;
end

scheme = lower(char(scheme));
switch scheme
    case 'tvb'
        K = normal(K_raw);
    case 'none'
        K = K_raw;
    case 'parkes'
        K = matrixNormalizationParkes(K_raw, parkesC);
    otherwise
        error('applyMouseNormalisation:UnknownScheme', ...
            'Unknown normalisation scheme "%s"', scheme);
end
end

%% ------------------------------------------------------------------
function A_norm = matrixNormalizationParkes(A, c)
%MATRIXNORMALIZATIONPARKES  MATLAB port of nctpy.utils.matrix_normalization
% (Parkes et al., Nature Protocols 2024). Returns
%
%       A_norm = A / (|lambda(A)|_max + c)
%
% i.e. the discrete-time form of the Parkes normalisation. The paper's
% continuous-time variant additionally subtracts I; in this codebase
% that subtraction is applied implicitly by con2cov (which solves
% dX = -X (I - C) dt + dW), so we deliberately stop at the rescaling
% step. Preserves the rank ordering of the eigenvalues of A and the
% eigenspaces of A.

if ~isnumeric(A) || size(A, 1) ~= size(A, 2)
    error('matrixNormalizationParkes:NotSquare', ...
        'A must be a square numeric matrix.');
end
if nargin < 2 || isempty(c)
    c = 1.0;
end
rho = max(abs(eig(A)));
if ~isfinite(rho) || rho <= 0
    error('matrixNormalizationParkes:DegenerateSpectrum', ...
        'rho(A) = %.3g is not a positive finite number; cannot normalise.', rho);
end
A_norm = A ./ (rho + c);
end

%% ------------------------------------------------------------------
function [x0_c, rho_c] = findCriticalX0Mouse(i, N, K, tau0, x0_base, x0_upper, ...
                                             x0_step, bisectTol, opt)
%FINDCRITICALX0MOUSE Coarse forward sweep + bisection on x0_i to locate
% the smallest x0_i at which rho(C) >= 1 (others held at x0_base).
% Mirrors findCriticalX0 in runEZ1_section45.m.

x0_vec_base = x0_base + zeros(N, 1);
Z0 = 3 + zeros(N, 1);

last_stable     = x0_base;
last_stable_rho = NaN;
first_unstable  = NaN;
first_unstable_rho = NaN;

coarse_grid = x0_base:x0_step:x0_upper;
prev_z = Z0;
for k = 1:length(coarse_grid)
    [stable, rho_val, z_out] = evalStabilityMouse(coarse_grid(k), i, x0_vec_base, ...
        K, tau0, prev_z, opt);
    if stable
        last_stable     = coarse_grid(k);
        last_stable_rho = rho_val;
        prev_z = z_out;
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
hi_rho = first_unstable_rho;
while (hi - lo) > bisectTol
    mid = 0.5 * (lo + hi);
    [stable, rho_val, ~] = evalStabilityMouse(mid, i, x0_vec_base, K, tau0, Z0, opt);
    if stable
        lo = mid;
    else
        hi = mid;
        hi_rho = rho_val;
    end
end
x0_c  = hi;
rho_c = hi_rho;
end

function [stable, rho_val, z_out] = evalStabilityMouse(x0_i, i, x0_vec_base, K, tau0, Z0, opt)
%EVALSTABILITYMOUSE Solve fixed point with x0(i) = x0_i and test rho(C) < 1.
% Returns stability flag, spectral radius of C, and fixed point.

x0_vec = x0_vec_base;
x0_vec(i) = x0_i;

try
    [z_out, ~, exitflag] = fsolve(@(z) oneDepileptor(z, x0_vec, K, tau0), Z0, opt);
catch
    stable = false; rho_val = Inf; z_out = NaN(size(Z0)); return;
end
if exitflag <= 0 || any(~isfinite(z_out))
    stable = false; rho_val = Inf; z_out = NaN(size(Z0)); return;
end

C_test = CouplingMatrix(z_out, K, tau0);
lambdas = eig(C_test);
if any(~isfinite(lambdas))
    stable = false; rho_val = Inf; return;
end
rho_val = max(abs(lambdas));
stable  = rho_val < 1;
end
