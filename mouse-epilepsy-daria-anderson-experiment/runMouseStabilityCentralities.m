function results = runMouseStabilityCentralities(mouseId, varargin)
%RUNMOUSESTABILITYCENTRALITIES  Stability centralities and critical x0 for one mouse.
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
%     6. Render the two-panel scatter and persist results.
%
%   Usage
%     results = runMouseStabilityCentralities('Anderson_1');
%     results = runMouseStabilityCentralities('Arnold_2', ...
%                  'Normalisation', 'tvb', ...
%                  'x0Upper', -1.0, ...
%                  'Plot', true, 'SaveResults', true);
%
%   Name-value options (all optional):
%     Normalisation : 'tvb' (default) | 'none' | 'parkes' | 'column'
%     ParkesC       : 1.0    c parameter in the Parkes normalisation.
%     ColScale      : 0.95   target column sum for the 'column' scheme.
%     x0Base        : -2.3   baseline (healthy) excitability (Epileptor schemes)
%     x0Upper       : -1.0   upper bound of per-node search    (Epileptor schemes)
%     x0Step        : 0.01   coarse sweep step                 (Epileptor schemes)
%     BisectTol     : 1e-4   bisection tolerance (in x0)        (Epileptor schemes)
%     MaxK          : 1e8    covariancesGaussianNet maxIterations
%     tau0          : 6667   1-D Epileptor slow timescale       (Epileptor schemes)
%     SaveResults   : true   write results .mat / figures to results/
%     Plot          : true   produce the figure
%     Verbose       : true   per-node fprintf
%     ResultsDir    : ''     output folder (default: results/ via setupMousePaths)

setupMousePaths();

p = inputParser;
addRequired(p,  'mouseId', @(s) ischar(s) || isstring(s));
addParameter(p, 'Normalisation', 'tvb', @(s) ischar(s) || isstring(s));
addParameter(p, 'ParkesC', 1.0,  @isscalar);
addParameter(p, 'ColScale', 0.95, @isscalar);
addParameter(p, 'x0Base',  -2.3,  @isscalar);
addParameter(p, 'x0Upper', -1.0,  @isscalar);
addParameter(p, 'x0Step',   0.01, @isscalar);
addParameter(p, 'BisectTol', 1e-4, @isscalar);
addParameter(p, 'MaxK',  100000000, @isscalar);
addParameter(p, 'tau0',  6667,   @isscalar);
addParameter(p, 'SaveResults', true, @islogical);
addParameter(p, 'Plot',        true, @islogical);
addParameter(p, 'Verbose',     true, @islogical);
addParameter(p, 'ResultsDir',  '', @(s) ischar(s) || isstring(s));
addParameter(p, 'RunParameters', [], @(x) isempty(x) || isstruct(x));
parse(p, mouseId, varargin{:});
opts = p.Results;
mouseId = char(opts.mouseId);

resultsDir = resolveMouseResultsDir(opts.ResultsDir);
stabilityDir = mouseResultsDir(resultsDir, 'stability');

%% Load + normalise
[K_raw, labels, info] = loadMouseConnectome(mouseId);
N = info.N;

if opts.Verbose
    fprintf('\n=== %s : N = %d, symmetryError = %.3g, trivial nodes = %d ===\n', ...
        mouseId, N, info.symmetryError, sum(info.isTrivial));
end

K = applyConnectomeNormalisation(K_raw, opts.Normalisation, opts.ParkesC, opts.ColScale);

schemeLower    = lower(char(opts.Normalisation));
isParkes       = strcmp(schemeLower, 'parkes');
isColumn       = ismember(schemeLower, {'column', 'col', 'colnorm'});
isLinearScheme = isParkes || isColumn;

%% Healthy fixed point + coupling matrix C
fopt = optimset('TolFun', 1e-14, 'TolX', 1e-14, 'Display', 'off');

if isLinearScheme
    z_fixed_healthy = NaN(N, 1);
    C_healthy       = K;
    rho_healthy     = max(abs(eig(C_healthy)));
    if opts.Verbose
        fprintf('%s normalisation -- skipping Epileptor solve. rho(C) = %.6f (target < 1)\n', ...
            schemeLower, rho_healthy);
    end
    if rho_healthy >= 1
        if isParkes
            hint = 'Try a larger ParkesC (default 1.0).';
        else
            hint = 'Try a smaller ColScale (default 0.95).';
        end
        error('runMouseStabilityCentralities:LinearSchemeUnstable', ...
            '%s-normalised %s has rho(C) = %.4f >= 1. %s', ...
            schemeLower, mouseId, rho_healthy, hint);
    end
else
    [z_fixed_healthy, C_healthy, rho_healthy] = healthyEpileptorCoupling( ...
        K, opts.x0Base, opts.tau0, 'runMouseStabilityCentralities');
    if opts.Verbose
        fprintf('Healthy rho(C) = %.6f\n', rho_healthy);
    end
end

%% D(-> i), D(k ->) from healthy C
stab = computeStabilityCentralities(C_healthy, opts.MaxK, 'runMouseStabilityCentralities');
D_susceptibility = stab.D_susceptibility;
D_influence      = stab.D_influence;
Omega            = stab.Omega;
OmegaT           = stab.OmegaTranspose;
err_fwd          = stab.err_fwd;
err_trans        = stab.err_trans;
D_st_healthy     = stab.D_st_healthy;
if opts.Verbose
    fprintf('D_st (healthy) = %.6f (via trace(Omega)/N)\n', D_st_healthy);
end

%% Classical graph-theoretic centralities (for comparison with D(->i), D(k->))
try
    classicalCentralities = computeNetworkCentralities(C_healthy);
    if opts.Verbose
        fprintf('Classical centralities computed on C_healthy.\n');
    end
catch ME
    warning('runMouseStabilityCentralities:Centralities', ...
        'computeNetworkCentralities failed: %s', ME.message);
    classicalCentralities = struct();
end

%% Per-node critical excitability x^c_{0,i}
x0_crit  = NaN(N, 1);
rho_crit = NaN(N, 1);
sweepTime = 0;
if isLinearScheme
    if opts.Verbose
        fprintf('%s normalisation -- skipping per-node x0^c sweep (no Epileptor model).\n', ...
            schemeLower);
    end
else
    tic;
    for i = 1:N
        if info.isTrivial(i)
            if opts.Verbose
                fprintf('node %3d/%d %-30s : skipped (trivial / no edges)\n', ...
                    i, N, labels{i});
            end
            continue;
        end
        [x0_crit(i), rho_crit(i)] = findCriticalX0(i, N, K, opts.tau0, ...
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
        fig = figure('Name', sprintf('Stability centralities — mouse %s (Parkes)', mouseId), ...
                     'Position', [100 100 800 520]);

        bar(D_susceptibility, 'FaceColor', [0.20 0.40 0.80]);
        xlabel('node index');
        ylabel('D(\rightarrow i)');
        xlim([0.5 N + 0.5]); grid on;

        sgtitle({sprintf('Mouse %s  --  Parkes normalisation (c = %.2f), D_{st} = %.4f', ...
                         strrep(mouseId, '_', '\_'), opts.ParkesC, D_st_healthy), ...
                 'D(\rightarrow i)  (susceptibility);  D(k\rightarrow) \equiv D(\rightarrow i) for symmetric K'}, ...
                'Interpreter', 'tex');
    elseif isColumn
        fig = figure('Name', sprintf('Stability centralities — mouse %s (column)', mouseId), ...
                     'Position', [100 100 1300 520]);

        subplot(1, 2, 1);
        bar(D_susceptibility, 'FaceColor', [0.20 0.40 0.80]);
        xlabel('node index'); ylabel('D(\rightarrow i)');
        xlim([0.5 N + 0.5]); grid on;
        title('Stability susceptibility centrality');

        subplot(1, 2, 2);
        bar(D_influence, 'FaceColor', [0.85 0.33 0.10]);
        xlabel('node index'); ylabel('D(k \rightarrow)');
        xlim([0.5 N + 0.5]); grid on;
        title('Stability influence centrality');

        sgtitle(sprintf(['Mouse %s  --  column normalisation (scale = %.2f), ' ...
                         'D_{st} = %.4f'], strrep(mouseId, '_', '\_'), ...
                         opts.ColScale, D_st_healthy), ...
                'Interpreter', 'tex');
    else
        subjectLabel = sprintf('Mouse %s, normalisation = %s', ...
            strrep(mouseId, '_', '\_'), opts.Normalisation);
        fig = plotStabilityScatterFigure(x0_crit, D_susceptibility, D_influence, ...
            subjectLabel, opts.x0Base, ...
            sprintf('Stability centralities — mouse %s', mouseId));
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
results.colScale         = opts.ColScale;
results.isParkes         = isParkes;
results.isColumn         = isColumn;
results.isLinearScheme   = isLinearScheme;
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

if isempty(opts.RunParameters)
    results.runParameters = mouseExperimentRunParameters('buildFromStabilityCentralities', opts, mouseId);
else
    results.runParameters = mouseExperimentRunParameters('merge', opts.RunParameters, ...
        struct('stabilityCentralities', struct('mouseId', mouseId)));
end

%% Persist
if opts.SaveResults
    prefix = mouseExperimentResultPrefix();
    baseName = sprintf('%s_%s_%s', prefix, mouseId, opts.Normalisation);
    save(fullfile(stabilityDir, [baseName '_results.mat']), 'results');
    if opts.Plot && ~isempty(fig)
        savefig(fig, fullfile(stabilityDir, [baseName '_figure.fig']));
        try
            exportgraphics(fig, fullfile(stabilityDir, [baseName '_figure.png']), 'Resolution', 200);
        catch
            saveas(fig, fullfile(stabilityDir, [baseName '_figure.png']));
        end
        close(fig);
    end
    if opts.Verbose
        fprintf('Saved results to %s\n', ...
            fullfile(stabilityDir, [baseName '_results.mat']));
    end
end
end
