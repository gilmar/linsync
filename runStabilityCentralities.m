function results = runStabilityCentralities(subjectId, varargin)
%RUNSTABILITYCENTRALITIES  Stability centralities for one subject's connectome.
%
%   results = runStabilityCentralities('Case_1', 'DataRoot', dataRoot)
%   results = runStabilityCentralities('Case_1', 'K', K, 'Labels', labels)
%   results = runStabilityCentralities(id, 'Normalisation', 'tvb', 'Plot', false)
%
%   This is the per-subject workhorse of the connectome pipeline. Given one
%   weighted structural connectome it:
%
%     1. Loads (or accepts) the raw NxN matrix K_raw and its node labels.
%     2. Normalises K_raw for the chosen scheme (see
%        APPLYCONNECTOMENORMALISATION) -- a raw connectome has arbitrary
%        units and must be put on a scale where it can act as the coupling
%        matrix of a stable linear system.
%     3. Forms the effective coupling matrix C:
%          * Epileptor schemes ('tvb', 'none') solve the healthy 1-D
%            Epileptor fixed point and linearise around it, so C is the
%            Jacobian of a neural-mass model driven by the connectome.
%          * Linear schemes ('parkes', 'column') use the normalised K
%            directly as C, with no neural model and no fixed point.
%     4. Computes the stationary covariance Omega of the linear
%        fluctuations around that state, and reads the stability
%        centralities off its diagonal:
%          D(-> i) = Omega_ii            how strongly node i absorbs
%                                        network-wide fluctuations
%          D(k ->) = (Omega_{C'})_kk     how strongly node k projects them
%        plus the network-level D_st = trace(Omega) / N.
%     5. Computes the classical graph centralities for comparison.
%     6. For Epileptor schemes, sweeps the per-node critical excitability
%        x0^c as a model-based epileptogenicity ground truth.
%     7. Plots and saves, embedding the full parameter record.
%
%   Inputs
%     subjectId : subject ID; also the folder name under 'DataRoot'
%
%   Data source (pick one)
%     'DataRoot'       : folder of per-subject subfolders (default: pwd/data)
%     'ConnectomeFile' : filename inside the subject folder ('connectome.csv')
%     'K'              : supply the raw matrix directly and skip loading
%     'Labels'         : node labels to go with 'K'
%     'IsTrivial'      : trivial-node mask to go with 'K'
%     'TrivialLabelSuffixes' : label suffixes marking placeholder rows
%
%   Model options
%     'Normalisation' : 'column' (default) | 'parkes' | 'tvb' | 'none'
%     'ParkesC'       : 1.0   c in A / (|lambda|_max + c)
%     'ColScale'      : 0.95  target column sum for column normalisation
%     'TvbPercentile' : 95    truncation percentile for the tvb scheme
%     'Tau0'          : 6667  1-D Epileptor slow timescale
%     'X0Base'        : -2.3  healthy excitability
%     'X0Upper'       : -1.0  upper bound of the x0^c search
%     'X0Step'        : 0.01  coarse sweep step
%     'BisectTol'     : 1e-4  bisection tolerance in x0
%     'MaxK'          : 1e8   covariancesGaussianNet iteration cap
%     'SweepCriticalX0' : true -- set false to skip the (slow) sweep
%
%   Output options
%     'SaveResults'   : true  write the .mat and figures
%     'Plot'          : true  build the figure
%     'Verbose'       : true  per-node progress
%     'ResultsDir'    : ''    run folder (default: pwd/results)
%     'RunParameters' : []    run-level provenance to embed
%
%   Returns the packed `results` struct, also saved as
%   <prefix>_<subjectId>_<scheme>_results.mat under the run's
%   stability_centralities/ folder.
%
%   Epileptor schemes need the Optimization Toolbox (fsolve).
%
%   See also RUNCOHORTSTABILITYCENTRALITIES, COMPUTESTABILITYCENTRALITIES,
%   SWEEPCRITICALEXCITABILITY, APPLYCONNECTOMENORMALISATION.

p = inputParser;
addRequired(p,  'subjectId', @(s) ischar(s) || isstring(s));
addParameter(p, 'DataRoot',       '', @(s) ischar(s) || isstring(s));
addParameter(p, 'ConnectomeFile', 'connectome.csv', @(s) ischar(s) || isstring(s));
addParameter(p, 'K',              [], @(x) isempty(x) || isnumeric(x));
addParameter(p, 'Labels',         {}, @(c) iscell(c) || isstring(c));
addParameter(p, 'IsTrivial',      [], @(x) isempty(x) || islogical(x) || isnumeric(x));
addParameter(p, 'TrivialLabelSuffixes', {}, @(c) iscell(c) || isstring(c) || ischar(c));
addParameter(p, 'Normalisation',  'column', @(s) ischar(s) || isstring(s));
addParameter(p, 'ParkesC',        1.0,  @isscalar);
addParameter(p, 'ColScale',       0.95, @isscalar);
addParameter(p, 'TvbPercentile',  95,   @isscalar);
addParameter(p, 'Tau0',           6667, @isscalar);
addParameter(p, 'X0Base',         -2.3, @isscalar);
addParameter(p, 'X0Upper',        -1.0, @isscalar);
addParameter(p, 'X0Step',         0.01, @isscalar);
addParameter(p, 'BisectTol',      1e-4, @isscalar);
addParameter(p, 'MaxK',           1e8,  @isscalar);
addParameter(p, 'SweepCriticalX0', true, @islogical);
addParameter(p, 'SaveResults',    true, @islogical);
addParameter(p, 'Plot',           true, @islogical);
addParameter(p, 'Verbose',        true, @islogical);
addParameter(p, 'ResultsDir',     '',   @(s) ischar(s) || isstring(s));
addParameter(p, 'RunParameters',  [],   @(x) isempty(x) || isstruct(x));
parse(p, subjectId, varargin{:});
opts = p.Results;
subjectId = char(opts.subjectId);

schemeInfo = normalisationSchemeInfo(opts.Normalisation);
scheme = schemeInfo.scheme;

resultsDir = resolveExperimentResultsDir(opts.ResultsDir);
stabilityDir = experimentResultsDir(resultsDir, 'stability');

%% ---- Load the connectome ------------------------------------------
if isempty(opts.K)
    dataRoot = char(opts.DataRoot);
    if isempty(dataRoot)
        dataRoot = fullfile(pwd, 'data');
    end
    [K_raw, labels, info] = loadConnectome(subjectId, ...
        'DataRoot', dataRoot, ...
        'FileName', opts.ConnectomeFile, ...
        'TrivialLabelSuffixes', opts.TrivialLabelSuffixes);
else
    K_raw = double(opts.K);
    info = struct();
    info.subjectId = subjectId;
    info.sourceFile = '';
    info.N = size(K_raw, 1);
    labels = opts.Labels;
    if isempty(labels)
        labels = arrayfun(@(i) sprintf('node_%03d', i), (1:info.N).', ...
            'UniformOutput', false);
    else
        labels = cellstr(labels);
        labels = labels(:);
    end
    info.labels = labels;
    if isempty(opts.IsTrivial)
        Koff = K_raw - diag(diag(K_raw));
        info.isTrivial = (sum(abs(Koff), 1).' + sum(abs(Koff), 2)) == 0;
    else
        info.isTrivial = logical(opts.IsTrivial(:));
    end
    info.symmetryError = max(abs(K_raw - K_raw.'), [], 'all');
    info.isSymmetric = info.symmetryError == 0;
end
N = info.N;

if opts.Verbose
    fprintf('\n=== %s : N = %d, symmetry error = %.3g, trivial nodes = %d ===\n', ...
        subjectId, N, info.symmetryError, sum(info.isTrivial));
end

%% ---- Normalise -----------------------------------------------------
K = applyConnectomeNormalisation(K_raw, scheme, ...
    opts.ParkesC, opts.ColScale, opts.TvbPercentile);

%% ---- Effective coupling matrix C ----------------------------------
if schemeInfo.isLinear
    z_fixed_healthy = NaN(N, 1);
    C_healthy = K;
    rho_healthy = max(abs(eig(C_healthy)));
    if opts.Verbose
        fprintf(['%s normalisation is linear -- using K directly as C ' ...
                 '(no Epileptor solve). rho(C) = %.6f (needs < 1)\n'], ...
            scheme, rho_healthy);
    end
    if rho_healthy >= 1
        error('runStabilityCentralities:UnstableCoupling', ...
            ['%s-normalised connectome for %s has rho(C) = %.4f >= 1, so the ' ...
             'stationary covariance does not exist. %s'], ...
            scheme, subjectId, rho_healthy, stabilityHint(scheme));
    end
else
    [z_fixed_healthy, C_healthy, rho_healthy] = healthyEpileptorCoupling( ...
        K, opts.X0Base, opts.Tau0, 'runStabilityCentralities');
    if opts.Verbose
        fprintf('Healthy fixed point solved; rho(C) = %.6f\n', rho_healthy);
    end
end

%% ---- Stability centralities ---------------------------------------
stab = computeStabilityCentralities(C_healthy, opts.MaxK, 'runStabilityCentralities');
if opts.Verbose
    fprintf('D_st (healthy) = %.6f  [trace(Omega)/N]\n', stab.D_st_healthy);
end

%% ---- Classical centralities for comparison ------------------------
try
    classicalCentralities = computeNetworkCentralities(C_healthy);
catch ME
    warning('runStabilityCentralities:Centralities', ...
        'computeNetworkCentralities failed: %s', ME.message);
    classicalCentralities = struct();
end

%% ---- Per-node critical excitability -------------------------------
x0Crit = NaN(N, 1);
rhoAtCrit = NaN(N, 1);
sweepTime = 0;
if schemeInfo.usesEpileptor && opts.SweepCriticalX0
    sweep = sweepCriticalExcitability(K, ...
        'Tau0',      opts.Tau0, ...
        'X0Base',    opts.X0Base, ...
        'X0Upper',   opts.X0Upper, ...
        'X0Step',    opts.X0Step, ...
        'BisectTol', opts.BisectTol, ...
        'IsTrivial', info.isTrivial, ...
        'Labels',    labels, ...
        'Verbose',   opts.Verbose);
    x0Crit = sweep.x0Crit;
    rhoAtCrit = sweep.rhoAtCrit;
    sweepTime = sweep.elapsedSec;
elseif opts.Verbose && schemeInfo.isLinear
    fprintf(['%s normalisation has no Epileptor model -- skipping the ' ...
             'critical-excitability sweep.\n'], scheme);
end

%% ---- Pack ----------------------------------------------------------
results = struct();
results.subjectId        = subjectId;
results.labels           = labels;
results.isTrivial        = info.isTrivial;
results.symmetryError    = info.symmetryError;
results.sourceFile       = info.sourceFile;
results.N                = N;
results.normalisation    = scheme;
results.schemeInfo       = schemeInfo;
results.parkesC          = opts.ParkesC;
results.colScale         = opts.ColScale;
results.tvbPercentile    = opts.TvbPercentile;
results.tau0             = opts.Tau0;
results.x0_base          = opts.X0Base;
results.x0_upper         = opts.X0Upper;
results.x0_step          = opts.X0Step;
results.bisectTol        = opts.BisectTol;
results.K_raw            = K_raw;
results.K                = K;
results.z_fixed_healthy  = z_fixed_healthy;
results.C_healthy        = C_healthy;
results.rho_healthy      = rho_healthy;
results.Omega            = stab.Omega;
results.OmegaTranspose   = stab.OmegaTranspose;
results.D_susceptibility = stab.D_susceptibility;
results.D_influence      = stab.D_influence;
results.D_st_healthy     = stab.D_st_healthy;
results.centralities     = classicalCentralities;
results.x0_crit          = x0Crit;
results.rho_at_crit      = rhoAtCrit;
results.err_fwd          = stab.err_fwd;
results.err_trans        = stab.err_trans;
results.sweepTime        = sweepTime;

if isempty(opts.RunParameters)
    results.runParameters = experimentRunParameters('fromOptions', opts, ...
        'Analysis', 'runStabilityCentralities', ...
        'ResultsDir', resultsDir, 'Scheme', scheme, 'Subjects', {subjectId});
else
    results.runParameters = experimentRunParameters('merge', opts.RunParameters, ...
        struct('stabilityCentralities', struct('subjectId', subjectId)));
end

%% ---- Plot and persist ---------------------------------------------
fig = [];
if opts.Plot
    fig = plotStabilityCentralitiesFigure(results);
end

if opts.SaveResults
    baseName = sprintf('%s_%s_%s', experimentResultPrefix(), subjectId, scheme);
    matPath = fullfile(stabilityDir, [baseName '_results.mat']);
    save(matPath, 'results');
    if ~isempty(fig) && isgraphics(fig, 'figure')
        saveFigureBoth(fig, fullfile(stabilityDir, [baseName '_figure']));
        close(fig);
    end
    if opts.Verbose
        fprintf('Saved %s\n', matPath);
    end
end
end

%% ------------------------------------------------------------------
function hint = stabilityHint(scheme)
switch scheme
    case 'parkes'
        hint = 'Increase ParkesC (default 1.0) to push the spectrum further inside the unit circle.';
    case 'column'
        hint = 'Decrease ColScale (default 0.95) to push the spectrum further inside the unit circle.';
    otherwise
        hint = 'Choose a normalisation that rescales the spectrum below 1.';
end
end
