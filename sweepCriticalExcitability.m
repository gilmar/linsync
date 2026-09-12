function out = sweepCriticalExcitability(K, varargin)
%SWEEPCRITICALEXCITABILITY  Per-node critical excitability x0^c for a connectome.
%
%   out = sweepCriticalExcitability(K, 'Tau0', 6667, 'X0Base', -2.3, ...
%                                   'X0Upper', -1.0, 'IsTrivial', isTrivial)
%
%   For each node i in turn: hold every other node at the healthy
%   excitability x0Base, raise x0_i toward the epileptogenic regime, and
%   record the smallest x0_i at which the linearised network loses
%   stability (spectral radius of the effective coupling matrix C reaches
%   1). That threshold, x0^c_i, is a model-based "how hard is it to make
%   this node trigger a seizure" score, and serves as the ground truth that
%   the stability centralities D(->i) and D(k->) are validated against.
%
%   A low x0^c means the node needs only a small excitability increase to
%   destabilise the whole network, i.e. it is highly epileptogenic.
%
%   Each node is searched with a coarse forward sweep at step X0Step,
%   followed by bisection to BisectTol between the last stable and first
%   unstable value -- cheap where the threshold is far away, precise where
%   it matters. See FINDCRITICALX0 for the per-node search itself.
%
%   Inputs
%     K : NxN normalised structural connectivity (coupling weights)
%
%   Name-value options
%     'Tau0'          : 6667   1-D Epileptor slow timescale
%     'X0Base'        : -2.3   healthy excitability held on all other nodes
%     'X0Upper'       : -1.0   upper bound of the search
%     'X0Step'        : 0.01   coarse sweep step
%     'BisectTol'     : 1e-4   bisection tolerance in x0
%     'IsTrivial'     : []     Nx1 logical; trivial nodes are skipped (NaN)
%     'Labels'        : {}     Nx1 cellstr, used only for verbose output
%     'Verbose'       : true   print one line per node
%     'FsolveOptions' : []     optimset struct for fsolve
%
%   Output struct
%     x0Crit     : Nx1 critical excitability; NaN for skipped nodes and for
%                  nodes that stay stable all the way to X0Upper
%     rhoAtCrit  : Nx1 spectral radius at the recorded threshold
%     nEvaluated : number of nodes actually searched
%     elapsedSec : wall-clock time of the sweep
%
%   Requires the Optimization Toolbox (fsolve) via findCriticalX0.
%
%   See also FINDCRITICALX0, RUNSTABILITYCENTRALITIES, HEALTHYEPILEPTORCOUPLING.

p = inputParser;
addParameter(p, 'Tau0',          6667, @isscalar);
addParameter(p, 'X0Base',        -2.3, @isscalar);
addParameter(p, 'X0Upper',       -1.0, @isscalar);
addParameter(p, 'X0Step',        0.01, @isscalar);
addParameter(p, 'BisectTol',     1e-4, @isscalar);
addParameter(p, 'IsTrivial',     [],   @(x) isempty(x) || islogical(x) || isnumeric(x));
addParameter(p, 'Labels',        {},   @(c) iscell(c) || isstring(c));
addParameter(p, 'Verbose',       true, @islogical);
addParameter(p, 'FsolveOptions', [],   @(x) isempty(x) || isstruct(x));
parse(p, varargin{:});
opts = p.Results;

N = size(K, 1);
if size(K, 2) ~= N
    error('sweepCriticalExcitability:NotSquare', ...
        'K must be square (got %dx%d).', size(K, 1), size(K, 2));
end
if opts.X0Upper <= opts.X0Base
    error('sweepCriticalExcitability:BadRange', ...
        'X0Upper (%g) must exceed X0Base (%g).', opts.X0Upper, opts.X0Base);
end

isTrivial = opts.IsTrivial;
if isempty(isTrivial)
    isTrivial = false(N, 1);
else
    isTrivial = logical(isTrivial(:));
end

labels = opts.Labels;
if isempty(labels)
    labels = arrayfun(@(i) sprintf('node %d', i), (1:N).', 'UniformOutput', false);
else
    labels = cellstr(labels);
end

fsolveOptions = opts.FsolveOptions;
if isempty(fsolveOptions)
    fsolveOptions = optimset('TolFun', 1e-14, 'TolX', 1e-14, 'Display', 'off');
end

x0Crit = NaN(N, 1);
rhoAtCrit = NaN(N, 1);
nEvaluated = 0;

t0 = tic;
for i = 1:N
    if isTrivial(i)
        if opts.Verbose
            fprintf('node %3d/%d %-32s : skipped (trivial / disconnected)\n', ...
                i, N, labels{i});
        end
        continue;
    end
    [x0Crit(i), rhoAtCrit(i)] = findCriticalX0(i, N, K, opts.Tau0, ...
        opts.X0Base, opts.X0Upper, opts.X0Step, opts.BisectTol, fsolveOptions);
    nEvaluated = nEvaluated + 1;
    if opts.Verbose
        fprintf('node %3d/%d %-32s : x0^c = %+.4f  rho(C)@crit = %.4f\n', ...
            i, N, labels{i}, x0Crit(i), rhoAtCrit(i));
    end
end
elapsedSec = toc(t0);

if opts.Verbose
    fprintf('Critical excitability sweep: %d node(s) in %.1f s.\n', nEvaluated, elapsedSec);
end

out = struct();
out.x0Crit     = x0Crit;
out.rhoAtCrit  = rhoAtCrit;
out.nEvaluated = nEvaluated;
out.elapsedSec = elapsedSec;
end
