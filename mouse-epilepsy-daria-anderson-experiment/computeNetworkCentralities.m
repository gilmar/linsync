function cents = computeNetworkCentralities(C, varargin)
%COMPUTENETWORKCENTRALITIES  Classical graph-theoretic centrality measures
% for a (directed, weighted) coupling matrix C, computed for comparison
% with the stability-based D(->i) and D(k->) measures from
% Liao & Lizier (2026), following the directed-network adaptations in
% their Appendix G (originally based on Oldham et al. 2019, PLoS ONE
% 14, e0220061).
%
% Convention: C(i,j) is the directed edge weight from node i (source) to
% node j (target), matching the linsync coupling-matrix convention.
%
% Returned struct fields (each is N-by-1):
%   in_strength        Weighted in-degree:  sum(C(:, i))         (in-DC)
%   out_strength       Weighted out-degree: sum(C(i, :))         (out-DC)
%   eigenvector_left   Leading left eigenvector of C (in-style)  (EC^L)
%   eigenvector_right  Leading right eigenvector of C (out-style)(EC^R)
%   pagerank           PageRank with damping AlphaPR             (PR)
%   katz               Katz centrality, attenuation AlphaKZ/rho  (KZ)
%   self_comm          diag(expm(C)) self-communicability        (SelfC)
%   betweenness        Betweenness on inverse-weight distances   (BC)
%   closeness_in       (N-1) / sum_j d(j -> i)                   (CC, in)
%   closeness_out      (N-1) / sum_j d(i -> j)                   (CC, out)
%
% Name-value options:
%   'AlphaPR'         0.85   PageRank damping factor (0,1)
%   'AlphaKZ'         0.5    Katz attenuation as fraction of 1/rho(C);
%                            must lie in (0, 1) so the Neumann series
%                            sum_k (alpha * C^T)^k converges.
%   'IgnoreSelfLoops' true   Zero out diag(C) before computing centralities
%
% Betweenness and closeness require the Brain Connectivity Toolbox
% (https://sites.google.com/site/bctnet/). If those BCT functions are not
% on the MATLAB path, the corresponding fields are filled with NaN and a
% one-time warning is issued.

p = inputParser;
addParameter(p, 'AlphaPR',         0.85, @(x) isscalar(x) && x > 0 && x < 1);
addParameter(p, 'AlphaKZ',         0.5,  @(x) isscalar(x) && x > 0 && x < 1);
addParameter(p, 'IgnoreSelfLoops', true, @islogical);
parse(p, varargin{:});
opts = p.Results;

A = full(C);
N = size(A, 1);
if size(A, 2) ~= N
    error('computeNetworkCentralities:NotSquare', ...
        'C must be a square matrix (got %dx%d).', size(A,1), size(A,2));
end
if opts.IgnoreSelfLoops
    A(1:N+1:end) = 0;
end

cents = struct();

%% Weighted in/out-degree (in-DC, out-DC)
cents.in_strength  = sum(A, 1).';
cents.out_strength = sum(A, 2);

%% Eigenvector centrality -- leading right & left eigenvectors of C
% Right (out-style):  C * v   = lambda * v   (paper Eq G6)
% Left  (in-style):   u' * C  = lambda * u'  (paper Eq G5),
%                  i.e. v_left = leading right eigenvector of C'
cents.eigenvector_right = leadingEigvec(A);
cents.eigenvector_left  = leadingEigvec(A.');

%% PageRank (paper Eq G7-G8): PR = beta*(I - alpha*A'*D^-1)^-1 * 1
cents.pagerank = pagerankCent(A, opts.AlphaPR);

%% Katz centrality (incoming): KZ = (I - alpha*A')^-1 * 1
cents.katz = katzCent(A, opts.AlphaKZ);

%% Self-communicability: diag(expm(C))
try
    cents.self_comm = real(diag(expm(A)));
catch ME
    warning('computeNetworkCentralities:Expm', ...
        'expm(C) failed (%s) -- self_comm set to NaN.', ME.message);
    cents.self_comm = NaN(N, 1);
end

%% Betweenness centrality (BCT)
if exist('betweenness_wei', 'file') == 2
    cents.betweenness = betweenness_wei(weightToLength(A));
else
    warnNoBCT('betweenness_wei');
    cents.betweenness = NaN(N, 1);
end

%% Closeness centrality (BCT distance_wei)
if exist('distance_wei', 'file') == 2
    Dmat    = distance_wei(weightToLength(A));
    Dmat(~isfinite(Dmat)) = 0;     % unreachable pairs contribute 0 to sum
    inDist  = sum(Dmat, 1).';      % sum_j d(j -> i)
    outDist = sum(Dmat, 2);        % sum_j d(i -> j)
    cents.closeness_in  = safeRatio(N - 1, inDist);
    cents.closeness_out = safeRatio(N - 1, outDist);
else
    warnNoBCT('distance_wei');
    cents.closeness_in  = NaN(N, 1);
    cents.closeness_out = NaN(N, 1);
end

end

%% ------------------------------------------------------------------
function v = leadingEigvec(M)
%LEADINGEIGVEC  Eigenvector of M for the eigenvalue with largest real
% part (Perron-Frobenius vector for non-negative matrices). Returned as
% a positive, unit-2-norm vector.
N = size(M, 1);
[V, D] = eig(M);
[~, idx] = max(real(diag(D)));
v = real(V(:, idx));
% Flip sign if dominantly negative so the vector is non-negative.
if sum(v) < 0
    v = -v;
end
v = abs(v);
nrm = norm(v);
if nrm > 0
    v = v / nrm;
else
    v = NaN(N, 1);
end
end

%% ------------------------------------------------------------------
function pr = pagerankCent(A, alpha)
%PAGERANKCENT  PageRank in {source-row, target-col} convention:
% PR = (1 - alpha) * (I - alpha * A' * D_out^-1)^-1 * 1, normalised so
% that sum(PR) = 1. Dangling nodes (out-degree 0) are treated by
% redistributing their mass uniformly via a teleportation column.
N = size(A, 1);
d_out = sum(A, 2);
dangling = (d_out == 0);
d_safe = d_out;
d_safe(dangling) = 1;
M = A.' * diag(1 ./ d_safe);
% Add uniform teleportation from dangling nodes
if any(dangling)
    M(:, dangling) = 1 / N;
end
try
    pr = (eye(N) - alpha * M) \ ((1 - alpha) / N * ones(N, 1));
catch ME
    warning('computeNetworkCentralities:PageRank', ...
        'PageRank linear solve failed (%s) -- pagerank set to NaN.', ME.message);
    pr = NaN(N, 1);
    return;
end
s = sum(pr);
if s > 0
    pr = pr / s;
end
end

%% ------------------------------------------------------------------
function k = katzCent(A, alphaFrac)
%KATZCENT  Katz centrality (incoming): KZ = (I - alpha*A')^-1 * 1, with
% alpha = alphaFrac / rho(A). Convergence requires alpha * rho(A) < 1.
N = size(A, 1);
ev = eig(A);
rho = max(abs(ev));
if ~isfinite(rho) || rho == 0
    k = NaN(N, 1);
    return;
end
alpha = alphaFrac / rho;
try
    k = (eye(N) - alpha * A.') \ ones(N, 1);
catch ME
    warning('computeNetworkCentralities:Katz', ...
        'Katz linear solve failed (%s) -- katz set to NaN.', ME.message);
    k = NaN(N, 1);
    return;
end
% Subtract the constant baseline so values reflect path-based centrality
k = k - 1;
end

%% ------------------------------------------------------------------
function L = weightToLength(A)
%WEIGHTTOLENGTH  Convert a weighted connection matrix to a connection-
% length matrix for shortest-path computations: L_ij = 1 / A_ij for
% positive weights, 0 otherwise. Self-loops are zeroed out.
L = zeros(size(A));
pos = A > 0;
L(pos) = 1 ./ A(pos);
L(1:size(A,1)+1:end) = 0;
end

%% ------------------------------------------------------------------
function r = safeRatio(num, den)
%SAFERATIO  Element-wise num ./ den, returning 0 where den == 0.
r = zeros(size(den));
ok = den > 0;
r(ok) = num ./ den(ok);
end

%% ------------------------------------------------------------------
function warnNoBCT(funcName)
%WARNNOBCT  Issue a single warning if a required BCT function is missing.
persistent warned;
if isempty(warned); warned = false; end
if warned; return; end
warning('computeNetworkCentralities:NoBCT', ...
    ['Brain Connectivity Toolbox function ''%s'' not found on the ' ...
     'MATLAB path -- the corresponding centralities will be NaN. ' ...
     'Download BCT from https://sites.google.com/site/bctnet/ and ' ...
     'addpath() its folder before running.'], funcName);
warned = true;
end
