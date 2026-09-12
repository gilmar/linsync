function A_norm = matrixNormalizationTvb(A, percentile)
%MATRIXNORMALIZATIONTVB  TVB-style truncate-and-rescale connectome normalisation.
%
%   A_norm = matrixNormalizationTvb(A)
%   A_norm = matrixNormalizationTvb(A, percentile)
%
%   Reproduces the normalisation used for structural connectomes in The
%   Virtual Brain (TVB) Epileptor workflows:
%
%     1. Take the strictly upper-triangular entries of A (i.e. one value
%        per undirected pair) and sort them.
%     2. Take the `percentile` (default 95) value of that sorted list as a
%        ceiling `uu`.
%     3. Zero the diagonal of A, clip every entry above `uu` down to `uu`.
%     4. Rescale by the maximum so that max(A_norm) == 1.
%
%   Truncating at the upper percentile stops a handful of very heavy
%   tracts from dominating the dynamics once the matrix is rescaled;
%   rescaling to [0, 1] puts every subject on a common coupling scale so
%   that a single global coupling constant is comparable across subjects.
%
%   Inputs
%     A          : NxN weighted adjacency matrix (assumed symmetric; only
%                  the upper triangle sets the ceiling).
%     percentile : ceiling percentile in (0, 100]. Default 95.
%
%   Output
%     A_norm     : NxN normalised matrix, zero diagonal, max == 1.
%
%   See also APPLYCONNECTOMENORMALISATION, MATRIXNORMALIZATIONPARKES,
%   MATRIXNORMALIZATIONCOLUMN.

if ~isnumeric(A) || size(A, 1) ~= size(A, 2)
    error('matrixNormalizationTvb:NotSquare', ...
        'A must be a square numeric matrix.');
end
if nargin < 2 || isempty(percentile)
    percentile = 95;
end
if ~isscalar(percentile) || ~isfinite(percentile) || percentile <= 0 || percentile > 100
    error('matrixNormalizationTvb:BadPercentile', ...
        'percentile must be a scalar in (0, 100].');
end

N = size(A, 1);
if N < 2
    error('matrixNormalizationTvb:TooSmall', ...
        'A must have at least 2 nodes (got %d).', N);
end

% Strictly upper-triangular entries, one per undirected pair.
u = sort(A(triu(true(N), 1)));

idx = floor(numel(u) * (percentile / 100));
idx = min(max(idx, 1), numel(u));
ceiling = u(idx);

A_norm = A - diag(diag(A));
A_norm(A_norm > ceiling) = ceiling;

peak = max(A_norm(:));
if ~isfinite(peak) || peak <= 0
    error('matrixNormalizationTvb:DegenerateWeights', ...
        'Maximum off-diagonal weight is %.3g; cannot rescale to [0, 1].', peak);
end
A_norm = A_norm / peak;
end
