function A_norm = matrixNormalizationColumn(A, scale)
%MATRIXNORMALIZATIONCOLUMN  Weighted in-degree (column) normalisation.
%
%   Each column is rescaled to sum to `scale` (default 0.95). Trivial nodes
%   (column sum = 0) stay all-zero.

if ~isnumeric(A) || size(A, 1) ~= size(A, 2)
    error('matrixNormalizationColumn:NotSquare', ...
        'A must be a square numeric matrix.');
end
if nargin < 2 || isempty(scale)
    scale = 0.95;
end
if ~isscalar(scale) || ~isfinite(scale) || scale <= 0
    error('matrixNormalizationColumn:BadScale', ...
        'scale must be a positive finite scalar.');
end
colSum = sum(A, 1);
A_norm = zeros(size(A), 'like', A);
nz = colSum > 0;
A_norm(:, nz) = A(:, nz) .* (scale ./ colSum(nz));
end
