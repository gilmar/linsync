function A_norm = matrixNormalizationParkes(A, c)
%MATRIXNORMALIZATIONPARKES  MATLAB port of nctpy.utils.matrix_normalization
% (Parkes et al., Nature Protocols 2024). Returns
%
%       A_norm = A / (|lambda(A)|_max + c)
%
% The continuous-time -I subtraction is applied implicitly by con2cov.

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
