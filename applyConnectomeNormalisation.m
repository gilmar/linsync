function K = applyConnectomeNormalisation(K_raw, scheme, parkesC, colScale)
%APPLYCONNECTOMENORMALISATION  Dispatch over connectome normalisation schemes.
%
%   'tvb'    -- normal(): 95th percentile truncation, zero diagonal, [0,1]
%   'none'   -- return K_raw unchanged
%   'parkes' -- matrixNormalizationParkes (Parkes et al. Nat Protoc 2024)
%   'column' -- matrixNormalizationColumn (Liao & Lizier 2026 Appendix H.1)

if nargin < 3 || isempty(parkesC)
    parkesC = 1.0;
end
if nargin < 4 || isempty(colScale)
    colScale = 0.95;
end

scheme = lower(char(scheme));
switch scheme
    case 'tvb'
        K = normal(K_raw);
    case 'none'
        K = K_raw;
    case 'parkes'
        K = matrixNormalizationParkes(K_raw, parkesC);
    case {'column', 'col', 'colnorm'}
        K = matrixNormalizationColumn(K_raw, colScale);
    otherwise
        error('applyConnectomeNormalisation:UnknownScheme', ...
            'Unknown normalisation scheme "%s"', scheme);
end
end
