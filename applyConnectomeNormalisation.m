function K = applyConnectomeNormalisation(K_raw, scheme, parkesC, colScale, tvbPercentile)
%APPLYCONNECTOMENORMALISATION  Dispatch over connectome normalisation schemes.
%
%   K = applyConnectomeNormalisation(K_raw, scheme)
%   K = applyConnectomeNormalisation(K_raw, scheme, parkesC, colScale, tvbPercentile)
%
%   A raw structural connectome carries arbitrary units (streamline counts,
%   fibre densities, ...), so it must be rescaled before it can be read as
%   the coupling matrix of a linear dynamical system. The scheme decides
%   both the scale and whether the result stays symmetric.
%
%   Schemes
%     'tvb'    -- matrixNormalizationTvb: truncate at the `tvbPercentile`
%                 percentile of the upper triangle, zero diagonal, rescale
%                 to [0, 1]. Symmetric. Intended to be fed through an
%                 Epileptor fixed point (see healthyEpileptorCoupling).
%     'none'   -- return K_raw unchanged (also for an Epileptor solve).
%     'parkes' -- matrixNormalizationParkes: A / (|lambda|_max + c), after
%                 Parkes et al., Nature Protocols 2024. Symmetric, so the
%                 resulting D(k->) is identical to D(->i).
%     'column' -- matrixNormalizationColumn: rescale each column to sum to
%                 `colScale`. Asymmetric in general, so D(k->) and D(->i)
%                 carry different information.
%
%   Defaults: parkesC = 1.0, colScale = 0.95, tvbPercentile = 95.
%
%   See also MATRIXNORMALIZATIONTVB, MATRIXNORMALIZATIONPARKES,
%   MATRIXNORMALIZATIONCOLUMN, NORMALISATIONSCHEMEINFO.

if nargin < 3 || isempty(parkesC)
    parkesC = 1.0;
end
if nargin < 4 || isempty(colScale)
    colScale = 0.95;
end
if nargin < 5 || isempty(tvbPercentile)
    tvbPercentile = 95;
end

scheme = normalisationSchemeInfo(scheme).scheme;
switch scheme
    case 'tvb'
        K = matrixNormalizationTvb(K_raw, tvbPercentile);
    case 'none'
        K = K_raw;
    case 'parkes'
        K = matrixNormalizationParkes(K_raw, parkesC);
    case 'column'
        K = matrixNormalizationColumn(K_raw, colScale);
    otherwise
        error('applyConnectomeNormalisation:UnknownScheme', ...
            'Unknown normalisation scheme "%s"', scheme);
end
end
