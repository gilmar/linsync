function LI = lateralityIndex(values, leftIdx, rightIdx, basis)
%LATERALITYINDEX  Left-vs-right asymmetry per homotopic node pair.
%
%   LI = lateralityIndex(values, leftIdx, rightIdx)            % 'norm'
%   LI = lateralityIndex(values, leftIdx, rightIdx, 'signed')
%
%   Inputs
%     values   : N x S matrix of a node-level metric (N nodes, S subjects)
%     leftIdx  : nPairs x 1 row indices of the left members  (pairHemisphereNodes)
%     rightIdx : nPairs x 1 row indices of the right members
%     basis    : 'norm' (default) or 'signed'
%
%   Bases
%     'norm'   LI = (L - R) / (L + R), bounded in [-1, +1]. Comparable
%              across metrics and subjects because it divides out the
%              overall magnitude, which makes it the right basis for
%              group-level tests. Undefined (NaN) where L + R <= 0.
%     'signed' LI = L - R, in the metric's own units. Keeps the effect
%              size interpretable ("how much more susceptible is the left
%              node"), which is what you want when reporting a specific
%              subject's outlying pairs.
%
%   Positive values mean the left-labelled node carries more of the metric.
%
%   See also PAIRHEMISPHERENODES, COMPAREHEMISPHERICASYMMETRY.

if nargin < 4 || isempty(basis)
    basis = 'norm';
end
basis = lower(strtrim(char(basis)));

leftIdx = leftIdx(:);
rightIdx = rightIdx(:);
if numel(leftIdx) ~= numel(rightIdx)
    error('lateralityIndex:PairMismatch', ...
        'leftIdx and rightIdx must have the same length (%d vs %d).', ...
        numel(leftIdx), numel(rightIdx));
end

L = values(leftIdx, :);
R = values(rightIdx, :);
LI = NaN(size(L));

switch basis
    case 'norm'
        denom = L + R;
        ok = isfinite(L) & isfinite(R) & denom > 0;
        LI(ok) = (L(ok) - R(ok)) ./ denom(ok);
    case 'signed'
        ok = isfinite(L) & isfinite(R);
        LI(ok) = L(ok) - R(ok);
    otherwise
        error('lateralityIndex:BadBasis', ...
            'basis must be ''norm'' or ''signed'' (got ''%s'').', basis);
end
end
