function out = zScoreOutliers(values, refMean, refStd, varargin)
%ZSCOREOUTLIERS  Per-node z-scores against a reference cohort, Bonferroni-corrected.
%
%   out = zScoreOutliers(values, refMean, refStd)
%   out = zScoreOutliers(values, refMean, refStd, 'Alpha', 0.05, 'Include', mask)
%
%   Standardises one subject's per-node values against the mean and SD of a
%   reference (control) cohort, converts them to two-tailed p-values under a
%   normal null, and flags the entries that survive Bonferroni correction.
%
%   Pass every metric you intend to test as a separate column so that the
%   correction family covers all of them at once -- the family size is the
%   number of eligible (finite, included) entries across the whole array,
%   which is what makes the correction honest when several metrics are
%   inspected per node.
%
%   Inputs
%     values  : NxM subject values (M metrics side by side)
%     refMean : NxM reference-cohort means (a Nx1 column is expanded)
%     refStd  : NxM reference-cohort standard deviations
%
%   Name-value options
%     'Alpha'   : 0.05 -- family-wise error rate
%     'Include' : logical mask, NxM or Nx1, marking entries eligible to be
%                 tested (e.g. ~isTrivial). Default: all true. Excluded
%                 entries get z / p values but never count toward the
%                 family size and are never flagged.
%
%   Output struct
%     z          : NxM z-scores; NaN where refStd is 0 or inputs are NaN
%     p          : NxM two-tailed p-values, uncorrected
%     pBonf      : NxM p-values multiplied by nTests and capped at 1
%     isOutlier  : NxM logical, pBonf < Alpha and entry included
%     nTests     : size of the correction family (>= 1)
%     zThreshold : |z| equivalent to the Bonferroni threshold, for display
%     alpha      : the Alpha used
%
%   Uses erfc / erfinv rather than normcdf / norminv so no Statistics and
%   Machine Learning Toolbox licence is needed.
%
%   See also COMPARECOHORTGROUPS, COMPAREHEMISPHERICASYMMETRY.

p = inputParser;
addParameter(p, 'Alpha', 0.05, @(x) isscalar(x) && x > 0 && x < 1);
addParameter(p, 'Include', [], @(x) isempty(x) || islogical(x) || isnumeric(x));
parse(p, varargin{:});
alpha = p.Results.Alpha;

values = double(values);
[N, M] = size(values);
refMean = expandToSize(double(refMean), N, M, 'refMean');
refStd  = expandToSize(double(refStd),  N, M, 'refStd');

include = p.Results.Include;
if isempty(include)
    include = true(N, M);
else
    include = expandToSize(logical(include), N, M, 'Include');
end

z = NaN(N, M);
ok = isfinite(values) & isfinite(refMean) & isfinite(refStd) & refStd > 0;
z(ok) = (values(ok) - refMean(ok)) ./ refStd(ok);

% Two-tailed normal p-value: 2 * (1 - Phi(|z|)) == erfc(|z| / sqrt(2)).
pVal = NaN(N, M);
pVal(isfinite(z)) = erfc(abs(z(isfinite(z))) / sqrt(2));

eligible = include & isfinite(z);
nTests = max(sum(eligible(:)), 1);

pBonf = min(1, pVal * nTests);
isOutlier = eligible & isfinite(pBonf) & pBonf < alpha;

% |z| at which an eligible test would clear Bonferroni, i.e.
% norminv(1 - alpha / (2 * nTests)) written with erfinv.
zThreshold = sqrt(2) * erfinv(1 - alpha / nTests);

out = struct( ...
    'z',          z, ...
    'p',          pVal, ...
    'pBonf',      pBonf, ...
    'isOutlier',  isOutlier, ...
    'nTests',     nTests, ...
    'zThreshold', zThreshold, ...
    'alpha',      alpha);
end

%% ------------------------------------------------------------------
function A = expandToSize(A, N, M, name)
if isvector(A) && numel(A) == N && M > 1 && size(A, 2) == 1
    A = repmat(A, 1, M);
end
if ~isequal(size(A), [N M])
    error('zScoreOutliers:SizeMismatch', ...
        '%s must be %dx%d (or %dx1); got %dx%d.', name, N, M, N, size(A, 1), size(A, 2));
end
end
