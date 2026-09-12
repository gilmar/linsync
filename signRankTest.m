function out = signRankTest(v, mu, varargin)
%SIGNRANKTEST  Wilcoxon signed-rank test against a constant, normal approximation.
%
%   out = signRankTest(v)        % tests median(v) == 0
%   out = signRankTest(v, mu)    % tests median(v) == mu
%   out = signRankTest(v, mu, 'ContinuityCorrection', true)
%
%   Answers "is this sample systematically shifted away from mu?" without
%   assuming normality -- the appropriate test for pooled laterality
%   indices, which are bounded and typically not normal.
%
%   Implemented with the tie-corrected normal approximation using only base
%   MATLAB, so no Statistics Toolbox licence is needed. The approximation is
%   the standard choice above roughly 20 observations; below that the exact
%   test would be preferable.
%
%   Zero differences are dropped (Wilcoxon's convention) and NaNs ignored.
%
%   'ContinuityCorrection' (default false) shrinks the test statistic by 0.5
%   toward its null mean before standardising. It is left off by default so
%   the p-value matches the Statistics Toolbox `signrank(v, mu, 'method',
%   'approximate')` exactly; turn it on for the slightly more conservative
%   value that R's wilcox.test reports by default.
%
%   Output struct
%     p        : two-tailed p-value (NaN when fewer than 2 non-zero values)
%     z        : normal-approximation z statistic
%     W        : sum of ranks of the positive differences
%     n        : number of non-zero differences used
%     median   : median of the input (NaN-ignoring)
%
%   See also WELCHTTEST, COMPAREHEMISPHERICASYMMETRY.

if nargin < 2 || isempty(mu)
    mu = 0;
end

p = inputParser;
addParameter(p, 'ContinuityCorrection', false, @islogical);
parse(p, varargin{:});
useContinuity = p.Results.ContinuityCorrection;

v = v(:);
v = v(isfinite(v));
out = struct('p', NaN, 'z', NaN, 'W', NaN, 'n', 0, 'median', NaN);
if isempty(v)
    return;
end
out.median = median(v);

d = v - mu;
d = d(d ~= 0);
n = numel(d);
out.n = n;
if n < 2
    return;
end

r = tiedRanks(abs(d));
W = sum(r(d > 0));
out.W = W;

muW = n * (n + 1) / 4;
varW = n * (n + 1) * (2 * n + 1) / 24;

% Tie correction: each group of t equal |d| values removes (t^3 - t)/48.
[sorted, ~] = sort(abs(d));
runStart = 1;
for k = 2:n + 1
    if k > n || sorted(k) ~= sorted(runStart)
        t = k - runStart;
        if t > 1
            varW = varW - (t^3 - t) / 48;
        end
        runStart = k;
    end
end

if varW <= 0
    return;
end

deviation = W - muW;
if useContinuity
    deviation = sign(deviation) * max(abs(deviation) - 0.5, 0);
end
out.z = deviation / sqrt(varW);
out.p = erfc(abs(out.z) / sqrt(2));
end

%% ------------------------------------------------------------------
function r = tiedRanks(x)
%TIEDRANKS  Ranks of x with ties replaced by their average rank.
n = numel(x);
[sorted, order] = sort(x(:));
r = zeros(n, 1);
k = 1;
while k <= n
    j = k;
    while j < n && sorted(j + 1) == sorted(k)
        j = j + 1;
    end
    r(order(k:j)) = (k + j) / 2;
    k = j + 1;
end
end
