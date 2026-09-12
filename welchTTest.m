function out = welchTTest(x, y)
%WELCHTTEST  Two-sample Welch t-test, without the Statistics Toolbox.
%
%   out = welchTTest(x, y)
%
%   Tests equality of the means of two independent samples without
%   assuming equal variances -- the right default for small cohorts, where
%   the equal-variance assumption is untestable and rarely true.
%
%   NaN entries are dropped. Fewer than two usable observations in either
%   sample yields NaN statistics rather than an error, so a cohort-wide
%   loop over nodes need not special-case sparse rows.
%
%   The two-tailed p-value uses the regularised incomplete beta function
%   (betainc), which is base MATLAB:
%       P(|T| > t) = I_{df / (df + t^2)}(df/2, 1/2)
%
%   Output struct
%     tStat  : Welch t statistic
%     df     : Welch-Satterthwaite degrees of freedom
%     p      : two-tailed p-value
%     meanX, meanY, nX, nY
%
%   See also SIGNRANKTEST, ZSCOREOUTLIERS.

x = x(:);
y = y(:);
x = x(isfinite(x));
y = y(isfinite(y));

out = struct('tStat', NaN, 'df', NaN, 'p', NaN, ...
    'meanX', NaN, 'meanY', NaN, 'nX', numel(x), 'nY', numel(y));

if numel(x) < 2 || numel(y) < 2
    if ~isempty(x), out.meanX = mean(x); end
    if ~isempty(y), out.meanY = mean(y); end
    return;
end

nx = numel(x);
ny = numel(y);
mx = mean(x);
my = mean(y);
vx = var(x);
vy = var(y);

out.meanX = mx;
out.meanY = my;

se2 = vx / nx + vy / ny;
if se2 <= 0
    return;
end

out.tStat = (mx - my) / sqrt(se2);
out.df = se2^2 / ((vx / nx)^2 / (nx - 1) + (vy / ny)^2 / (ny - 1));

if isfinite(out.tStat) && isfinite(out.df) && out.df > 0
    out.p = betainc(out.df / (out.df + out.tStat^2), out.df / 2, 0.5);
end
end
