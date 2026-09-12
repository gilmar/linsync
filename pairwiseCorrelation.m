function C = pairwiseCorrelation(X, corrType)
%PAIRWISECORRELATION  Correlation matrix of the columns of X, NaN-tolerant.
%
%   C = pairwiseCorrelation(X)                % Pearson
%   C = pairwiseCorrelation(X, 'Spearman')
%
%   Computes the M x M correlation matrix of the columns of an N x M
%   matrix, using only the rows where both columns of a pair are finite
%   (pairwise deletion). Columns that are constant or have fewer than two
%   usable rows give NaN rather than an error, so a measure that is
%   unavailable for one subject does not sink the whole matrix.
%
%   Uses the Statistics Toolbox `corr` when it is available, and otherwise
%   falls back to an equivalent base-MATLAB implementation: Spearman and
%   Kendall-free rank handling is done by ranking each column (ties
%   averaged) before a Pearson correlation, which is the definition of
%   Spearman's rho.
%
%   corrType : 'Pearson' (default) | 'Spearman'
%
%   See also COMPARECENTRALITYMEASURES.

if nargin < 2 || isempty(corrType)
    corrType = 'Pearson';
end
corrType = char(corrType);

X = double(X);
M = size(X, 2);
C = NaN(M, M);

usable = false(1, M);
for j = 1:M
    usable(j) = sum(isfinite(X(:, j))) >= 2;
end
if ~any(usable)
    return;
end

if exist('corr', 'file') == 2
    try
        sub = X(:, usable);
        C(usable, usable) = corr(sub, 'Type', corrType, 'Rows', 'pairwise');
        return;
    catch ME
        warning('pairwiseCorrelation:ToolboxCorr', ...
            'corr() failed (%s); falling back to the base-MATLAB implementation.', ...
            ME.message);
    end
end

switch lower(corrType)
    case 'pearson'
        R = X;
    case 'spearman'
        R = NaN(size(X));
        for j = 1:M
            R(:, j) = rankColumn(X(:, j));
        end
    otherwise
        error('pairwiseCorrelation:UnsupportedType', ...
            ['Without the Statistics Toolbox only Pearson and Spearman are ' ...
             'supported (got "%s").'], corrType);
end

idx = find(usable);
for a = 1:numel(idx)
    for b = a:numel(idx)
        ja = idx(a);
        jb = idx(b);
        ok = isfinite(R(:, ja)) & isfinite(R(:, jb));
        if sum(ok) < 2
            continue;
        end
        u = R(ok, ja);
        v = R(ok, jb);
        if strcmpi(corrType, 'spearman')
            % Re-rank within the retained rows so ties and dropped rows are
            % handled exactly as pairwise Spearman requires.
            u = rankColumn(u);
            v = rankColumn(v);
        end
        r = pearson(u, v);
        C(ja, jb) = r;
        C(jb, ja) = r;
    end
end
end

%% ------------------------------------------------------------------
function r = pearson(u, v)
u = u - mean(u);
v = v - mean(v);
denom = sqrt(sum(u.^2) * sum(v.^2));
if denom <= 0
    r = NaN;
else
    r = sum(u .* v) / denom;
end
end

%% ------------------------------------------------------------------
function r = rankColumn(x)
%RANKCOLUMN  Average ranks of the finite entries; NaN stays NaN.
r = NaN(size(x));
ok = isfinite(x);
v = x(ok);
n = numel(v);
if n == 0
    return;
end
[sorted, order] = sort(v);
ranks = zeros(n, 1);
k = 1;
while k <= n
    j = k;
    while j < n && sorted(j + 1) == sorted(k)
        j = j + 1;
    end
    ranks(order(k:j)) = (k + j) / 2;
    k = j + 1;
end
r(ok) = ranks;
end
