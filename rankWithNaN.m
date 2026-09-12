function r = rankWithNaN(values, direction)
%RANKWITHNAN  Rank a vector, sending NaN entries to the bottom.
%
%   r = rankWithNaN(values)              % 'descend' (rank 1 = largest)
%   r = rankWithNaN(values, 'ascend')    % rank 1 = smallest
%
%   Finite entries are ranked 1..nFinite in the requested direction; NaN
%   entries all receive rank numel(values) + 1 so they sort last and never
%   pull a mean rank toward the top.
%
%   Used to aggregate per-subject rankings across a cohort where some
%   nodes are missing a value in some subjects.

if nargin < 2 || isempty(direction)
    direction = 'descend';
end

v = values(:);
N = numel(v);
r = NaN(N, 1);
finite = ~isnan(v);
nFinite = sum(finite);

[~, ord] = sort(v(finite), direction);
ranks = NaN(nFinite, 1);
ranks(ord) = 1:nFinite;
r(finite) = ranks;
r(~finite) = N + 1;
end
