function v = extractMetricVector(results, fieldPath, N)
%EXTRACTMETRICVECTOR  Pull one node-level metric out of a results struct.
%
%   v = extractMetricVector(results, 'D_susceptibility', N)
%   v = extractMetricVector(results, 'centralities.betweenness', N)
%
%   Walks the dotted field path into the per-subject results struct and
%   returns an Nx1 column. Anything missing -- an absent field, an absent
%   nested struct, a vector of the wrong length -- returns all-NaN rather
%   than erroring, so a cohort mixing older result files (before a metric
%   existed) with newer ones still aggregates: NaN-aware means simply drop
%   the subjects that lack the metric.
%
%   See also CONNECTOMEMETRICCATALOG, COMPARECOHORTGROUPS.

v = NaN(N, 1);
if ~isstruct(results)
    return;
end

parts = strsplit(char(fieldPath), '.');
node = results;
for k = 1:numel(parts)
    if ~isstruct(node) || ~isfield(node, parts{k})
        return;
    end
    node = node.(parts{k});
end

if ~isnumeric(node) && ~islogical(node)
    return;
end
node = double(node(:));
if numel(node) ~= N
    return;
end
v = node;
end
