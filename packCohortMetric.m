function M = packCohortMetric(cohort, fieldPath, mask)
%PACKCOHORTMETRIC  Gather one metric across a cohort into an N x nSubjects matrix.
%
%   M = packCohortMetric(cohort, 'D_susceptibility')
%   M = packCohortMetric(cohort, 'centralities.betweenness', cohort.isControl)
%
%   Columns follow the order of cohort.subjects (restricted to `mask` when
%   given). Subjects missing the metric contribute an all-NaN column, which
%   the NaN-aware statistics downstream simply skip.
%
%   See also LOADCOHORTRESULTS, EXTRACTMETRICVECTOR, CONNECTOMEMETRICCATALOG.

if nargin < 3 || isempty(mask)
    mask = true(numel(cohort.subjects), 1);
end
mask = logical(mask(:));

idx = find(mask);
M = NaN(cohort.N, numel(idx));
for k = 1:numel(idx)
    M(:, k) = extractMetricVector(cohort.results{idx(k)}, fieldPath, cohort.N);
end
end
