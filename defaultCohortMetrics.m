function fields = defaultCohortMetrics(cohort, callerId)
%DEFAULTCOHORTMETRICS  The node metrics a cohort comparison compares by default.
%
%   fields = defaultCohortMetrics(cohort)
%
%   Returns the dotted field paths that every cohort-level comparison uses
%   unless the caller names its own 'Metrics':
%
%     D_susceptibility            always
%     D_influence                 only when the coupling matrix is
%                                 asymmetric -- for a symmetric K it is
%                                 identical to D(->i) by construction, so
%                                 including it would double the
%                                 multiple-comparison family for no
%                                 additional information
%     centralities.betweenness    when the Brain Connectivity Toolbox was
%                                 available during the run
%
%   Betweenness is the classical centrality worth carrying alongside the
%   stability measures: on normalised structural connectomes the others
%   (PageRank, Katz, eigenvector, degree, self-communicability) correlate
%   almost perfectly with D(->i), whereas betweenness reflects path
%   structure and stays partly independent.
%
%   See also CONNECTOMEMETRICCATALOG, COMPARECOHORTGROUPS,
%   COMPAREHEMISPHERICASYMMETRY.

if nargin < 2 || isempty(callerId)
    callerId = 'defaultCohortMetrics';
end

fields = {'D_susceptibility'};

schemeInfo = normalisationSchemeInfo(cohort.scheme);
if ~schemeInfo.isSymmetric
    fields{end+1} = 'D_influence';
end

bc = packCohortMetric(cohort, 'centralities.betweenness');
if any(isfinite(bc(:)))
    fields{end+1} = 'centralities.betweenness';
else
    warning('%s:NoBetweenness', callerId, ...
        ['Betweenness centrality is unavailable in these results (the Brain ' ...
         'Connectivity Toolbox was not on the MATLAB path when the stability ' ...
         'run executed), so it is omitted from the comparison.']);
end
end
