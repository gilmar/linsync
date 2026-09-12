function catalog = connectomeMetricCatalog(fieldNames)
%CONNECTOMEMETRICCATALOG  Names and labels for the node-level metrics of a run.
%
%   catalog = connectomeMetricCatalog()
%   catalog = connectomeMetricCatalog({'D_susceptibility', 'centralities.betweenness'})
%
%   Every per-node quantity a stability run produces is addressed by the
%   dotted field path where it lives inside the results struct. This
%   catalog maps those paths to the short key, TeX label and filename
%   suffix used in figures, CSV columns and output filenames, so that
%   adding a metric to a comparison is a one-line change rather than a
%   hunt through plotting code.
%
%   With no argument, returns the full catalog. With a cellstr of field
%   paths, returns those entries in that order; an unknown path still
%   produces an entry, with labels derived from the path itself, so a
%   user-supplied metric works without editing this file.
%
%   Entry fields
%     field       : dotted path into the results struct (e.g. 'centralities.katz')
%     key         : short identifier, valid as a struct field name
%     label       : TeX label for axes and titles
%     fileSuffix  : filename-safe suffix for per-metric outputs
%     description : one-line explanation
%
%   See also EXTRACTMETRICVECTOR, COMPARECOHORTGROUPS, COMPARECENTRALITYMEASURES.

known = { ...
    'D_susceptibility',           'D_susc',   'D(\rightarrow i)',  'D_to_i',    'stability susceptibility centrality: Omega_ii'; ...
    'D_influence',                'D_infl',   'D(k \rightarrow)',  'D_k_to',    'stability influence centrality: diagonal of Omega for C'''; ...
    'x0_crit',                    'x0_crit',  'x^c_{0,i}',         'x0_crit',   'critical excitability (lower = more epileptogenic)'; ...
    'centralities.betweenness',   'BC',       'BC',                'BC',        'betweenness centrality on inverse-weight distances'; ...
    'centralities.in_strength',   'in_DC',    'in-DC',             'in_DC',     'weighted in-degree'; ...
    'centralities.out_strength',  'out_DC',   'out-DC',            'out_DC',    'weighted out-degree'; ...
    'centralities.eigenvector_left',  'EC_L', 'EC^{L}',            'EC_L',      'leading left eigenvector centrality'; ...
    'centralities.eigenvector_right', 'EC_R', 'EC^{R}',            'EC_R',      'leading right eigenvector centrality'; ...
    'centralities.pagerank',      'PR',       'PR',                'PR',        'PageRank'; ...
    'centralities.katz',          'KZ',       'KZ',                'KZ',        'Katz centrality'; ...
    'centralities.self_comm',     'SelfC',    'SelfC',             'SelfC',     'self-communicability, diag(expm(C))'; ...
    'centralities.closeness_in',  'CC_in',    'CC_{in}',           'CC_in',     'closeness over incoming shortest paths'; ...
    'centralities.closeness_out', 'CC_out',   'CC_{out}',          'CC_out',    'closeness over outgoing shortest paths'};

full = repmat(emptyEntry(), size(known, 1), 1);
for k = 1:size(known, 1)
    full(k) = makeEntry(known{k, 1}, known{k, 2}, known{k, 3}, known{k, 4}, known{k, 5});
end

if nargin < 1 || isempty(fieldNames)
    catalog = full;
    return;
end

fieldNames = cellstr(fieldNames);
catalog = repmat(emptyEntry(), numel(fieldNames), 1);
knownPaths = {full.field};
for k = 1:numel(fieldNames)
    idx = find(strcmp(knownPaths, fieldNames{k}), 1);
    if isempty(idx)
        catalog(k) = derivedEntry(fieldNames{k});
    else
        catalog(k) = full(idx);
    end
end
end

%% ------------------------------------------------------------------
function e = makeEntry(field, key, label, fileSuffix, description)
e = struct('field', field, 'key', key, 'label', label, ...
    'fileSuffix', fileSuffix, 'description', description);
end

%% ------------------------------------------------------------------
function e = emptyEntry()
e = makeEntry('', '', '', '', '');
end

%% ------------------------------------------------------------------
function e = derivedEntry(field)
parts = strsplit(field, '.');
leaf = parts{end};
key = matlab.lang.makeValidName(leaf);
e = makeEntry(field, key, strrep(leaf, '_', '\_'), ...
    matlab.lang.makeValidName(leaf), ...
    sprintf('user-supplied metric "%s"', field));
end
