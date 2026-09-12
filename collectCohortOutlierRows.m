function rows = collectCohortOutlierRows(summary, metricField)
%COLLECTCOHORTOUTLIERROWS  Significant per-node findings from a cohort comparison.
%
%   rows = collectCohortOutlierRows(summary)
%   rows = collectCohortOutlierRows(summary, 'D_susceptibility')
%
%   Flattens the nested per-subject deviation structs of a
%   COMPARECOHORTGROUPS summary into one row per (case subject, flagged
%   node) for the chosen metric, ready for PLOTCOHORTOUTLIERSUMMARY.
%
%   Reading the findings from the in-memory summary rather than from the
%   written CSVs keeps the figures tied to exactly the numbers the
%   statistics produced, with no re-parsing step to drift out of sync.
%
%   Row fields
%     subjectId    : case subject the finding belongs to
%     region       : full node label
%     baseName     : label with the hemisphere prefix stripped
%     hemisphere   : 'L', 'R' or ''
%     value        : the case subject's value
%     controlMean  : control-cohort mean for that node
%     effectPct    : percentage excess over the control mean
%     z            : z-score against the control cohort
%     pBonf        : Bonferroni-corrected p-value
%
%   See also COMPARECOHORTGROUPS, PLOTCOHORTOUTLIERSUMMARY,
%   COLLECTLATERALITYOUTLIERROWS.

if nargin < 2 || isempty(metricField)
    metricField = 'D_susceptibility';
end

rows = emptyRows();
if ~isstruct(summary) || ~isfield(summary, 'deviations')
    return;
end

metricIdx = find(strcmp({summary.metrics.field}, metricField), 1);
if isempty(metricIdx)
    warning('collectCohortOutlierRows:MetricMissing', ...
        'Metric "%s" was not part of this comparison.', metricField);
    return;
end

labels = summary.labels;
controlMean = summary.controlMean(:, metricIdx);
caseIds = summary.caseIds;

for a = 1:numel(caseIds)
    fld = matlab.lang.makeValidName(caseIds{a});
    if ~isfield(summary.deviations, fld)
        continue;
    end
    dev = summary.deviations.(fld);

    flagged = find(dev.isOutlier(:, metricIdx));
    for k = 1:numel(flagged)
        i = flagged(k);
        [baseName, hemisphere] = stripHemispherePrefix(labels{i});
        value = dev.values(i, metricIdx);
        mu = controlMean(i);
        if ~isfinite(mu) || mu == 0
            effectPct = NaN;
        else
            effectPct = 100 * (value - mu) / mu;
        end

        rows(end+1) = struct( ...
            'subjectId',   caseIds{a}, ...
            'region',      char(labels{i}), ...
            'baseName',    baseName, ...
            'hemisphere',  hemisphere, ...
            'value',       value, ...
            'controlMean', mu, ...
            'effectPct',   effectPct, ...
            'z',           dev.z(i, metricIdx), ...
            'pBonf',       dev.pBonf(i, metricIdx)); %#ok<AGROW>
    end
end
end

%% ------------------------------------------------------------------
function rows = emptyRows()
rows = struct('subjectId', {}, 'region', {}, 'baseName', {}, 'hemisphere', {}, ...
    'value', {}, 'controlMean', {}, 'effectPct', {}, 'z', {}, 'pBonf', {});
end
