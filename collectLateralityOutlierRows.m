function rows = collectLateralityOutlierRows(summary, metricField)
%COLLECTLATERALITYOUTLIERROWS  Significant laterality findings from an asymmetry run.
%
%   rows = collectLateralityOutlierRows(summary)
%   rows = collectLateralityOutlierRows(summary, 'D_influence')
%
%   Flattens the per-subject results of a COMPAREHEMISPHERICASYMMETRY
%   summary into one row per (case subject, flagged region pair) for the
%   chosen metric, in the same row shape that PLOTCOHORTOUTLIERSUMMARY
%   consumes.
%
%   The effect is reported in the metric's own units as the signed
%   difference L - R, together with the bounded LI_norm, so the figure can
%   show magnitude while the table keeps the scale-free version.
%
%   Row fields
%     subjectId, region (the region-pair base name), baseName,
%     hemisphere ('' -- a pair spans both), value (LI_signed),
%     controlMean (control LI_signed mean), effectPct (signed difference
%     relative to the control mean magnitude), z, pBonf, liNorm
%
%   See also COMPAREHEMISPHERICASYMMETRY, PLOTCOHORTOUTLIERSUMMARY.

if nargin < 2 || isempty(metricField)
    metricField = 'D_influence';
end

rows = emptyRows();
if ~isstruct(summary) || ~isfield(summary, 'perSubject')
    return;
end

metricIdx = find(strcmp({summary.metrics.field}, metricField), 1);
if isempty(metricIdx)
    % Fall back to the first available metric rather than returning nothing:
    % a symmetric scheme has no D(k->) and the caller should still get a figure.
    if isempty(summary.metrics)
        return;
    end
    metricIdx = 1;
    warning('collectLateralityOutlierRows:MetricMissing', ...
        'Metric "%s" was not part of this analysis; using "%s" instead.', ...
        metricField, summary.metrics(1).field);
end

pairLabels = summary.pairLabels;
ctrlSigned = summary.controlMeanSigned(:, metricIdx);
caseIds = summary.caseIds;

for a = 1:numel(caseIds)
    fld = matlab.lang.makeValidName(caseIds{a});
    if ~isfield(summary.perSubject, fld)
        continue;
    end
    dev = summary.perSubject.(fld);
    subjectCol = find(strcmp(summary.subjects, caseIds{a}), 1);

    flagged = find(dev.isOutlier(:, metricIdx));
    for k = 1:numel(flagged)
        i = flagged(k);
        signedLI = summary.LI_signed{metricIdx}(i, subjectCol);
        mu = ctrlSigned(i);
        if ~isfinite(mu) || mu == 0
            effectPct = NaN;
        else
            effectPct = 100 * (signedLI - mu) / abs(mu);
        end

        rows(end+1) = struct( ...
            'subjectId',   caseIds{a}, ...
            'region',      char(pairLabels{i}), ...
            'baseName',    char(pairLabels{i}), ...
            'hemisphere',  '', ...
            'value',       signedLI, ...
            'controlMean', mu, ...
            'effectPct',   effectPct, ...
            'z',           dev.z(i, metricIdx), ...
            'pBonf',       dev.pBonf(i, metricIdx), ...
            'liNorm',      summary.LI_norm{metricIdx}(i, subjectCol)); %#ok<AGROW>
    end
end
end

%% ------------------------------------------------------------------
function rows = emptyRows()
rows = struct('subjectId', {}, 'region', {}, 'baseName', {}, 'hemisphere', {}, ...
    'value', {}, 'controlMean', {}, 'effectPct', {}, 'z', {}, 'pBonf', {}, ...
    'liNorm', {});
end
