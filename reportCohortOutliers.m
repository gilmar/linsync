function report = reportCohortOutliers(varargin)
%REPORTCOHORTOUTLIERS  Print the flagged nodes of a case-vs-control comparison.
%
%   reportCohortOutliers('ResultsDir', resultsDir, 'Normalisation', 'column')
%   report = reportCohortOutliers(..., 'ZThreshold', 2.5)
%
%   Reads the summary written by COMPARECOHORTGROUPS and prints, per case
%   subject, every node that survived Bonferroni correction, with its
%   z-score and corrected p-value.
%
%   'ZThreshold' filters the printed list only. Whether a node counts as an
%   outlier was decided by the Bonferroni correction at run time; raising
%   the threshold here just shortens the console output for a quick scan,
%   and lowering it cannot resurrect nodes that failed the correction.
%
%   Name-value options
%     'ResultsDir'    : ''
%     'Normalisation' : 'column'
%     'ZThreshold'    : 2.0 -- minimum |z| to display
%
%   Returns a table of the printed rows (empty when nothing is flagged).
%
%   See also COMPARECOHORTGROUPS, REPORTTOPNODES.

p = inputParser;
addParameter(p, 'ResultsDir',    '', @(s) ischar(s) || isstring(s));
addParameter(p, 'Normalisation', 'column', @(s) ischar(s) || isstring(s));
addParameter(p, 'ZThreshold',    2.0, @isscalar);
parse(p, varargin{:});
opts = p.Results;

scheme = normalisationSchemeInfo(opts.Normalisation).scheme;
resultsDir = resolveExperimentResultsDir(opts.ResultsDir);

matFile = locateExperimentResultsFile(resultsDir, 'cohort_comparison', ...
    sprintf('compare_cohort_groups_%s.mat', scheme));
if isempty(matFile)
    error('reportCohortOutliers:NoSummary', ...
        ['No compare_cohort_groups_%s.mat under %s. Run compareCohortGroups ' ...
         'for this scheme first.'], scheme, resultsDir);
end
S = load(matFile);
if ~isfield(S, 'summary')
    error('reportCohortOutliers:BadSummary', ...
        '%s does not contain a summary struct.', matFile);
end
summary = S.summary;
info = cohortGroupInfo(summary);

rowsOut = cell(0, 1);
fprintf('\nCase-vs-control outliers  --  scheme = %s\n', scheme);
fprintf('%d %s subject(s) vs %d %s subject(s); Bonferroni alpha = %.3g; display |z| >= %.2f\n', ...
    info.nCase, info.caseLabel, info.nControl, info.controlLabel, ...
    summary.alpha, opts.ZThreshold);

for a = 1:numel(summary.caseIds)
    caseId = summary.caseIds{a};
    fld = matlab.lang.makeValidName(caseId);
    if ~isfield(summary.deviations, fld)
        continue;
    end
    dev = summary.deviations.(fld);

    fprintf('\n--- %s  (family size m = %d, Bonferroni |z| >= %.2f) ---\n', ...
        caseId, dev.nTests, dev.zThreshold);

    printed = 0;
    for m = 1:numel(summary.metrics)
        flagged = find(dev.isOutlier(:, m) & abs(dev.z(:, m)) >= opts.ZThreshold);
        if isempty(flagged)
            continue;
        end
        [~, order] = sort(dev.pBonf(flagged, m), 'ascend');
        flagged = flagged(order);

        fprintf('  %s:\n', summary.metrics(m).field);
        for k = 1:numel(flagged)
            i = flagged(k);
            fprintf('    %-32s  value = %10.4g  control mean = %10.4g  z = %+6.2f  p_bonf = %.3g\n', ...
                summary.labels{i}, dev.values(i, m), summary.controlMean(i, m), ...
                dev.z(i, m), dev.pBonf(i, m));
            rowsOut{end+1, 1} = {caseId, summary.metrics(m).field, ...
                summary.labels{i}, dev.values(i, m), summary.controlMean(i, m), ...
                dev.z(i, m), dev.pBonf(i, m)}; %#ok<AGROW>
            printed = printed + 1;
        end
    end

    if printed == 0
        fprintf('  (nothing flagged at |z| >= %.2f)\n', opts.ZThreshold);
    end
end
fprintf('\n');

if isempty(rowsOut)
    report = table();
else
    report = cell2table(vertcat(rowsOut{:}), 'VariableNames', ...
        {'subject', 'metric', 'region', 'value', 'control_mean', 'z', 'p_bonf'});
end
end
