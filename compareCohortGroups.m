function summary = compareCohortGroups(varargin)
%COMPARECOHORTGROUPS  Per-node case-vs-control comparison of node metrics.
%
%   summary = compareCohortGroups('ResultsDir', resultsDir, ...
%                'Normalisation', 'column', 'CasePrefix', 'Case', ...
%                'ControlPrefix', 'Ctrl')
%
%   For each case subject, compares every requested node-level metric
%   against the control cohort's per-node mean and SD, and flags the nodes
%   whose deviation survives Bonferroni correction across
%   (non-trivial nodes x metrics) for that subject.
%
%   Why Bonferroni rather than a bare |z| >= 2 rule: with N nodes and
%   several metrics there are hundreds of implicit tests per subject, so
%   |z| >= 2 would be expected to "find" a dozen regions in pure noise.
%   Correction is what makes a flagged region worth reporting. The
%   'ZThreshold' option therefore only controls what downstream reports
%   display -- it never decides what counts as an outlier.
%
%   Reads the per-subject files written by RUNCOHORTSTABILITYCENTRALITIES.
%
%   Name-value options
%     'ResultsDir'    : ''       -- run folder to read and write
%     'Normalisation' : 'column' -- which scheme's results to load
%     'CasePrefix'    : 'Case'
%     'ControlPrefix' : 'Control'
%     'Metrics'       : {}       -- dotted field paths; default is
%                       D(->i), D(k->) and betweenness when available
%     'Alpha'         : 0.05     -- family-wise level for the correction
%     'ZThreshold'    : 2.0      -- display-only threshold for reports
%     'SaveResults'   : true
%     'Plot'          : true
%     'RunParameters' : []
%
%   Outputs (under <ResultsDir>/cohort_comparison/)
%     compare_cohort_groups_<subject>_<scheme>_<metric>.{fig,png}
%     compare_cohort_groups_<subject>_<scheme>_outliers.csv  (when any)
%     compare_cohort_groups_<scheme>.mat
%
%   See also RUNCOHORTSTABILITYCENTRALITIES, COMPAREHEMISPHERICASYMMETRY,
%   ZSCOREOUTLIERS, REPORTCOHORTOUTLIERS.

p = inputParser;
addParameter(p, 'ResultsDir',    '',        @(s) ischar(s) || isstring(s));
addParameter(p, 'Normalisation', 'column',  @(s) ischar(s) || isstring(s));
addParameter(p, 'CasePrefix',    'Case',    @(s) ischar(s) || isstring(s));
addParameter(p, 'ControlPrefix', 'Control', @(s) ischar(s) || isstring(s));
addParameter(p, 'Metrics',       {},   @(c) iscell(c) || isstring(c));
addParameter(p, 'Alpha',         0.05, @(x) isscalar(x) && x > 0 && x < 1);
addParameter(p, 'ZThreshold',    2.0,  @isscalar);
addParameter(p, 'SaveResults',   true, @islogical);
addParameter(p, 'Plot',          true, @islogical);
addParameter(p, 'RunParameters', [],   @(x) isempty(x) || isstruct(x));
parse(p, varargin{:});
opts = p.Results;

scheme = normalisationSchemeInfo(opts.Normalisation).scheme;
resultsDir = resolveExperimentResultsDir(opts.ResultsDir);
outDir = experimentResultsDir(resultsDir, 'cohort_comparison');

cohort = loadCohortResults(resultsDir, scheme, ...
    'CasePrefix', opts.CasePrefix, ...
    'ControlPrefix', opts.ControlPrefix, ...
    'RequireGroups', true);

%% ---- Decide which metrics to compare -------------------------------
metricFields = cellstr(opts.Metrics);
if isempty(metricFields)
    metricFields = defaultCohortMetrics(cohort, 'compareCohortGroups');
end
metrics = connectomeMetricCatalog(metricFields);
nMetrics = numel(metrics);

caseIds = cohort.caseIds;
controlIds = cohort.controlIds;
nControl = numel(controlIds);

fprintf('compareCohortGroups: %d case / %d control subject(s), %d metric(s), scheme %s.\n', ...
    numel(caseIds), nControl, nMetrics, scheme);
fprintf('  case   : %s\n', strjoin(caseIds, ', '));
fprintf('  control: %s\n', strjoin(controlIds, ', '));

%% ---- Pack metrics and control statistics ---------------------------
caseValues    = cell(nMetrics, 1);
controlValues = cell(nMetrics, 1);
controlMean   = NaN(cohort.N, nMetrics);
controlStd    = NaN(cohort.N, nMetrics);
metricUsable  = false(nMetrics, 1);
for m = 1:nMetrics
    caseValues{m}    = packCohortMetric(cohort, metrics(m).field, cohort.isCase);
    controlValues{m} = packCohortMetric(cohort, metrics(m).field, cohort.isControl);
    controlMean(:, m) = mean(controlValues{m}, 2, 'omitnan');
    controlStd(:, m)  = std(controlValues{m}, 0, 2, 'omitnan');
    metricUsable(m) = any(isfinite(controlValues{m}(:))) && any(isfinite(caseValues{m}(:)));
    if ~metricUsable(m)
        warning('compareCohortGroups:MetricUnavailable', ...
            ['Metric "%s" is missing or all-NaN in the loaded results and will ' ...
             'be skipped.'], metrics(m).field);
    end
end

keepMetric = metricUsable;
if ~any(keepMetric)
    error('compareCohortGroups:NoUsableMetrics', ...
        'None of the requested metrics are present in the results under %s.', resultsDir);
end
metrics = metrics(keepMetric);
caseValues = caseValues(keepMetric);
controlValues = controlValues(keepMetric);
controlMean = controlMean(:, keepMetric);
controlStd = controlStd(:, keepMetric);
nMetrics = numel(metrics);

% A node that is trivial in any control subject has a meaningless control
% mean, so it is excluded from the testing family for every case subject.
includeMask = repmat(~cohort.isTrivialControl(:), 1, nMetrics);

%% ---- One case subject at a time ------------------------------------
deviations = struct();
for a = 1:numel(caseIds)
    caseId = caseIds{a};
    caseTex = strrep(caseId, '_', '\_');

    values = NaN(cohort.N, nMetrics);
    for m = 1:nMetrics
        values(:, m) = caseValues{m}(:, a);
    end

    stats = zScoreOutliers(values, controlMean, controlStd, ...
        'Alpha', opts.Alpha, 'Include', includeMask);

    baseName = sprintf('compare_cohort_groups_%s_%s', caseId, scheme);
    figureTitle = sprintf('Per-node metrics: %s vs %s cohort  --  normalisation = %s', ...
        caseTex, texEscape(cohort.controlPrefix), scheme);

    if opts.Plot
        for m = 1:nMetrics
            spec = struct();
            spec.caseValues    = values(:, m);
            spec.controlValues = controlValues{m};
            spec.controlMean   = controlMean(:, m);
            spec.controlStd    = controlStd(:, m);
            spec.z             = stats.z(:, m);
            spec.isOutlier     = stats.isOutlier(:, m);
            spec.labels        = cohort.labels;
            spec.isTrivial     = cohort.isTrivialControl;
            spec.yLabel        = metrics(m).label;
            spec.caseLabel     = caseId;
            spec.controlLabel  = cohort.controlPrefix;
            spec.figureTitle   = figureTitle;
            spec.figureName    = sprintf('%s vs %s - %s', ...
                caseId, cohort.controlPrefix, metrics(m).key);
            spec.panelTitle    = sprintf( ...
                '%s: %s vs %s mean \\pm SD (n=%d)  [Bonferroni |z|\\geq%.2f, \\alpha=%.3f, m=%d]', ...
                metrics(m).label, caseTex, texEscape(cohort.controlPrefix), ...
                nControl, stats.zThreshold, opts.Alpha, stats.nTests);

            fig = plotCohortComparisonFigure(spec);
            if opts.SaveResults
                saveFigureBoth(fig, fullfile(outDir, ...
                    sprintf('%s_%s', baseName, metrics(m).fileSuffix)));
            end
            close(fig);
        end
    end

    dev = struct();
    dev.subjectId = caseId;
    dev.metrics = {metrics.field};
    dev.values = values;
    dev.z = stats.z;
    dev.p = stats.p;
    dev.pBonf = stats.pBonf;
    dev.isOutlier = stats.isOutlier;
    dev.nTests = stats.nTests;
    dev.zThreshold = stats.zThreshold;
    deviations.(matlab.lang.makeValidName(caseId)) = dev;

    T = buildOutlierTable(cohort, metrics, values, controlMean, controlStd, stats);
    deviations.(matlab.lang.makeValidName(caseId)).outlierTable = T;

    if isempty(T)
        fprintf('  %s: no outlier nodes at Bonferroni alpha=%.3f (|z|>=%.2f, m=%d).\n', ...
            caseId, opts.Alpha, stats.zThreshold, stats.nTests);
    else
        fprintf('  %s: %d outlier node(s) (Bonferroni alpha=%.3f, |z|>=%.2f, m=%d).\n', ...
            caseId, height(T), opts.Alpha, stats.zThreshold, stats.nTests);
        if opts.SaveResults
            csvPath = fullfile(outDir, [baseName '_outliers.csv']);
            writetable(T, csvPath);
            fprintf('    -> %s\n', csvPath);
        end
    end
end

%% ---- Pack the summary ----------------------------------------------
summary = struct();
summary.scheme           = scheme;
summary.resultsDir       = resultsDir;
summary.subjects         = cohort.subjects;
summary.labels           = cohort.labels;
summary.isCase           = cohort.isCase;
summary.isControl        = cohort.isControl;
summary.caseIds          = caseIds;
summary.controlIds       = controlIds;
summary.casePrefix       = cohort.casePrefix;
summary.controlPrefix    = cohort.controlPrefix;
summary.metrics          = metrics;
summary.caseValues       = caseValues;
summary.controlValues    = controlValues;
summary.controlMean      = controlMean;
summary.controlStd       = controlStd;
summary.isTrivialControl = cohort.isTrivialControl;
summary.deviations       = deviations;
summary.alpha            = opts.Alpha;
summary.zThreshold       = opts.ZThreshold;

patch = experimentRunParameters('fromOptions', opts, ...
    'Analysis', 'compareCohortGroups', 'ResultsDir', resultsDir, ...
    'Scheme', scheme, 'Subjects', cohort.subjects);
summary.runParameters = experimentRunParameters('merge', opts.RunParameters, patch);

if opts.SaveResults
    save(fullfile(outDir, sprintf('compare_cohort_groups_%s.mat', scheme)), ...
        'summary');
end
end

%% ------------------------------------------------------------------
function T = buildOutlierTable(cohort, metrics, values, controlMean, controlStd, stats)
%BUILDOUTLIERTABLE  One row per flagged node, columns per metric.
outIdx = find(any(stats.isOutlier, 2));
if isempty(outIdx)
    T = [];
    return;
end

T = table(outIdx, string(cohort.labels(outIdx)), ...
    'VariableNames', {'NodeIdx', 'Region'});
for m = 1:numel(metrics)
    key = metrics(m).key;
    T.(['case_' key])         = values(outIdx, m);
    T.(['control_mean_' key]) = controlMean(outIdx, m);
    T.(['control_sd_' key])   = controlStd(outIdx, m);
    T.(['z_' key])            = stats.z(outIdx, m);
    T.(['p_' key])            = stats.p(outIdx, m);
    T.(['p_bonf_' key])       = stats.pBonf(outIdx, m);
end

% Most significant first, across whichever metric flagged the node.
minBonf = min(stats.pBonf(outIdx, :), [], 2, 'omitnan');
[~, order] = sort(minBonf, 'ascend');
T = T(order, :);
end

%% ------------------------------------------------------------------
function s = texEscape(str)
s = strrep(char(str), '_', '\_');
end
