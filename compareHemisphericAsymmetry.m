function summary = compareHemisphericAsymmetry(varargin)
%COMPAREHEMISPHERICASYMMETRY  Left-right laterality of node metrics, case vs control.
%
%   summary = compareHemisphericAsymmetry('ResultsDir', resultsDir, ...
%                'Normalisation', 'column', 'CasePrefix', 'Case', ...
%                'ControlPrefix', 'Ctrl')
%
%   Pairs each left-hemisphere node with its right-hemisphere twin and asks
%   whether the case cohort is more lateralised than the control cohort.
%   A unilateral pathology should show up as asymmetry between homotopic
%   regions even when the absolute values vary between subjects, which is
%   what makes laterality a more sensitive readout than raw per-node values.
%
%   Three complementary views are produced, because they answer different
%   questions and can disagree:
%
%     1. Per case subject -- z-score of each pair's laterality index
%        against the control cohort mean +- SD, Bonferroni-corrected over
%        (pairs x metrics). Answers "which regions are unusual in THIS
%        subject".
%     2. Group level -- Welch t-test of the laterality index per pair,
%        case group vs control group, Bonferroni-corrected. Answers "which
%        regions differ between the groups".
%     3. Systematic direction -- pooled laterality across all pairs, tested
%        against zero with a signed-rank test, per group. Answers "does one
%        hemisphere dominate overall".
%
%   Laterality bases
%     LI_norm   = (L - R) / (L + R)   bounded, scale-free, comparable
%                                     across metrics and subjects
%     LI_signed = L - R               in the metric's own units, so the
%                                     effect size stays interpretable
%
%   'LateralityBasis' selects which one drives the per-subject z-scores and
%   flagging. Views 2 and 3 always use LI_norm: pooling or comparing across
%   subjects requires the scale-free version.
%
%   Name-value options
%     'ResultsDir'      : ''
%     'Normalisation'   : 'column'
%     'CasePrefix'      : 'Case'
%     'ControlPrefix'   : 'Control'
%     'LeftPrefix'      : 'L-'     -- atlas prefix marking the left hemisphere
%     'RightPrefix'     : 'R-'
%     'Metrics'         : {}       -- default: see DEFAULTCOHORTMETRICS
%     'LateralityBasis' : 'signed' | 'norm'
%     'Alpha'           : 0.05
%     'ZThreshold'      : 2.0      -- display only
%     'SaveResults'     : true
%     'Plot'            : true
%     'RunParameters'   : []
%
%   Outputs (under <ResultsDir>/hemispheric_asymmetry/)
%     LR_asymmetry_<subject>_<scheme>.{fig,png}, _outliers.csv
%     LR_asymmetry_groupTest_<scheme>.{csv,fig,png}
%     LR_asymmetry_systematic_<scheme>.{csv,fig,png}
%     LR_asymmetry_<scheme>.mat
%
%   See also PAIRHEMISPHERENODES, LATERALITYINDEX, COMPARECOHORTGROUPS,
%   WELCHTTEST, SIGNRANKTEST.

p = inputParser;
addParameter(p, 'ResultsDir',      '',        @(s) ischar(s) || isstring(s));
addParameter(p, 'Normalisation',   'column',  @(s) ischar(s) || isstring(s));
addParameter(p, 'CasePrefix',      'Case',    @(s) ischar(s) || isstring(s));
addParameter(p, 'ControlPrefix',   'Control', @(s) ischar(s) || isstring(s));
addParameter(p, 'LeftPrefix',      'L-',      @(s) ischar(s) || isstring(s));
addParameter(p, 'RightPrefix',     'R-',      @(s) ischar(s) || isstring(s));
addParameter(p, 'Metrics',         {},        @(c) iscell(c) || isstring(c));
addParameter(p, 'LateralityBasis', 'signed',  @(s) ischar(s) || isstring(s));
addParameter(p, 'Alpha',           0.05, @(x) isscalar(x) && x > 0 && x < 1);
addParameter(p, 'ZThreshold',      2.0,  @isscalar);
addParameter(p, 'SaveResults',     true, @islogical);
addParameter(p, 'Plot',            true, @islogical);
addParameter(p, 'RunParameters',   [],   @(x) isempty(x) || isstruct(x));
parse(p, varargin{:});
opts = p.Results;

basis = lower(strtrim(char(opts.LateralityBasis)));
if ~ismember(basis, {'norm', 'signed'})
    error('compareHemisphericAsymmetry:BadBasis', ...
        'LateralityBasis must be ''signed'' or ''norm'' (got ''%s'').', basis);
end

scheme = normalisationSchemeInfo(opts.Normalisation).scheme;
resultsDir = resolveExperimentResultsDir(opts.ResultsDir);
outDir = experimentResultsDir(resultsDir, 'hemispheric_asymmetry');

cohort = loadCohortResults(resultsDir, scheme, ...
    'CasePrefix', opts.CasePrefix, ...
    'ControlPrefix', opts.ControlPrefix, ...
    'RequireGroups', true);

%% ---- Pair the hemispheres -------------------------------------------
pairs = pairHemisphereNodes(cohort.labels, ...
    'LeftPrefix', opts.LeftPrefix, 'RightPrefix', opts.RightPrefix, ...
    'WarnUnpaired', false);
if pairs.nPairs == 0
    error('compareHemisphericAsymmetry:NoPairs', ...
        ['No left/right node pairs found with prefixes "%s" and "%s". ' ...
         'Set LeftPrefix and RightPrefix to match this atlas.'], ...
        char(opts.LeftPrefix), char(opts.RightPrefix));
end

% A pair whose left or right node is trivial in any control subject cannot
% be compared against a control mean, so it is dropped up front.
pairExcluded = cohort.isTrivialControl(pairs.leftIdx) | ...
               cohort.isTrivialControl(pairs.rightIdx);
pairLeftIdx  = pairs.leftIdx(~pairExcluded);
pairRightIdx = pairs.rightIdx(~pairExcluded);
pairLabels   = pairs.baseLabels(~pairExcluded);
nPairs = numel(pairLeftIdx);
if nPairs == 0
    error('compareHemisphericAsymmetry:AllPairsExcluded', ...
        'Every left/right pair is trivial in at least one control subject.');
end

%% ---- Metrics and laterality indices ---------------------------------
metricFields = cellstr(opts.Metrics);
if isempty(metricFields)
    metricFields = defaultCohortMetrics(cohort, 'compareHemisphericAsymmetry');
end
metrics = connectomeMetricCatalog(metricFields);

keep = false(numel(metrics), 1);
nodeValues = cell(numel(metrics), 1);
for m = 1:numel(metrics)
    nodeValues{m} = packCohortMetric(cohort, metrics(m).field);
    keep(m) = any(isfinite(nodeValues{m}(:)));
    if ~keep(m)
        warning('compareHemisphericAsymmetry:MetricUnavailable', ...
            'Metric "%s" is missing or all-NaN; skipping it.', metrics(m).field);
    end
end
metrics = metrics(keep);
nodeValues = nodeValues(keep);
nMetrics = numel(metrics);
if nMetrics == 0
    error('compareHemisphericAsymmetry:NoUsableMetrics', ...
        'None of the requested metrics are present in the results under %s.', resultsDir);
end

LI = struct('norm', {cell(nMetrics, 1)}, 'signed', {cell(nMetrics, 1)});
for m = 1:nMetrics
    LI.norm{m}   = lateralityIndex(nodeValues{m}, pairLeftIdx, pairRightIdx, 'norm');
    LI.signed{m} = lateralityIndex(nodeValues{m}, pairLeftIdx, pairRightIdx, 'signed');
end

caseIds = cohort.caseIds;
controlIds = cohort.controlIds;
nControl = numel(controlIds);

fprintf('compareHemisphericAsymmetry: %d case / %d control subject(s), %d pair(s), %d metric(s).\n', ...
    numel(caseIds), nControl, nPairs, nMetrics);
fprintf('  per-subject basis: LI_%s\n', basis);

% Control statistics for both bases, so the outlier tables can always
% report the raw values regardless of which basis drove the flagging.
ctrlMean = struct('norm', {NaN(nPairs, nMetrics)}, 'signed', {NaN(nPairs, nMetrics)});
ctrlStd  = struct('norm', {NaN(nPairs, nMetrics)}, 'signed', {NaN(nPairs, nMetrics)});
for m = 1:nMetrics
    for b = {'norm', 'signed'}
        vals = LI.(b{1}){m}(:, cohort.isControl);
        ctrlMean.(b{1})(:, m) = mean(vals, 2, 'omitnan');
        ctrlStd.(b{1})(:, m)  = std(vals, 0, 2, 'omitnan');
    end
end

if strcmp(basis, 'signed')
    basisTex = 'LI_{signed}';
    zeroLineLabel = 'L - R = 0';
    clampToUnit = false;
else
    basisTex = 'LI_{norm}';
    zeroLineLabel = 'LI_{norm} = 0';
    clampToUnit = true;
end

%% ---- View 1: per case subject ---------------------------------------
perSubject = struct();
outlierTables = struct();
for a = 1:numel(caseIds)
    caseId = caseIds{a};
    subjectCol = find(strcmp(cohort.subjects, caseId), 1);
    caseTex = strrep(caseId, '_', '\_');

    values = NaN(nPairs, nMetrics);
    controlBasis = cell(nMetrics, 1);
    for m = 1:nMetrics
        values(:, m) = LI.(basis){m}(:, subjectCol);
        controlBasis{m} = LI.(basis){m}(:, cohort.isControl);
    end

    stats = zScoreOutliers(values, ctrlMean.(basis), ctrlStd.(basis), 'Alpha', opts.Alpha);

    if opts.Plot
        fig = figure('Name', sprintf('L-R asymmetry: %s vs %s', caseId, cohort.controlPrefix), ...
            'Position', [60 60 1700 380 * nMetrics + 80]);
        for m = 1:nMetrics
            panelTitle = sprintf('%s (%s): %s vs %s (n=%d)  [Bonferroni |z|\\geq%.2f]', ...
                basisTex, metrics(m).label, caseTex, ...
                strrep(cohort.controlPrefix, '_', '\_'), nControl, stats.zThreshold);
            plotLateralityPanel(nMetrics, m, values(:, m), controlBasis{m}, ...
                ctrlMean.(basis)(:, m), ctrlStd.(basis)(:, m), stats.z(:, m), ...
                pairLabels, stats.isOutlier(:, m), panelTitle, ...
                sprintf('%s (%s)', basisTex, metrics(m).label), ...
                clampToUnit, zeroLineLabel, caseId, cohort.controlPrefix);
        end
        sgtitle(sprintf('Left-right laterality: %s vs %s  --  %s  (basis = %s)', ...
            caseTex, strrep(cohort.controlPrefix, '_', '\_'), scheme, basis), ...
            'Interpreter', 'tex');

        if opts.SaveResults
            saveFigureBoth(fig, fullfile(outDir, ...
                sprintf('LR_asymmetry_%s_%s', caseId, scheme)));
        end
        close(fig);
    end

    dev = struct();
    dev.subjectId = caseId;
    dev.basis = basis;
    dev.metrics = {metrics.field};
    dev.LI = values;
    dev.z = stats.z;
    dev.p = stats.p;
    dev.pBonf = stats.pBonf;
    dev.isOutlier = stats.isOutlier;
    dev.nTests = stats.nTests;
    dev.zThreshold = stats.zThreshold;
    perSubject.(matlab.lang.makeValidName(caseId)) = dev;

    T = buildLateralityOutlierTable(pairLabels, pairLeftIdx, pairRightIdx, ...
        metrics, LI, ctrlMean, ctrlStd, subjectCol, stats, basis);
    outlierTables.(matlab.lang.makeValidName(caseId)) = T;

    if isempty(T)
        fprintf('  %s: no Bonferroni-significant pairs at alpha=%.3f.\n', caseId, opts.Alpha);
    else
        fprintf('  %s: %d significant pair(s) (Bonferroni alpha=%.3f, m=%d).\n', ...
            caseId, height(T), opts.Alpha, stats.nTests);
        if opts.SaveResults
            writetable(T, fullfile(outDir, ...
                sprintf('LR_asymmetry_%s_%s_outliers.csv', caseId, scheme)));
        end
    end
end

%% ---- View 2: group-level Welch tests --------------------------------
groupTest = struct();
groupRows = cell(0, 1);
pAll = [];
for m = 1:nMetrics
    caseLI = LI.norm{m}(:, cohort.isCase);
    ctrlLI = LI.norm{m}(:, cohort.isControl);

    tStat = NaN(nPairs, 1);
    dfVal = NaN(nPairs, 1);
    pVal  = NaN(nPairs, 1);
    for q = 1:nPairs
        res = welchTTest(caseLI(q, :), ctrlLI(q, :));
        tStat(q) = res.tStat;
        dfVal(q) = res.df;
        pVal(q)  = res.p;
    end

    groupTest.(metrics(m).key) = struct( ...
        'metricField', metrics(m).field, ...
        'tStat', tStat, 'df', dfVal, 'p', pVal, ...
        'LI_case_mean', mean(caseLI, 2, 'omitnan'), ...
        'LI_control_mean', mean(ctrlLI, 2, 'omitnan'));
    pAll = [pAll; pVal(:)]; %#ok<AGROW>
end

nGroupTests = max(sum(isfinite(pAll)), 1);
for m = 1:nMetrics
    key = metrics(m).key;
    gt = groupTest.(key);
    gt.p_bonf = min(1, gt.p * nGroupTests);
    gt.bonf_significant = isfinite(gt.p_bonf) & gt.p_bonf < opts.Alpha;
    groupTest.(key) = gt;
    for q = 1:nPairs
        groupRows{end+1, 1} = {pairLabels{q}, key, ...
            gt.LI_case_mean(q), gt.LI_control_mean(q), ...
            gt.tStat(q), gt.df(q), gt.p(q), gt.p_bonf(q), gt.bonf_significant(q)}; %#ok<AGROW>
    end
end
Tgroup = cell2table(vertcat(groupRows{:}), 'VariableNames', ...
    {'Region', 'Metric', 'LI_case_mean', 'LI_control_mean', ...
     't_stat', 'df', 'p', 'p_bonf', 'bonf_significant'});

%% ---- View 3: systematic direction -----------------------------------
systematic = struct();
sysRows = cell(0, 1);
groupTags = {'case', 'control'};
groupMasks = {cohort.isCase, cohort.isControl};
groupNames = {cohort.casePrefix, cohort.controlPrefix};
for m = 1:nMetrics
    key = metrics(m).key;
    for g = 1:numel(groupTags)
        vals = LI.norm{m}(:, groupMasks{g});
        v = vals(:);
        v = v(isfinite(v));
        res = signRankTest(v, 0);
        meanLI = mean(v, 'omitnan');
        sdLI = std(v, 0, 'omitnan');
        side = inferDominantSide(meanLI, res.p, opts.Alpha);
        systematic.(key).(groupTags{g}) = struct( ...
            'groupName', groupNames{g}, ...
            'nObs', numel(v), 'meanLI', meanLI, 'sdLI', sdLI, ...
            'signrankP', res.p, 'signrankZ', res.z, 'dominantSide', side);
        sysRows{end+1, 1} = {groupTags{g}, groupNames{g}, key, numel(v), ...
            meanLI, sdLI, res.p, side}; %#ok<AGROW>
    end
end
Tsys = cell2table(vertcat(sysRows{:}), 'VariableNames', ...
    {'group', 'group_name', 'metric', 'n_obs', 'mean_LI', 'sd_LI', ...
     'signrank_p', 'inferred_dominant_side'});

%% ---- Group-level figures --------------------------------------------
if opts.Plot
    figGroup = plotGroupTestFigure(scheme, pairLabels, metrics, LI.norm, ...
        cohort, groupTest, opts.Alpha);
    figSys = plotSystematicFigure(scheme, metrics, systematic, groupTags, groupNames);
    if opts.SaveResults
        saveFigureBoth(figGroup, fullfile(outDir, sprintf('LR_asymmetry_groupTest_%s', scheme)));
        saveFigureBoth(figSys, fullfile(outDir, sprintf('LR_asymmetry_systematic_%s', scheme)));
    end
    close(figGroup);
    close(figSys);
end

if opts.SaveResults
    writetable(Tgroup, fullfile(outDir, sprintf('LR_asymmetry_groupTest_%s.csv', scheme)));
    writetable(Tsys,   fullfile(outDir, sprintf('LR_asymmetry_systematic_%s.csv', scheme)));
end

%% ---- Summary ---------------------------------------------------------
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
summary.leftPrefix       = char(opts.LeftPrefix);
summary.rightPrefix      = char(opts.RightPrefix);
summary.metrics          = metrics;
summary.pairLabels       = pairLabels;
summary.pairLeftIdx      = pairLeftIdx;
summary.pairRightIdx     = pairRightIdx;
summary.pairExcluded     = pairExcluded;
summary.LI_norm          = LI.norm;
summary.LI_signed        = LI.signed;
summary.controlMeanNorm  = ctrlMean.norm;
summary.controlStdNorm   = ctrlStd.norm;
summary.controlMeanSigned = ctrlMean.signed;
summary.controlStdSigned  = ctrlStd.signed;
summary.lateralityBasis  = basis;
summary.perSubject       = perSubject;
summary.outlierTables    = outlierTables;
summary.groupTest        = groupTest;
summary.groupTestTable   = Tgroup;
summary.systematic       = systematic;
summary.systematicTable  = Tsys;
summary.alpha            = opts.Alpha;
summary.zThreshold       = opts.ZThreshold;

patch = experimentRunParameters('fromOptions', opts, ...
    'Analysis', 'compareHemisphericAsymmetry', 'ResultsDir', resultsDir, ...
    'Scheme', scheme, 'Subjects', cohort.subjects);
summary.runParameters = experimentRunParameters('merge', opts.RunParameters, patch);

if opts.SaveResults
    save(fullfile(outDir, sprintf('LR_asymmetry_%s.mat', scheme)), 'summary');
end
end

%% ==================================================================
function plotLateralityPanel(nPanels, panelIdx, caseVals, controlVals, ...
    ctrlMean, ctrlStd, zScore, pairLabels, outlierMask, panelTitle, ylab, ...
    clampToUnit, zeroLineLabel, caseLabel, controlLabel)
%PLOTLATERALITYPANEL  One metric's laterality panel: case against the control band.

ax = subplot(nPanels, 1, panelIdx);
hold(ax, 'on');
nP = numel(caseVals);
x = 1:nP;

yl = [min([caseVals(:); ctrlMean - ctrlStd; controlVals(:)], [], 'omitnan'), ...
      max([caseVals(:); ctrlMean + ctrlStd; controlVals(:)], [], 'omitnan')];
if ~all(isfinite(yl)) || diff(yl) <= 0
    yl = [-1 1];
end
pad = 0.05 * max(diff(yl), 0.1);
yl = yl + [-pad pad];
if clampToUnit
    % LI_norm is bounded in [-1, 1]; clamping keeps the axis honest.
    yl(1) = max(yl(1), -1.05);
    yl(2) = min(yl(2),  1.05);
end

yline(ax, 0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, ...
    'DisplayName', zeroLineLabel);

valid = isfinite(ctrlMean) & isfinite(ctrlStd);
if any(valid)
    xv = x(valid);
    yLo = ctrlMean - ctrlStd;
    yHi = ctrlMean + ctrlStd;
    fill(ax, [xv fliplr(xv)], [yLo(valid).' fliplr(yHi(valid).')], ...
        [0.80 0.84 0.95], 'EdgeColor', 'none', 'FaceAlpha', 0.65, ...
        'DisplayName', sprintf('%s mean \\pm SD', strrep(controlLabel, '_', '\_')));
end

plot(ax, x, ctrlMean, '-', 'Color', [0.20 0.30 0.70], 'LineWidth', 1.2, ...
    'DisplayName', sprintf('%s mean', strrep(controlLabel, '_', '\_')));

for j = 1:size(controlVals, 2)
    plot(ax, x, controlVals(:, j), '.', 'Color', [0.55 0.60 0.85], ...
        'MarkerSize', 6, 'HandleVisibility', 'off');
end

plot(ax, x, caseVals, 'o', 'Color', [0.65 0.10 0.10], ...
    'MarkerFaceColor', [0.85 0.33 0.10], 'MarkerSize', 5, ...
    'DisplayName', strrep(caseLabel, '_', '\_'));

outlierIdx = find(outlierMask & isfinite(zScore));
if ~isempty(outlierIdx)
    plot(ax, x(outlierIdx), caseVals(outlierIdx), 'o', 'Color', 'k', ...
        'MarkerSize', 10, 'LineWidth', 1.5, 'HandleVisibility', 'off');
    for k = 1:numel(outlierIdx)
        i = outlierIdx(k);
        text(ax, i, caseVals(i), sprintf('  %s (z=%+.1f)', pairLabels{i}, zScore(i)), ...
            'Interpreter', 'none', 'FontSize', 7, 'Rotation', 25, ...
            'VerticalAlignment', 'bottom');
    end
end

xlim(ax, [0.5 nP + 0.5]);
ylim(ax, yl);
set(ax, 'XTick', 1:nP, 'XTickLabel', pairLabels, ...
    'TickLabelInterpreter', 'none', 'FontSize', 6);
xtickangle(ax, 45);
grid(ax, 'on');
box(ax, 'on');
ylabel(ax, ylab, 'Interpreter', 'tex');
title(ax, panelTitle, 'Interpreter', 'tex');
legend(ax, 'Location', 'best');
end

%% ------------------------------------------------------------------
function fig = plotGroupTestFigure(scheme, pairLabels, metrics, LInorm, ...
    cohort, groupTest, alpha)
%PLOTGROUPTESTFIGURE  Case vs control mean laterality per pair, one panel per metric.
% Points on the diagonal mean the two groups lateralise identically; the
% distance from the diagonal is the group effect.

nMet = numel(metrics);
fig = figure('Name', sprintf('L-R group test (%s)', scheme), ...
    'Position', [80 80 520 * nMet 480]);
nPairs = numel(pairLabels);

for m = 1:nMet
    key = metrics(m).key;
    ax = subplot(1, nMet, m);
    hold(ax, 'on');
    gt = groupTest.(key);

    caseLI = LInorm{m}(:, cohort.isCase);
    ctrlLI = LInorm{m}(:, cohort.isControl);
    all = [caseLI(:); ctrlLI(:)];
    all = all(isfinite(all));
    if isempty(all)
        lims = [-1 1];
    else
        pad = 0.1 * max(max(all) - min(all), 0.1);
        lims = [min(all) - pad, max(all) + pad];
    end

    for q = 1:nPairs
        xA = gt.LI_case_mean(q);
        xR = gt.LI_control_mean(q);
        if ~isfinite(xA) || ~isfinite(xR)
            continue;
        end
        if gt.bonf_significant(q)
            col = [0.85 0.20 0.10];
            mk = 8;
        else
            col = [0.40 0.40 0.40];
            mk = 4;
        end
        plot(ax, xR, xA, 'o', 'Color', col, 'MarkerFaceColor', col, ...
            'MarkerSize', mk, 'HandleVisibility', 'off');
    end

    plot(ax, lims, lims, 'k--', 'HandleVisibility', 'off');
    xlabel(ax, sprintf('%s mean LI_{norm} (%s)', ...
        strrep(cohort.controlPrefix, '_', '\_'), metrics(m).label), 'Interpreter', 'tex');
    ylabel(ax, sprintf('%s mean LI_{norm} (%s)', ...
        strrep(cohort.casePrefix, '_', '\_'), metrics(m).label), 'Interpreter', 'tex');
    title(ax, sprintf('%s: pair means (red = Bonferroni p<%.2f)', metrics(m).label, alpha), ...
        'Interpreter', 'tex');
    axis(ax, 'equal');
    xlim(ax, lims);
    ylim(ax, lims);
    grid(ax, 'on');
end

sgtitle(sprintf('Group-level L-R laterality: %s vs %s  --  %s', ...
    strrep(cohort.casePrefix, '_', '\_'), strrep(cohort.controlPrefix, '_', '\_'), scheme), ...
    'Interpreter', 'tex');
end

%% ------------------------------------------------------------------
function fig = plotSystematicFigure(scheme, metrics, systematic, groupTags, groupNames)
%PLOTSYSTEMATICFIGURE  Pooled laterality per group, with sign-rank p-values.

nMet = numel(metrics);
fig = figure('Name', sprintf('L-R systematic bias (%s)', scheme), ...
    'Position', [100 100 760 440]);
ax = axes(fig);
hold(ax, 'on');

colors = [0.85 0.33 0.10; 0.20 0.30 0.70];
xBase = 1:nMet;
offset = 0.35;

for g = 1:numel(groupTags)
    means = NaN(nMet, 1);
    sds = NaN(nMet, 1);
    for m = 1:nMet
        key = metrics(m).key;
        if isfield(systematic, key) && isfield(systematic.(key), groupTags{g})
            means(m) = systematic.(key).(groupTags{g}).meanLI;
            sds(m)   = systematic.(key).(groupTags{g}).sdLI;
        end
    end
    xPos = xBase + (g - 1.5) * offset;
    errorbar(ax, xPos, means, sds, 'o', 'Color', colors(g, :), ...
        'MarkerFaceColor', colors(g, :), 'LineWidth', 1.2, 'CapSize', 8, ...
        'DisplayName', strrep(groupNames{g}, '_', '\_'));

    for m = 1:nMet
        key = metrics(m).key;
        if ~isfield(systematic, key) || ~isfield(systematic.(key), groupTags{g})
            continue;
        end
        pv = systematic.(key).(groupTags{g}).signrankP;
        if isfinite(pv) && isfinite(means(m))
            text(ax, xPos(m), means(m) + abs(sds(m)) + 0.02, sprintf('p=%.3f', pv), ...
                'FontSize', 7, 'HorizontalAlignment', 'center');
        end
    end
end

yline(ax, 0, '--', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
set(ax, 'XTick', xBase, 'XTickLabel', {metrics.label}, 'TickLabelInterpreter', 'tex');
xlim(ax, [0.5, nMet + 0.5]);
ylabel(ax, 'Mean LI_{norm}  (all pairs \times subjects)', 'Interpreter', 'tex');
title(ax, sprintf('Systematic L-R bias by group  --  %s', scheme), 'Interpreter', 'tex');
legend(ax, 'Location', 'best');
grid(ax, 'on');
end

%% ------------------------------------------------------------------
function side = inferDominantSide(meanLI, pValue, alpha)
%INFERDOMINANTSIDE  Which hemisphere dominates, or 'none' when not significant.
if ~isfinite(meanLI) || ~isfinite(pValue) || pValue >= alpha
    side = 'none';
elseif meanLI > 0
    side = 'L';
else
    side = 'R';
end
end

%% ------------------------------------------------------------------
function T = buildLateralityOutlierTable(pairLabels, pairLeftIdx, pairRightIdx, ...
    metrics, LI, ctrlMean, ctrlStd, subjectCol, stats, basis)
%BUILDLATERALITYOUTLIERTABLE  One row per flagged pair, both bases reported.
% The z and p columns are named after the basis that produced them, while
% the raw LI_norm and LI_signed values are always emitted, so a reader can
% see the effect in bounded and in native units regardless of the setting.

outIdx = find(any(stats.isOutlier, 2));
if isempty(outIdx)
    T = [];
    return;
end

if strcmp(basis, 'signed')
    zPrefix = 'z_signed';
    pPrefix = 'p_bonf_signed';
else
    zPrefix = 'z';
    pPrefix = 'p_bonf';
end

T = table(outIdx, string(pairLabels(outIdx)), ...
    pairLeftIdx(outIdx), pairRightIdx(outIdx), ...
    'VariableNames', {'PairIdx', 'Region', 'L_NodeIdx', 'R_NodeIdx'});

for m = 1:numel(metrics)
    key = metrics(m).key;
    T.(['LI_norm_case_' key])            = LI.norm{m}(outIdx, subjectCol);
    T.(['LI_norm_control_mean_' key])    = ctrlMean.norm(outIdx, m);
    T.(['LI_norm_control_sd_' key])      = ctrlStd.norm(outIdx, m);
    T.(['LI_signed_case_' key])          = LI.signed{m}(outIdx, subjectCol);
    T.(['LI_signed_control_mean_' key])  = ctrlMean.signed(outIdx, m);
    T.(['LI_signed_control_sd_' key])    = ctrlStd.signed(outIdx, m);
    T.([zPrefix '_' key])                = stats.z(outIdx, m);
    T.([pPrefix '_' key])                = stats.pBonf(outIdx, m);
end

minP = min(stats.pBonf(outIdx, :), [], 2, 'omitnan');
[~, order] = sort(minP, 'ascend');
T = T(order, :);
end
