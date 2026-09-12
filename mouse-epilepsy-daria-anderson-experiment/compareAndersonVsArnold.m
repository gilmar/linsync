function summary = compareAndersonVsArnold(varargin)
%COMPAREANDERSONVSARNOLD  Per-node sanity check: compare D(->i),
% D(k->), and Betweenness Centrality (BC) for every Anderson mouse
% against the per-node average across the Arnold cohort, to flag
% regions whose centrality changes dramatically between the two strains.
%
% BC is included alongside the stability centralities because, on the
% Parkes-normalised symmetric mouse connectomes, the other classical
% centralities (PR, KZ, EC, SelfC, in/out-DC) collapse onto D(->i)
% (Spearman rho >= 0.96), while BC carries genuinely independent
% information (rho ~ 0.78). See results/centrality_comparison_analysis.md.
%
% Reads the per-mouse results that runAllMiceStabilityCentralities writes under
% results/ (stabilityCentralities_<mouseId>_<scheme>_results.mat). For each Anderson
% mouse, produces one figure per metric (D(->i), D(k->), BC when available):
%
%   * shaded band  -- Arnold mean +- SD per node
%   * blue line    -- Arnold mean
%   * grey dots    -- individual Arnold values per node
%   * red markers  -- the Anderson mouse's per-node value
%   * black ring + label  -- nodes with |z| >= ZThreshold
%
% Outputs (under results/<experiment>/anderson_vs_arnold/):
%   compare_anderson_vs_arnold_<mouseId>_<scheme>_D_to_i.{fig,png}
%   compare_anderson_vs_arnold_<mouseId>_<scheme>_D_k_to.{fig,png}
%   compare_anderson_vs_arnold_<mouseId>_<scheme>_BC.{fig,png}  (when BC present)
%   compare_anderson_vs_arnold_<mouseId>_<scheme>_outliers.csv
%   compare_anderson_vs_arnold_<scheme>.mat
%
% Usage:
%   compareAndersonVsArnold;
%   compareAndersonVsArnold('Normalisation', 'tvb', 'ZThreshold', 2.0);
%
% Name-value options:
%   Normalisation : 'tvb' (default)  -- which scheme to load
%   ZThreshold    : 2.0              -- |z| above which a node is flagged
%                                       (kept for display; flagging uses Alpha)
%   Alpha         : 0.05             -- family-wise significance level for
%                                       Bonferroni correction across all
%                                       non-trivial nodes x metrics per mouse
%   SaveResults   : true             -- write figures + CSVs
%   CasePrefix    : 'Anderson'        -- mouseId prefix for case cohort
%   ControlPrefix : 'Arnold'         -- mouseId prefix for control cohort

setupMousePaths();

p = inputParser;
addParameter(p, 'Normalisation', 'column', @(s) ischar(s) || isstring(s));
addParameter(p, 'ZThreshold',    2.0,  @isscalar);
addParameter(p, 'Alpha',         0.05, @(x) isscalar(x) && x > 0 && x < 1);
addParameter(p, 'SaveResults',   true, @islogical);
addParameter(p, 'ResultsDir',    '',   @(s) ischar(s) || isstring(s));
addParameter(p, 'RunParameters', [],   @(x) isempty(x) || isstruct(x));
addParameter(p, 'CasePrefix',    'Anderson', @(s) ischar(s) || isstring(s));
addParameter(p, 'ControlPrefix', 'Arnold',   @(s) ischar(s) || isstring(s));
parse(p, varargin{:});
opts = p.Results;
scheme = char(opts.Normalisation);
casePrefix = char(opts.CasePrefix);
controlPrefix = char(opts.ControlPrefix);

resultsDir = resolveMouseResultsDir(opts.ResultsDir);
outDir = mouseResultsDir(resultsDir, 'anderson_vs_arnold');

%% Discover and load per-mouse result files
files = listPerMouseResultFiles(resultsDir, scheme);
if isempty(files)
    error('compareAndersonVsArnold:NoResults', ...
        'No stabilityCentralities_*_%s_results.mat found under %s. Run runAllMiceStabilityCentralities first.', ...
        scheme, resultsDir);
end

mice    = cell(0, 1);
loaded  = cell(0, 1);
for k = 1:numel(files)
    s = load(fullfile(files(k).folder, files(k).name));
    if ~isfield(s, 'results'); continue; end
    r = s.results;
    if ~ischar(r.mouseId) && ~isstring(r.mouseId); continue; end
    mice{end+1, 1}   = char(r.mouseId); %#ok<AGROW>
    loaded{end+1, 1} = r;                %#ok<AGROW>
end
if isempty(mice)
    error('compareAndersonVsArnold:EmptyResults', ...
        'Loaded files contained no usable results structs.');
end

[isCase, isControl] = mouseCohortGroupMask(mice, casePrefix, controlPrefix);
isAnderson = isCase;
isArnold   = isControl;

if ~any(isCase)
    error('compareAndersonVsArnold:NoCaseMice', ...
        'No mice with case prefix "%s" found in %s (mice: %s).', ...
        casePrefix, resultsDir, strjoin(mice, ', '));
end
if ~any(isControl)
    error('compareAndersonVsArnold:NoControlMice', ...
        'No mice with control prefix "%s" found in %s (mice: %s).', ...
        controlPrefix, resultsDir, strjoin(mice, ', '));
end

labels = loaded{1}.labels;
N      = numel(labels);

% Cross-check that all mice share the same parcellation
for k = 2:numel(loaded)
    if numel(loaded{k}.labels) ~= N || ...
       ~all(string(loaded{k}.labels) == string(labels))
        warning('compareAndersonVsArnold:LabelMismatch', ...
            'Mouse %s labels differ from %s -- assuming index alignment.', ...
            mice{k}, mice{1});
    end
end

%% Pack per-mouse vectors into NxM matrices
D_susc_anderson = packField(loaded, isAnderson, 'D_susceptibility', N);
D_infl_anderson = packField(loaded, isAnderson, 'D_influence',      N);
D_susc_arnold   = packField(loaded, isArnold,   'D_susceptibility', N);
D_infl_arnold   = packField(loaded, isArnold,   'D_influence',      N);

% Betweenness centrality from the classical-centralities struct on each
% per-mouse result. Older results (pre-centralities) get NaN and BC is
% silently skipped.
BC_anderson = packCentrality(loaded, isAnderson, 'betweenness', N);
BC_arnold   = packCentrality(loaded, isArnold,   'betweenness', N);
hasBC = any(~isnan(BC_arnold(:))) && any(~isnan(BC_anderson(:)));
if ~hasBC
    warning('compareAndersonVsArnold:NoBC', ...
        ['Betweenness centrality not available on any loaded result ' ...
         '(centralities.betweenness missing or all-NaN). BC figure will ' ...
         'be skipped. Re-run runAllMiceStabilityCentralities with BCT on the path.']);
end

D_susc_arnold_mean = mean(D_susc_arnold, 2, 'omitnan');
D_susc_arnold_std  = std(D_susc_arnold,  0, 2, 'omitnan');
D_infl_arnold_mean = mean(D_infl_arnold, 2, 'omitnan');
D_infl_arnold_std  = std(D_infl_arnold,  0, 2, 'omitnan');
BC_arnold_mean     = mean(BC_arnold,     2, 'omitnan');
BC_arnold_std      = std(BC_arnold,      0, 2, 'omitnan');

% A node is "trivial" for the comparison if it is trivial in any Arnold
% run (Arnold mean would be uninformative there).
trivialAcrossArnold = false(N, 1);
arnoldIdx = find(isArnold);
for ii = 1:numel(arnoldIdx)
    trivialAcrossArnold = trivialAcrossArnold | loaded{arnoldIdx(ii)}.isTrivial(:);
end

andersonNames = mice(isAnderson);
arnoldNames   = mice(isArnold);
nArnold       = numel(arnoldNames);

fprintf('compareAndersonVsArnold: %d Anderson mouse(/mice), %d Arnold mouse(/mice).\n', ...
    numel(andersonNames), nArnold);
fprintf('  Anderson: %s\n', strjoin(andersonNames, ', '));
fprintf('  Arnold  : %s\n', strjoin(arnoldNames,   ', '));

%% One figure per Anderson mouse x metric
deviations = struct();
for a = 1:numel(andersonNames)
    aname = andersonNames{a};
    aSusc = D_susc_anderson(:, a);
    aInfl = D_infl_anderson(:, a);
    aBC   = BC_anderson(:, a);

    zSusc = safeZ(aSusc, D_susc_arnold_mean, D_susc_arnold_std);
    zInfl = safeZ(aInfl, D_infl_arnold_mean, D_infl_arnold_std);
    zBC   = safeZ(aBC,   BC_arnold_mean,     BC_arnold_std);

    % Bonferroni correction: one two-tailed test per non-trivial node per
    % included metric. BC tests are added to the family only when BC data
    % is present, so the threshold tightens accordingly.
    nTests = sum(~trivialAcrossArnold & ~isnan(zSusc)) + ...
             sum(~trivialAcrossArnold & ~isnan(zInfl));
    if hasBC
        nTests = nTests + sum(~trivialAcrossArnold & ~isnan(zBC));
    end
    if nTests == 0; nTests = 1; end
    pSusc      = 2 * (1 - normcdf(abs(zSusc)));
    pInfl      = 2 * (1 - normcdf(abs(zInfl)));
    pBC        = 2 * (1 - normcdf(abs(zBC)));
    pSusc_bonf = min(1, pSusc * nTests);
    pInfl_bonf = min(1, pInfl * nTests);
    pBC_bonf   = min(1, pBC   * nTests);
    % Corrected z-threshold for display (equivalent to Bonferroni at Alpha)
    zThresh_bonf = norminv(1 - opts.Alpha / (2 * nTests));
    outlierMaskSusc = ~trivialAcrossArnold & ~isnan(pSusc_bonf) & pSusc_bonf < opts.Alpha;
    outlierMaskInfl = ~trivialAcrossArnold & ~isnan(pInfl_bonf) & pInfl_bonf < opts.Alpha;
    outlierMaskBC   = ~trivialAcrossArnold & ~isnan(pBC_bonf)   & pBC_bonf   < opts.Alpha;

    anameTex = strrep(aname, '_', '\_');
    panelTitleFmt = '%s: %s vs Arnold mean \\pm SD (n=%d)  [Bonferroni |z|\\geq%.2f, \\alpha=%.3f]';
    cohortTitle = sprintf('Per-node centralities: %s vs Arnold cohort  --  normalisation = %s', ...
        anameTex, scheme);

    metricPanels = { ...
        struct('suffix', 'D_to_i', 'short', 'D(\rightarrow i)', ...
            'vals', aSusc, 'arnold', D_susc_arnold, 'mean', D_susc_arnold_mean, ...
            'std', D_susc_arnold_std, 'z', zSusc, 'outliers', outlierMaskSusc, ...
            'title', sprintf(panelTitleFmt, 'D(\rightarrow i)', anameTex, nArnold, zThresh_bonf, opts.Alpha)); ...
        struct('suffix', 'D_k_to', 'short', 'D(k \rightarrow)', ...
            'vals', aInfl, 'arnold', D_infl_arnold, 'mean', D_infl_arnold_mean, ...
            'std', D_infl_arnold_std, 'z', zInfl, 'outliers', outlierMaskInfl, ...
            'title', sprintf(panelTitleFmt, 'D(k \rightarrow)', anameTex, nArnold, zThresh_bonf, opts.Alpha)) ...
        };
    if hasBC
        metricPanels{end+1} = struct('suffix', 'BC', 'short', 'BC', ...
            'vals', aBC, 'arnold', BC_arnold, 'mean', BC_arnold_mean, ...
            'std', BC_arnold_std, 'z', zBC, 'outliers', outlierMaskBC, ...
            'title', sprintf(panelTitleFmt, 'BC', anameTex, nArnold, zThresh_bonf, opts.Alpha));
    end

    baseName = sprintf('compare_anderson_vs_arnold_%s_%s', aname, scheme);
    for p = 1:numel(metricPanels)
        mp = metricPanels{p};
        fig = plotComparisonFigure(mp.vals, mp.arnold, mp.mean, mp.std, mp.z, ...
            labels, trivialAcrossArnold, mp.outliers, mp.title, mp.short, cohortTitle, ...
            sprintf('%s vs Arnold — %s', aname, mp.short));
        if opts.SaveResults
            saveComparisonFigure(fig, outDir, sprintf('%s_%s', baseName, mp.suffix));
        else
            close(fig);
        end
    end

    devStruct = struct( ...
        'D_susc_anderson', aSusc,      'D_infl_anderson', aInfl, ...
        'BC_anderson',     aBC, ...
        'z_susc',          zSusc,      'z_infl',          zInfl,    'z_BC',        zBC, ...
        'p_susc',          pSusc,      'p_infl',          pInfl,    'p_BC',        pBC, ...
        'p_bonf_susc',     pSusc_bonf, 'p_bonf_infl',     pInfl_bonf,'p_bonf_BC',  pBC_bonf, ...
        'nTests',          nTests,     'zThresh_bonf',    zThresh_bonf, ...
        'hasBC',           hasBC);
    deviations.(matlab.lang.makeValidName(aname)) = devStruct;

    if opts.SaveResults
        outMask = outlierMaskSusc | outlierMaskInfl | outlierMaskBC;
        outIdx = find(outMask);
        if ~isempty(outIdx)
            T = table(outIdx, string(labels(outIdx)), ...
                aSusc(outIdx), D_susc_arnold_mean(outIdx), ...
                D_susc_arnold_std(outIdx), zSusc(outIdx), ...
                pSusc(outIdx), pSusc_bonf(outIdx), ...
                aInfl(outIdx), D_infl_arnold_mean(outIdx), ...
                D_infl_arnold_std(outIdx), zInfl(outIdx), ...
                pInfl(outIdx), pInfl_bonf(outIdx), ...
                aBC(outIdx),   BC_arnold_mean(outIdx), ...
                BC_arnold_std(outIdx), zBC(outIdx), ...
                pBC(outIdx),   pBC_bonf(outIdx), ...
                'VariableNames', {'NodeIdx','Region', ...
                    'D_susc_anderson','D_susc_arnold_mean','D_susc_arnold_std', ...
                    'z_susc','p_susc','p_bonf_susc', ...
                    'D_infl_anderson','D_infl_arnold_mean','D_infl_arnold_std', ...
                    'z_infl','p_infl','p_bonf_infl', ...
                    'BC_anderson','BC_arnold_mean','BC_arnold_std', ...
                    'z_BC','p_BC','p_bonf_BC'});
            % Sort by smallest Bonferroni p-value across the three metrics
            minBonfP = min([T.p_bonf_susc, T.p_bonf_infl, T.p_bonf_BC], [], 2);
            [~, ord] = sort(minBonfP, 'ascend');
            T = T(ord, :);
            writetable(T, fullfile(outDir, [baseName '_outliers.csv']));
            fprintf('  %s: %d outlier node(s) (Bonferroni p<%.3f, |z|>=%.2f, m=%d) -> %s\n', ...
                aname, numel(outIdx), opts.Alpha, zThresh_bonf, nTests, [baseName '_outliers.csv']);
        else
            fprintf('  %s: no outlier nodes at Bonferroni p<%.3f (|z|>=%.2f, m=%d).\n', ...
                aname, opts.Alpha, zThresh_bonf, nTests);
        end
    end
end

%% Pack overall summary
summary = struct();
summary.scheme              = scheme;
summary.mice                = mice;
summary.casePrefix          = casePrefix;
summary.controlPrefix       = controlPrefix;
summary.isCase              = isCase;
summary.isControl           = isControl;
summary.isAnderson          = isAnderson;
summary.isArnold            = isArnold;
summary.labels              = labels;
summary.D_susc_anderson     = D_susc_anderson;
summary.D_infl_anderson     = D_infl_anderson;
summary.D_susc_arnold       = D_susc_arnold;
summary.D_infl_arnold       = D_infl_arnold;
summary.D_susc_arnold_mean  = D_susc_arnold_mean;
summary.D_susc_arnold_std   = D_susc_arnold_std;
summary.D_infl_arnold_mean  = D_infl_arnold_mean;
summary.D_infl_arnold_std   = D_infl_arnold_std;
summary.BC_anderson         = BC_anderson;
summary.BC_arnold           = BC_arnold;
summary.BC_arnold_mean      = BC_arnold_mean;
summary.BC_arnold_std       = BC_arnold_std;
summary.hasBC               = hasBC;
summary.trivialAcrossArnold = trivialAcrossArnold;
summary.deviations          = deviations;
summary.zThreshold          = opts.ZThreshold;
summary.alpha               = opts.Alpha;

if isempty(opts.RunParameters)
    summary.runParameters = mouseExperimentRunParameters('buildFromCompareAnderson', ...
        opts, scheme, resultsDir);
else
    summary.runParameters = mouseExperimentRunParameters('merge', opts.RunParameters, ...
        mouseExperimentRunParameters('buildFromCompareAnderson', opts, scheme, resultsDir));
end

if opts.SaveResults
    save(fullfile(outDir, ...
        sprintf('compare_anderson_vs_arnold_%s.mat', scheme)), ...
        '-struct', 'summary');
end
end

%% ------------------------------------------------------------------
function M = packField(loaded, mask, fieldName, N)
%PACKFIELD Pack the chosen field from each selected result struct into
% an N-by-K matrix (K = nnz(mask)).
idx = find(mask);
M = NaN(N, numel(idx));
for ii = 1:numel(idx)
    v = loaded{idx(ii)}.(fieldName);
    M(:, ii) = v(:);
end
end

%% ------------------------------------------------------------------
function M = packCentrality(loaded, mask, fieldName, N)
%PACKCENTRALITY Pull a named field from loaded{i}.centralities into an
% N-by-K matrix. Mice without a centralities struct, or without the named
% field, contribute an all-NaN column so they are silently excluded by
% downstream NaN-aware aggregation.
idx = find(mask);
M = NaN(N, numel(idx));
for ii = 1:numel(idx)
    r = loaded{idx(ii)};
    if ~isfield(r, 'centralities') || ~isstruct(r.centralities)
        continue;
    end
    if ~isfield(r.centralities, fieldName)
        continue;
    end
    v = r.centralities.(fieldName);
    if numel(v) ~= N
        continue;
    end
    M(:, ii) = v(:);
end
end

%% ------------------------------------------------------------------
function z = safeZ(x, mu, sigma)
%SAFEZ Per-node z-score with NaN where sigma == 0 or sigma is NaN
% (i.e. the cohort contains a single non-NaN value or all-NaN node).
z = NaN(size(x));
ok = ~isnan(sigma) & sigma > 0 & ~isnan(mu) & ~isnan(x);
z(ok) = (x(ok) - mu(ok)) ./ sigma(ok);
end

%% ------------------------------------------------------------------
function fig = plotComparisonFigure(anderVals, arnoldVals, ...
    arnMean, arnStd, zScore, labels, trivial, outlierMask, panelTitle, ylab, ...
    sgTitleText, figName)
%PLOTCOMPARISONFIGURE Render one centrality comparison figure.
% outlierMask: logical N-vector, pre-computed from Bonferroni correction.
% Atlas placeholders (BACKGROUND, *_MASK) flagged in trivial are omitted
% from the x-axis; they are already excluded from Bonferroni testing.
% Outliers are marked on the x-axis (bold, coloured tick labels) and with
% horizontal z-score labels above the Anderson marker.

plotMask = ~trivial(:);
if any(plotMask)
    anderVals   = anderVals(plotMask);
    arnoldVals  = arnoldVals(plotMask, :);
    arnMean     = arnMean(plotMask);
    arnStd      = arnStd(plotMask);
    zScore      = zScore(plotMask);
    outlierMask = outlierMask(plotMask);
    labels = labels(plotMask);
end

figW = max(1600, min(2200, 1100 + 14 * numel(anderVals)));
figH = max(640, min(920, 500 + 5 * numel(anderVals)));
fig = figure('Name', figName, 'Position', [60 60 figW figH], 'Color', 'w');
ax = axes(fig); hold(ax, 'on');
N = numel(anderVals);
x = 1:N;

yl = [min([anderVals; arnMean - arnStd; arnoldVals(:)], [], 'omitnan') ...
      max([anderVals; arnMean + arnStd; arnoldVals(:)], [], 'omitnan')];
if ~all(isfinite(yl))
    yl = [0 1];
end
ySpan = max(diff(yl), eps);
ylPad = 0.06 * ySpan;
yl = yl + [-1 1] * ylPad;
zHeadroom = 0.12 * ySpan;

% Arnold mean +- SD shaded band
valid = ~isnan(arnMean) & ~isnan(arnStd);
if any(valid)
    xv = x(valid);
    yLo = (arnMean - arnStd); yLo = yLo(valid);
    yHi = (arnMean + arnStd); yHi = yHi(valid);
    fill(ax, [xv fliplr(xv)], [yLo' fliplr(yHi')], ...
        [0.80 0.84 0.95], 'EdgeColor', 'none', 'FaceAlpha', 0.65, ...
        'DisplayName', 'Arnold mean \pm SD');
end

% Arnold mean line
plot(ax, x, arnMean, '-', 'Color', [0.20 0.30 0.70], ...
    'LineWidth', 1.2, 'DisplayName', 'Arnold mean');

% Individual Arnold values
for j = 1:size(arnoldVals, 2)
    plot(ax, x, arnoldVals(:, j), '.', 'Color', [0.55 0.60 0.85], ...
        'MarkerSize', 6, 'HandleVisibility', 'off');
end

% Anderson values
plot(ax, x, anderVals, 'o', 'Color', [0.65 0.10 0.10], ...
    'MarkerFaceColor', [0.85 0.33 0.10], 'MarkerSize', 5, ...
    'DisplayName', 'Anderson');

outlierIdx = find(outlierMask & ~isnan(zScore));
if ~isempty(outlierIdx)
    peakY = max([anderVals(outlierIdx); ...
        arnMean(outlierIdx) + arnStd(outlierIdx)], [], 'omitnan');
    yl(2) = max(yl(2), peakY + zHeadroom);
end
xlim(ax, [0.5 N + 0.5]);
ylim(ax, yl);

zLabelDy = 0.025 * ySpan;
if ~isempty(outlierIdx)
    plot(ax, x(outlierIdx), anderVals(outlierIdx), 'o', ...
        'Color', 'k', 'MarkerSize', 10, 'LineWidth', 1.5, ...
        'HandleVisibility', 'off');
    for k = 1:numel(outlierIdx)
        i = outlierIdx(k);
        yi = anderVals(i);
        if ~isfinite(yi)
            continue;
        end
        text(ax, i, yi + zLabelDy, sprintf('z=%+.1f', zScore(i)), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', ...
            'Rotation', 0, ...
            'FontSize', 8, 'Color', [0.15 0.15 0.15], ...
            'Interpreter', 'none', 'Clipping', 'on');
    end
end

tickStyle = styleComparisonXAxis(ax, labels, outlierIdx);
grid(ax, 'on'); box(ax, 'on');
ylabel(ax, ylab, 'Interpreter', 'tex', 'FontSize', 10);
title(ax, panelTitle, 'Interpreter', 'tex', 'FontSize', 10);
lg = legend(ax, 'Location', 'northeastoutside', 'FontSize', 8, 'Box', 'off');
sg = sgtitle(fig, sgTitleText, 'Interpreter', 'tex', 'FontSize', 11);
layoutComparisonAxes(fig, ax, N, lg, sg, tickStyle);
finalizeComparisonXLabel(fig, ax);
end

%% ------------------------------------------------------------------
function tickStyle = styleComparisonXAxis(ax, labels, outlierIdx)
%STYLECOMPARISONXAXIS  Region names on x-axis; highlight outlier ticks.
% Uses hand-placed text labels so FontWeight/Color work on all MATLAB versions.
% Returns tickStyle for layout and x-axis title placement.
N = numel(labels);
tickLabs = cell(N, 1);
for i = 1:N
    tickLabs{i} = labelToChar(labels, i);
end
outlierMask = false(N, 1);
outlierMask(outlierIdx) = true;

if N > 40
    tickAngle = 90;
    tickFont = 6;
else
    tickAngle = 45;
    tickFont = 7;
end

highlightColor = [0.78 0.10 0.05];
normalColor    = [0.32 0.32 0.32];
yl = ax.YLim;

xticks(ax, 1:N);
xticklabels(ax, repmat({''}, N, 1));
ax.TickLabelInterpreter = 'none';
if isprop(ax, 'XAxis')
    ax.XAxis.TickLabelGapMultiplier = 1.8;
end

for i = 1:N
    if outlierMask(i)
        fw = 'bold';
        fc = highlightColor;
    else
        fw = 'normal';
        fc = normalColor;
    end
    text(ax, i, yl(1), tickLabs{i}, ...
        'Rotation', tickAngle, ...
        'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'top', ...
        'FontSize', tickFont, ...
        'FontWeight', fw, ...
        'Color', fc, ...
        'Interpreter', 'none', ...
        'Clipping', 'off', ...
        'Tag', 'comparisonXTickLabel');
end

tickStyle = struct('tickAngle', tickAngle, 'tickFont', tickFont, 'nRegions', N);
end

%% ------------------------------------------------------------------
function s = labelToChar(labels, idx)
%LABELTOCHAR  Extract one region label as char (cell / string / char vector).
if iscell(labels)
    s = char(labels{idx});
elseif isstring(labels)
    s = char(labels(idx));
else
    s = char(labels(idx));
end
end

%% ------------------------------------------------------------------
function layoutComparisonAxes(fig, ax, nRegions, legendObj, sgTitleObj, tickStyle)
%LAYOUTCOMPARISONAXES  Resize axes so tick labels, titles, and legend clear the data.
if nargin < 4
    legendObj = [];
end
if nargin < 5
    sgTitleObj = [];
end
if nargin < 6
    tickStyle = struct('tickAngle', 90, 'nRegions', nRegions);
end
drawnow;
ax.Units = 'normalized';
ti = ax.TightInset;
left   = max(ti(1), 0.08);
bottom = max(ti(2), min(0.42, 0.14 + 0.004 * nRegions));
if tickStyle.tickAngle == 90
    customPad = min(0.16, 0.04 + 0.0020 * tickStyle.nRegions);
else
    customPad = min(0.12, 0.03 + 0.0014 * tickStyle.nRegions);
end
bottom = max(bottom, ti(2) + customPad + 0.038);  % custom ticks + "Region" caption
top    = max(ti(4), 0.10);
if ~isempty(sgTitleObj) && isgraphics(sgTitleObj)
    top = top + 0.05;
end
right = max(ti(3), 0.06);
if ~isempty(legendObj) && isgraphics(legendObj)
    legendObj.Location = 'northeastoutside';
    drawnow;
    ti = ax.TightInset;
    right = max(right, ti(3) + 0.12);
end
w = max(0.45, 1 - left - right);
h = max(0.30, 1 - bottom - top);
ax.Position = [left, bottom, w, h];
drawnow;
end

%% ------------------------------------------------------------------
function finalizeComparisonXLabel(fig, ax)
%FINALIZECOMPARISONXLABEL  Place "Region" below custom tick labels (figure coords).
% Axes xlabel overlaps hand-drawn tick text; measure label extent and use a
% figure-level textbox centred under the lowest tick label.
drawnow;
fig.Units = 'normalized';
ax.Units = 'normalized';

if isprop(ax, 'XLabel')
    ax.XLabel.String = '';
end
delete(findobj(fig, 'Tag', 'comparisonXLabel'));

gap = 0.016;
labelH = 0.026;
minYFloor = 0.012;
minTickY = lowestComparisonTickLabelY(fig, ax);
yPos = minTickY - gap - labelH;

if yPos < minYFloor
    shift = minYFloor - yPos;
    pos = ax.Position;
    ax.Position = [pos(1), pos(2) + shift, pos(3), max(pos(4) - 0.5 * shift, 0.28)];
    drawnow;
    minTickY = lowestComparisonTickLabelY(fig, ax);
    yPos = max(minYFloor, minTickY - gap - labelH);
end

ap = ax.Position;
annotation(fig, 'textbox', [ap(1), yPos, ap(3), labelH], ...
    'String', 'Region', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'FontSize', 10, ...
    'Interpreter', 'none', ...
    'Tag', 'comparisonXLabel');
end

%% ------------------------------------------------------------------
function yMin = lowestComparisonTickLabelY(fig, ax)
%LOWESTCOMPARISONTICKLABELY  Lowest figure-normalized Y of custom tick labels.
yMin = ax.Position(2);
ht = findobj(ax, 'Tag', 'comparisonXTickLabel', 'Type', 'text');
if isempty(ht)
    return;
end
fig.Units = 'normalized';
ax.Units = 'normalized';
ap = ax.Position;
yLim = ax.YLim;
ySpan = diff(yLim);
if ySpan <= 0
    return;
end
for k = 1:numel(ht)
    ext = ht(k).Extent;
    cornersY = [ext(2), ext(2) + ext(4)];
    for yd = cornersY
        yNorm = (yd - yLim(1)) / ySpan;
        yFig = ap(2) + yNorm * ap(4);
        yMin = min(yMin, yFig);
    end
end
end

%% ------------------------------------------------------------------
function saveComparisonFigure(fig, outDir, baseName)
%SAVECOMPARISONFIGURE Write .fig and .png for one comparison figure.
pngPath = fullfile(outDir, [baseName '.png']);
savefig(fig, fullfile(outDir, [baseName '.fig']));
try
    exportgraphics(fig, pngPath, 'Resolution', 200, 'BackgroundColor', 'white');
catch
    try
        exportgraphics(fig, pngPath, 'Resolution', 200, ...
            'BackgroundColor', 'white', 'Padding', 12);
    catch
        saveas(fig, pngPath);
    end
end
close(fig);
end
