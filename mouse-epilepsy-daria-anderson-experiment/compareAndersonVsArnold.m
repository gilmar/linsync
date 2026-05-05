function summary = compareAndersonVsArnold(varargin)
%COMPAREANDERSONVSARNOLD  Per-node sanity check: compare D(->i) and
% D(k->) for every Anderson mouse against the per-node average across
% the Arnold cohort, to flag regions whose centrality changes
% dramatically between the two strains.
%
% Reads the per-mouse results that runAllMiceSection45 writes under
% results/ (section45_<mouseId>_<scheme>_results.mat). For each Anderson
% mouse, produces a 2x1 figure (D(->i) on top, D(k->) on bottom):
%
%   * shaded band  -- Arnold mean +- SD per node
%   * blue line    -- Arnold mean
%   * grey dots    -- individual Arnold values per node
%   * red markers  -- the Anderson mouse's per-node value
%   * black ring + label  -- nodes with |z| >= ZThreshold
%
% Outputs (under results/):
%   compare_anderson_vs_arnold_<mouseId>_<scheme>.{fig,png}
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

setupMousePaths();

p = inputParser;
addParameter(p, 'Normalisation', 'parkes', @(s) ischar(s) || isstring(s));
addParameter(p, 'ZThreshold',    2.0,  @isscalar);
addParameter(p, 'Alpha',         0.05, @(x) isscalar(x) && x > 0 && x < 1);
addParameter(p, 'SaveResults',   true, @islogical);
parse(p, varargin{:});
opts = p.Results;
scheme = char(opts.Normalisation);

resultsDir = setupMousePaths();

%% Discover and load per-mouse result files
pattern = sprintf('section45_*_%s_results.mat', scheme);
files = dir(fullfile(resultsDir, pattern));
if isempty(files)
    error('compareAndersonVsArnold:NoResults', ...
        'No files matching %s found under %s. Run runAllMiceSection45 first.', ...
        pattern, resultsDir);
end

mice    = cell(0, 1);
loaded  = cell(0, 1);
for k = 1:numel(files)
    s = load(fullfile(resultsDir, files(k).name));
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

isAnderson = startsWith(string(mice), 'Anderson');
isArnold   = startsWith(string(mice), 'Arnold');

if ~any(isAnderson)
    error('compareAndersonVsArnold:NoAnderson', ...
        'No Anderson mice found in %s (mice: %s).', resultsDir, strjoin(mice, ', '));
end
if ~any(isArnold)
    error('compareAndersonVsArnold:NoArnold', ...
        'No Arnold mice found in %s (mice: %s).', resultsDir, strjoin(mice, ', '));
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

D_susc_arnold_mean = mean(D_susc_arnold, 2, 'omitnan');
D_susc_arnold_std  = std(D_susc_arnold,  0, 2, 'omitnan');
D_infl_arnold_mean = mean(D_infl_arnold, 2, 'omitnan');
D_infl_arnold_std  = std(D_infl_arnold,  0, 2, 'omitnan');

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

%% One figure per Anderson mouse
deviations = struct();
for a = 1:numel(andersonNames)
    aname = andersonNames{a};
    aSusc = D_susc_anderson(:, a);
    aInfl = D_infl_anderson(:, a);

    zSusc = safeZ(aSusc, D_susc_arnold_mean, D_susc_arnold_std);
    zInfl = safeZ(aInfl, D_infl_arnold_mean, D_infl_arnold_std);

    % Bonferroni correction: one two-tailed test per non-trivial node per metric
    nTests = sum(~trivialAcrossArnold & ~isnan(zSusc)) + ...
             sum(~trivialAcrossArnold & ~isnan(zInfl));
    if nTests == 0; nTests = 1; end
    pSusc      = 2 * (1 - normcdf(abs(zSusc)));
    pInfl      = 2 * (1 - normcdf(abs(zInfl)));
    pSusc_bonf = min(1, pSusc * nTests);
    pInfl_bonf = min(1, pInfl * nTests);
    % Corrected z-threshold for display (equivalent to Bonferroni at Alpha)
    zThresh_bonf = norminv(1 - opts.Alpha / (2 * nTests));
    outlierMaskSusc = ~trivialAcrossArnold & ~isnan(pSusc_bonf) & pSusc_bonf < opts.Alpha;
    outlierMaskInfl = ~trivialAcrossArnold & ~isnan(pInfl_bonf) & pInfl_bonf < opts.Alpha;

    fig = figure('Name', sprintf('%s vs Arnold per-node centralities', aname), ...
                 'Position', [60 60 1700 820]);

    plotComparisonPanel(1, aSusc, D_susc_arnold, ...
        D_susc_arnold_mean, D_susc_arnold_std, zSusc, ...
        labels, trivialAcrossArnold, outlierMaskSusc, ...
        sprintf('D(\\rightarrow i): %s vs Arnold mean \\pm SD (n=%d)  [Bonferroni |z|\\geq%.2f, \\alpha=%.3f]', ...
            strrep(aname, '_', '\_'), nArnold, zThresh_bonf, opts.Alpha), ...
        'D(\rightarrow i)');

    plotComparisonPanel(2, aInfl, D_infl_arnold, ...
        D_infl_arnold_mean, D_infl_arnold_std, zInfl, ...
        labels, trivialAcrossArnold, outlierMaskInfl, ...
        sprintf('D(k \\rightarrow): %s vs Arnold mean \\pm SD (n=%d)  [Bonferroni |z|\\geq%.2f, \\alpha=%.3f]', ...
            strrep(aname, '_', '\_'), nArnold, zThresh_bonf, opts.Alpha), ...
        'D(k \rightarrow)');

    sgtitle(sprintf(['Per-node stability centralities: %s vs Arnold ' ...
        'cohort  --  normalisation = %s'], strrep(aname, '_', '\_'), scheme), ...
        'Interpreter', 'tex');

    deviations.(matlab.lang.makeValidName(aname)) = struct( ...
        'D_susc_anderson', aSusc,      'D_infl_anderson', aInfl, ...
        'z_susc',          zSusc,      'z_infl',          zInfl, ...
        'p_susc',          pSusc,      'p_infl',          pInfl, ...
        'p_bonf_susc',     pSusc_bonf, 'p_bonf_infl',     pInfl_bonf, ...
        'nTests',          nTests,     'zThresh_bonf',    zThresh_bonf);

    if opts.SaveResults
        baseName = sprintf('compare_anderson_vs_arnold_%s_%s', aname, scheme);
        savefig(fig, fullfile(resultsDir, [baseName '.fig']));
        try
            exportgraphics(fig, fullfile(resultsDir, [baseName '.png']), 'Resolution', 200);
        catch
            saveas(fig, fullfile(resultsDir, [baseName '.png']));
        end

        outMask = outlierMaskSusc | outlierMaskInfl;
        outIdx = find(outMask);
        if ~isempty(outIdx)
            T = table(outIdx, string(labels(outIdx)), ...
                aSusc(outIdx), D_susc_arnold_mean(outIdx), ...
                D_susc_arnold_std(outIdx), zSusc(outIdx), ...
                pSusc(outIdx), pSusc_bonf(outIdx), ...
                aInfl(outIdx), D_infl_arnold_mean(outIdx), ...
                D_infl_arnold_std(outIdx), zInfl(outIdx), ...
                pInfl(outIdx), pInfl_bonf(outIdx), ...
                'VariableNames', {'NodeIdx','Region', ...
                    'D_susc_anderson','D_susc_arnold_mean','D_susc_arnold_std', ...
                    'z_susc','p_susc','p_bonf_susc', ...
                    'D_infl_anderson','D_infl_arnold_mean','D_infl_arnold_std', ...
                    'z_infl','p_infl','p_bonf_infl'});
            % Sort by smallest Bonferroni p-value in either centrality
            minBonfP = min(T.p_bonf_susc, T.p_bonf_infl);
            [~, ord] = sort(minBonfP, 'ascend');
            T = T(ord, :);
            writetable(T, fullfile(resultsDir, [baseName '_outliers.csv']));
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
summary.trivialAcrossArnold = trivialAcrossArnold;
summary.deviations          = deviations;
summary.zThreshold          = opts.ZThreshold;
summary.alpha               = opts.Alpha;

if opts.SaveResults
    save(fullfile(resultsDir, ...
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
function z = safeZ(x, mu, sigma)
%SAFEZ Per-node z-score with NaN where sigma == 0 or sigma is NaN
% (i.e. the cohort contains a single non-NaN value or all-NaN node).
z = NaN(size(x));
ok = ~isnan(sigma) & sigma > 0 & ~isnan(mu) & ~isnan(x);
z(ok) = (x(ok) - mu(ok)) ./ sigma(ok);
end

%% ------------------------------------------------------------------
function plotComparisonPanel(panelIdx, anderVals, arnoldVals, ...
    arnMean, arnStd, zScore, labels, trivial, outlierMask, panelTitle, ylab)
%PLOTCOMPARISONPANEL Render one centrality comparison panel into the
% panelIdx-th subplot of a 2x1 layout.
% outlierMask: logical N-vector, pre-computed from Bonferroni correction.

ax = subplot(2, 1, panelIdx); hold(ax, 'on');
N = numel(anderVals);
x = 1:N;

% Greyed background bands for trivial nodes
trivialIdx = find(trivial);
yl = [min([anderVals; arnMean - arnStd; arnoldVals(:)], [], 'omitnan') ...
      max([anderVals; arnMean + arnStd; arnoldVals(:)], [], 'omitnan')];
if ~all(isfinite(yl))
    yl = [0 1];
end
yl = yl + [-1 1] * 0.05 * max(diff(yl), eps);
for k = 1:numel(trivialIdx)
    patch(ax, trivialIdx(k) + 0.5 * [-1 1 1 -1], yl([1 1 2 2]), ...
        [0.94 0.94 0.94], 'EdgeColor', 'none', 'HandleVisibility', 'off');
end

% Arnold mean +- SD shaded band (only over connected nodes for clarity)
valid = ~isnan(arnMean) & ~isnan(arnStd) & ~trivial;
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

% Highlight Bonferroni-significant outliers
outlierIdx = find(outlierMask & ~trivial & ~isnan(zScore));
if ~isempty(outlierIdx)
    plot(ax, x(outlierIdx), anderVals(outlierIdx), 'o', ...
        'Color', 'k', 'MarkerSize', 10, 'LineWidth', 1.5, ...
        'HandleVisibility', 'off');
    for k = 1:numel(outlierIdx)
        i = outlierIdx(k);
        text(ax, i, anderVals(i), ...
            sprintf('  %s (z=%+.1f)', labels{i}, zScore(i)), ...
            'Interpreter', 'none', 'FontSize', 7, 'Rotation', 45, ...
            'VerticalAlignment', 'bottom');
    end
end

xlim(ax, [0.5 N + 0.5]);
ylim(ax, yl);
set(ax, 'XTick', 1:N, 'XTickLabel', labels, ...
    'TickLabelInterpreter', 'none', 'FontSize', 6);
xtickangle(ax, 75);
grid(ax, 'on'); box(ax, 'on');
ylabel(ax, ylab, 'Interpreter', 'tex');
title(ax, panelTitle, 'Interpreter', 'tex');
legend(ax, 'Location', 'best');
end
