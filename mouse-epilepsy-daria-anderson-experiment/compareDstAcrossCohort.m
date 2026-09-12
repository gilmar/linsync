function summary = compareDstAcrossCohort(varargin)
%COMPAREDSTACROSSCOHORT  Tabulate and chart the network-level Deviation
% from Stability D_st = trace(Omega)/N for every mouse in the cohort.
%
% D_st is the scalar early-warning signal for instability: with C
% Parkes-normalised so that rho(C) is pinned just below 1, D_st reports
% how much the full eigenvalue spectrum (not just the leading mode)
% amplifies the stationary covariance over its uncoupled baseline.
%
% Reads stabilityCentralities_<mouseId>_<scheme>_results.mat under results/, pulls
% results.D_st_healthy from each, and produces:
%   * a printed table to the console
%   * a CSV  (D_st_cohort_<scheme>.csv)
%   * a per-mouse bar chart with strain means as horizontal lines
%     (D_st_cohort_<scheme>.{fig,png})
%
% Usage:
%   compareDstAcrossCohort;
%   compareDstAcrossCohort('Normalisation', 'tvb', 'SaveResults', false);
%
% Name-value options:
%   Normalisation : 'parkes' (default) -- which scheme to load
%   SaveResults   : true               -- write CSV and figure

setupMousePaths();

p = inputParser;
addParameter(p, 'Normalisation', 'column', @(s) ischar(s) || isstring(s));
addParameter(p, 'SaveResults',   true,     @islogical);
addParameter(p, 'ResultsDir',    '',       @(s) ischar(s) || isstring(s));
addParameter(p, 'RunParameters', [],       @(x) isempty(x) || isstruct(x));
parse(p, varargin{:});
opts   = p.Results;
scheme = char(opts.Normalisation);

resultsDir = resolveMouseResultsDir(opts.ResultsDir);
outDir = mouseResultsDir(resultsDir, 'dst_cohort');

files = listPerMouseResultFiles(resultsDir, scheme);
if isempty(files)
    error('compareDstAcrossCohort:NoResults', ...
        'No stabilityCentralities_*_%s_results.mat under %s.', scheme, resultsDir);
end

mice  = strings(0, 1);
Dst   = nan(0, 1);
rhoC  = nan(0, 1);
for k = 1:numel(files)
    s = load(fullfile(files(k).folder, files(k).name));
    if ~isfield(s, 'results'); continue; end
    r = s.results;
    if ~isfield(r, 'D_st_healthy') || ~isfield(r, 'C_healthy')
        continue;
    end
    mice(end+1, 1) = string(r.mouseId);   %#ok<AGROW>
    Dst(end+1, 1)  = r.D_st_healthy;       %#ok<AGROW>
    rhoC(end+1, 1) = max(abs(eig(r.C_healthy))); %#ok<AGROW>
end

isAnderson = startsWith(mice, "Anderson");
isArnold   = startsWith(mice, "Arnold");

[~, ord] = sort(mice);
mice = mice(ord); Dst = Dst(ord); rhoC = rhoC(ord);
isAnderson = isAnderson(ord); isArnold = isArnold(ord);

%% Per-mouse table
T = table(mice, Dst, rhoC, 'VariableNames', {'mouse', 'D_st', 'rho_C'});
disp(T);

mAnd = mean(Dst(isAnderson), 'omitnan');
sAnd = std (Dst(isAnderson), 0, 'omitnan');
mArn = mean(Dst(isArnold),   'omitnan');
sArn = std (Dst(isArnold),   0, 'omitnan');

fprintf('Anderson : n=%d  mean=%.3f  SD=%.3f  range=[%.3f, %.3f]\n', ...
    nnz(isAnderson), mAnd, sAnd, min(Dst(isAnderson)), max(Dst(isAnderson)));
fprintf('Arnold   : n=%d  mean=%.3f  SD=%.3f  range=[%.3f, %.3f]\n', ...
    nnz(isArnold),   mArn, sArn, min(Dst(isArnold)),   max(Dst(isArnold)));

%% Welch t-test if Statistics Toolbox available
pVal = NaN; tStat = NaN;
if exist('ttest2', 'file') == 2 && nnz(isAnderson) > 1 && nnz(isArnold) > 1
    [~, pVal, ~, st] = ttest2(Dst(isAnderson), Dst(isArnold), ...
        'Vartype', 'unequal');
    tStat = st.tstat;
    fprintf('Welch t-test (Anderson vs Arnold): t=%.3f, p=%.3f\n', tStat, pVal);
else
    fprintf(['Skipping Welch t-test (need n>=2 per group and ' ...
             'ttest2 from the Statistics Toolbox).\n']);
end

%% Chart
fig = figure('Name', sprintf('D_st across cohort (%s)', scheme), ...
             'Position', [80 80 880 520]);
ax = axes(fig); hold(ax, 'on');

barColors = repmat([0.55 0.65 0.85], numel(mice), 1);
barColors(isAnderson, :) = repmat([0.85 0.33 0.10], nnz(isAnderson), 1);
b = bar(ax, 1:numel(mice), Dst, 0.7, 'FaceColor', 'flat', ...
    'EdgeColor', [0.2 0.2 0.2]);
b.CData = barColors;

% Cohort mean lines
xlimX = [0.4 numel(mice) + 0.6];
plot(ax, xlimX, [mAnd mAnd], '--', 'Color', [0.85 0.33 0.10], ...
    'LineWidth', 1.2, 'DisplayName', sprintf('Anderson mean (%.3f)', mAnd));
plot(ax, xlimX, [mArn mArn], '--', 'Color', [0.20 0.30 0.70], ...
    'LineWidth', 1.2, 'DisplayName', sprintf('Arnold mean (%.3f)',   mArn));

% Value labels
for i = 1:numel(mice)
    text(ax, i, Dst(i), sprintf('  %.3f', Dst(i)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
        'FontSize', 9);
end

set(ax, 'XTick', 1:numel(mice), 'XTickLabel', mice, ...
    'TickLabelInterpreter', 'none');
xtickangle(ax, 30);
xlim(ax, xlimX);
ylim(ax, [0, max(Dst) * 1.15]);
ylabel(ax, 'D_{st} = trace(\Omega) / N');
ttlExtra = '';
if isfinite(pVal)
    ttlExtra = sprintf(' (Welch t = %.2f, p = %.3f)', tStat, pVal);
end
title(ax, sprintf(['Network-level D_{st} per mouse  --  scheme = %s' ...
    '%s'], scheme, ttlExtra), 'Interpreter', 'tex');
grid(ax, 'on'); box(ax, 'on');
legend(ax, 'Location', 'best');

if isempty(opts.RunParameters)
    runParameters = mouseExperimentRunParameters('buildFromCompareDst', ...
        opts, scheme, resultsDir);
else
    runParameters = mouseExperimentRunParameters('merge', opts.RunParameters, ...
        mouseExperimentRunParameters('buildFromCompareDst', opts, scheme, resultsDir));
end

summary = struct( ...
    'scheme',     scheme, ...
    'mice',       mice, ...
    'D_st',       Dst, ...
    'rho_C',      rhoC, ...
    'isAnderson', isAnderson, ...
    'isArnold',   isArnold, ...
    'mean_Anderson', mAnd, 'std_Anderson', sAnd, ...
    'mean_Arnold',   mArn, 'std_Arnold',   sArn, ...
    'tStat',      tStat, 'pValue', pVal, ...
    'runParameters', runParameters);

if opts.SaveResults
    csvPath = fullfile(outDir, sprintf('D_st_cohort_%s.csv', scheme));
    writetable(T, csvPath);
    figBase = fullfile(outDir, sprintf('D_st_cohort_%s', scheme));
    savefig(fig, [figBase '.fig']);
    try
        exportgraphics(fig, [figBase '.png'], 'Resolution', 200);
    catch
        saveas(fig, [figBase '.png']);
    end
    save(fullfile(outDir, sprintf('D_st_cohort_%s.mat', scheme)), ...
        'summary', 'runParameters');
    fprintf('Wrote %s, %s.{fig,png}, D_st_cohort_%s.mat\n', csvPath, figBase, scheme);
    close(fig);
end
end
