function summary = compareDstAcrossCohort(varargin)
%COMPAREDSTACROSSCOHORT  Network-level deviation from stability across a cohort.
%
%   summary = compareDstAcrossCohort('ResultsDir', resultsDir, ...
%                'Normalisation', 'column', 'CasePrefix', 'Case', ...
%                'ControlPrefix', 'Ctrl')
%
%   Tabulates and charts D_st = trace(Omega)/N for every subject in a run,
%   and tests the case group against the control group.
%
%   Where the per-node centralities say WHICH region is fragile, D_st is
%   the single number saying how fragile the network is as a whole: the
%   average stationary variance a unit noise input produces once the whole
%   eigenvalue spectrum -- not just the leading mode -- has amplified it.
%   Two networks with the same spectral radius can have very different
%   D_st, which is exactly why it is worth reporting alongside rho(C).
%
%   The group difference uses a Welch t-test (WELCHTTEST), which needs no
%   Statistics Toolbox and assumes no equality of variances.
%
%   Name-value options
%     'ResultsDir'    : ''
%     'Normalisation' : 'column'
%     'CasePrefix'    : 'Case'
%     'ControlPrefix' : 'Control'
%     'SaveResults'   : true
%     'Plot'          : true
%     'RunParameters' : []
%
%   Outputs (under <ResultsDir>/dst_cohort/)
%     D_st_cohort_<scheme>.{csv,fig,png,mat}
%
%   See also COMPUTESTABILITYCENTRALITIES, COVARIANCESGAUSSIANNET, WELCHTTEST.

p = inputParser;
addParameter(p, 'ResultsDir',    '',        @(s) ischar(s) || isstring(s));
addParameter(p, 'Normalisation', 'column',  @(s) ischar(s) || isstring(s));
addParameter(p, 'CasePrefix',    'Case',    @(s) ischar(s) || isstring(s));
addParameter(p, 'ControlPrefix', 'Control', @(s) ischar(s) || isstring(s));
addParameter(p, 'SaveResults',   true, @islogical);
addParameter(p, 'Plot',          true, @islogical);
addParameter(p, 'RunParameters', [],   @(x) isempty(x) || isstruct(x));
parse(p, varargin{:});
opts = p.Results;

scheme = normalisationSchemeInfo(opts.Normalisation).scheme;
resultsDir = resolveExperimentResultsDir(opts.ResultsDir);
outDir = experimentResultsDir(resultsDir, 'dst_cohort');

cohort = loadCohortResults(resultsDir, scheme, ...
    'CasePrefix', opts.CasePrefix, 'ControlPrefix', opts.ControlPrefix);

nSubjects = numel(cohort.subjects);
Dst = NaN(nSubjects, 1);
rhoC = NaN(nSubjects, 1);
for k = 1:nSubjects
    r = cohort.results{k};
    if isfield(r, 'D_st_healthy')
        Dst(k) = r.D_st_healthy;
    end
    if isfield(r, 'rho_healthy')
        rhoC(k) = r.rho_healthy;
    elseif isfield(r, 'C_healthy')
        rhoC(k) = max(abs(eig(r.C_healthy)));
    end
end

T = table(string(cohort.subjects), Dst, rhoC, ...
    groupTags(cohort), 'VariableNames', {'subject', 'D_st', 'rho_C', 'group'});
disp(T);

caseD = Dst(cohort.isCase);
ctrlD = Dst(cohort.isControl);
welch = welchTTest(caseD, ctrlD);

fprintf('%-10s : n=%d  mean=%.4f  SD=%.4f\n', cohort.casePrefix, ...
    numel(caseD), mean(caseD, 'omitnan'), std(caseD, 0, 'omitnan'));
fprintf('%-10s : n=%d  mean=%.4f  SD=%.4f\n', cohort.controlPrefix, ...
    numel(ctrlD), mean(ctrlD, 'omitnan'), std(ctrlD, 0, 'omitnan'));
if isfinite(welch.p)
    fprintf('Welch t-test (%s vs %s): t=%.3f, df=%.1f, p=%.4f\n', ...
        cohort.casePrefix, cohort.controlPrefix, welch.tStat, welch.df, welch.p);
else
    fprintf('Welch t-test skipped (needs at least 2 subjects in each group).\n');
end

%% ---- Chart -----------------------------------------------------------
fig = [];
if opts.Plot
    fig = plotDstFigure(cohort, Dst, welch, scheme);
end

%% ---- Summary ---------------------------------------------------------
summary = struct();
summary.scheme        = scheme;
summary.resultsDir    = resultsDir;
summary.subjects      = cohort.subjects;
summary.isCase        = cohort.isCase;
summary.isControl     = cohort.isControl;
summary.casePrefix    = cohort.casePrefix;
summary.controlPrefix = cohort.controlPrefix;
summary.D_st          = Dst;
summary.rho_C         = rhoC;
summary.table         = T;
summary.meanCase      = mean(caseD, 'omitnan');
summary.stdCase       = std(caseD, 0, 'omitnan');
summary.meanControl   = mean(ctrlD, 'omitnan');
summary.stdControl    = std(ctrlD, 0, 'omitnan');
summary.welch         = welch;

patch = experimentRunParameters('fromOptions', opts, ...
    'Analysis', 'compareDstAcrossCohort', 'ResultsDir', resultsDir, ...
    'Scheme', scheme, 'Subjects', cohort.subjects);
summary.runParameters = experimentRunParameters('merge', opts.RunParameters, patch);

if opts.SaveResults
    csvPath = fullfile(outDir, sprintf('D_st_cohort_%s.csv', scheme));
    writetable(T, csvPath);
    if ~isempty(fig) && isgraphics(fig, 'figure')
        saveFigureBoth(fig, fullfile(outDir, sprintf('D_st_cohort_%s', scheme)));
    end
    save(fullfile(outDir, sprintf('D_st_cohort_%s.mat', scheme)), 'summary');
    fprintf('Wrote %s and D_st_cohort_%s.{fig,png,mat}\n', csvPath, scheme);
end

if ~isempty(fig) && isgraphics(fig, 'figure')
    close(fig);
end
end

%% ------------------------------------------------------------------
function tags = groupTags(cohort)
tags = strings(numel(cohort.subjects), 1);
tags(:) = "unassigned";
tags(cohort.isCase) = string(cohort.casePrefix);
tags(cohort.isControl) = string(cohort.controlPrefix);
end

%% ------------------------------------------------------------------
function fig = plotDstFigure(cohort, Dst, welch, scheme)
%PLOTDSTFIGURE  Per-subject D_st bars with the two group means overlaid.
n = numel(cohort.subjects);
fig = figure('Name', sprintf('D_st across cohort (%s)', scheme), ...
    'Position', [80 80 max(880, 90 * n) 520], 'Color', 'w');
ax = axes(fig);
hold(ax, 'on');

caseColor = [0.85 0.33 0.10];
ctrlColor = [0.55 0.65 0.85];
barColors = repmat([0.70 0.70 0.70], n, 1);
barColors(cohort.isControl, :) = repmat(ctrlColor, nnz(cohort.isControl), 1);
barColors(cohort.isCase, :) = repmat(caseColor, nnz(cohort.isCase), 1);

b = bar(ax, 1:n, Dst, 0.7, 'FaceColor', 'flat', 'EdgeColor', [0.2 0.2 0.2]);
b.CData = barColors;

xlims = [0.4, n + 0.6];
meanCase = mean(Dst(cohort.isCase), 'omitnan');
meanCtrl = mean(Dst(cohort.isControl), 'omitnan');
if isfinite(meanCase)
    plot(ax, xlims, [meanCase meanCase], '--', 'Color', caseColor, 'LineWidth', 1.2, ...
        'DisplayName', sprintf('%s mean (%.4f)', strrep(cohort.casePrefix, '_', '\_'), meanCase));
end
if isfinite(meanCtrl)
    plot(ax, xlims, [meanCtrl meanCtrl], '--', 'Color', [0.20 0.30 0.70], 'LineWidth', 1.2, ...
        'DisplayName', sprintf('%s mean (%.4f)', strrep(cohort.controlPrefix, '_', '\_'), meanCtrl));
end

for i = 1:n
    if isfinite(Dst(i))
        text(ax, i, Dst(i), sprintf('  %.3f', Dst(i)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 9);
    end
end

set(ax, 'XTick', 1:n, 'XTickLabel', cohort.subjects, 'TickLabelInterpreter', 'none');
xtickangle(ax, 30);
xlim(ax, xlims);
maxD = max(Dst, [], 'omitnan');
if isfinite(maxD) && maxD > 0
    ylim(ax, [0, maxD * 1.15]);
end
ylabel(ax, 'D_{st} = trace(\Omega) / N', 'Interpreter', 'tex');

titleText = sprintf('Network-level D_{st} per subject  --  scheme = %s', scheme);
if isfinite(welch.p)
    titleText = sprintf('%s  (Welch t = %.2f, p = %.4f)', titleText, welch.tStat, welch.p);
end
title(ax, titleText, 'Interpreter', 'tex');
grid(ax, 'on');
box(ax, 'on');
legend(ax, 'Location', 'best', 'Interpreter', 'tex');
end
