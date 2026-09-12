function summary = runCohortStabilityCentralities(varargin)
%RUNCOHORTSTABILITYCENTRALITIES  Stability centralities for every subject in a cohort.
%
%   summary = runCohortStabilityCentralities('DataRoot', dataRoot, ...
%                'ResultsDir', resultsDir, 'Normalisation', 'column')
%
%   Runs RUNSTABILITYCENTRALITIES over every subject discovered under
%   DataRoot (or the explicit 'Subjects' list), then builds the cross-subject
%   summary that the comparison steps consume:
%
%     * per-node values packed into N x nSubjects matrices
%     * per-subject ranks, and mean rank across the cohort, so that regions
%       are comparable even when the metrics themselves are on different
%       scales between subjects
%     * a top-nodes CSV, an overview heatmap, and a packed .mat
%
%   The cohort is treated as all-or-nothing: if any subject fails, the run
%   errors instead of writing a partial summary. A summary silently missing
%   a subject would quietly bias every cohort mean computed from it.
%
%   Name-value options
%     'DataRoot'       : '' -- subject folders (default: pwd/data)
%     'ConnectomeFile' : 'connectome.csv'
%     'Subjects'       : {} -- explicit ID list; default is discovery
%     'Normalisation'  : 'column' | 'parkes' | 'tvb' | 'none'
%     'TopK'           : 10 -- rows listed per subject in the console report
%     'ParkesC', 'ColScale', 'TvbPercentile', 'Tau0',
%     'X0Base', 'X0Upper', 'X0Step', 'BisectTol', 'MaxK',
%     'TrivialLabelSuffixes', 'SweepCriticalX0'
%                      : forwarded to runStabilityCentralities
%     'Plot', 'Verbose', 'SaveResults', 'ResultsDir', 'RunParameters'
%
%   Outputs (under <ResultsDir>/stability_centralities/)
%     <prefix>_<subjectId>_<scheme>_results.mat  + _figure.{fig,png}
%     <prefix>_summary_topnodes_<scheme>.csv
%     <prefix>_summary_overview_<scheme>.{fig,png}
%     <prefix>_summary_<scheme>.mat
%
%   See also RUNSTABILITYCENTRALITIES, COMPARECOHORTGROUPS, REPORTTOPNODES.

p = inputParser;
addParameter(p, 'DataRoot',       '', @(s) ischar(s) || isstring(s));
addParameter(p, 'ConnectomeFile', 'connectome.csv', @(s) ischar(s) || isstring(s));
addParameter(p, 'Subjects',       {}, @(c) iscell(c) || isstring(c));
addParameter(p, 'Normalisation',  'column', @(s) ischar(s) || isstring(s));
addParameter(p, 'TopK',           10,   @isscalar);
addParameter(p, 'ParkesC',        1.0,  @isscalar);
addParameter(p, 'ColScale',       0.95, @isscalar);
addParameter(p, 'TvbPercentile',  95,   @isscalar);
addParameter(p, 'Tau0',           6667, @isscalar);
addParameter(p, 'X0Base',         -2.3, @isscalar);
addParameter(p, 'X0Upper',        -1.0, @isscalar);
addParameter(p, 'X0Step',         0.01, @isscalar);
addParameter(p, 'BisectTol',      1e-4, @isscalar);
addParameter(p, 'MaxK',           1e8,  @isscalar);
addParameter(p, 'SweepCriticalX0', true, @islogical);
addParameter(p, 'TrivialLabelSuffixes', {}, @(c) iscell(c) || isstring(c) || ischar(c));
addParameter(p, 'Plot',           true, @islogical);
addParameter(p, 'Verbose',        true, @islogical);
addParameter(p, 'SaveResults',    true, @islogical);
addParameter(p, 'ResultsDir',     '',   @(s) ischar(s) || isstring(s));
addParameter(p, 'RunParameters',  [],   @(x) isempty(x) || isstruct(x));
parse(p, varargin{:});
opts = p.Results;

scheme = normalisationSchemeInfo(opts.Normalisation).scheme;
topK = opts.TopK;

resultsDir = resolveExperimentResultsDir(opts.ResultsDir);
stabilityDir = experimentResultsDir(resultsDir, 'stability');
prefix = experimentResultPrefix();

dataRoot = char(opts.DataRoot);
if isempty(dataRoot)
    dataRoot = fullfile(pwd, 'data');
end

subjects = opts.Subjects;
if isempty(subjects)
    subjects = listConnectomeSubjects(dataRoot, opts.ConnectomeFile);
else
    subjects = cellstr(subjects);
    subjects = subjects(:);
end
if isempty(subjects)
    error('runCohortStabilityCentralities:NoSubjects', ...
        'No subjects with a "%s" file found under %s.', opts.ConnectomeFile, dataRoot);
end

%% ---- Per-subject runs ---------------------------------------------
M = numel(subjects);
allResults = cell(M, 1);
failures = cell(0, 1);
totalTic = tic;
for m = 1:M
    subjectId = subjects{m};
    fprintf('\n========== %s (%d/%d) ==========\n', subjectId, m, M);
    try
        allResults{m} = runStabilityCentralities(subjectId, ...
            'DataRoot',        dataRoot, ...
            'ConnectomeFile',  opts.ConnectomeFile, ...
            'TrivialLabelSuffixes', opts.TrivialLabelSuffixes, ...
            'Normalisation',   scheme, ...
            'ParkesC',         opts.ParkesC, ...
            'ColScale',        opts.ColScale, ...
            'TvbPercentile',   opts.TvbPercentile, ...
            'Tau0',            opts.Tau0, ...
            'X0Base',          opts.X0Base, ...
            'X0Upper',         opts.X0Upper, ...
            'X0Step',          opts.X0Step, ...
            'BisectTol',       opts.BisectTol, ...
            'MaxK',            opts.MaxK, ...
            'SweepCriticalX0', opts.SweepCriticalX0, ...
            'ResultsDir',      resultsDir, ...
            'RunParameters',   opts.RunParameters, ...
            'SaveResults',     opts.SaveResults, ...
            'Plot',            opts.Plot, ...
            'Verbose',         opts.Verbose);
    catch ME
        warning('runCohortStabilityCentralities:SubjectFailed', ...
            '%s failed: %s', subjectId, ME.message);
        allResults{m} = [];
        failures{end+1, 1} = sprintf('%s (%s)', subjectId, ME.message); %#ok<AGROW>
    end
end
fprintf('\nCohort processed in %.1f s.\n', toc(totalTic));

if ~isempty(failures)
    error('runCohortStabilityCentralities:IncompleteCohort', ...
        ['%d of %d subject(s) failed, so no cohort summary was written ' ...
         '(a partial summary would bias every cohort mean):\n  %s'], ...
        numel(failures), M, strjoin(failures, sprintf('\n  ')));
end

%% ---- Pack the cohort ----------------------------------------------
canonicalLabels = allResults{1}.labels;
N = numel(canonicalLabels);
for m = 2:M
    if numel(allResults{m}.labels) ~= N || ...
            ~all(string(allResults{m}.labels) == string(canonicalLabels))
        warning('runCohortStabilityCentralities:LabelMismatch', ...
            ['Node labels for %s differ from %s. The summary assumes the ' ...
             'subjects share one parcellation and aligns them by index.'], ...
            subjects{m}, subjects{1});
    end
end

D_susc_all = NaN(N, M);
D_infl_all = NaN(N, M);
x0_crit_all = NaN(N, M);
D_st_all = NaN(1, M);
isTrivialAny = false(N, 1);
for m = 1:M
    r = allResults{m};
    D_susc_all(:, m)  = r.D_susceptibility;
    D_infl_all(:, m)  = r.D_influence;
    x0_crit_all(:, m) = r.x0_crit;
    D_st_all(m)       = r.D_st_healthy;
    isTrivialAny = isTrivialAny | logical(r.isTrivial(:));
end

rank_susc = NaN(N, M);
rank_infl = NaN(N, M);
rank_x0   = NaN(N, M);
for m = 1:M
    rank_susc(:, m) = rankWithNaN(D_susc_all(:, m), 'descend');
    rank_infl(:, m) = rankWithNaN(D_infl_all(:, m), 'descend');
    rank_x0(:, m)   = rankWithNaN(x0_crit_all(:, m), 'ascend');
end
mean_rank_susc = mean(rank_susc, 2, 'omitnan');
mean_rank_infl = mean(rank_infl, 2, 'omitnan');
mean_rank_x0   = mean(rank_x0,   2, 'omitnan');

hasCriticalX0 = any(isfinite(x0_crit_all(:)));
if hasCriticalX0
    sortKey = mean_rank_x0;
    sortDesc = 'smallest mean rank of x^c_0 (most epileptogenic)';
else
    % No sweep was run, so every x0 rank is the NaN sentinel. Report it as
    % missing rather than as a meaningless constant.
    mean_rank_x0(:) = NaN;
    sortKey = mean_rank_susc;
    sortDesc = 'smallest mean rank of D(->i) (most susceptible)';
end
[~, aggOrder] = sort(sortKey, 'ascend');

%% ---- Console report -----------------------------------------------
fprintf('\n----- Aggregate top-%d nodes: %s -----\n', min(topK, N), sortDesc);
for kk = 1:min(topK, N)
    i = aggOrder(kk);
    if hasCriticalX0
        fprintf('  %2d. %-32s  mean rank: x0^c = %5.1f | D(->i) = %5.1f | D(k->) = %5.1f%s\n', ...
            kk, canonicalLabels{i}, mean_rank_x0(i), mean_rank_susc(i), ...
            mean_rank_infl(i), trivialTag(isTrivialAny(i)));
    else
        fprintf('  %2d. %-32s  mean rank: D(->i) = %5.1f | D(k->) = %5.1f%s\n', ...
            kk, canonicalLabels{i}, mean_rank_susc(i), mean_rank_infl(i), ...
            trivialTag(isTrivialAny(i)));
    end
end

%% ---- Summary CSV ---------------------------------------------------
header = {'rank', 'region', 'is_trivial', 'mean_rank_x0', 'mean_rank_susc', 'mean_rank_infl'};
for m = 1:M, header{end+1} = sprintf('x0_crit__%s', subjects{m}); end %#ok<AGROW>
for m = 1:M, header{end+1} = sprintf('D_susc__%s',  subjects{m}); end %#ok<AGROW>
for m = 1:M, header{end+1} = sprintf('D_infl__%s',  subjects{m}); end %#ok<AGROW>

rows = cell(N, numel(header));
for kk = 1:N
    i = aggOrder(kk);
    rows{kk, 1} = kk;
    rows{kk, 2} = canonicalLabels{i};
    rows{kk, 3} = double(isTrivialAny(i));
    rows{kk, 4} = mean_rank_x0(i);
    rows{kk, 5} = mean_rank_susc(i);
    rows{kk, 6} = mean_rank_infl(i);
    col = 7;
    for m = 1:M, rows{kk, col} = x0_crit_all(i, m);  col = col + 1; end
    for m = 1:M, rows{kk, col} = D_susc_all(i, m);   col = col + 1; end
    for m = 1:M, rows{kk, col} = D_infl_all(i, m);   col = col + 1; end
end
T = cell2table(rows, 'VariableNames', header);

%% ---- Overview figure -----------------------------------------------
overviewFig = [];
if opts.Plot
    overviewFig = buildOverviewFigure(scheme, subjects, canonicalLabels, ...
        D_susc_all, D_infl_all, x0_crit_all, hasCriticalX0);
end

%% ---- Provenance and save -------------------------------------------
if isempty(opts.RunParameters)
    runParameters = experimentRunParameters('fromOptions', opts, ...
        'Analysis', 'runCohortStabilityCentralities', ...
        'ResultsDir', resultsDir, 'Scheme', scheme, 'Subjects', subjects);
else
    runParameters = opts.RunParameters;
end

summary = struct();
summary.scheme          = scheme;
summary.subjects        = subjects;
summary.labels          = canonicalLabels;
summary.D_susc          = D_susc_all;
summary.D_infl          = D_infl_all;
summary.x0_crit         = x0_crit_all;
summary.D_st            = D_st_all;
summary.meanRankSusc    = mean_rank_susc;
summary.meanRankInfl    = mean_rank_infl;
summary.meanRankX0      = mean_rank_x0;
summary.isTrivialAny    = isTrivialAny;
summary.topNodesTable   = T;
summary.runParameters   = runParameters;

if opts.SaveResults
    csvFile = fullfile(stabilityDir, sprintf('%s_summary_topnodes_%s.csv', prefix, scheme));
    writetable(T, csvFile);
    fprintf('\nWrote %s\n', csvFile);

    if ~isempty(overviewFig) && isgraphics(overviewFig, 'figure')
        saveFigureBoth(overviewFig, ...
            fullfile(stabilityDir, sprintf('%s_summary_overview_%s', prefix, scheme)));
        close(overviewFig);
    end

    % Named copies so the .mat is readable field-by-field without
    % unpacking the summary struct.
    subjectIds = subjects;
    labels = canonicalLabels;
    normalisation = scheme;
    save(fullfile(stabilityDir, sprintf('%s_summary_%s.mat', prefix, scheme)), ...
        'summary', 'subjectIds', 'labels', 'normalisation', ...
        'D_susc_all', 'D_infl_all', 'x0_crit_all', ...
        'mean_rank_susc', 'mean_rank_infl', 'mean_rank_x0', ...
        'isTrivialAny', 'runParameters');
elseif ~isempty(overviewFig) && isgraphics(overviewFig, 'figure')
    close(overviewFig);
end
end

%% ------------------------------------------------------------------
function tag = trivialTag(isTrivial)
if isTrivial
    tag = '   [trivial]';
else
    tag = '';
end
end

%% ------------------------------------------------------------------
function fig = buildOverviewFigure(scheme, subjects, labels, ...
    D_susc_all, D_infl_all, x0_crit_all, hasCriticalX0)
%BUILDOVERVIEWFIGURE  Region x subject heatmaps, one panel per informative metric.
schemeInfo = normalisationSchemeInfo(scheme);

panels = {};
if hasCriticalX0
    panels{end+1} = struct('data', x0_crit_all, ...
        'title', 'x^c_{0,i}  (lower = more epileptogenic)');
end
panels{end+1} = struct('data', D_susc_all, 'title', 'D(\rightarrow i)  (susceptibility)');
if ~schemeInfo.isSymmetric
    panels{end+1} = struct('data', D_infl_all, 'title', 'D(k \rightarrow)  (influence)');
end

nPanels = numel(panels);
N = numel(labels);
M = numel(subjects);

fig = figure('Name', sprintf('Cohort stability centralities (%s)', scheme), ...
             'Position', [80 80 min(1600, 560 * nPanels) 700]);
for k = 1:nPanels
    subplot(1, nPanels, k);
    imagesc(panels{k}.data);
    yticks(1:N); yticklabels(labels);
    set(gca, 'FontSize', 6, 'YDir', 'normal', 'TickLabelInterpreter', 'none');
    xticks(1:M); xticklabels(subjects); xtickangle(45);
    colorbar; colormap(gca, parula);
    title(panels{k}.title, 'Interpreter', 'tex');
end

if schemeInfo.isSymmetric
    sgtitle(sprintf(['Cross-subject stability centralities  --  normalisation = %s' ...
                     '  [D(k\\rightarrow) \\equiv D(\\rightarrow i) for symmetric K]'], scheme));
else
    sgtitle(sprintf('Cross-subject stability centralities  --  normalisation = %s', scheme));
end
end
