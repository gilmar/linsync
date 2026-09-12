function report = reportTopNodes(varargin)
%REPORTTOPNODES  Print the most notable regions of a cohort stability run.
%
%   reportTopNodes('ResultsDir', resultsDir, 'Normalisation', 'column')
%   report = reportTopNodes(..., 'N', 15, 'SortBy', 'susceptibility')
%
%   Reads the cohort summary written by RUNCOHORTSTABILITYCENTRALITIES and
%   prints a ranked table to the console (and returns it as a table), so a
%   run can be sanity-checked without loading anything into the workspace.
%
%   Ranking uses the mean rank across subjects rather than the mean value:
%   the metrics are on different scales in different subjects, so averaging
%   raw values would let one subject with large values decide the order.
%
%   Name-value options
%     'ResultsDir'    : ''
%     'Normalisation' : 'column'
%     'N'             : 10 -- rows to print
%     'SortBy'        : 'auto' -- 'criticalX0' when a sweep was run, else
%                       'susceptibility'; also accepts 'influence'
%     'IncludeTrivial': false -- trivial placeholder nodes are hidden by
%                       default because normalisation can push them to the
%                       top of a ranking despite carrying no anatomy
%
%   See also RUNCOHORTSTABILITYCENTRALITIES, REPORTCOHORTOUTLIERS.

p = inputParser;
addParameter(p, 'ResultsDir',     '', @(s) ischar(s) || isstring(s));
addParameter(p, 'Normalisation',  'column', @(s) ischar(s) || isstring(s));
addParameter(p, 'N',              10, @isscalar);
addParameter(p, 'SortBy',         'auto', @(s) ischar(s) || isstring(s));
addParameter(p, 'IncludeTrivial', false, @islogical);
parse(p, varargin{:});
opts = p.Results;

scheme = normalisationSchemeInfo(opts.Normalisation).scheme;
resultsDir = resolveExperimentResultsDir(opts.ResultsDir);

matFile = locateExperimentResultsFile(resultsDir, 'stability', ...
    sprintf('%s_summary_%s.mat', experimentResultPrefix(), scheme));
if isempty(matFile)
    error('reportTopNodes:NoSummary', ...
        ['No %s_summary_%s.mat under %s. Run runCohortStabilityCentralities ' ...
         'for this scheme first.'], experimentResultPrefix(), scheme, resultsDir);
end
S = load(matFile);
if ~isfield(S, 'summary')
    error('reportTopNodes:BadSummary', '%s does not contain a summary struct.', matFile);
end
summary = S.summary;

labels = summary.labels(:);
isTrivial = logical(summary.isTrivialAny(:));

sortBy = lower(char(opts.SortBy));
if strcmp(sortBy, 'auto')
    if any(isfinite(summary.meanRankX0))
        sortBy = 'criticalx0';
    else
        sortBy = 'susceptibility';
    end
end

switch sortBy
    case {'criticalx0', 'x0', 'x0crit'}
        key = summary.meanRankX0;
        keyName = 'mean rank of x0^c (lower = more epileptogenic)';
    case {'susceptibility', 'susc', 'd_to_i'}
        key = summary.meanRankSusc;
        keyName = 'mean rank of D(->i) (lower = more susceptible)';
    case {'influence', 'infl', 'd_k_to'}
        key = summary.meanRankInfl;
        keyName = 'mean rank of D(k->) (lower = more influential)';
    otherwise
        error('reportTopNodes:BadSortBy', ...
            'SortBy must be auto, criticalX0, susceptibility or influence.');
end

keep = true(numel(labels), 1);
if ~opts.IncludeTrivial
    keep = ~isTrivial;
end
idx = find(keep);
[~, order] = sort(key(idx), 'ascend');
idx = idx(order);
nShow = min(opts.N, numel(idx));
idx = idx(1:nShow);

meanSusc = mean(summary.D_susc, 2, 'omitnan');
meanInfl = mean(summary.D_infl, 2, 'omitnan');
meanX0 = mean(summary.x0_crit, 2, 'omitnan');

report = table((1:nShow).', string(labels(idx)), key(idx), ...
    meanX0(idx), meanSusc(idx), meanInfl(idx), isTrivial(idx), ...
    'VariableNames', {'rank', 'region', 'sort_key', ...
    'mean_x0_crit', 'mean_D_to_i', 'mean_D_k_to', 'is_trivial'});

fprintf('\nTop %d regions  --  scheme = %s, %d subject(s)\n', ...
    nShow, scheme, numel(summary.subjects));
fprintf('Ranked by %s\n', keyName);
if ~opts.IncludeTrivial && any(isTrivial)
    fprintf('(%d trivial placeholder node(s) hidden; pass IncludeTrivial=true to show them)\n', ...
        sum(isTrivial));
end
fprintf('%-5s %-34s %10s %12s %12s %12s\n', ...
    'rank', 'region', 'sortKey', 'x0^c', 'D(->i)', 'D(k->)');
for k = 1:nShow
    i = idx(k);
    fprintf('%-5d %-34s %10.2f %12s %12s %12s\n', k, labels{i}, key(i), ...
        numOrDash(meanX0(i)), numOrDash(meanSusc(i)), numOrDash(meanInfl(i)));
end
fprintf('\n');
end

%% ------------------------------------------------------------------
function s = numOrDash(v)
if isfinite(v)
    s = sprintf('%.4g', v);
else
    s = '-';
end
end
