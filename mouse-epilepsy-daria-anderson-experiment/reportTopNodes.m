function reportTopNodes(varargin)
%REPORTTOPNODES  Print top-N structurally susceptible (or influential, or
% epileptogenic) regions from a stabilityCentralities summary CSV produced by
% runAllMiceStabilityCentralities.
%
% The CSV is sorted by the chosen centrality's mean rank across mice
% (ascending = best / most notable first).  Under the 'parkes' scheme,
% 'x0' sorting is meaningless (all NaN) and a warning is issued.
%
% Usage:
%   reportTopNodes                                      % top 10, parkes, D(->i)
%   reportTopNodes('N', 20)
%   reportTopNodes('Normalisation', 'tvb', 'SortBy', 'x0')
%   reportTopNodes('Normalisation', 'tvb', 'SortBy', 'infl')
%   reportTopNodes('IncludeTrivial', true)              % keep BACKGROUND / *_MASK rows
%
% Parameters:
%   'N'              – number of top nodes to print  (default 10)
%   'Normalisation'  – 'parkes', 'tvb', or 'column'  (default 'parkes')
%   'SortBy'         – 'susc' | 'infl' | 'x0'        (default 'susc')
%   'IncludeTrivial' – include placeholder/mask rows (BACKGROUND,
%                      *_MASK) in the ranking. Defaults to false:
%                      trivial nodes are excluded so the top-K reflects
%                      actual brain regions. Requires the summary CSV
%                      to carry an 'is_trivial' column (added in
%                      runAllMiceStabilityCentralities); legacy CSVs
%                      without that column behave as if true.
%   'ResultsDir'     – output folder (default results/)

p = inputParser();
p.addParameter('N',             10,       @(x) isnumeric(x) && x > 0 && isscalar(x));
p.addParameter('Normalisation', 'parkes', @(x) ischar(x) || isstring(x));
p.addParameter('SortBy',        'susc',   @(x) ischar(x) || isstring(x));
p.addParameter('IncludeTrivial', false,   @(x) islogical(x) && isscalar(x));
p.addParameter('ResultsDir',    '',       @(x) ischar(x) || isstring(x));
p.parse(varargin{:});
opts = p.Results;

normStr  = char(opts.Normalisation);
sortStr  = lower(char(opts.SortBy));
isParkes = strcmpi(normStr, 'parkes');
isColumn = ismember(lower(normStr), {'column', 'col', 'colnorm'});

if (isParkes || isColumn) && strcmp(sortStr, 'x0')
    warning('reportTopNodes:LinearX0', ...
        'x0^c is not computed under %s normalisation. Falling back to ''susc''.', normStr);
    sortStr = 'susc';
end

resultsDir = resolveMouseResultsDir(opts.ResultsDir);
prefix = mouseExperimentResultPrefix();
csvFile = fullfile(resultsDir, sprintf('%s_summary_topnodes_%s.csv', prefix, normStr));
if ~exist(csvFile, 'file')
    error('reportTopNodes: file not found:\n  %s\nRun runAllMiceStabilityCentralities first.', csvFile);
end

T = readtable(csvFile, 'TextType', 'string');

% Choose sort column and display columns
switch sortStr
    case 'susc'
        rankCol  = 'mean_rank_susc';
        dataPrefix = 'D_susc__';
        colLabel = 'D(\rightarrow i)';
    case 'infl'
        rankCol  = 'mean_rank_infl';
        dataPrefix = 'D_infl__';
        colLabel = 'D(k\rightarrow)';
    case 'x0'
        rankCol  = 'mean_rank_x0';
        dataPrefix = 'x0_crit__';
        colLabel = 'x0^c';
    otherwise
        error('reportTopNodes: SortBy must be ''susc'', ''infl'', or ''x0''.');
end

T = sortrows(T, rankCol, 'ascend');

% Optionally drop trivial / mask rows (BACKGROUND, *_MASK) from the ranking.
% Trivial nodes have isolated diagonal entries in C that inflate Omega_ii and
% otherwise dominate the susc / infl top-K under tvb. We keep them in the CSV
% so consumers can still find them, but exclude them from the printed top-K
% unless the caller explicitly opts in.
nTrivialDropped = 0;
if ~opts.IncludeTrivial
    if ismember('is_trivial', T.Properties.VariableNames)
        triv = logical(T.is_trivial);
        nTrivialDropped = sum(triv);
        T = T(~triv, :);
    else
        warning('reportTopNodes:NoTrivialColumn', ...
            ['Summary CSV %s has no ''is_trivial'' column (likely a legacy ' ...
             'run). Showing all rows including trivial placeholders. ' ...
             'Re-run runAllMiceStabilityCentralities to regenerate.'], csvFile);
    end
end

% Identify per-mouse data columns
allVars  = T.Properties.VariableNames;
isCols   = strncmp(allVars, dataPrefix, length(dataPrefix));
dataCols = allVars(isCols);
if isempty(dataCols)
    error('reportTopNodes: no columns with prefix ''%s'' found in %s.', dataPrefix, csvFile);
end
mouseNames = strrep(dataCols, dataPrefix, '');

N = min(opts.N, height(T));

% ---- Print ----
sepLen = 42 + 10 * numel(mouseNames);

fprintf('\nTop-%d regions  |  normalisation=%s  |  sorted by %s\n', N, normStr, rankCol);
if isParkes
    fprintf('Note: D(k->)=D(->i) for symmetric K; x0^c not computed.\n');
elseif isColumn
    fprintf('Note: x0^c not computed under column normalisation.\n');
end
if nTrivialDropped > 0
    fprintf('Note: %d trivial / mask region(s) (BACKGROUND, *_MASK) excluded ', nTrivialDropped);
    fprintf('(pass IncludeTrivial=true to keep them).\n');
end
valueFmt = centralityValueFormat(normStr, sortStr);
fprintf('%s\n', repmat('-', 1, sepLen));

% Header
fprintf('  %2s  %-32s  %8s', '#', 'Region', 'mean_rank');
for m = 1:numel(mouseNames)
    fprintf('  %9s', mouseNames{m});
end
fprintf('\n%s\n', repmat('-', 1, sepLen));

for k = 1:N
    row      = T(k, :);
    rankVal  = row.(rankCol);
    fprintf('  %2d. %-32s  %8.1f', k, row.region, rankVal);
    for m = 1:numel(dataCols)
        v = row.(dataCols{m});
        if isnan(v)
            fprintf('  %9s', 'NaN');
        else
            fprintf(['  ' valueFmt], v);
        end
    end
    fprintf('\n');
end
fprintf('%s\n', repmat('-', 1, sepLen));
fprintf('Columns: %s per mouse.\n\n', colLabel);
end

%% ------------------------------------------------------------------
function fmt = centralityValueFormat(normStr, sortStr)
%CENTRALITYVALUEFORMAT  printf format for per-mouse centrality magnitudes.
if strcmp(sortStr, 'x0')
    fmt = '%+9.4f';
    return;
end
switch lower(normStr)
    case {'column', 'col', 'colnorm'}
        fmt = '%9.4f';
    case 'parkes'
        fmt = '%10.1f';
    case 'tvb'
        fmt = '%9.2f';
    otherwise
        fmt = '%9.4g';
end
end
