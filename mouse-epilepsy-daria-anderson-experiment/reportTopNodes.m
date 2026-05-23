function reportTopNodes(varargin)
%REPORTTOPNODES  Print top-N structurally susceptible (or influential, or
% epileptogenic) regions from a section45 summary CSV produced by
% runAllMiceSection45.
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
%
% Parameters:
%   'N'              – number of top nodes to print  (default 10)
%   'Normalisation'  – 'parkes', 'tvb', or 'column'  (default 'parkes')
%   'SortBy'         – 'susc' | 'infl' | 'x0'        (default 'susc')
%   'ResultsDir'     – output folder (default results/)

p = inputParser();
p.addParameter('N',             10,       @(x) isnumeric(x) && x > 0 && isscalar(x));
p.addParameter('Normalisation', 'parkes', @(x) ischar(x) || isstring(x));
p.addParameter('SortBy',        'susc',   @(x) ischar(x) || isstring(x));
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
csvFile = fullfile(resultsDir, sprintf('section45_summary_topnodes_%s.csv', normStr));
if ~exist(csvFile, 'file')
    error('reportTopNodes: file not found:\n  %s\nRun runAllMiceSection45 first.', csvFile);
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
end
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
        elseif strcmp(sortStr, 'x0')
            fprintf('  %+9.4f', v);
        else
            fprintf('  %9.0f', v);
        end
    end
    fprintf('\n');
end
fprintf('%s\n', repmat('-', 1, sepLen));
fprintf('Columns: %s per mouse.\n\n', colLabel);
end
