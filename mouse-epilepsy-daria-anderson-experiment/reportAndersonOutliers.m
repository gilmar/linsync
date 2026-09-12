function reportAndersonOutliers(varargin)
%REPORTANDERSONOUTLIERS  Print per-Anderson outlier tables from the CSVs
% produced by compareAndersonVsArnold, then summarise which nodes are
% flagged across multiple Anderson mice.
%
% Usage:
%   reportAndersonOutliers                          % parkes, |z|>=2, all
%   reportAndersonOutliers('Normalisation', 'tvb')
%   reportAndersonOutliers('ZThreshold', 1.5)       % stricter filter
%   reportAndersonOutliers('N', 10)                 % cap per-mouse table
%   reportAndersonOutliers('SortBy', 'infl')        % sort by |z_infl|
%
% Parameters:
%   'Normalisation'  – 'parkes', 'tvb', or 'column' (default 'parkes')
%   'ZThreshold'     – minimum |z| to show  (default 2)
%   'N'              – max rows per mouse   (default Inf = show all)
%   'SortBy'         – 'susc' or 'infl'     (default 'susc')
%   'ResultsDir'     – output folder (default results/)

p = inputParser();
p.addParameter('Normalisation', 'parkes', @(x) ischar(x) || isstring(x));
p.addParameter('ZThreshold',    2,        @(x) isnumeric(x) && isscalar(x) && x >= 0);
p.addParameter('N',             Inf,      @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('SortBy',        'susc',   @(x) ischar(x) || isstring(x));
p.addParameter('ResultsDir',    '',       @(x) ischar(x) || isstring(x));
p.parse(varargin{:});
opts = p.Results;

normStr = char(opts.Normalisation);
sortStr = lower(char(opts.SortBy));

switch sortStr
    case 'susc',  zCol = 'z_susc'; valCol = 'D_susc_anderson'; meanCol = 'D_susc_arnold_mean'; sdCol = 'D_susc_arnold_std';
    case 'infl',  zCol = 'z_infl'; valCol = 'D_infl_anderson'; meanCol = 'D_infl_arnold_mean'; sdCol = 'D_infl_arnold_std';
    otherwise,    error('reportAndersonOutliers: SortBy must be ''susc'' or ''infl''.');
end

resultsDir = resolveMouseResultsDir(opts.ResultsDir);

% Discover Anderson outlier CSVs for this scheme (subfolder or legacy flat root)
pattern = sprintf('compare_anderson_vs_arnold_Anderson_*_%s_outliers.csv', normStr);
filePaths = {};
searchDirs = mouseResultsSearchDirs(resultsDir, 'anderson_vs_arnold');
for d = 1:numel(searchDirs)
    found = dir(fullfile(searchDirs{d}, pattern));
    for k = 1:numel(found)
        filePaths{end+1, 1} = fullfile(found(k).folder, found(k).name); %#ok<AGROW>
    end
end
if isempty(filePaths)
    error('reportAndersonOutliers: no outlier CSVs found matching:\n  %s\nRun compareAndersonVsArnold first.', pattern);
end
[filePaths, ord] = sort(filePaths);
filePaths = filePaths(ord);

% Extract mouse IDs
mouseIds = cell(numel(filePaths), 1);
for f = 1:numel(filePaths)
    [~, fname] = fileparts(filePaths{f});
    tok = regexp(fname, 'compare_anderson_vs_arnold_(Anderson_\w+)_', 'tokens', 'once');
    mouseIds{f} = tok{1};
end

sepLen = 74;
[valFmt, meanFmt, sdFmt] = andersonOutlierValueFormats(normStr);
fprintf('\nAnderson vs Arnold outliers  |  normalisation=%s  |  |z|>=%.1f  |  sorted by |%s|\n', ...
        normStr, opts.ZThreshold, zCol);
fprintf('%s\n', repmat('=', 1, sepLen));

% Per-mouse tables; also accumulate region->mice map for shared summary
sharedMap = containers.Map('KeyType', 'char', 'ValueType', 'any');

for f = 1:numel(filePaths)
    mouseId = mouseIds{f};
    T = readtable(filePaths{f}, 'TextType', 'string');

    % Filter and sort by |z|
    absZ   = abs(T.(zCol));
    T      = T(absZ >= opts.ZThreshold, :);
    absZ   = abs(T.(zCol));
    [~, ord] = sort(absZ, 'descend');
    T = T(ord, :);

    nShow = min(opts.N, height(T));
    fprintf('\n%s  --  %d outlier(s) at |z|>=%.1f\n', mouseId, height(T), opts.ZThreshold);
    fprintf('%s\n', repmat('-', 1, sepLen));
    fprintf('  %-32s  %8s  %8s  %8s  %8s  %s\n', ...
            'Region', 'Anderson', 'Arnold mu', 'Arnold sd', 'z', 'dir');
    fprintf('%s\n', repmat('-', 1, sepLen));

    for k = 1:nShow
        row  = T(k, :);
        zval = row.(zCol);
        if isnan(zval)
            zStr = '    NaN';
        else
            if zval > 0, dir_ = 'HIGH'; else, dir_ = 'LOW'; end
            zStr = sprintf('%+7.2f', zval);
        end
        if isnan(zval)
            dir_ = '-';
        end
        fprintf(['  %-32s  ' valFmt '  ' meanFmt '  ' sdFmt '  %s  %s\n'], ...
                row.Region, row.(valCol), row.(meanCol), row.(sdCol), zStr, dir_);
    end
    if height(T) > nShow
        fprintf('  ... (%d more not shown)\n', height(T) - nShow);
    end
    fprintf('%s\n', repmat('-', 1, sepLen));

    % Accumulate for shared-outlier summary
    for k = 1:height(T)
        reg = char(T.Region(k));
        zval = T.(zCol)(k);
        if sharedMap.isKey(reg)
            entry = sharedMap(reg);
            entry.mice{end+1}  = mouseId;
            entry.zvals(end+1) = zval;
            sharedMap(reg) = entry;
        else
            sharedMap(reg) = struct('mice', {{mouseId}}, 'zvals', zval);
        end
    end
end

% Shared outliers section
regKeys = keys(sharedMap);
sharedRegs = regKeys(cellfun(@(k) numel(sharedMap(k).mice) > 1, regKeys));

fprintf('\nShared outliers (flagged in >1 Anderson mouse):\n');
if isempty(sharedRegs)
    fprintf('  (none at |z|>=%.1f)\n', opts.ZThreshold);
else
    fprintf('%s\n', repmat('-', 1, sepLen));
    fprintf('  %-32s  %-16s  %s\n', 'Region', 'Mice', 'z per mouse');
    fprintf('%s\n', repmat('-', 1, sepLen));
    % Sort shared regions by number of mice flagged (desc), then alphabetically
    nMice = cellfun(@(k) numel(sharedMap(k).mice), sharedRegs);
    [~, ord] = sort(nMice, 'descend');
    sharedRegs = sharedRegs(ord);
    for s = 1:numel(sharedRegs)
        reg   = sharedRegs{s};
        entry = sharedMap(reg);
        miceStr = strjoin(entry.mice, ', ');
        zStr  = strjoin(arrayfun(@(z) sprintf('%+.2f', z), entry.zvals, ...
                        'UniformOutput', false), '  ');
        fprintf('  %-32s  %-16s  %s\n', reg, miceStr, zStr);
    end
    fprintf('%s\n', repmat('-', 1, sepLen));
end
fprintf('\n');
end

%% ------------------------------------------------------------------
function [valFmt, meanFmt, sdFmt] = andersonOutlierValueFormats(normStr)
%ANDERSONOUTLIERVALUEFORMATS  printf formats matched to scheme magnitude.
switch lower(normStr)
    case {'column', 'col', 'colnorm'}
        valFmt  = '%8.4f';
        meanFmt = '%8.4f';
        sdFmt   = '%8.4f';
    case 'parkes'
        valFmt  = '%10.1f';
        meanFmt = '%10.1f';
        sdFmt   = '%8.2f';
    case 'tvb'
        valFmt  = '%9.2f';
        meanFmt = '%9.2f';
        sdFmt   = '%7.3f';
    otherwise
        valFmt  = '%9.4g';
        meanFmt = '%9.4g';
        sdFmt   = '%8.4g';
end
end
