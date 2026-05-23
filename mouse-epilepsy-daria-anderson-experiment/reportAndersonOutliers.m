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

% Discover Anderson outlier CSVs for this scheme
pattern = fullfile(resultsDir, sprintf('compare_anderson_vs_arnold_Anderson_*_%s_outliers.csv', normStr));
files = dir(pattern);
if isempty(files)
    error('reportAndersonOutliers: no outlier CSVs found matching:\n  %s\nRun compareAndersonVsArnold first.', pattern);
end
files = sort({files.name});   % alphabetical = Anderson_1, Anderson_2, ...

% Extract mouse IDs
mouseIds = cell(numel(files), 1);
for f = 1:numel(files)
    tok = regexp(files{f}, 'compare_anderson_vs_arnold_(Anderson_\w+)_', 'tokens', 'once');
    mouseIds{f} = tok{1};
end

sepLen = 74;
fprintf('\nAnderson vs Arnold outliers  |  normalisation=%s  |  |z|>=%.1f  |  sorted by |%s|\n', ...
        normStr, opts.ZThreshold, zCol);
fprintf('%s\n', repmat('=', 1, sepLen));

% Per-mouse tables; also accumulate region->mice map for shared summary
sharedMap = containers.Map('KeyType', 'char', 'ValueType', 'any');

for f = 1:numel(files)
    mouseId = mouseIds{f};
    T = readtable(fullfile(resultsDir, files{f}), 'TextType', 'string');

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
        if zval > 0, dir_ = 'HIGH'; else, dir_ = 'LOW'; end
        fprintf('  %-32s  %8.1f  %8.1f  %8.1f  %+7.2f  %s\n', ...
                row.Region, row.(valCol), row.(meanCol), row.(sdCol), zval, dir_);
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
