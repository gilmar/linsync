function files = listPerMouseResultFiles(resultsDir, scheme)
%LISTPERMOUSERESULTFILES  Dir listing for per-mouse *_<scheme>_results.mat files.
%
%   Primary prefix: stabilityCentralities_* (current pipeline).
%   Also accepts legacy section45_* (e.g. mouse_experiment_reference/).

if nargin < 2
    error('listPerMouseResultFiles:Args', 'resultsDir and scheme are required.');
end

scheme = char(scheme);
allFiles = [];
for prefix = {'stabilityCentralities', 'section45'}
    found = dir(fullfile(resultsDir, sprintf('%s_*_%s_results.mat', prefix{1}, scheme)));
    allFiles = [allFiles; found]; %#ok<AGROW>
end

if isempty(allFiles)
    files = allFiles;
    return;
end

byMouse = containers.Map('KeyType', 'char', 'ValueType', 'any');
for k = 1:numel(allFiles)
    mouseId = extractMouseIdFromResultFile(allFiles(k).name);
    if isempty(mouseId)
        continue;
    end
    if ~isKey(byMouse, mouseId)
        byMouse(mouseId) = allFiles(k);
    elseif startsWith(allFiles(k).name, 'stabilityCentralities_')
        byMouse(mouseId) = allFiles(k);
    end
end

vals = byMouse.values;
if isempty(vals)
    files = allFiles(1:0);
else
    files = vertcat(vals{:});
    files = files(:);
end
end

function mouseId = extractMouseIdFromResultFile(fname)
[~, base, ext] = fileparts(fname);
if ~strcmp(ext, '.mat') || ~endsWith(base, '_results')
    mouseId = '';
    return;
end
base = regexprep(base, '^(stabilityCentralities|section45)_', '');
base = regexprep(base, '_results$', '');
parts = strsplit(base, '_');
if numel(parts) < 2
    mouseId = '';
    return;
end
scheme = parts{end};
if ~ismember(scheme, {'tvb', 'parkes', 'column', 'none', 'col', 'colnorm'})
    mouseId = '';
    return;
end
mouseId = strjoin(parts(1:end-1), '_');
end
