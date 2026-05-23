function propsPath = newMouseExperiment(scheme)
%NEWMOUSEEXPERIMENT  Create a timestamped experiment .properties file.
%
%   propsPath = newMouseExperiment('column')
%   propsPath = newMouseExperiment()   % default scheme: column
%
%   Writes configs/<yyyy-MM-dd_HH-mm-ss>_<scheme>.properties from the
%   template and prints how to run it.

if nargin < 1 || isempty(scheme)
    scheme = 'column';
end
scheme = normalizeScheme(scheme);

setupMousePaths();

experimentRoot = fileparts(mfilename('fullpath'));
configsDir = fullfile(experimentRoot, 'configs');
templateFile = fullfile(configsDir, 'experiment.template.properties');
if ~isfile(templateFile)
    error('newMouseExperiment:NoTemplate', 'Missing %s', templateFile);
end

expName = [datestr(now, 'yyyy-mm-dd_HH-MM-SS') '_' scheme];
propsPath = fullfile(configsDir, [expName '.properties']);

copyfile(templateFile, propsPath);

text = fileread(propsPath);
text = strrep(text, 'experiment.name=REPLACE_ME', ['experiment.name=' expName]);
text = strrep(text, 'experiment.description=', ...
    ['experiment.description=Created ' datestr(now, 'yyyy-mm-dd HH:MM:SS')]);
text = regexprep(text, 'normalisation=\w+', ['normalisation=' scheme], 'once');
fid = fopen(propsPath, 'w');
if fid < 0
    error('newMouseExperiment:WriteFailed', 'Could not write %s', propsPath);
end
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s', text);

fprintf('Created experiment config:\n  %s\n', propsPath);
fprintf('Edit parameters if needed, then run:\n');
fprintf('  runMouseExperiment(''%s'')\n', propsPath);
end

function scheme = normalizeScheme(scheme)
scheme = lower(strtrim(char(scheme)));
switch scheme
    case {'col', 'colnorm'}, scheme = 'column';
    case {'tvb', 'none', 'parkes', 'column'}
    otherwise
        error('newMouseExperiment:BadScheme', ...
            'scheme must be tvb, none, parkes, or column.');
end
end
