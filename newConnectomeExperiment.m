function propsPath = newConnectomeExperiment(varargin)
%NEWCONNECTOMEEXPERIMENT  Create a new experiment config from the template.
%
%   propsPath = newConnectomeExperiment()
%   propsPath = newConnectomeExperiment('Name', 'pilot', 'Normalisation', 'tvb')
%   propsPath = newConnectomeExperiment('ProjectRoot', studyRoot)
%
%   Copies configs/experiment.template.properties into the project's
%   configs/ folder under a new name, fills in the experiment name and
%   normalisation scheme, and prints the command to run it.
%
%   Starting every study from the annotated template means new configs
%   carry the documentation of each key with them, instead of being a bare
%   list of values copied from someone else's run.
%
%   Name-value options
%     'Name'          : experiment name. Default: <timestamp>_<scheme>
%     'Normalisation' : 'column' (default) | 'parkes' | 'tvb' | 'none'
%     'ProjectRoot'   : study folder (default: current folder)
%     'Template'      : explicit template path
%
%   See also RUNCONNECTOMEEXPERIMENT, LOADEXPERIMENTCONFIG.

p = inputParser;
addParameter(p, 'Name',          '', @(s) ischar(s) || isstring(s));
addParameter(p, 'Normalisation', 'column', @(s) ischar(s) || isstring(s));
addParameter(p, 'ProjectRoot',   '', @(s) ischar(s) || isstring(s));
addParameter(p, 'Template',      '', @(s) ischar(s) || isstring(s));
parse(p, varargin{:});
opts = p.Results;

scheme = normalisationSchemeInfo(opts.Normalisation).scheme;

projectRoot = char(opts.ProjectRoot);
if isempty(projectRoot)
    projectRoot = pwd;
end
configsDir = fullfile(projectRoot, 'configs');
if ~isfolder(configsDir)
    mkdir(configsDir);
end

templateFile = char(opts.Template);
if isempty(templateFile)
    toolkitRoot = fileparts(mfilename('fullpath'));
    candidates = { ...
        fullfile(configsDir, 'experiment.template.properties'), ...
        fullfile(toolkitRoot, 'configs', 'experiment.template.properties')};
    templateFile = '';
    for k = 1:numel(candidates)
        if isfile(candidates{k})
            templateFile = candidates{k};
            break;
        end
    end
end
if isempty(templateFile) || ~isfile(templateFile)
    error('newConnectomeExperiment:NoTemplate', ...
        'Could not find experiment.template.properties (looked in %s and the toolkit configs/).', ...
        configsDir);
end

name = char(opts.Name);
if isempty(name)
    name = sprintf('%s_%s', datestr(now, 'yyyy-mm-dd_HH-MM-SS'), scheme); %#ok<DATST,TNOW1>
end

propsPath = fullfile(configsDir, [name '.properties']);
if isfile(propsPath)
    error('newConnectomeExperiment:Exists', ...
        'A config named %s already exists. Choose another Name.', propsPath);
end

% Rewrite the three identity keys in place, leaving every comment and every
% other key exactly as the template has them.
%
% Note the [^\r\n]* rather than .* : MATLAB's regexp treats '.' as matching
% newlines by default, so a greedy '.*$' here would swallow the rest of the
% file and the new config would be truncated to a single line.
text = fileread(templateFile);
text = replaceKey(text, 'experiment.name', name);
text = replaceKey(text, 'experiment.description', ...
    ['Created ' datestr(now, 'yyyy-mm-dd HH:MM:SS')]); %#ok<DATST,TNOW1>
text = replaceKey(text, 'normalisation', scheme);

fid = fopen(propsPath, 'w');
if fid < 0
    error('newConnectomeExperiment:WriteFailed', 'Could not write %s', propsPath);
end
closer = onCleanup(@() fclose(fid));
fprintf(fid, '%s', text);

fprintf('Created experiment config:\n  %s\n', propsPath);
fprintf('Edit the data.root / data.file and cohort prefixes to match your study, then run:\n');
fprintf('  runConnectomeExperiment(''%s'')\n', propsPath);
end

%% ------------------------------------------------------------------
function text = replaceKey(text, key, value)
%REPLACEKEY  Set one key=value line in a .properties file body.
%
%   Appends the key when the template does not already carry it, so a
%   trimmed-down template still produces a usable config.

pattern = ['(?m)^' regexptranslate('escape', key) '=[^\r\n]*$'];
% In a replacement string only '\' and '$' are special (token references).
replacement = strrep(strrep([key '=' value], '\', '\\'), '$', '\$');
if isempty(regexp(text, pattern, 'once'))
    if ~isempty(text) && text(end) ~= newline
        text = [text newline];
    end
    text = [text key '=' value newline];
else
    text = regexprep(text, pattern, replacement, 'once');
end
end
