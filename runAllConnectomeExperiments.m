function summaries = runAllConnectomeExperiments(configsDir, varargin)
%RUNALLCONNECTOMEEXPERIMENTS  Run every .properties config in a folder.
%
%   summaries = runAllConnectomeExperiments()
%   summaries = runAllConnectomeExperiments('configs')
%   summaries = runAllConnectomeExperiments(configsDir, 'ProjectRoot', studyRoot)
%
%   Runs the configs in alphabetical order, skipping the template. Useful
%   for sweeping the same cohort through several normalisation schemes:
%   each config writes to its own timestamped results folder, so the runs
%   are directly comparable afterwards and none can overwrite another.
%
%   A config that fails is recorded and the batch continues, so one broken
%   parameter set does not cost the whole sweep.
%
%   Returns a struct array: configFile, experimentName, allSucceeded,
%   errorMessage.
%
%   See also RUNCONNECTOMEEXPERIMENT, NEWCONNECTOMEEXPERIMENT.

if nargin < 1 || isempty(configsDir)
    configsDir = fullfile(pwd, 'configs');
end
configsDir = char(configsDir);

p = inputParser;
addParameter(p, 'ProjectRoot', '', @(s) ischar(s) || isstring(s));
parse(p, varargin{:});

if ~isfolder(configsDir)
    error('runAllConnectomeExperiments:NoConfigsDir', ...
        'Config folder does not exist: %s', configsDir);
end

listing = dir(fullfile(configsDir, '*.properties'));
names = sort({listing.name});
names = names(~strcmpi(names, 'experiment.template.properties'));
if isempty(names)
    warning('runAllConnectomeExperiments:NoConfigs', ...
        'No .properties files (other than the template) in %s', configsDir);
    summaries = struct('configFile', {}, 'experimentName', {}, ...
        'allSucceeded', {}, 'errorMessage', {});
    return;
end

summaries = struct('configFile', cell(numel(names), 1), ...
    'experimentName', '', 'allSucceeded', false, 'errorMessage', '');

fprintf('\n=== Running %d experiment config(s) from %s ===\n', numel(names), configsDir);

for k = 1:numel(names)
    propsPath = fullfile(configsDir, names{k});
    fprintf('\n---------- [%d/%d] %s ----------\n', k, numel(names), names{k});
    summaries(k).configFile = propsPath;
    try
        manifest = runConnectomeExperiment(propsPath, ...
            'ProjectRoot', p.Results.ProjectRoot);
        summaries(k).experimentName = manifest.experimentName;
        summaries(k).allSucceeded = manifest.allSucceeded;
        if ~manifest.allSucceeded
            failed = manifest.steps(~[manifest.steps.success]);
            summaries(k).errorMessage = strjoin( ...
                strcat({failed.name}, ': ', {failed.errorMessage}), ' | ');
        end
    catch ME
        summaries(k).allSucceeded = false;
        summaries(k).errorMessage = ME.message;
        warning('runAllConnectomeExperiments:ConfigFailed', ...
            '%s: %s', names{k}, ME.message);
    end
end

nOk = sum([summaries.allSucceeded]);
fprintf('\n=== Batch complete: %d/%d config(s) succeeded ===\n', nOk, numel(names));
for k = 1:numel(summaries)
    if ~summaries(k).allSucceeded
        fprintf('  FAILED %s: %s\n', names{k}, summaries(k).errorMessage);
    end
end
end
