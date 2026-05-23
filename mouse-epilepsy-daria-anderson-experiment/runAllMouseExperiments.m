function summaries = runAllMouseExperiments(configsDir)
%RUNALLMOUSEEXPERIMENTS  Run every .properties file in configs/ (except template).
%
%   summaries = runAllMouseExperiments()
%   summaries = runAllMouseExperiments('configs')
%
%   Typical first use: runs initial_column, initial_parkes, and initial_tvb.

if nargin < 1 || isempty(configsDir)
    configsDir = fullfile(fileparts(mfilename('fullpath')), 'configs');
end

setupMousePaths();

listing = dir(fullfile(configsDir, '*.properties'));
names = {listing.name};
names = names(~strcmp(names, 'experiment.template.properties'));
if isempty(names)
    warning('runAllMouseExperiments:NoConfigs', 'No .properties files in %s', configsDir);
    summaries = struct([]);
    return;
end

names = sort(names);
summaries = struct('configFile', cell(numel(names), 1), ...
    'experimentName', '', 'allSucceeded', false);

fprintf('\n=== Running %d experiment config(s) from %s ===\n', numel(names), configsDir);

for k = 1:numel(names)
    propsPath = fullfile(configsDir, names{k});
    fprintf('\n---------- [%d/%d] %s ----------\n', k, numel(names), names{k});
    summaries(k).configFile = propsPath;
    try
        manifest = runMouseExperiment(propsPath);
        summaries(k).experimentName = manifest.experimentName;
        summaries(k).allSucceeded = manifest.allSucceeded;
    catch ME
        summaries(k).experimentName = '';
        summaries(k).allSucceeded = false;
        summaries(k).errorMessage = ME.message;
        warning('runAllMouseExperiments:Failed', '%s: %s', names{k}, ME.message);
    end
end

nOk = sum([summaries.allSucceeded]);
fprintf('\n=== Batch complete: %d/%d succeeded ===\n', nOk, numel(names));
end
