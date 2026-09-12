function fpath = experimentProvenanceFile(resultsDir, experimentName, kind)
%EXPERIMENTPROVENANCEFILE  Path of one provenance artefact for a run.
%
%   fpath = experimentProvenanceFile(resultsDir, experimentName, 'parametersJson')
%
%   The run name is embedded in the filename, not only in the directory
%   name, so that a provenance file stays identifiable after someone copies
%   it out of its folder to share.
%
%   experimentName is the run folder name (config name + stamp), e.g.
%   'baseline_column_2026-05-23_1430'.
%
%   kind
%     'properties'     -> experiment_<name>.properties
%     'parametersMat'  -> experiment_parameters_<name>.mat
%     'parametersJson' -> experiment_parameters_<name>.json
%     'reports'        -> reports_<name>.log
%     'manifest'       -> run_manifest.mat  (name-independent)
%
%   See also EXPERIMENTRUNPARAMETERS, EXPERIMENTFOLDERNAME.

experimentName = char(experimentName);
switch lower(char(kind))
    case 'properties'
        fname = sprintf('experiment_%s.properties', experimentName);
    case 'parametersmat'
        fname = sprintf('experiment_parameters_%s.mat', experimentName);
    case 'parametersjson'
        fname = sprintf('experiment_parameters_%s.json', experimentName);
    case 'reports'
        fname = sprintf('reports_%s.log', experimentName);
    case 'manifest'
        fname = 'run_manifest.mat';
    otherwise
        error('experimentProvenanceFile:BadKind', ...
            ['kind must be properties, parametersMat, parametersJson, ' ...
             'reports or manifest (got "%s").'], kind);
end
fpath = fullfile(experimentResultsDir(resultsDir, 'provenance'), fname);
end
