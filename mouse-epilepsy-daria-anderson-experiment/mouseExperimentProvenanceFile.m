function fpath = mouseExperimentProvenanceFile(resultsDir, experimentName, kind)
%MOUSEEXPERIMENTPROVENANCEFILE  Suffixed provenance artefact path for one run.
%
%   fpath = mouseExperimentProvenanceFile(resultsDir, experimentName, kind)
%
%   experimentName is the results subfolder name (config + run stamp), e.g.
%   'initial_column_2026-05-23_1430'. Filenames embed that suffix so copies
%   outside the folder stay identifiable.
%
%   kind:
%     'properties'      -> experiment_<name>.properties
%     'parametersMat'   -> experiment_parameters_<name>.mat
%     'parametersJson'  -> experiment_parameters_<name>.json
%     'reports'         -> reports_<name>.log

experimentName = char(experimentName);
switch lower(kind)
    case 'properties'
        fname = sprintf('experiment_%s.properties', experimentName);
    case 'parametersmat'
        fname = sprintf('experiment_parameters_%s.mat', experimentName);
    case 'parametersjson'
        fname = sprintf('experiment_parameters_%s.json', experimentName);
    case 'reports'
        fname = sprintf('reports_%s.log', experimentName);
    otherwise
        error('mouseExperimentProvenanceFile:BadKind', ...
            'kind must be properties, parametersMat, parametersJson, or reports.');
end
fpath = fullfile(mouseResultsDir(resultsDir, 'provenance'), fname);
end
