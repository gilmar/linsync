function varargout = experimentRunParameters(action, varargin)
%EXPERIMENTRUNPARAMETERS  Build, merge and save the provenance record of a run.
%
%   runParams = experimentRunParameters('fromConfig', cfg, 'ResultsDir', d, 'Subjects', s)
%   runParams = experimentRunParameters('fromOptions', opts, 'Analysis', name)
%   runParams = experimentRunParameters('merge', runParams, patch)
%   experimentRunParameters('save', resultsDir, runParams)
%
%   A figure is only as trustworthy as the record of what produced it. This
%   function builds one struct carrying every parameter, the software
%   environment and the cohort, which the orchestrator writes to disk twice
%   (at run start and at run end) and which every analysis embeds inside
%   its own .mat output. Copy a result file anywhere and it still says how
%   it was made.
%
%   Actions
%     'fromConfig'  Build a full record from a parsed .properties config.
%                   Name-value: 'ResultsDir', 'Subjects'.
%     'fromOptions' Build a partial record from an inputParser Results
%                   struct, for analyses called directly rather than
%                   through the orchestrator. Name-value: 'Analysis',
%                   'ResultsDir', 'Scheme', 'Subjects'.
%     'merge'       Recursively merge a patch struct into a record, so an
%                   analysis can add its own options to the run-level
%                   record without discarding it.
%     'save'        Write experiment_parameters_<name>.{mat,json} under the
%                   run's provenance/ folder. JSON is written alongside the
%                   .mat so the record is diffable and readable without
%                   MATLAB.
%
%   See also RUNCONNECTOMEEXPERIMENT, EXPERIMENTPROVENANCEFILE.

switch lower(char(action))
    case {'fromconfig', 'buildfromconfig'}
        varargout{1} = buildFromConfig(varargin{:});
    case {'fromoptions', 'buildfromoptions'}
        varargout{1} = buildFromOptions(varargin{:});
    case 'merge'
        varargout{1} = mergeParams(varargin{:});
    case 'save'
        saveParams(varargin{:});
    otherwise
        error('experimentRunParameters:BadAction', ...
            'Unknown action "%s". Use fromConfig, fromOptions, merge or save.', action);
end
end

%% ------------------------------------------------------------------
function runParams = buildFromConfig(cfg, varargin)
p = inputParser;
addParameter(p, 'ResultsDir', '', @(s) ischar(s) || isstring(s));
addParameter(p, 'Subjects',   {}, @(c) iscell(c) || isstring(c));
parse(p, varargin{:});

runParams = struct();
runParams.schemaVersion = 3;
runParams.recordedAt = timestamp();
runParams.experimentName = cfg.experimentName;
runParams.experimentDescription = cfg.experimentDescription;
runParams.configFile = cfg.propsFile;
runParams.resultsDir = char(p.Results.ResultsDir);
if isfield(cfg, 'configExperimentName')
    runParams.configExperimentName = cfg.configExperimentName;
end

runParams.environment = environmentBlock();

runParams.project = struct( ...
    'projectRoot',    cfg.projectRoot, ...
    'dataRoot',       cfg.dataRoot, ...
    'resultsRoot',    cfg.resultsRoot, ...
    'connectomeFile', cfg.connectomeFile);

subjects = p.Results.Subjects;
if isempty(subjects)
    subjects = cfg.subjects;
end
runParams.cohort = struct( ...
    'subjects',      {cellstr(subjects)}, ...
    'nSubjects',     numel(subjects), ...
    'casePrefix',    cfg.casePrefix, ...
    'controlPrefix', cfg.controlPrefix);

runParams.normalisation = cfg.normalisation;
runParams.resultFilePrefix = experimentResultPrefix();

runParams.stabilityCentralities = struct( ...
    'normalisation', cfg.normalisation, ...
    'parkesC',       cfg.parkesC, ...
    'colScale',      cfg.colScale, ...
    'tvbPercentile', cfg.tvbPercentile, ...
    'x0Base',        cfg.x0Base, ...
    'x0Upper',       cfg.x0Upper, ...
    'x0Step',        cfg.x0Step, ...
    'bisectTol',     cfg.bisectTol, ...
    'maxK',          cfg.maxK, ...
    'tau0',          cfg.tau0, ...
    'topK',          cfg.topK, ...
    'discreteTime',  false, ...
    'fsolveTolFun',  1e-14, ...
    'fsolveTolX',    1e-14, ...
    'classicalCentralities', struct( ...
        'alphaPR', 0.85, 'alphaKZ', 0.5, 'ignoreSelfLoops', true));

runParams.atlas = struct( ...
    'leftPrefix',           cfg.leftPrefix, ...
    'rightPrefix',          cfg.rightPrefix, ...
    'trivialLabelSuffixes', {cfg.trivialLabelSuffixes});

runParams.comparison = struct( ...
    'corrType',        cfg.corrType, ...
    'alpha',           cfg.alpha, ...
    'zThreshold',      cfg.zThreshold, ...
    'lateralityBasis', cfg.lateralityBasis);

runParams.pipeline = struct( ...
    'runQc',                    cfg.pipelineRunQc, ...
    'runStabilityCentralities', cfg.pipelineRunStabilityCentralities, ...
    'runCentralityCorr',        cfg.pipelineRunCentralityCorr, ...
    'runCohortComparison',      cfg.pipelineRunCohortComparison, ...
    'runHemisphericAsymmetry',  cfg.pipelineRunHemisphericAsymmetry, ...
    'runDstCohort',             cfg.pipelineRunDstCohort, ...
    'runSummaryFigures',        cfg.pipelineRunSummaryFigures, ...
    'runReports',               cfg.pipelineRunReports, ...
    'stopOnError',              cfg.pipelineStopOnError);

runParams.output = struct( ...
    'saveResults', cfg.saveResults, ...
    'plot',        cfg.plot, ...
    'verbose',     cfg.verbose);

runParams.config = stripRawConfig(cfg);
end

%% ------------------------------------------------------------------
function runParams = buildFromOptions(opts, varargin)
p = inputParser;
addParameter(p, 'Analysis',   '', @(s) ischar(s) || isstring(s));
addParameter(p, 'ResultsDir', '', @(s) ischar(s) || isstring(s));
addParameter(p, 'Scheme',     '', @(s) ischar(s) || isstring(s));
addParameter(p, 'Subjects',   {}, @(c) iscell(c) || isstring(c));
parse(p, varargin{:});

runParams = struct();
runParams.schemaVersion = 3;
runParams.recordedAt = timestamp();
runParams.experimentName = '';
runParams.configFile = '';
runParams.analysis = char(p.Results.Analysis);

resultsDir = char(p.Results.ResultsDir);
if isempty(resultsDir) && isstruct(opts) && isfield(opts, 'ResultsDir')
    resultsDir = char(opts.ResultsDir);
end
runParams.resultsDir = resultsDir;

scheme = char(p.Results.Scheme);
if isempty(scheme) && isstruct(opts) && isfield(opts, 'Normalisation')
    scheme = char(opts.Normalisation);
end
runParams.normalisation = scheme;

runParams.environment = environmentBlock();
runParams.cohort = struct( ...
    'subjects',  {cellstr(p.Results.Subjects)}, ...
    'nSubjects', numel(p.Results.Subjects));

% Record the analysis options verbatim, minus the handles and the
% already-recorded run parameters, so nothing is silently lost.
runParams.options = scrubOptions(opts);
end

%% ------------------------------------------------------------------
function runParams = mergeParams(runParams, patch)
if nargin < 2 || isempty(patch)
    return;
end
if isempty(runParams)
    runParams = patch;
    return;
end
runParams = mergeStruct(runParams, patch);
end

%% ------------------------------------------------------------------
function saveParams(resultsDir, runParams)
if nargin < 2
    error('experimentRunParameters:SaveArgs', 'resultsDir and runParams are required.');
end
if ~isfield(runParams, 'experimentName') || isempty(runParams.experimentName)
    error('experimentRunParameters:NoExperimentName', ...
        'runParams.experimentName is required to name the provenance files.');
end

runParams.recordedAt = timestamp();
expName = char(runParams.experimentName);

matPath = experimentProvenanceFile(resultsDir, expName, 'parametersMat');
save(matPath, 'runParams');

jsonPath = experimentProvenanceFile(resultsDir, expName, 'parametersJson');
try
    jsonText = jsonencode(runParams);
    fid = fopen(jsonPath, 'w');
    if fid < 0
        warning('experimentRunParameters:JsonWrite', 'Could not write %s', jsonPath);
    else
        closer = onCleanup(@() fclose(fid));
        fprintf(fid, '%s', prettifyJson(jsonText));
    end
catch ME
    warning('experimentRunParameters:JsonEncode', ...
        'JSON export skipped: %s', ME.message);
end
end

%% ------------------------------------------------------------------
function env = environmentBlock()
env = struct();
env.matlabVersion = version;
env.computer = computer('arch');
env.hostname = '';
try
    if usejava('jvm')
        env.hostname = char(java.net.InetAddress.getLocalHost.getHostName);
    end
catch
    env.hostname = '';
end
if ispc
    env.user = getenv('USERNAME');
else
    env.user = getenv('USER');
end
env.toolkitRoot = fileparts(mfilename('fullpath'));
end

%% ------------------------------------------------------------------
function s = timestamp()
s = datestr(now, 'yyyy-mm-dd HH:MM:SS'); %#ok<DATST,TNOW1>
end

%% ------------------------------------------------------------------
function cfg = stripRawConfig(cfg)
% The raw key/value cell arrays are already captured by the .properties
% snapshot; keeping them here would only bloat every embedded copy.
if isfield(cfg, 'raw')
    cfg = rmfield(cfg, 'raw');
end
end

%% ------------------------------------------------------------------
function opts = scrubOptions(opts)
if ~isstruct(opts)
    opts = struct();
    return;
end
drop = {'RunParameters'};
for k = 1:numel(drop)
    if isfield(opts, drop{k})
        opts = rmfield(opts, drop{k});
    end
end
names = fieldnames(opts);
for k = 1:numel(names)
    if isa(opts.(names{k}), 'function_handle')
        opts.(names{k}) = func2str(opts.(names{k}));
    end
end
end

%% ------------------------------------------------------------------
function s = mergeStruct(s, t)
if ~isstruct(t)
    return;
end
names = fieldnames(t);
for i = 1:numel(names)
    fn = names{i};
    if isstruct(t.(fn)) && isfield(s, fn) && isstruct(s.(fn)) && ...
            isscalar(t.(fn)) && isscalar(s.(fn))
        s.(fn) = mergeStruct(s.(fn), t.(fn));
    else
        s.(fn) = t.(fn);
    end
end
end

%% ------------------------------------------------------------------
function txt = prettifyJson(jsonText)
% Break after commas that start a new key so the file diffs line by line.
txt = regexprep(jsonText, ',"', sprintf(',\n"'));
end
