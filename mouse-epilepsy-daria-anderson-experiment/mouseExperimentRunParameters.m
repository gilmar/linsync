function varargout = mouseExperimentRunParameters(action, varargin)
%MOUSEEXPERIMENTRUNPARAMETERS  Build and save provenance for experiment runs.
%
%   runParams = mouseExperimentRunParameters('buildFromConfig', cfg, ...)
%   runParams = mouseExperimentRunParameters('buildFromStabilityCentralities', opts, mouseId)
%   runParams = mouseExperimentRunParameters('merge', runParams, patchStruct)
%   mouseExperimentRunParameters('save', resultsDir, runParams)
%
%   Writes experiment_parameters_<experimentName>.{mat,json} under resultsDir.
%   The struct is also stored in run_manifest.mat and embedded
%   in per-mouse / summary / comparison .mat files when passed through the
%   pipeline.

switch lower(action)
    case 'buildfromconfig'
        varargout{1} = buildFromConfig(varargin{:});
    case 'buildfromstabilitycentralities'
        varargout{1} = buildFromStabilityCentralities(varargin{:});
    case 'buildfromcompareanderson'
        varargout{1} = buildFromCompareAnderson(varargin{:});
    case 'buildfromcomparecentrality'
        varargout{1} = buildFromCompareCentrality(varargin{:});
    case 'buildfromcomparedst'
        varargout{1} = buildFromCompareDst(varargin{:});
    case 'buildfromcomparelr'
        varargout{1} = buildFromCompareLR(varargin{:});
    case 'merge'
        varargout{1} = mergeParams(varargin{:});
    case 'save'
        saveParams(varargin{:});
    otherwise
        error('mouseExperimentRunParameters:BadAction', ...
            'Unknown action "%s".', action);
end
end

%% ------------------------------------------------------------------
function runParams = buildFromConfig(cfg, varargin)
p = inputParser;
addParameter(p, 'ResultsDir', '', @(s) ischar(s) || isstring(s));
addParameter(p, 'Mice', {{}}, @(c) iscell(c) || isstring(c));
parse(p, varargin{:});

runParams = struct();
runParams.schemaVersion = 2;
runParams.recordedAt = datestr(now, 'yyyy-mm-dd HH:MM:SS');
runParams.experimentName = cfg.experimentName;
if isfield(cfg, 'configExperimentName')
    runParams.configExperimentName = cfg.configExperimentName;
end
runParams.experimentDescription = cfg.experimentDescription;
runParams.configFile = cfg.propsFile;
runParams.resultsDir = char(p.Results.ResultsDir);

runParams.environment = environmentBlock();

runParams.cohort = struct();
if isempty(p.Results.Mice)
    runParams.cohort.mice = listAvailableMice();
else
    runParams.cohort.mice = cellstr(p.Results.Mice);
end
runParams.cohort.nMice = numel(runParams.cohort.mice);

runParams.normalisation = cfg.normalisation;
runParams.resultFilePrefix = mouseExperimentResultPrefix();

runParams.stabilityCentralities = struct( ...
    'normalisation', cfg.normalisation, ...
    'parkesC', cfg.parkesC, ...
    'colScale', cfg.colScale, ...
    'x0Base', cfg.x0Base, ...
    'x0Upper', cfg.x0Upper, ...
    'x0Step', cfg.x0Step, ...
    'bisectTol', cfg.bisectTol, ...
    'maxK', cfg.maxK, ...
    'tau0', cfg.tau0, ...
    'topK', cfg.topK, ...
    'discreteTime', false, ...
    'fsolveTolFun', 1e-14, ...
    'fsolveTolX', 1e-14, ...
    'classicalCentralities', struct('alphaPR', 0.85, 'alphaKZ', 0.5, 'ignoreSelfLoops', true));

runParams.comparison = struct( ...
    'corrType', cfg.corrType, ...
    'alpha', cfg.alpha, ...
    'zThreshold', cfg.zThreshold, ...
    'lateralityBasis', cfg.lateralityBasis, ...
    'compareCentralityReorder', 'cluster', ...
    'compareCentralityPerMouse', true);

runParams.pipeline = struct( ...
    'runHeatmap', cfg.pipelineRunHeatmap, ...
    'runStabilityCentralities', cfg.pipelineRunStabilityCentralities, ...
    'runCentralityCorr', cfg.pipelineRunCentralityCorr, ...
    'runAndersonVsArnold', cfg.pipelineRunAndersonVsArnold, ...
    'runLRAsymmetry', cfg.pipelineRunLRAsymmetry, ...
    'runDstCohort', cfg.pipelineRunDstCohort, ...
    'runReports', cfg.pipelineRunReports, ...
    'stopOnError', cfg.pipelineStopOnError);

runParams.output = struct( ...
    'saveResults', cfg.saveResults, ...
    'plot', cfg.plot, ...
    'verbose', cfg.verbose);

runParams.config = cfg;
end

%% ------------------------------------------------------------------
function runParams = buildFromStabilityCentralities(opts, mouseId)
% opts: inputParser Results from runMouseStabilityCentralities
runParams = struct();
runParams.schemaVersion = 2;
runParams.recordedAt = datestr(now, 'yyyy-mm-dd HH:MM:SS');
runParams.experimentName = '';
runParams.configFile = '';
runParams.resultsDir = '';
if isfield(opts, 'ResultsDir') && ~isempty(opts.ResultsDir)
    runParams.resultsDir = char(opts.ResultsDir);
end
runParams.environment = environmentBlock();
runParams.cohort = struct('mice', {mouseId}, 'nMice', 1);
runParams.normalisation = char(opts.Normalisation);
runParams.stabilityCentralities = struct( ...
    'normalisation', char(opts.Normalisation), ...
    'parkesC', opts.ParkesC, ...
    'colScale', opts.ColScale, ...
    'x0Base', opts.x0Base, ...
    'x0Upper', opts.x0Upper, ...
    'x0Step', opts.x0Step, ...
    'bisectTol', opts.BisectTol, ...
    'maxK', opts.MaxK, ...
    'tau0', opts.tau0, ...
    'topK', localTopK(opts), ...
    'discreteTime', false, ...
    'fsolveTolFun', 1e-14, ...
    'fsolveTolX', 1e-14, ...
    'mouseId', char(mouseId));
runParams.comparison = struct();
runParams.pipeline = struct();
runParams.output = struct( ...
    'saveResults', opts.SaveResults, ...
    'plot', opts.Plot, ...
    'verbose', opts.Verbose);
end

%% ------------------------------------------------------------------
function runParams = buildFromCompareAnderson(opts, scheme, resultsDir)
runParams = struct('schemaVersion', 1, 'recordedAt', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
runParams.normalisation = char(scheme);
runParams.resultsDir = char(resultsDir);
runParams.comparison = struct( ...
    'normalisation', char(scheme), ...
    'alpha', opts.Alpha, ...
    'zThreshold', opts.ZThreshold, ...
    'saveResults', opts.SaveResults);
end

%% ------------------------------------------------------------------
function runParams = buildFromCompareCentrality(opts, scheme, resultsDir)
runParams = struct('schemaVersion', 1, 'recordedAt', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
runParams.normalisation = char(scheme);
runParams.resultsDir = char(resultsDir);
runParams.comparison = struct( ...
    'normalisation', char(scheme), ...
    'corrType', char(opts.CorrType), ...
    'perMouse', opts.PerMouse, ...
    'reorder', char(opts.Reorder), ...
    'saveResults', opts.SaveResults);
end

%% ------------------------------------------------------------------
function runParams = buildFromCompareDst(opts, scheme, resultsDir)
runParams = struct('schemaVersion', 1, 'recordedAt', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
runParams.normalisation = char(scheme);
runParams.resultsDir = char(resultsDir);
runParams.comparison = struct( ...
    'normalisation', char(scheme), ...
    'saveResults', opts.SaveResults);
end

%% ------------------------------------------------------------------
function runParams = buildFromCompareLR(opts, scheme, resultsDir)
runParams = struct('schemaVersion', 1, 'recordedAt', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
runParams.normalisation = char(scheme);
runParams.resultsDir = char(resultsDir);
basis = 'signed';
if isfield(opts, 'LateralityBasis') && ~isempty(opts.LateralityBasis)
    basis = lower(strtrim(char(opts.LateralityBasis)));
end
runParams.comparison = struct( ...
    'normalisation', char(scheme), ...
    'alpha', opts.Alpha, ...
    'zThreshold', opts.ZThreshold, ...
    'lateralityBasis', basis, ...
    'saveResults', opts.SaveResults, ...
    'analysis', 'leftRightAsymmetry');
end

%% ------------------------------------------------------------------
function runParams = mergeParams(runParams, patch)
if nargin < 2 || isempty(patch)
    return;
end
runParams = mergeStruct(runParams, patch);
end

%% ------------------------------------------------------------------
function saveParams(resultsDir, runParams)
if nargin < 2
    error('mouseExperimentRunParameters:SaveArgs', 'resultsDir and runParams required.');
end
runParams.recordedAt = datestr(now, 'yyyy-mm-dd HH:MM:SS');
if ~isfield(runParams, 'experimentName') || isempty(runParams.experimentName)
    error('mouseExperimentRunParameters:NoExperimentName', ...
        'runParams.experimentName required to name provenance files.');
end
expName = char(runParams.experimentName);
matPath = mouseExperimentProvenanceFile(resultsDir, expName, 'parametersMat');
save(matPath, 'runParams');

jsonPath = mouseExperimentProvenanceFile(resultsDir, expName, 'parametersJson');
try
    jsonText = jsonencode(runParams);
    fid = fopen(jsonPath, 'w');
    if fid < 0
        warning('mouseExperimentRunParameters:JsonWrite', 'Could not write %s', jsonPath);
    else
        c = onCleanup(@() fclose(fid));
        fprintf(fid, '%s', prettifyJson(jsonText));
    end
catch ME
    warning('mouseExperimentRunParameters:JsonEncode', ...
        'JSON export skipped: %s', ME.message);
end
end

%% ------------------------------------------------------------------
function env = environmentBlock()
env = struct();
env.matlabVersion = version;
env.computer = computer('arch');
try
    env.hostname = char(java.net.InetAddress.getLocalHost.getHostName);
catch
    env.hostname = '';
end
if ispc
    env.user = getenv('USERNAME');
else
    env.user = getenv('USER');
end
env.toolkitRoot = fileparts(fileparts(mfilename('fullpath')));
env.experimentRoot = fileparts(mfilename('fullpath'));
end

%% ------------------------------------------------------------------
function s = mergeStruct(s, t)
if ~isstruct(t)
    return;
end
names = fieldnames(t);
for i = 1:numel(names)
    fn = names{i};
    if isstruct(t.(fn)) && isfield(s, fn) && isstruct(s.(fn))
        s.(fn) = mergeStruct(s.(fn), t.(fn));
    else
        s.(fn) = t.(fn);
    end
end
end

%% ------------------------------------------------------------------
function txt = prettifyJson(jsonText)
% Insert newlines after commas between top-level keys for readability.
txt = regexprep(jsonText, ',"', sprintf(',\n"'));
end

function v = localTopK(opts)
v = NaN;
if isfield(opts, 'TopK')
    v = opts.TopK;
end
end
