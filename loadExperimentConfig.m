function cfg = loadExperimentConfig(propsFile, varargin)
%LOADEXPERIMENTCONFIG  Parse a connectome-experiment .properties file.
%
%   cfg = loadExperimentConfig('configs/baseline_column.properties')
%   cfg = loadExperimentConfig(propsFile, 'ProjectRoot', '/path/to/study')
%
%   One experiment = one config file = one normalisation scheme + one
%   parameter set + one results folder. Driving runs from a plain text file
%   (rather than edited-in-place script constants) is what makes a run
%   reproducible: the file is snapshotted into the results folder, so the
%   exact inputs that produced a figure stay recoverable.
%
%   Syntax: Java-style `key=value` (or `key: value`), one per line, with
%   `#` or `!` starting a comment. Parsed in pure MATLAB, so it works under
%   -nojvm. Unknown keys are kept in cfg.raw and otherwise ignored.
%
%   The project root (which holds configs/, data/ and results/) defaults to
%   the parent of the config file's folder when the config sits in a
%   configs/ directory, and to the config file's own folder otherwise.
%
%   See configs/experiment.template.properties for the documented key list,
%   and RUNCONNECTOMEEXPERIMENT for what consumes each key.

if nargin < 1 || isempty(propsFile)
    error('loadExperimentConfig:NoFile', 'propsFile is required.');
end

p = inputParser;
addParameter(p, 'ProjectRoot', '', @(s) ischar(s) || isstring(s));
parse(p, varargin{:});

propsFile = resolveConfigPath(char(propsFile));
raw = readPropertiesFile(propsFile);

cfg = struct();
cfg.propsFile = propsFile;
cfg.raw = raw;

% ---- Project layout -------------------------------------------------
projectRoot = char(p.Results.ProjectRoot);
if isempty(projectRoot)
    projectRoot = getProp(raw, {'project.root'}, '');
end
if isempty(projectRoot)
    projectRoot = defaultProjectRoot(propsFile);
end
cfg.projectRoot    = resolveRelative(projectRoot, fileparts(propsFile));
cfg.dataRoot       = resolveRelative(getProp(raw, {'data.root'}, 'data'), cfg.projectRoot);
cfg.resultsRoot    = resolveRelative(getProp(raw, {'results.root'}, 'results'), cfg.projectRoot);
cfg.connectomeFile = getProp(raw, {'data.file', 'connectome.file'}, 'connectome.csv');
cfg.subjects       = parseList(getProp(raw, {'subjects', 'cohort.subjects'}, ''));

% ---- Identity -------------------------------------------------------
cfg.experimentName = getProp(raw, {'experiment.name'}, '');
cfg.experimentDescription = getProp(raw, {'experiment.description'}, '');
if isempty(cfg.experimentName)
    [~, cfg.experimentName] = fileparts(propsFile);
end

% ---- Model and normalisation ---------------------------------------
schemeInfo = normalisationSchemeInfo(getProp(raw, {'normalisation'}, 'column'));
cfg.normalisation = schemeInfo.scheme;
cfg.schemeInfo    = schemeInfo;
cfg.parkesC       = parseNum(getProp(raw, {'parkes.c'}, '1.0'), 'parkes.c');
cfg.colScale      = parseNum(getProp(raw, {'col.scale'}, '0.95'), 'col.scale');
cfg.tvbPercentile = parseNum(getProp(raw, {'tvb.percentile'}, '95'), 'tvb.percentile');
cfg.x0Base        = parseNum(getProp(raw, {'x0.base'}, '-2.3'), 'x0.base');
cfg.x0Upper       = parseNum(getProp(raw, {'x0.upper'}, '-1.0'), 'x0.upper');
cfg.x0Step        = parseNum(getProp(raw, {'x0.step'}, '0.01'), 'x0.step');
cfg.bisectTol     = parseNum(getProp(raw, {'bisect.tol'}, '1e-4'), 'bisect.tol');
cfg.maxK          = parseNum(getProp(raw, {'max.k'}, '1e8'), 'max.k');
cfg.tau0          = parseNum(getProp(raw, {'tau0'}, '6667'), 'tau0');
cfg.topK          = parseNum(getProp(raw, {'top.k'}, '10'), 'top.k');

% ---- Atlas conventions ---------------------------------------------
cfg.leftPrefix  = getProp(raw, {'hemisphere.left.prefix'}, 'L-');
cfg.rightPrefix = getProp(raw, {'hemisphere.right.prefix'}, 'R-');
cfg.trivialLabelSuffixes = parseList(getProp(raw, {'trivial.label.suffixes'}, ''));

% ---- Cohort ---------------------------------------------------------
cfg.casePrefix    = getProp(raw, {'cohort.case.prefix', 'case.prefix'}, 'Case');
cfg.controlPrefix = getProp(raw, {'cohort.control.prefix', 'control.prefix'}, 'Control');

% ---- Statistics -----------------------------------------------------
cfg.corrType   = getProp(raw, {'corr.type'}, 'Spearman');
cfg.alpha      = parseNum(getProp(raw, {'alpha'}, '0.05'), 'alpha');
cfg.zThreshold = parseNum(getProp(raw, {'z.threshold'}, '2'), 'z.threshold');
cfg.lateralityBasis = normaliseLateralityBasis(getProp(raw, {'laterality.basis'}, 'signed'));

% ---- Pipeline -------------------------------------------------------
cfg.pipelineRunQc = parseBool(getProp(raw, ...
    {'pipeline.runQc', 'pipeline.runHeatmap'}, 'false'));
cfg.pipelineRunStabilityCentralities = parseBool(getProp(raw, ...
    {'pipeline.runStabilityCentralities'}, 'true'));
cfg.pipelineRunCentralityCorr = parseBool(getProp(raw, ...
    {'pipeline.runCentralityCorr'}, 'true'));
cfg.pipelineRunCohortComparison = parseBool(getProp(raw, ...
    {'pipeline.runCohortComparison', 'pipeline.runCaseVsControl'}, 'true'));
cfg.pipelineRunHemisphericAsymmetry = parseBool(getProp(raw, ...
    {'pipeline.runHemisphericAsymmetry', 'pipeline.runLRAsymmetry'}, 'true'));
cfg.pipelineRunDstCohort = parseBool(getProp(raw, ...
    {'pipeline.runDstCohort'}, 'true'));
cfg.pipelineRunSummaryFigures = parseBool(getProp(raw, ...
    {'pipeline.runSummaryFigures'}, 'true'));
cfg.pipelineRunReports = parseBool(getProp(raw, ...
    {'pipeline.runReports'}, 'true'));
cfg.pipelineStopOnError = parseBool(getProp(raw, ...
    {'pipeline.stopOnError'}, 'false'));

% ---- Output ---------------------------------------------------------
cfg.saveResults = parseBool(getProp(raw, {'save.results'}, 'true'));
cfg.plot        = parseBool(getProp(raw, {'plot'}, 'true'));
cfg.verbose     = parseBool(getProp(raw, {'verbose'}, 'true'));

validateConfig(cfg);
end

%% ------------------------------------------------------------------
function propsFile = resolveConfigPath(propsFile)
if isfile(propsFile)
    propsFile = fullPathOf(propsFile);
    return;
end
candidate = fullfile(pwd, propsFile);
if isfile(candidate)
    propsFile = fullPathOf(candidate);
    return;
end
error('loadExperimentConfig:FileNotFound', 'Properties file not found: %s', propsFile);
end

function p = fullPathOf(f)
d = dir(f);
p = fullfile(d(1).folder, d(1).name);
end

%% ------------------------------------------------------------------
function root = defaultProjectRoot(propsFile)
configDir = fileparts(propsFile);
[parent, leaf] = fileparts(configDir);
if strcmpi(leaf, 'configs') || strcmpi(leaf, 'config')
    root = parent;
else
    root = configDir;
end
end

%% ------------------------------------------------------------------
function out = resolveRelative(pathStr, baseDir)
pathStr = char(pathStr);
if isempty(pathStr)
    out = char(baseDir);
elseif isAbsolutePath(pathStr)
    out = pathStr;
else
    out = fullfile(char(baseDir), pathStr);
end
end

%% ------------------------------------------------------------------
function raw = readPropertiesFile(propsFile)
%READPROPERTIESFILE  Minimal Java-properties reader (no JVM required).
raw = struct('keys', {{}}, 'values', {{}});
text = fileread(propsFile);
lines = regexp(text, '\r\n|\r|\n', 'split');
for k = 1:numel(lines)
    line = strtrim(lines{k});
    if isempty(line) || line(1) == '#' || line(1) == '!'
        continue;
    end
    sep = regexp(line, '[=:]', 'once');
    if isempty(sep)
        continue;
    end
    key = strtrim(line(1:sep - 1));
    value = strtrim(line(sep + 1:end));
    if isempty(key)
        continue;
    end
    existing = find(strcmp(raw.keys, key), 1);
    if isempty(existing)
        raw.keys{end+1} = key;
        raw.values{end+1} = value;
    else
        raw.values{existing} = value;
    end
end
end

%% ------------------------------------------------------------------
function v = getProp(raw, keyNames, defaultValue)
v = defaultValue;
for i = 1:numel(keyNames)
    idx = find(strcmpi(raw.keys, keyNames{i}), 1);
    if ~isempty(idx) && ~isempty(strtrim(raw.values{idx}))
        v = strtrim(raw.values{idx});
        return;
    end
end
end

%% ------------------------------------------------------------------
function tf = parseBool(s)
tf = ismember(lower(strtrim(char(s))), {'true', '1', 'yes', 'on'});
end

%% ------------------------------------------------------------------
function v = parseNum(s, keyName)
v = str2double(strtrim(char(s)));
if isnan(v)
    error('loadExperimentConfig:BadNumeric', ...
        'Key "%s" must be numeric (got "%s").', keyName, s);
end
end

%% ------------------------------------------------------------------
function items = parseList(s)
s = strtrim(char(s));
if isempty(s)
    items = {};
    return;
end
parts = strtrim(strsplit(s, {',', ';'}));
items = parts(~cellfun(@isempty, parts));
items = items(:)';
end

%% ------------------------------------------------------------------
function basis = normaliseLateralityBasis(basis)
basis = lower(strtrim(char(basis)));
switch basis
    case {'signed', 'raw', 'diff'}
        basis = 'signed';
    case {'norm', 'normalised', 'normalized'}
        basis = 'norm';
    otherwise
        error('loadExperimentConfig:BadLateralityBasis', ...
            'laterality.basis must be "signed" or "norm" (got "%s").', basis);
end
end

%% ------------------------------------------------------------------
function validateConfig(cfg)
if cfg.alpha <= 0 || cfg.alpha >= 1
    error('loadExperimentConfig:BadAlpha', 'alpha must be in (0, 1).');
end
if cfg.colScale <= 0
    error('loadExperimentConfig:BadColScale', 'col.scale must be positive.');
end
if cfg.tvbPercentile <= 0 || cfg.tvbPercentile > 100
    error('loadExperimentConfig:BadTvbPercentile', 'tvb.percentile must be in (0, 100].');
end
if cfg.x0Upper <= cfg.x0Base
    error('loadExperimentConfig:BadX0Range', ...
        'x0.upper (%g) must exceed x0.base (%g).', cfg.x0Upper, cfg.x0Base);
end
if isempty(cfg.casePrefix) || isempty(cfg.controlPrefix)
    error('loadExperimentConfig:BadCohortPrefix', ...
        'cohort.case.prefix and cohort.control.prefix must be non-empty.');
end
end
