function cfg = loadMouseExperimentConfig(propsFile)
%LOADMOUSEEXPERIMENTCONFIG  Parse a mouse experiment .properties file.
%
%   cfg = loadMouseExperimentConfig('configs/initial_column.properties')
%
%   Returns a struct with typed fields used by runMouseExperiment.

if nargin < 1 || isempty(propsFile)
    error('loadMouseExperimentConfig:NoFile', 'propsFile is required.');
end

experimentRoot = fileparts(mfilename('fullpath'));
if ~isfile(propsFile)
    candidate = fullfile(experimentRoot, propsFile);
    if isfile(candidate)
        propsFile = candidate;
    else
        error('loadMouseExperimentConfig:FileNotFound', ...
            'Properties file not found: %s', propsFile);
    end
end
propsFile = char(java.io.File(propsFile).getCanonicalPath());

props = java.util.Properties();
fis = java.io.FileInputStream(propsFile);
try
    props.load(fis);
finally
    fis.close();
end

keys = cell(props.stringPropertyNames.toArray);
raw = struct();
for k = 1:numel(keys)
    key = char(keys{k});
    raw.(matlab.lang.makeValidName(strrep(key, '.', '_'))) = char(props.getProperty(keys{k}));
end

cfg = struct();
cfg.propsFile = propsFile;
cfg.experimentName = getProp(raw, {'experiment_name'}, '');
cfg.experimentDescription = getProp(raw, {'experiment_description'}, '');
cfg.normalisation = normalizeMouseNormalisation( ...
    getProp(raw, {'normalisation'}, 'column'));
cfg.parkesC    = str2double(getProp(raw, {'parkes_c'}, '1.0'));
cfg.colScale   = str2double(getProp(raw, {'col_scale'}, '0.95'));
cfg.x0Base     = str2double(getProp(raw, {'x0_base'}, '-2.3'));
cfg.x0Upper    = str2double(getProp(raw, {'x0_upper'}, '-1.0'));
cfg.x0Step     = str2double(getProp(raw, {'x0_step'}, '0.01'));
cfg.bisectTol  = str2double(getProp(raw, {'bisect_tol'}, '1e-4'));
cfg.maxK       = str2double(getProp(raw, {'max_k'}, '1e8'));
cfg.tau0       = str2double(getProp(raw, {'tau0'}, '6667'));
cfg.topK       = str2double(getProp(raw, {'top_k'}, '10'));
cfg.corrType   = getProp(raw, {'corr_type'}, 'Spearman');
cfg.alpha      = str2double(getProp(raw, {'alpha'}, '0.05'));
cfg.zThreshold = str2double(getProp(raw, {'z_threshold'}, '2'));
cfg.lateralityBasis = normalizeLateralityBasis( ...
    getProp(raw, {'laterality_basis'}, 'signed'));

cfg.pipelineRunHeatmap          = parseBool(getProp(raw, {'pipeline_runHeatmap'}, 'false'));
cfg.pipelineRunStabilityCentralities = parseBool(getProp(raw, ...
    {'pipeline_runStabilityCentralities', 'pipeline_runSection45'}, 'true'));
cfg.pipelineRunCentralityCorr   = parseBool(getProp(raw, {'pipeline_runCentralityCorr'}, 'true'));
cfg.pipelineRunAndersonVsArnold = parseBool(getProp(raw, {'pipeline_runAndersonVsArnold'}, 'true'));
cfg.pipelineRunLRAsymmetry      = parseBool(getProp(raw, {'pipeline_runLRAsymmetry'}, 'true'));
cfg.pipelineRunDstCohort        = parseBool(getProp(raw, {'pipeline_runDstCohort'}, 'true'));
cfg.pipelineRunReports          = parseBool(getProp(raw, {'pipeline_runReports'}, 'true'));
cfg.pipelineStopOnError         = parseBool(getProp(raw, {'pipeline_stopOnError'}, 'false'));

cfg.saveResults = parseBool(getProp(raw, {'save_results'}, 'true'));
cfg.plot        = parseBool(getProp(raw, {'plot'}, 'true'));
cfg.verbose     = parseBool(getProp(raw, {'verbose'}, 'true'));

if isempty(cfg.experimentName)
    [~, cfg.experimentName] = fileparts(propsFile);
end

validateConfig(cfg);
end

%% ------------------------------------------------------------------
function v = getProp(raw, fieldNames, defaultVal)
v = defaultVal;
for i = 1:numel(fieldNames)
    fn = fieldNames{i};
    if isfield(raw, fn) && ~isempty(strtrim(raw.(fn)))
        v = strtrim(raw.(fn));
        return;
    end
end
end

function tf = parseBool(s)
s = lower(strtrim(char(s)));
tf = ismember(s, {'true', '1', 'yes', 'on'});
end

function scheme = normalizeMouseNormalisation(scheme)
scheme = lower(strtrim(char(scheme)));
switch scheme
    case {'col', 'colnorm'}
        scheme = 'column';
    case {'tvb', 'none', 'parkes', 'column'}
        % ok
    otherwise
        error('loadMouseExperimentConfig:BadScheme', ...
            'Unknown normalisation "%s". Use tvb, none, parkes, or column.', scheme);
end
end

function basis = normalizeLateralityBasis(basis)
basis = lower(strtrim(char(basis)));
switch basis
    case {'signed', 'raw', 'diff'}
        basis = 'signed';
    case {'norm', 'normalised', 'normalized'}
        basis = 'norm';
    otherwise
        error('loadMouseExperimentConfig:BadLateralityBasis', ...
            'Unknown laterality.basis "%s". Use signed or norm.', basis);
end
end

function validateConfig(cfg)
if isnan(cfg.parkesC) || isnan(cfg.colScale)
    error('loadMouseExperimentConfig:BadNumeric', 'Invalid parkes.c or col.scale.');
end
if cfg.alpha <= 0 || cfg.alpha >= 1
    error('loadMouseExperimentConfig:BadAlpha', 'alpha must be in (0, 1).');
end
end
