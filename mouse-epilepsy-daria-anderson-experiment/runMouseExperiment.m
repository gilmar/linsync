function manifest = runMouseExperiment(propsFile)
%RUNMOUSEEXPERIMENT  Run the full mouse experiment pipeline from a .properties file.
%
%   manifest = runMouseExperiment('configs/initial_column.properties')
%
%   Writes all outputs under results/<experiment.name>/ and saves
%   run_manifest.mat, experiment.properties, and experiment_parameters.{mat,json}.

if nargin < 1 || isempty(propsFile)
    error('runMouseExperiment:NoFile', ...
        'Provide a properties file, e.g. runMouseExperiment(''configs/initial_column.properties'')');
end

setupMousePaths();
cfg = loadMouseExperimentConfig(propsFile);
scheme = cfg.normalisation;

resultsDir = setupMousePaths('ExperimentName', cfg.experimentName);
copyfile(cfg.propsFile, fullfile(resultsDir, 'experiment.properties'));

manifest = struct();
manifest.experimentName = cfg.experimentName;
manifest.normalisation = scheme;
manifest.startedAt = datestr(now, 'yyyy-mm-dd HH:MM:SS');
manifest.matlabVersion = version;
manifest.configFile = cfg.propsFile;
manifest.steps = struct([]);

fprintf('\n=== Mouse experiment: %s  |  scheme=%s ===\n', cfg.experimentName, scheme);
fprintf('Results directory: %s\n\n', resultsDir);

mice = listAvailableMice();
if isempty(mice)
    error('runMouseExperiment:NoMice', 'No mouse CSVs found under data/.');
end

runParams = mouseExperimentRunParameters('buildFromConfig', cfg, ...
    'ResultsDir', resultsDir, 'Mice', mice);
mouseExperimentRunParameters('save', resultsDir, runParams);
manifest.runParameters = runParams;

commonArgs = {'ResultsDir', resultsDir, ...
    'RunParameters', runParams, ...
    'Normalisation', scheme, ...
    'ParkesC', cfg.parkesC, ...
    'ColScale', cfg.colScale, ...
    'x0Base', cfg.x0Base, ...
    'x0Upper', cfg.x0Upper, ...
    'x0Step', cfg.x0Step, ...
    'BisectTol', cfg.bisectTol, ...
    'MaxK', cfg.maxK, ...
    'tau0', cfg.tau0, ...
    'TopK', cfg.topK, ...
    'Plot', cfg.plot, ...
    'Verbose', cfg.verbose, ...
    'SaveResults', cfg.saveResults};

if cfg.pipelineRunHeatmap
    manifest.steps(end+1) = runPipelineStep(@() compareMouseHeatmap( ...
        'ResultsDir', resultsDir), 'compareMouseHeatmap', cfg.pipelineStopOnError);
end

if cfg.pipelineRunSection45
    manifest.steps(end+1) = runPipelineStep(@() runAllMiceSection45(commonArgs{:}), ...
        'runAllMiceSection45', cfg.pipelineStopOnError);
    manifest.steps(end+1) = runPipelineStep(@() assertCohortResultsComplete(resultsDir, scheme, mice), ...
        'postflight_section45', cfg.pipelineStopOnError);
end

if cfg.pipelineRunCentralityCorr
    manifest.steps(end+1) = runPipelineStep(@() compareCentralityMeasures( ...
        'ResultsDir', resultsDir, 'Normalisation', scheme, ...
        'CorrType', cfg.corrType, 'SaveResults', cfg.saveResults, ...
        'RunParameters', runParams), ...
        'compareCentralityMeasures', cfg.pipelineStopOnError);
end

if cfg.pipelineRunAndersonVsArnold
    manifest.steps(end+1) = runPipelineStep(@() compareAndersonVsArnold( ...
        'ResultsDir', resultsDir, 'Normalisation', scheme, ...
        'Alpha', cfg.alpha, 'ZThreshold', cfg.zThreshold, ...
        'SaveResults', cfg.saveResults, 'RunParameters', runParams), ...
        'compareAndersonVsArnold', cfg.pipelineStopOnError);
end

if cfg.pipelineRunDstCohort
    manifest.steps(end+1) = runPipelineStep(@() compareDstAcrossCohort( ...
        'ResultsDir', resultsDir, 'Normalisation', scheme, ...
        'SaveResults', cfg.saveResults, 'RunParameters', runParams), ...
        'compareDstAcrossCohort', cfg.pipelineStopOnError);
end

if cfg.pipelineRunReports
    logFile = fullfile(resultsDir, 'reports.log');
    manifest.steps(end+1) = runPipelineStep(@() runReports(logFile, resultsDir, scheme, cfg), ...
        'reports', cfg.pipelineStopOnError);
end

manifest.finishedAt = datestr(now, 'yyyy-mm-dd HH:MM:SS');
manifest.allSucceeded = all([manifest.steps.success]);
runParams.finishedAt = manifest.finishedAt;
runParams.pipelineSteps = manifest.steps;
manifest.runParameters = runParams;
mouseExperimentRunParameters('save', resultsDir, runParams);
save(fullfile(resultsDir, 'run_manifest.mat'), 'manifest', 'cfg', 'runParams');

if manifest.allSucceeded
    fprintf('\nExperiment %s finished successfully.\n', cfg.experimentName);
else
    fprintf('\nExperiment %s finished with failures (see run_manifest.mat).\n', cfg.experimentName);
end
end

%% ------------------------------------------------------------------
function step = runPipelineStep(fn, stepName, stopOnError)
step = struct('name', stepName, 'success', false, ...
    'durationSec', NaN, 'errorMessage', '');
t0 = tic;
try
    fn();
    step.success = true;
catch ME
    step.errorMessage = ME.message;
    if stopOnError
        rethrow(ME);
    else
        warning('runMouseExperiment:StepFailed', '%s failed: %s', stepName, ME.message);
    end
end
step.durationSec = toc(t0);
fprintf('  [%s] %s (%.1f s)\n', iif(step.success, 'OK', 'FAIL'), stepName, step.durationSec);
end

function s = iif(cond, a, b)
if cond, s = a; else, s = b; end
end

function assertCohortResultsComplete(resultsDir, scheme, mice)
pattern = sprintf('section45_*_%s_results.mat', scheme);
files = dir(fullfile(resultsDir, pattern));
if numel(files) < numel(mice)
    error('runMouseExperiment:IncompleteCohort', ...
        'Expected %d section45 results for scheme "%s", found %d in %s.', ...
        numel(mice), scheme, numel(files), resultsDir);
end
end

function runReports(logFile, resultsDir, scheme, cfg)
diary(logFile);
cleanup = onCleanup(@() diary('off'));
fprintf('=== reportTopNodes ===\n');
reportTopNodes('ResultsDir', resultsDir, 'Normalisation', scheme, 'N', cfg.topK);
fprintf('\n=== reportAndersonOutliers ===\n');
reportAndersonOutliers('ResultsDir', resultsDir, 'Normalisation', scheme, ...
    'ZThreshold', cfg.zThreshold);
end
