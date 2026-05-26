function manifest = runMouseExperiment(propsFile)
%RUNMOUSEEXPERIMENT  Run the full mouse experiment pipeline from a .properties file.
%
%   manifest = runMouseExperiment('configs/initial_column.properties')
%
%   Writes all outputs under results/<experiment.name>_<yyyy-mm-dd_HHMM>/
%   (config experiment.name + run timestamp) and saves run_manifest.mat,
%   experiment_<name>.properties, experiment_parameters_<name>.{mat,json},
%   and reports_<name>.log (name = results subfolder).

if nargin < 1 || isempty(propsFile)
    error('runMouseExperiment:NoFile', ...
        'Provide a properties file, e.g. runMouseExperiment(''configs/initial_column.properties'')');
end

setupMousePaths();
cfg = loadMouseExperimentConfig(propsFile);
scheme = cfg.normalisation;

configExperimentName = cfg.experimentName;
[cfg.experimentName, runStamp] = mouseExperimentFolderName(configExperimentName);

resultsDir = setupMousePaths('ExperimentName', cfg.experimentName);
propsFile = mouseExperimentProvenanceFile(resultsDir, cfg.experimentName, 'properties');
copyfile(cfg.propsFile, propsFile);
appendExperimentResultsFolderNote(propsFile, cfg.experimentName, runStamp);

manifest = struct();
manifest.configExperimentName = configExperimentName;
manifest.experimentName = cfg.experimentName;
manifest.runStamp = runStamp;
manifest.normalisation = scheme;
manifest.startedAt = datestr(now, 'yyyy-mm-dd HH:MM:SS');
manifest.matlabVersion = version;
manifest.configFile = cfg.propsFile;
manifest.steps = emptyPipelineSteps();

fprintf('\n=== Mouse experiment: %s  |  scheme=%s ===\n', cfg.experimentName, scheme);
fprintf('Results directory: %s\n\n', resultsDir);

mice = listAvailableMice();
if isempty(mice)
    error('runMouseExperiment:NoMice', 'No mouse CSVs found under data/.');
end
fprintf('Cohort: %d mouse(s) with connectome CSVs: %s\n', numel(mice), strjoin(mice, ', '));
if numel(mice) < 5
    warning('runMouseExperiment:PartialCohort', ...
        ['Only %d mice with connectome CSVs (README lists 5: Anderson_1–2, Arnold_3–5; ' ...
         'Arnold_2 is skipped; Arnold_1 is optional if present).'], numel(mice));
end
hasAnderson = any(startsWith(string(mice), 'Anderson'));
hasArnold   = any(startsWith(string(mice), 'Arnold'));
canCompareStrains = hasAnderson && hasArnold;

cfg.configExperimentName = configExperimentName;
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

cohortReady = ~cfg.pipelineRunStabilityCentralities;
if cfg.pipelineRunStabilityCentralities
    stepSc = runPipelineStep(@() runAllMiceStabilityCentralities(commonArgs{:}), ...
        'runAllMiceStabilityCentralities', cfg.pipelineStopOnError);
    manifest.steps(end+1) = stepSc;
    stepPf = runPipelineStep(@() assertCohortResultsComplete(resultsDir, scheme, mice), ...
        'postflight_stabilityCentralities', cfg.pipelineStopOnError);
    manifest.steps(end+1) = stepPf;
    cohortReady = stepSc.success && stepPf.success;
end

if cfg.pipelineRunCentralityCorr
if cohortReady
    manifest.steps(end+1) = runPipelineStep(@() compareCentralityMeasures( ...
        'ResultsDir', resultsDir, 'Normalisation', scheme, ...
        'CorrType', cfg.corrType, 'SaveResults', cfg.saveResults, ...
        'RunParameters', runParams), ...
        'compareCentralityMeasures', cfg.pipelineStopOnError);
else
    manifest.steps(end+1) = skippedPipelineStep('compareCentralityMeasures', ...
        'per-mouse cohort incomplete');
end
end

andersonOut = {[]};
if cfg.pipelineRunAndersonVsArnold
if cohortReady && canCompareStrains
    manifest.steps(end+1) = runPipelineStep(@() captureAndersonSummary(), ...
        'compareAndersonVsArnold', cfg.pipelineStopOnError);
elseif ~cohortReady
    manifest.steps(end+1) = skippedPipelineStep('compareAndersonVsArnold', ...
        'per-mouse cohort incomplete');
else
    manifest.steps(end+1) = skippedPipelineStep('compareAndersonVsArnold', ...
        'no Anderson and/or Arnold mice in data/', true);
end
end

lrOut = {[]};
if cfg.pipelineRunLRAsymmetry
if cohortReady && canCompareStrains
    manifest.steps(end+1) = runPipelineStep(@() captureLRSummary(), ...
        'compareLeftRightAsymmetry', cfg.pipelineStopOnError);
elseif ~cohortReady
    manifest.steps(end+1) = skippedPipelineStep('compareLeftRightAsymmetry', ...
        'per-mouse cohort incomplete');
else
    manifest.steps(end+1) = skippedPipelineStep('compareLeftRightAsymmetry', ...
        'no Anderson and/or Arnold mice in data/', true);
end
end

andersonSummary = andersonOut{1};
lrSummary = lrOut{1};
if strcmp(scheme, 'column') && ~isempty(andersonSummary) && ~isempty(lrSummary)
    manifest.steps(end+1) = runPipelineStep(@() renderComparisonFigures( ...
        andersonSummary, lrSummary, 'OutputDir', fullfile(resultsDir, 'comparison_figures')), ...
        'renderComparisonFigures', cfg.pipelineStopOnError);
end

if cfg.pipelineRunDstCohort
if cohortReady
    manifest.steps(end+1) = runPipelineStep(@() compareDstAcrossCohort( ...
        'ResultsDir', resultsDir, 'Normalisation', scheme, ...
        'SaveResults', cfg.saveResults, 'RunParameters', runParams), ...
        'compareDstAcrossCohort', cfg.pipelineStopOnError);
else
    manifest.steps(end+1) = skippedPipelineStep('compareDstAcrossCohort', ...
        'per-mouse cohort incomplete');
end
end

if cfg.pipelineRunReports
if cohortReady
    logFile = mouseExperimentProvenanceFile(resultsDir, cfg.experimentName, 'reports');
    manifest.steps(end+1) = runPipelineStep(@() runReports(logFile, resultsDir, scheme, cfg, mice), ...
        'reports', cfg.pipelineStopOnError);
else
    manifest.steps(end+1) = skippedPipelineStep('reports', 'per-mouse cohort incomplete');
end
end

manifest.finishedAt = datestr(now, 'yyyy-mm-dd HH:MM:SS');
runParams.finishedAt = manifest.finishedAt;
runParams.pipelineSteps = manifest.steps;
manifest.runParameters = runParams;
mouseExperimentRunParameters('save', resultsDir, runParams);
save(fullfile(resultsDir, 'run_manifest.mat'), 'manifest', 'cfg', 'runParams');

if cohortReady
    manifest.steps(end+1) = runPipelineStep(@() assertMouseExperimentOutputs( ...
        resultsDir, scheme, cfg, mice), 'assertExpectedOutputs', cfg.pipelineStopOnError);
    manifest.allSucceeded = all([manifest.steps.success]);
    runParams.pipelineSteps = manifest.steps;
    save(fullfile(resultsDir, 'run_manifest.mat'), 'manifest', 'cfg', 'runParams');
else
    manifest.allSucceeded = all([manifest.steps.success]);
end

if manifest.allSucceeded
    fprintf('\nExperiment %s finished successfully.\n', cfg.experimentName);
else
    fprintf('\nExperiment %s finished with failures (see run_manifest.mat).\n', cfg.experimentName);
end

    function captureAndersonSummary()
        andersonOut{1} = compareAndersonVsArnold( ...
            'ResultsDir', resultsDir, 'Normalisation', scheme, ...
            'Alpha', cfg.alpha, 'ZThreshold', cfg.zThreshold, ...
            'SaveResults', cfg.saveResults, 'RunParameters', runParams);
    end

    function captureLRSummary()
        lrOut{1} = compareLeftRightAsymmetry( ...
            'ResultsDir', resultsDir, 'Normalisation', scheme, ...
            'Alpha', cfg.alpha, 'ZThreshold', cfg.zThreshold, ...
            'LateralityBasis', cfg.lateralityBasis, ...
            'SaveResults', cfg.saveResults, 'RunParameters', runParams);
    end

end

%% ------------------------------------------------------------------
function step = skippedPipelineStep(stepName, reason, countAsSuccess)
if nargin < 3
    countAsSuccess = false;
end
step = struct('name', stepName, 'success', logical(countAsSuccess), ...
    'durationSec', 0, 'errorMessage', ['skipped: ' reason]);
tag = iif(countAsSuccess, 'OK', 'SKIP');
fprintf('  [%s] %s (%s)\n', tag, stepName, reason);
end

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
files = listPerMouseResultFiles(resultsDir, scheme);
if numel(files) < numel(mice)
    error('runMouseExperiment:IncompleteCohort', ...
        'Expected %d stabilityCentralities results for scheme "%s", found %d in %s.', ...
        numel(mice), scheme, numel(files), resultsDir);
end
end

function appendExperimentResultsFolderNote(noteFile, folderName, runStamp)
%APPENDEXPERIMENTRESULTSFOLDERNOTE  Record actual results folder in snapshot.
fid = fopen(noteFile, 'a');
if fid < 0
    return;
end
c = onCleanup(@() fclose(fid));
fprintf(fid, '\n# Appended at run start\n');
fprintf(fid, 'experiment.resultsFolder=%s\n', folderName);
fprintf(fid, 'experiment.runStamp=%s\n', runStamp);
end

function runReports(logFile, resultsDir, scheme, cfg, mice)
diary(logFile);
cleanup = onCleanup(@() diary('off'));
fprintf('=== reportTopNodes ===\n');
reportTopNodes('ResultsDir', resultsDir, 'Normalisation', scheme, 'N', cfg.topK);
if any(startsWith(string(mice), 'Anderson'))
    fprintf('\n=== reportAndersonOutliers ===\n');
    try
        reportAndersonOutliers('ResultsDir', resultsDir, 'Normalisation', scheme, ...
            'ZThreshold', cfg.zThreshold);
    catch ME
        fprintf('reportAndersonOutliers skipped: %s\n', ME.message);
    end
else
    fprintf('\n=== reportAndersonOutliers === skipped (no Anderson mice in cohort)\n');
end
end

function steps = emptyPipelineSteps()
steps = struct('name', {}, 'success', {}, 'durationSec', {}, 'errorMessage', {});
end
