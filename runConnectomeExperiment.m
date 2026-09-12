function manifest = runConnectomeExperiment(propsFile, varargin)
%RUNCONNECTOMEEXPERIMENT  Run the whole connectome stability pipeline from a config.
%
%   manifest = runConnectomeExperiment('configs/baseline_column.properties')
%   manifest = runConnectomeExperiment(propsFile, 'ProjectRoot', studyRoot)
%
%   One call takes a study from connectome files to figures, statistics and
%   a provenance record. Everything lands under
%       <results root>/<experiment.name>_<yyyy-mm-dd_HHMM>/
%   so re-running a config never overwrites an earlier run.
%
%   Steps, in order (each can be switched off with a pipeline.run* key):
%     1. QC              -- plotConnectomeHeatmaps over the raw inputs
%     2. Stability       -- runCohortStabilityCentralities for every subject
%     3. Postflight      -- assert every subject produced a result
%     4. Centrality corr -- compareCentralityMeasures
%     5. Cohort compare  -- compareCohortGroups (case vs control)
%     6. Laterality      -- compareHemisphericAsymmetry
%     7. Summary figures -- renderCohortComparisonFigures
%     8. Network D_st    -- compareDstAcrossCohort
%     9. Reports         -- reportTopNodes and reportCohortOutliers, tee'd
%                           into reports_<run>.log
%
%   A failing step is recorded and the run continues, unless
%   pipeline.stopOnError is set. The comparison steps are skipped outright
%   when the per-subject stage did not produce a complete cohort: comparing
%   against a control mean computed from a partial cohort would produce
%   confident, wrong numbers. Every skip and failure is written to
%   run_manifest.mat with its reason.
%
%   Name-value options
%     'ProjectRoot' : override the project root from the config
%
%   Returns the manifest: per-step name, success, duration and any error.
%
%   See also LOADEXPERIMENTCONFIG, NEWCONNECTOMEEXPERIMENT,
%   RUNALLCONNECTOMEEXPERIMENTS, EXPERIMENTRUNPARAMETERS.

if nargin < 1 || isempty(propsFile)
    error('runConnectomeExperiment:NoFile', ...
        ['Provide a properties file, e.g. ' ...
         'runConnectomeExperiment(''configs/baseline_column.properties'')']);
end

p = inputParser;
addParameter(p, 'ProjectRoot', '', @(s) ischar(s) || isstring(s));
parse(p, varargin{:});

cfg = loadExperimentConfig(propsFile, 'ProjectRoot', p.Results.ProjectRoot);
setupConnectomePaths('ProjectRoot', cfg.projectRoot, 'Quiet', true);

scheme = cfg.normalisation;
configExperimentName = cfg.experimentName;
[cfg.experimentName, runStamp] = experimentFolderName(configExperimentName);
cfg.configExperimentName = configExperimentName;

resultsDir = fullfile(cfg.resultsRoot, cfg.experimentName);
experimentResultsLayout(resultsDir, true);

% Snapshot the config as run, and note where the outputs went.
snapshotPath = experimentProvenanceFile(resultsDir, cfg.experimentName, 'properties');
copyfile(cfg.propsFile, snapshotPath);
appendRunNote(snapshotPath, cfg.experimentName, runStamp);

manifest = struct();
manifest.configExperimentName = configExperimentName;
manifest.experimentName = cfg.experimentName;
manifest.runStamp = runStamp;
manifest.normalisation = scheme;
manifest.startedAt = timestamp();
manifest.matlabVersion = version;
manifest.configFile = cfg.propsFile;
manifest.resultsDir = resultsDir;
manifest.steps = emptySteps();

fprintf('\n=== Connectome experiment: %s  |  scheme = %s ===\n', cfg.experimentName, scheme);
fprintf('Project : %s\n', cfg.projectRoot);
fprintf('Data    : %s\n', cfg.dataRoot);
fprintf('Results : %s\n\n', resultsDir);

%% ---- Cohort discovery ------------------------------------------------
subjects = cfg.subjects;
if isempty(subjects)
    subjects = listConnectomeSubjects(cfg.dataRoot, cfg.connectomeFile);
end
if isempty(subjects)
    error('runConnectomeExperiment:NoSubjects', ...
        'No subjects with a "%s" file under %s.', cfg.connectomeFile, cfg.dataRoot);
end
subjects = cellstr(subjects);
fprintf('Cohort: %d subject(s): %s\n', numel(subjects), strjoin(subjects, ', '));

[isCase, isControl] = cohortGroupMask(subjects, cfg.casePrefix, cfg.controlPrefix);
canCompareGroups = any(isCase) && any(isControl);
if ~canCompareGroups
    fprintf(['Note: case prefix "%s" and control prefix "%s" do not both match a ' ...
             'subject, so the group comparisons will be skipped.\n'], ...
        cfg.casePrefix, cfg.controlPrefix);
end

%% ---- Provenance, written before any analysis runs --------------------
runParams = experimentRunParameters('fromConfig', cfg, ...
    'ResultsDir', resultsDir, 'Subjects', subjects);
experimentRunParameters('save', resultsDir, runParams);
manifest.runParameters = runParams;

stabilityArgs = { ...
    'DataRoot',        cfg.dataRoot, ...
    'ConnectomeFile',  cfg.connectomeFile, ...
    'Subjects',        subjects, ...
    'TrivialLabelSuffixes', cfg.trivialLabelSuffixes, ...
    'ResultsDir',      resultsDir, ...
    'RunParameters',   runParams, ...
    'Normalisation',   scheme, ...
    'ParkesC',         cfg.parkesC, ...
    'ColScale',        cfg.colScale, ...
    'TvbPercentile',   cfg.tvbPercentile, ...
    'X0Base',          cfg.x0Base, ...
    'X0Upper',         cfg.x0Upper, ...
    'X0Step',          cfg.x0Step, ...
    'BisectTol',       cfg.bisectTol, ...
    'MaxK',            cfg.maxK, ...
    'Tau0',            cfg.tau0, ...
    'TopK',            cfg.topK, ...
    'Plot',            cfg.plot, ...
    'Verbose',         cfg.verbose, ...
    'SaveResults',     cfg.saveResults};

groupArgs = { ...
    'ResultsDir',    resultsDir, ...
    'Normalisation', scheme, ...
    'CasePrefix',    cfg.casePrefix, ...
    'ControlPrefix', cfg.controlPrefix, ...
    'RunParameters', runParams, ...
    'SaveResults',   cfg.saveResults, ...
    'Plot',          cfg.plot};

%% ---- 1. QC ------------------------------------------------------------
if cfg.pipelineRunQc
    manifest.steps(end+1) = runStep(@() plotConnectomeHeatmaps( ...
        'DataRoot', cfg.dataRoot, 'ConnectomeFile', cfg.connectomeFile, ...
        'Subjects', subjects, 'Layout', 'hemisphere', ...
        'LeftPrefix', cfg.leftPrefix, 'RightPrefix', cfg.rightPrefix, ...
        'TrivialLabelSuffixes', cfg.trivialLabelSuffixes, ...
        'ResultsDir', resultsDir, 'SaveResults', cfg.saveResults), ...
        'plotConnectomeHeatmaps', cfg.pipelineStopOnError);
end

%% ---- 2-3. Stability centralities and postflight -----------------------
cohortReady = ~cfg.pipelineRunStabilityCentralities;
if cfg.pipelineRunStabilityCentralities
    stepRun = runStep(@() runCohortStabilityCentralities(stabilityArgs{:}), ...
        'runCohortStabilityCentralities', cfg.pipelineStopOnError);
    manifest.steps(end+1) = stepRun;

    stepCheck = runStep(@() assertCohortComplete(resultsDir, scheme, subjects), ...
        'postflight_cohortComplete', cfg.pipelineStopOnError);
    manifest.steps(end+1) = stepCheck;

    cohortReady = stepRun.success && stepCheck.success;
end

%% ---- 4. Centrality correlations ---------------------------------------
if cfg.pipelineRunCentralityCorr
    if cohortReady
        manifest.steps(end+1) = runStep(@() compareCentralityMeasures( ...
            groupArgs{:}, 'CorrType', cfg.corrType), ...
            'compareCentralityMeasures', cfg.pipelineStopOnError);
    else
        manifest.steps(end+1) = skipStep('compareCentralityMeasures', ...
            'per-subject cohort incomplete');
    end
end

%% ---- 5. Case vs control ------------------------------------------------
cohortSummary = {[]};
if cfg.pipelineRunCohortComparison
    if cohortReady && canCompareGroups
        manifest.steps(end+1) = runStep(@() captureCohortSummary(), ...
            'compareCohortGroups', cfg.pipelineStopOnError);
    elseif ~cohortReady
        manifest.steps(end+1) = skipStep('compareCohortGroups', ...
            'per-subject cohort incomplete');
    else
        manifest.steps(end+1) = skipStep('compareCohortGroups', ...
            sprintf('no subjects matching both "%s" and "%s"', ...
                cfg.casePrefix, cfg.controlPrefix), true);
    end
end

%% ---- 6. Hemispheric asymmetry -------------------------------------------
lateralitySummary = {[]};
if cfg.pipelineRunHemisphericAsymmetry
    if cohortReady && canCompareGroups
        manifest.steps(end+1) = runStep(@() captureLateralitySummary(), ...
            'compareHemisphericAsymmetry', cfg.pipelineStopOnError);
    elseif ~cohortReady
        manifest.steps(end+1) = skipStep('compareHemisphericAsymmetry', ...
            'per-subject cohort incomplete');
    else
        manifest.steps(end+1) = skipStep('compareHemisphericAsymmetry', ...
            sprintf('no subjects matching both "%s" and "%s"', ...
                cfg.casePrefix, cfg.controlPrefix), true);
    end
end

%% ---- 7. Summary figures --------------------------------------------------
if cfg.pipelineRunSummaryFigures
    if ~isempty(cohortSummary{1}) || ~isempty(lateralitySummary{1})
        manifest.steps(end+1) = runStep(@() renderCohortComparisonFigures( ...
            cohortSummary{1}, lateralitySummary{1}, ...
            'OutputDir', experimentResultsDir(resultsDir, 'comparison_figures')), ...
            'renderCohortComparisonFigures', cfg.pipelineStopOnError);
    else
        manifest.steps(end+1) = skipStep('renderCohortComparisonFigures', ...
            'no comparison summaries available', true);
    end
end

%% ---- 8. Network-level D_st -----------------------------------------------
if cfg.pipelineRunDstCohort
    if cohortReady
        manifest.steps(end+1) = runStep(@() compareDstAcrossCohort(groupArgs{:}), ...
            'compareDstAcrossCohort', cfg.pipelineStopOnError);
    else
        manifest.steps(end+1) = skipStep('compareDstAcrossCohort', ...
            'per-subject cohort incomplete');
    end
end

%% ---- 9. Reports -----------------------------------------------------------
if cfg.pipelineRunReports
    if cohortReady
        logFile = experimentProvenanceFile(resultsDir, cfg.experimentName, 'reports');
        manifest.steps(end+1) = runStep( ...
            @() writeReports(logFile, resultsDir, scheme, cfg, canCompareGroups), ...
            'reports', cfg.pipelineStopOnError);
    else
        manifest.steps(end+1) = skipStep('reports', 'per-subject cohort incomplete');
    end
end

%% ---- Finish ---------------------------------------------------------------
manifest.finishedAt = timestamp();
manifest.allSucceeded = all([manifest.steps.success]);

runParams.finishedAt = manifest.finishedAt;
runParams.pipelineSteps = manifest.steps;
manifest.runParameters = runParams;
experimentRunParameters('save', resultsDir, runParams);
save(experimentProvenanceFile(resultsDir, cfg.experimentName, 'manifest'), ...
    'manifest', 'cfg', 'runParams');

if manifest.allSucceeded
    fprintf('\nExperiment %s finished successfully.\n', cfg.experimentName);
else
    failed = manifest.steps(~[manifest.steps.success]);
    fprintf('\nExperiment %s finished with %d failed/skipped step(s): %s\n', ...
        cfg.experimentName, numel(failed), strjoin({failed.name}, ', '));
    fprintf('See %s\n', experimentProvenanceFile(resultsDir, cfg.experimentName, 'manifest'));
end

    function captureCohortSummary()
        cohortSummary{1} = compareCohortGroups(groupArgs{:}, ...
            'Alpha', cfg.alpha, 'ZThreshold', cfg.zThreshold);
    end

    function captureLateralitySummary()
        lateralitySummary{1} = compareHemisphericAsymmetry(groupArgs{:}, ...
            'Alpha', cfg.alpha, 'ZThreshold', cfg.zThreshold, ...
            'LateralityBasis', cfg.lateralityBasis, ...
            'LeftPrefix', cfg.leftPrefix, 'RightPrefix', cfg.rightPrefix);
    end
end

%% ------------------------------------------------------------------
function steps = emptySteps()
steps = struct('name', {}, 'success', {}, 'durationSec', {}, 'errorMessage', {});
end

%% ------------------------------------------------------------------
function step = runStep(fn, stepName, stopOnError)
step = struct('name', stepName, 'success', false, 'durationSec', NaN, 'errorMessage', '');
t0 = tic;
try
    fn();
    step.success = true;
catch ME
    step.errorMessage = ME.message;
    if stopOnError
        rethrow(ME);
    end
    warning('runConnectomeExperiment:StepFailed', '%s failed: %s', stepName, ME.message);
end
step.durationSec = toc(t0);
fprintf('  [%s] %s (%.1f s)\n', tag(step.success), stepName, step.durationSec);
end

%% ------------------------------------------------------------------
function step = skipStep(stepName, reason, countAsSuccess)
if nargin < 3
    countAsSuccess = false;
end
step = struct('name', stepName, 'success', logical(countAsSuccess), ...
    'durationSec', 0, 'errorMessage', ['skipped: ' reason]);
if countAsSuccess
    label = 'OK';
else
    label = 'SKIP';
end
fprintf('  [%s] %s (%s)\n', label, stepName, reason);
end

%% ------------------------------------------------------------------
function s = tag(success)
if success
    s = 'OK';
else
    s = 'FAIL';
end
end

%% ------------------------------------------------------------------
function assertCohortComplete(resultsDir, scheme, subjects)
files = listSubjectResultFiles(resultsDir, scheme);
if numel(files) < numel(subjects)
    error('runConnectomeExperiment:IncompleteCohort', ...
        ['Expected %d per-subject result files for scheme "%s" but found %d in %s. ' ...
         'Comparison steps will be skipped.'], ...
        numel(subjects), scheme, numel(files), resultsDir);
end
end

%% ------------------------------------------------------------------
function appendRunNote(noteFile, folderName, runStamp)
%APPENDRUNNOTE  Record where this run's outputs went, inside the snapshot.
fid = fopen(noteFile, 'a');
if fid < 0
    return;
end
closer = onCleanup(@() fclose(fid));
fprintf(fid, '\n# Appended by runConnectomeExperiment at run start\n');
fprintf(fid, 'experiment.resultsFolder=%s\n', folderName);
fprintf(fid, 'experiment.runStamp=%s\n', runStamp);
end

%% ------------------------------------------------------------------
function writeReports(logFile, resultsDir, scheme, cfg, canCompareGroups)
%WRITEREPORTS  Tee the console reports into the run's provenance log.
diary(logFile);
closer = onCleanup(@() diary('off'));

fprintf('=== reportTopNodes ===\n');
reportTopNodes('ResultsDir', resultsDir, 'Normalisation', scheme, 'N', cfg.topK);

if canCompareGroups
    fprintf('\n=== reportCohortOutliers ===\n');
    try
        reportCohortOutliers('ResultsDir', resultsDir, 'Normalisation', scheme, ...
            'ZThreshold', cfg.zThreshold);
    catch ME
        fprintf('reportCohortOutliers skipped: %s\n', ME.message);
    end
else
    fprintf('\n=== reportCohortOutliers === skipped (no case/control split)\n');
end
end

%% ------------------------------------------------------------------
function s = timestamp()
s = datestr(now, 'yyyy-mm-dd HH:MM:SS'); %#ok<DATST,TNOW1>
end
