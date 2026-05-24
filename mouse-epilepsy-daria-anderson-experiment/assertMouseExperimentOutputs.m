function report = assertMouseExperimentOutputs(resultsDir, scheme, cfg, mice)
%ASSERTMOUSEEXPERIMENTOUTPUTS  Verify expected artefacts for one experiment run.
%
%   report = assertMouseExperimentOutputs(resultsDir, scheme, cfg, mice)
%
%   Checks filenames produced by runMouseExperiment (stabilityCentralities_*,
%   centrality_corr_*, compare_anderson_vs_arnold_*, D_st_cohort_*, etc.)
%   against the artefact types in results/mouse_experiment_reference/ (which
%   used the legacy section45_* prefix for the same pipeline steps).
%
%   Throws assertMouseExperimentOutputs:MissingFiles if anything required
%   is absent.

if nargin < 4 || isempty(mice)
    mice = listAvailableMice();
end
scheme = char(scheme);
resultsDir = char(resultsDir);
prefix = mouseExperimentResultPrefix();

report = struct();
report.scheme = scheme;
report.resultsDir = resultsDir;
report.prefix = prefix;
report.mice = mice;
report.required = {};
report.missing = {};
report.present = {};

req = expectedOutputPaths(resultsDir, scheme, cfg, mice, prefix);
report.required = req;

for k = 1:numel(req)
    fpath = req{k};
    if isfile(fpath)
        report.present{end+1, 1} = fpath; %#ok<AGROW>
    else
        report.missing{end+1, 1} = fpath; %#ok<AGROW>
    end
end

fprintf('\nOutput check (%s): %d/%d required files present.\n', ...
    scheme, numel(report.present), numel(req));

if ~isempty(report.missing)
    fprintf('Missing:\n');
    for k = 1:numel(report.missing)
        fprintf('  %s\n', report.missing{k});
    end
    error('assertMouseExperimentOutputs:MissingFiles', ...
        '%d required file(s) missing under %s (scheme=%s). See list above.', ...
        numel(report.missing), resultsDir, scheme);
end
end

%% ------------------------------------------------------------------
function paths = expectedOutputPaths(resultsDir, scheme, cfg, mice, prefix)
paths = {};

% Orchestrator provenance (always)
paths = [paths; provenancePaths(resultsDir)];

if cfg.pipelineRunHeatmap
    paths = [paths; { ...
        fullfile(resultsDir, 'mouse_heatmaps_overview_reference.fig')
        fullfile(resultsDir, 'mouse_heatmaps_overview_reference.png')
        }];
end

if cfg.pipelineRunStabilityCentralities
    for m = 1:numel(mice)
        mid = mice{m};
        paths = [paths; perMouseStabilityPaths(resultsDir, prefix, mid, scheme, cfg)]; %#ok<AGROW>
    end
    paths = [paths; summaryStabilityPaths(resultsDir, prefix, scheme, cfg)];
end

if cfg.pipelineRunCentralityCorr
    paths = [paths; centralityCorrPaths(resultsDir, scheme, mice, cfg)];
end

if cfg.pipelineRunAndersonVsArnold
    andersonMice = mice(startsWith(mice, 'Anderson'));
    if ~isempty(andersonMice) && any(startsWith(mice, 'Arnold'))
        paths = [paths; andersonPaths(resultsDir, scheme, mice, cfg)];
    end
end

if cfg.pipelineRunLRAsymmetry
    andersonMice = mice(startsWith(mice, 'Anderson'));
    if ~isempty(andersonMice) && any(startsWith(mice, 'Arnold'))
        paths = [paths; lrAsymmetryPaths(resultsDir, scheme, mice, cfg)];
    end
end

if cfg.pipelineRunDstCohort
    paths = [paths; dstPaths(resultsDir, scheme, cfg)];
end

if cfg.pipelineRunReports
    paths{end+1, 1} = fullfile(resultsDir, 'reports.log');
end

paths = paths(:);
end

function paths = provenancePaths(resultsDir)
paths = {
    fullfile(resultsDir, 'experiment.properties')
    fullfile(resultsDir, 'experiment_parameters.mat')
    fullfile(resultsDir, 'experiment_parameters.json')
    fullfile(resultsDir, 'run_manifest.mat')
    };
end

function paths = perMouseStabilityPaths(resultsDir, prefix, mouseId, scheme, cfg)
paths = {fullfile(resultsDir, sprintf('%s_%s_%s_results.mat', prefix, mouseId, scheme))};
if cfg.saveResults
    paths = [paths; {
        fullfile(resultsDir, sprintf('%s_%s_%s_figure.fig', prefix, mouseId, scheme))
        fullfile(resultsDir, sprintf('%s_%s_%s_figure.png', prefix, mouseId, scheme))
        }];
end
end

function paths = summaryStabilityPaths(resultsDir, prefix, scheme, cfg)
paths = {};
if ~cfg.saveResults
    return;
end
paths = {
    fullfile(resultsDir, sprintf('%s_summary_topnodes_%s.csv', prefix, scheme))
    fullfile(resultsDir, sprintf('%s_summary_overview_%s.fig', prefix, scheme))
    fullfile(resultsDir, sprintf('%s_summary_overview_%s.png', prefix, scheme))
    fullfile(resultsDir, sprintf('%s_summary_%s.mat', prefix, scheme))
    };
end

function paths = centralityCorrPaths(resultsDir, scheme, mice, cfg)
paths = {};
if ~cfg.saveResults
    return;
end
for m = 1:numel(mice)
    mid = mice{m};
    paths = [paths; {
        fullfile(resultsDir, sprintf('centrality_corr_%s_%s.fig', mid, scheme))
        fullfile(resultsDir, sprintf('centrality_corr_%s_%s.png', mid, scheme))
        }]; %#ok<AGROW>
end
paths = [paths; {
    fullfile(resultsDir, sprintf('centrality_corr_mean_%s.fig', scheme))
    fullfile(resultsDir, sprintf('centrality_corr_mean_%s.png', scheme))
    fullfile(resultsDir, sprintf('centrality_corr_bars_%s.fig', scheme))
    fullfile(resultsDir, sprintf('centrality_corr_bars_%s.png', scheme))
    fullfile(resultsDir, sprintf('centrality_corr_%s.mat', scheme))
    }];
end

function paths = andersonPaths(resultsDir, scheme, mice, cfg)
paths = {};
if ~cfg.saveResults
    return;
end
andersonMice = mice(startsWith(mice, 'Anderson'));
for m = 1:numel(andersonMice)
    aname = andersonMice{m};
    paths = [paths; {
        fullfile(resultsDir, sprintf('compare_anderson_vs_arnold_%s_%s.fig', aname, scheme))
        fullfile(resultsDir, sprintf('compare_anderson_vs_arnold_%s_%s.png', aname, scheme))
        }]; %#ok<AGROW>
    % Outlier CSV is written only when outliers exist; do not require.
end
paths{end+1, 1} = fullfile(resultsDir, sprintf('compare_anderson_vs_arnold_%s.mat', scheme));
end

function paths = lrAsymmetryPaths(resultsDir, scheme, mice, cfg)
paths = {};
if ~cfg.saveResults
    return;
end
andersonMice = mice(startsWith(mice, 'Anderson'));
for m = 1:numel(andersonMice)
    aname = andersonMice{m};
    paths = [paths; {
        fullfile(resultsDir, sprintf('compare_LR_asymmetry_%s_%s.fig', aname, scheme))
        fullfile(resultsDir, sprintf('compare_LR_asymmetry_%s_%s.png', aname, scheme))
        }]; %#ok<AGROW>
end
paths = [paths; {
    fullfile(resultsDir, sprintf('LR_asymmetry_groupTest_%s.fig', scheme))
    fullfile(resultsDir, sprintf('LR_asymmetry_groupTest_%s.png', scheme))
    fullfile(resultsDir, sprintf('LR_asymmetry_groupTest_%s.csv', scheme))
    fullfile(resultsDir, sprintf('LR_asymmetry_systematic_%s.fig', scheme))
    fullfile(resultsDir, sprintf('LR_asymmetry_systematic_%s.png', scheme))
    fullfile(resultsDir, sprintf('LR_asymmetry_systematic_%s.csv', scheme))
    fullfile(resultsDir, sprintf('LR_asymmetry_%s.mat', scheme))
    }];
end

function paths = dstPaths(resultsDir, scheme, cfg)
paths = {};
if ~cfg.saveResults
    return;
end
paths = {
    fullfile(resultsDir, sprintf('D_st_cohort_%s.csv', scheme))
    fullfile(resultsDir, sprintf('D_st_cohort_%s.fig', scheme))
    fullfile(resultsDir, sprintf('D_st_cohort_%s.png', scheme))
    fullfile(resultsDir, sprintf('D_st_cohort_%s.mat', scheme))
    };
end
