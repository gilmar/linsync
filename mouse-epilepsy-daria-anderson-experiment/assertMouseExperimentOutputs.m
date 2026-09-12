function report = assertMouseExperimentOutputs(resultsDir, scheme, cfg, mice)
%ASSERTMOUSEEXPERIMENTOUTPUTS  Verify expected artefacts for one experiment run.
%
%   report = assertMouseExperimentOutputs(resultsDir, scheme, cfg, mice)
%
%   Checks filenames produced by runMouseExperiment (stabilityCentralities_*,
%   centrality_corr_*, compare_anderson_vs_arnold_*, D_st_cohort_*, etc.)
%   under the category subfolders defined by mouseExperimentResultsLayout.
%   Legacy flat layouts under the experiment root are also accepted when
%   reading (see locateResultsFile).
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

expName = cfg.experimentName;
req = expectedOutputPaths(resultsDir, scheme, cfg, mice, prefix, expName);
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
function paths = expectedOutputPaths(resultsDir, scheme, cfg, mice, prefix, expName)
paths = {};

paths = [paths; provenancePaths(resultsDir, expName)];

if cfg.pipelineRunHeatmap
    paths = [paths; qcPaths(resultsDir, cfg)];
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
    paths{end+1, 1} = mouseExperimentProvenanceFile(resultsDir, expName, 'reports');
end

paths = paths(:);
end

function paths = provenancePaths(resultsDir, expName)
paths = {
    mouseExperimentProvenanceFile(resultsDir, expName, 'properties')
    mouseExperimentProvenanceFile(resultsDir, expName, 'parametersMat')
    mouseExperimentProvenanceFile(resultsDir, expName, 'parametersJson')
    fullfile(mouseResultsDir(resultsDir, 'provenance'), 'run_manifest.mat')
    };
end

function paths = qcPaths(resultsDir, cfg)
paths = {};
if ~cfg.saveResults
    return;
end
d = mouseResultsDir(resultsDir, 'qc');
paths = {
    fullfile(d, 'mouse_heatmaps_overview_reference.fig')
    fullfile(d, 'mouse_heatmaps_overview_reference.png')
    };
end

function paths = perMouseStabilityPaths(resultsDir, prefix, mouseId, scheme, cfg)
d = mouseResultsDir(resultsDir, 'stability');
paths = {fullfile(d, sprintf('%s_%s_%s_results.mat', prefix, mouseId, scheme))};
if cfg.saveResults
    paths = [paths; {
        fullfile(d, sprintf('%s_%s_%s_figure.fig', prefix, mouseId, scheme))
        fullfile(d, sprintf('%s_%s_%s_figure.png', prefix, mouseId, scheme))
        }];
end
end

function paths = summaryStabilityPaths(resultsDir, prefix, scheme, cfg)
paths = {};
if ~cfg.saveResults
    return;
end
d = mouseResultsDir(resultsDir, 'stability');
paths = {
    fullfile(d, sprintf('%s_summary_topnodes_%s.csv', prefix, scheme))
    fullfile(d, sprintf('%s_summary_overview_%s.fig', prefix, scheme))
    fullfile(d, sprintf('%s_summary_overview_%s.png', prefix, scheme))
    fullfile(d, sprintf('%s_summary_%s.mat', prefix, scheme))
    };
end

function paths = centralityCorrPaths(resultsDir, scheme, mice, cfg)
paths = {};
if ~cfg.saveResults
    return;
end
d = mouseResultsDir(resultsDir, 'centrality_correlation');
for m = 1:numel(mice)
    mid = mice{m};
    paths = [paths; {
        fullfile(d, sprintf('centrality_corr_%s_%s.fig', mid, scheme))
        fullfile(d, sprintf('centrality_corr_%s_%s.png', mid, scheme))
        }]; %#ok<AGROW>
end
paths = [paths; {
    fullfile(d, sprintf('centrality_corr_mean_%s.fig', scheme))
    fullfile(d, sprintf('centrality_corr_mean_%s.png', scheme))
    fullfile(d, sprintf('centrality_corr_bars_%s.fig', scheme))
    fullfile(d, sprintf('centrality_corr_bars_%s.png', scheme))
    fullfile(d, sprintf('centrality_corr_%s.mat', scheme))
    }];
end

function paths = andersonPaths(resultsDir, scheme, mice, cfg)
paths = {};
if ~cfg.saveResults
    return;
end
d = mouseResultsDir(resultsDir, 'anderson_vs_arnold');
hasBC = true;
matPath = locateResultsFile(resultsDir, 'anderson_vs_arnold', ...
    sprintf('compare_anderson_vs_arnold_%s.mat', scheme));
if ~isempty(matPath)
    s = load(matPath, 'hasBC');
    if isfield(s, 'hasBC')
        hasBC = s.hasBC;
    end
end
metricSuffixes = {'D_to_i', 'D_k_to'};
if hasBC
    metricSuffixes{end+1} = 'BC';
end
andersonMice = mice(startsWith(mice, 'Anderson'));
for m = 1:numel(andersonMice)
    aname = andersonMice{m};
    for s = 1:numel(metricSuffixes)
        suffix = metricSuffixes{s};
        paths = [paths; {
            fullfile(d, sprintf('compare_anderson_vs_arnold_%s_%s_%s.fig', aname, scheme, suffix))
            fullfile(d, sprintf('compare_anderson_vs_arnold_%s_%s_%s.png', aname, scheme, suffix))
            }]; %#ok<AGROW>
    end
end
paths{end+1, 1} = fullfile(d, sprintf('compare_anderson_vs_arnold_%s.mat', scheme));
end

function paths = lrAsymmetryPaths(resultsDir, scheme, mice, cfg)
paths = {};
if ~cfg.saveResults
    return;
end
d = mouseResultsDir(resultsDir, 'lr_asymmetry');
andersonMice = mice(startsWith(mice, 'Anderson'));
for m = 1:numel(andersonMice)
    aname = andersonMice{m};
    paths = [paths; {
        fullfile(d, sprintf('compare_LR_asymmetry_%s_%s.fig', aname, scheme))
        fullfile(d, sprintf('compare_LR_asymmetry_%s_%s.png', aname, scheme))
        }]; %#ok<AGROW>
end
paths = [paths; {
    fullfile(d, sprintf('LR_asymmetry_groupTest_%s.fig', scheme))
    fullfile(d, sprintf('LR_asymmetry_groupTest_%s.png', scheme))
    fullfile(d, sprintf('LR_asymmetry_groupTest_%s.csv', scheme))
    fullfile(d, sprintf('LR_asymmetry_systematic_%s.fig', scheme))
    fullfile(d, sprintf('LR_asymmetry_systematic_%s.png', scheme))
    fullfile(d, sprintf('LR_asymmetry_systematic_%s.csv', scheme))
    fullfile(d, sprintf('LR_asymmetry_%s.mat', scheme))
    }];
end

function paths = dstPaths(resultsDir, scheme, cfg)
paths = {};
if ~cfg.saveResults
    return;
end
d = mouseResultsDir(resultsDir, 'dst_cohort');
paths = {
    fullfile(d, sprintf('D_st_cohort_%s.csv', scheme))
    fullfile(d, sprintf('D_st_cohort_%s.fig', scheme))
    fullfile(d, sprintf('D_st_cohort_%s.png', scheme))
    fullfile(d, sprintf('D_st_cohort_%s.mat', scheme))
    };
end
