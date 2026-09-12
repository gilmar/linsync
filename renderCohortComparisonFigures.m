function paths = renderCohortComparisonFigures(cohortSummary, lateralitySummary, varargin)
%RENDERCOHORTCOMPARISONFIGURES  Publication-style summary figures for one run.
%
%   paths = renderCohortComparisonFigures(cohortSummary, lateralitySummary)
%   paths = renderCohortComparisonFigures(cohortSummary, lateralitySummary, ...
%               'OutputDir', experimentResultsDir(resultsDir, 'comparison_figures'))
%
%   Turns the two cohort analyses into the two figures a reader actually
%   wants at the end of a run:
%
%     region_susceptibility_outliers -- which regions are significantly
%       more susceptible, D(->i), in the case subjects than in the controls
%     laterality_influence_outliers  -- which homotopic pairs show a
%       significantly abnormal left-right difference in influence, D(k->)
%
%   Both are built from the in-memory summary structs returned by
%   COMPARECOHORTGROUPS and COMPAREHEMISPHERICASYMMETRY -- not by
%   re-reading the CSVs those functions wrote -- so the figures cannot
%   drift away from the statistics that produced them. Titles, legends and
%   cohort sizes are derived from the summaries, so nothing about the group
%   names is hard-coded.
%
%   Name-value options
%     'OutputDir'    : '' -- default: comparison_figures/ under the run
%     'SusceptibilityMetric' : 'D_susceptibility'
%     'InfluenceMetric'      : 'D_influence'
%     'CloseFigures' : true
%
%   Returns a struct of the paths written; a figure with no significant
%   findings is skipped and its field is left empty.
%
%   See also COMPARECOHORTGROUPS, COMPAREHEMISPHERICASYMMETRY,
%   PLOTCOHORTOUTLIERSUMMARY.

p = inputParser;
addParameter(p, 'OutputDir', '', @(s) ischar(s) || isstring(s));
addParameter(p, 'SusceptibilityMetric', 'D_susceptibility', @(s) ischar(s) || isstring(s));
addParameter(p, 'InfluenceMetric',      'D_influence',      @(s) ischar(s) || isstring(s));
addParameter(p, 'CloseFigures', true, @islogical);
parse(p, varargin{:});
opts = p.Results;

if nargin < 2
    lateralitySummary = struct();
end

outDir = char(opts.OutputDir);
if isempty(outDir)
    outDir = defaultOutputDir(cohortSummary);
end
if ~isfolder(outDir)
    mkdir(outDir);
end

paths = struct('regionSusceptibility', '', 'lateralityInfluence', '');

%% ---- Susceptibility outliers ----------------------------------------
if isstruct(cohortSummary) && isfield(cohortSummary, 'deviations')
    info = cohortGroupInfo(cohortSummary);
    rows = collectCohortOutlierRows(cohortSummary, opts.SusceptibilityMetric);
    if isempty(rows)
        fprintf('renderCohortComparisonFigures: no significant D(->i) findings; figure skipped.\n');
    else
        fig = plotCohortOutlierSummary(rows, ...
            'Title', sprintf('Regions more susceptible in %s than in %s', ...
                info.caseLabel, info.controlLabel), ...
            'Subtitle', cohortCountLine(info), ...
            'Footnote', sprintf(['Bonferroni-corrected at alpha = %.3g across ' ...
                'non-trivial nodes x metrics, per case subject.'], info.alpha), ...
            'XLabel', 'excess over control mean D(\rightarrow i)   [%]', ...
            'ValueField', 'effectPct', 'ValueFormat', 'percent');
        if ~isempty(fig)
            base = fullfile(outDir, 'region_susceptibility_outliers');
            saveCohortSummaryFigure(fig, base);
            paths.regionSusceptibility = [base '.png'];
            if opts.CloseFigures
                close(fig);
            end
        end
    end
end

%% ---- Laterality outliers --------------------------------------------
if isstruct(lateralitySummary) && isfield(lateralitySummary, 'perSubject')
    info = cohortGroupInfo(lateralitySummary);
    rows = collectLateralityOutlierRows(lateralitySummary, opts.InfluenceMetric);
    if isempty(rows)
        fprintf('renderCohortComparisonFigures: no significant laterality findings; figure skipped.\n');
    else
        fig = plotCohortOutlierSummary(rows, ...
            'Title', sprintf('Left-right influence asymmetry: %s vs %s', ...
                info.caseLabel, info.controlLabel), ...
            'Subtitle', cohortCountLine(info), ...
            'Footnote', sprintf(['LI_signed = L - R in D(k->) units; Bonferroni-corrected ' ...
                'at alpha = %.3g across region pairs x metrics, per case subject.'], info.alpha), ...
            'XLabel', 'L - R difference in D(k \rightarrow)', ...
            'ValueField', 'value', 'ValueFormat', 'numeric');
        if ~isempty(fig)
            base = fullfile(outDir, 'laterality_influence_outliers');
            saveCohortSummaryFigure(fig, base);
            paths.lateralityInfluence = [base '.png'];
            if opts.CloseFigures
                close(fig);
            end
        end
    end
end

fprintf('Comparison figures -> %s\n', outDir);
end

%% ------------------------------------------------------------------
function outDir = defaultOutputDir(summary)
if isstruct(summary) && isfield(summary, 'resultsDir') && ~isempty(summary.resultsDir)
    outDir = experimentResultsDir(summary.resultsDir, 'comparison_figures');
else
    outDir = fullfile(pwd, 'comparison_figures');
end
end

%% ------------------------------------------------------------------
function line = cohortCountLine(info)
line = sprintf('%d %s subject(s) vs %d %s subject(s)', ...
    info.nCase, info.caseLabel, info.nControl, info.controlLabel);
end

%% ------------------------------------------------------------------
function saveCohortSummaryFigure(fig, basePath)
%SAVECOHORTSUMMARYFIGURE  Fixed 7.5 x 4.94 in at 200 dpi -> 1500 x 988 px.
% Pinning the print size keeps every summary figure in a study identical in
% dimensions, so they can be dropped into a document side by side.
wIn = 7.5;
hIn = 4.94;
dpi = 200;
fig.Units = 'inches';
fig.Position = [1 1 wIn hIn];
fig.PaperUnits = 'inches';
fig.PaperSize = [wIn hIn];
fig.PaperPosition = [0 0 wIn hIn];
fig.PaperPositionMode = 'manual';

drawnow;
savefig(fig, [basePath '.fig']);
drawnow;

% Force the vector renderer so text and thin rules stay crisp in the raster
% output. '-vector' is the current spelling; '-painters' is the same switch
% on releases before R2023a.
try
    print(fig, [basePath '.png'], '-dpng', sprintf('-r%d', dpi), '-vector');
catch
    print(fig, [basePath '.png'], '-dpng', sprintf('-r%d', dpi), '-painters');
end
end
