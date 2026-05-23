function runAllMiceSection45(varargin)
%RUNALLMICESECTION45  Apply the Liao §4.5 workflow to every Daria/Anderson
% mouse coarse connectome under data/, then summarise which nodes "stand
% out" (high stability susceptibility / influence centrality, equivalently
% low critical excitability x^c_{0,i}) and on which mice.
%
% Outputs (under results/):
%   section45_<mouseId>_<scheme>_results.mat   per-mouse results struct
%   section45_<mouseId>_<scheme>_figure.{fig,png}
%   section45_summary_topnodes_<scheme>.csv    per-mouse top-K rankings
%   section45_summary_overview_<scheme>.{fig,png}
%
% Name-value options:
%   Normalisation : 'tvb' (default) | 'parkes' | 'column' | 'none'
%   ParkesC       : 1.0     c parameter for the 'parkes' scheme
%   ColScale      : 0.95    target column sum for the 'column' scheme
%   TopK          : 10      number of top nodes to print/save per mouse
%   x0Base, x0Upper, x0Step, BisectTol, MaxK, tau0  -- forwarded to runMouseSection45
%   Plot, Verbose, SaveResults, ResultsDir

p = inputParser;
addParameter(p, 'Normalisation', 'tvb', @(s) ischar(s) || isstring(s));
addParameter(p, 'ParkesC',       1.0,      @isscalar);
addParameter(p, 'ColScale',      0.95,     @isscalar);
addParameter(p, 'TopK',          10,       @isscalar);
addParameter(p, 'x0Base',        -2.3,     @isscalar);
addParameter(p, 'x0Upper',       -1.0,     @isscalar);
addParameter(p, 'x0Step',         0.01,    @isscalar);
addParameter(p, 'BisectTol',      1e-4,    @isscalar);
addParameter(p, 'MaxK',           1e8,     @isscalar);
addParameter(p, 'tau0',           6667,    @isscalar);
addParameter(p, 'Plot',           true,    @islogical);
addParameter(p, 'Verbose',        true,    @islogical);
addParameter(p, 'SaveResults',    true,    @islogical);
addParameter(p, 'ResultsDir',     '',      @(s) ischar(s) || isstring(s));
addParameter(p, 'RunParameters',  [],      @(x) isempty(x) || isstruct(x));
parse(p, varargin{:});
opts          = p.Results;
normalisation = char(opts.Normalisation);
topK          = opts.TopK;

resultsDir = resolveMouseResultsDir(opts.ResultsDir);

mice = listAvailableMice();
if isempty(mice)
    error('runAllMiceSection45: no mouse CSVs found under data/.');
end

allResults = cell(numel(mice), 1);
totalTic = tic;
for m = 1:numel(mice)
    mouseId = mice{m};
    fprintf('\n========== Running %s (%d/%d) ==========\n', mouseId, m, numel(mice));
    try
        allResults{m} = runMouseSection45(mouseId, ...
            'Normalisation', normalisation, ...
            'ParkesC',       opts.ParkesC, ...
            'ColScale',      opts.ColScale, ...
            'x0Base',        opts.x0Base, ...
            'x0Upper',       opts.x0Upper, ...
            'x0Step',        opts.x0Step, ...
            'BisectTol',     opts.BisectTol, ...
            'MaxK',          opts.MaxK, ...
            'tau0',          opts.tau0, ...
            'ResultsDir',    resultsDir, ...
            'RunParameters', opts.RunParameters, ...
            'SaveResults',   opts.SaveResults, ...
            'Plot',          opts.Plot, ...
            'Verbose',       opts.Verbose);
    catch ME
        warning('runAllMiceSection45:MouseFailed', ...
            'Mouse %s failed: %s', mouseId, ME.message);
        allResults{m} = [];
    end
end
fprintf('\nAll mice processed in %.1f s.\n', toc(totalTic));

%% Cross-mouse summary -----------------------------------------------
% Use the first non-empty result to build the canonical label list (all
% mice share the same parcellation).
canonicalLabels = [];
for m = 1:numel(mice)
    if ~isempty(allResults{m})
        canonicalLabels = allResults{m}.labels;
        break;
    end
end
if isempty(canonicalLabels)
    error('runAllMiceSection45: no successful mouse runs.');
end
N = numel(canonicalLabels);

% Pack per-mouse vectors into NxM matrices (NaN-padded).
M = numel(mice);
D_susc_all = NaN(N, M);
D_infl_all = NaN(N, M);
x0_crit_all = NaN(N, M);
for m = 1:M
    if isempty(allResults{m}); continue; end
    r = allResults{m};
    D_susc_all(:, m)  = r.D_susceptibility;
    D_infl_all(:, m)  = r.D_influence;
    x0_crit_all(:, m) = r.x0_crit;
end

% Aggregate ranking (mean rank across mice, where a high-D / low-x0 node
% is rank 1). Trivial / NaN entries are excluded from the ranking on
% that mouse but counted with the worst-rank fallback for the average.
[rank_susc, rank_infl, rank_x0] = deal(NaN(N, M));
for m = 1:M
    if isempty(allResults{m}); continue; end
    rank_susc(:, m) = rankWithNaN(D_susc_all(:, m), 'descend');
    rank_infl(:, m) = rankWithNaN(D_infl_all(:, m), 'descend');
    rank_x0(:, m)   = rankWithNaN(x0_crit_all(:, m), 'ascend');
end
mean_rank_susc = mean(rank_susc, 2, 'omitnan');
mean_rank_infl = mean(rank_infl, 2, 'omitnan');
mean_rank_x0   = mean(rank_x0,   2, 'omitnan');

%% Print per-mouse top-K and aggregate top-K
fprintf('\n----- Per-mouse top-%d nodes (lowest x^c_{0,i}, i.e. most susceptible) -----\n', topK);
for m = 1:M
    if isempty(allResults{m})
        fprintf('%s : (no result)\n', mice{m});
        continue;
    end
    [~, ord] = sort(x0_crit_all(:, m), 'ascend');  % NaN goes to end
    fprintf('\n%s top %d:\n', mice{m}, topK);
    for kk = 1:min(topK, N)
        i = ord(kk);
        fprintf('  %2d. %-32s  x0^c = %+.4f  D(->i) = %8.2f  D(k->) = %8.2f\n', ...
            kk, canonicalLabels{i}, x0_crit_all(i, m), ...
            D_susc_all(i, m), D_infl_all(i, m));
    end
end

fprintf('\n----- Aggregate top-%d (smallest mean rank of x^c_{0,i} across mice) -----\n', topK);
[~, agg_ord] = sort(mean_rank_x0, 'ascend');
for kk = 1:min(topK, N)
    i = agg_ord(kk);
    fprintf('  %2d. %-32s  mean rank x0^c = %5.1f  | susc = %5.1f  | infl = %5.1f\n', ...
        kk, canonicalLabels{i}, mean_rank_x0(i), ...
        mean_rank_susc(i), mean_rank_infl(i));
end

%% Save summary CSV
summaryRows = cell(0, 1);
header = {'rank', 'region', 'mean_rank_x0', 'mean_rank_susc', 'mean_rank_infl'};
for m = 1:M
    header{end+1} = sprintf('x0_crit__%s', mice{m}); %#ok<SAGROW>
end
for m = 1:M
    header{end+1} = sprintf('D_susc__%s', mice{m}); %#ok<SAGROW>
end
for m = 1:M
    header{end+1} = sprintf('D_infl__%s', mice{m}); %#ok<SAGROW>
end

summaryTable = cell(N, length(header));
for kk = 1:N
    i = agg_ord(kk);
    summaryTable{kk, 1} = kk;
    summaryTable{kk, 2} = canonicalLabels{i};
    summaryTable{kk, 3} = mean_rank_x0(i);
    summaryTable{kk, 4} = mean_rank_susc(i);
    summaryTable{kk, 5} = mean_rank_infl(i);
    col = 6;
    for m = 1:M
        summaryTable{kk, col} = x0_crit_all(i, m); col = col + 1;
    end
    for m = 1:M
        summaryTable{kk, col} = D_susc_all(i, m); col = col + 1;
    end
    for m = 1:M
        summaryTable{kk, col} = D_infl_all(i, m); col = col + 1;
    end
end

T = cell2table(summaryTable, 'VariableNames', header);
csvFile = fullfile(resultsDir, sprintf('section45_summary_topnodes_%s.csv', normalisation));
writetable(T, csvFile);
fprintf('\nWrote summary table %s\n', csvFile);

%% Cross-mouse overview figure
isParkes  = strcmpi(normalisation, 'parkes');
isColumn  = ismember(lower(normalisation), {'column', 'col', 'colnorm'});
isLinear  = isParkes || isColumn;
% Parkes on symmetric K: D(k->)=D(->i); only show D(->i) (1 panel).
% Column: D(k->) != D(->i), no x0; show both D panels (2 panels).
% TVB / non-linear: x0^c + both centrality panels (3 panels).
overviewFig = figure('Name', sprintf('§4.5 mouse overview (%s)', normalisation), ...
                     'Position', [80 80 1600 700]);

if isParkes
    subplot(1, 1, 1);
    imagesc(D_susc_all);
    yticks(1:N); yticklabels(canonicalLabels);
    set(gca, 'FontSize', 6, 'YDir', 'normal');
    xticks(1:M); xticklabels(mice); xtickangle(45);
    colorbar; colormap(gca, parula);
    title('D(\rightarrow i)  (susceptibility)');
elseif isColumn
    subplot(1, 2, 1);
    imagesc(D_susc_all);
    yticks(1:N); yticklabels(canonicalLabels);
    set(gca, 'FontSize', 6, 'YDir', 'normal');
    xticks(1:M); xticklabels(mice); xtickangle(45);
    colorbar; colormap(gca, parula);
    title('D(\rightarrow i)  (susceptibility)');

    subplot(1, 2, 2);
    imagesc(D_infl_all);
    yticks(1:N); yticklabels(canonicalLabels);
    set(gca, 'FontSize', 6, 'YDir', 'normal');
    xticks(1:M); xticklabels(mice); xtickangle(45);
    colorbar; colormap(gca, parula);
    title('D(k \rightarrow)  (influence)');
else
    subplot(1, 3, 1);
    imagesc(x0_crit_all);
    yticks(1:N); yticklabels(canonicalLabels);
    set(gca, 'FontSize', 6, 'YDir', 'normal');
    xticks(1:M); xticklabels(mice); xtickangle(45);
    colorbar; colormap(gca, parula);
    title('x^c_{0,i}  (lower = more epileptogenic)');

    subplot(1, 3, 2);
    imagesc(D_susc_all);
    yticks(1:N); yticklabels(canonicalLabels);
    set(gca, 'FontSize', 6, 'YDir', 'normal');
    xticks(1:M); xticklabels(mice); xtickangle(45);
    colorbar; colormap(gca, parula);
    title('D(\rightarrow i)  (susceptibility)');

    subplot(1, 3, 3);
    imagesc(D_infl_all);
    yticks(1:N); yticklabels(canonicalLabels);
    set(gca, 'FontSize', 6, 'YDir', 'normal');
    xticks(1:M); xticklabels(mice); xtickangle(45);
    colorbar; colormap(gca, parula);
    title('D(k \rightarrow)  (influence)');
end

if isParkes
    % Avoid sprintf here: \r in \rightarrow would be mis-parsed as carriage return.
    sgtitle(['Cross-mouse §4.5 summary  --  normalisation = ' normalisation ...
             '  [D(k\rightarrow) \equiv D(\rightarrow i) for symmetric K]']);
else
    sgtitle(sprintf('Cross-mouse §4.5 summary  --  normalisation = %s', normalisation));
end

savefig(overviewFig, fullfile(resultsDir, sprintf('section45_summary_overview_%s.fig', normalisation)));
try
    exportgraphics(overviewFig, fullfile(resultsDir, sprintf('section45_summary_overview_%s.png', normalisation)), 'Resolution', 200);
catch
    saveas(overviewFig, fullfile(resultsDir, sprintf('section45_summary_overview_%s.png', normalisation)));
end

if isempty(opts.RunParameters)
    runParameters = mouseExperimentRunParameters('buildFromSection45', opts, 'cohort');
    runParameters.cohort.mice = mice;
    runParameters.cohort.nMice = numel(mice);
    runParameters.section45.topK = topK;
else
    runParameters = opts.RunParameters;
end
save(fullfile(resultsDir, sprintf('section45_summary_%s.mat', normalisation)), ...
    'mice', 'canonicalLabels', 'D_susc_all', 'D_infl_all', 'x0_crit_all', ...
    'mean_rank_susc', 'mean_rank_infl', 'mean_rank_x0', 'normalisation', ...
    'runParameters');

fprintf('\nDone.\n');
end

%% ------------------------------------------------------------------
function r = rankWithNaN(values, dir)
%RANKWITHNAN  Return rank vector for values; NaN entries get rank N+1.
v = values(:);
N = numel(v);
r = NaN(N, 1);
finite = ~isnan(v);
[~, ord] = sort(v(finite), dir);
ranks = NaN(sum(finite), 1);
ranks(ord) = 1:sum(finite);
r(finite) = ranks;
r(~finite) = N + 1;  % NaN -> worst rank (won't bias mean towards top)
end
