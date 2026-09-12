function summary = compareLeftRightAsymmetry(varargin)
%COMPARELEFTRIGHTASYMMETRY  Left-right laterality: Anderson vs Arnold.
%
% Pairs each L-<region> node with its R-<region> counterpart, computes
% laterality indices LI_norm = (L-R)/(L+R) and LI_signed = L-R for
% D(->i), D(k->), and BC, then contrasts epileptic (Anderson) vs control
% (Arnold) cohorts in three ways:
%
%   1. Per-Anderson mouse: z-score of LI_<basis> vs Arnold mean+-SD per
%      pair, where <basis> is set by 'LateralityBasis' (default 'signed',
%      i.e. the raw difference L-R; pass 'norm' to use LI_norm instead).
%   2. Group-level: Welch t-test of LI_norm per region pair x metric
%      (always on LI_norm, basis-independent).
%   3. Systematic direction: cohort-wide mean LI_norm and sign-rank vs
%      zero (always on LI_norm, basis-independent).
%
% Reads stabilityCentralities_*_<scheme>_results.mat from results/.
%
% Usage:
%   compareLeftRightAsymmetry;
%   compareLeftRightAsymmetry('Normalisation', 'column', 'Alpha', 0.05);
%   compareLeftRightAsymmetry('LateralityBasis', 'norm');
%
% Name-value options:
%   Normalisation   : 'column' (default)
%   ZThreshold      : 2.0   -- display-only (flagging uses Alpha + Bonferroni)
%   Alpha           : 0.05
%   LateralityBasis : 'signed' (default) | 'norm'
%                     'signed' -> per-Anderson z-score and flagging on
%                                 LI_signed = L - R (raw difference in
%                                 metric units). Matches the spec used in
%                                 the manuscript.
%                     'norm'   -> per-Anderson z-score and flagging on
%                                 LI_norm = (L - R)/(L + R).
%   SaveResults     : true
%   ResultsDir      : ''
%   RunParameters   : []
%   CasePrefix      : 'Anderson'
%   ControlPrefix   : 'Arnold'

setupMousePaths();

p = inputParser;
addParameter(p, 'Normalisation',   'column', @(s) ischar(s) || isstring(s));
addParameter(p, 'ZThreshold',      2.0,  @isscalar);
addParameter(p, 'Alpha',           0.05, @(x) isscalar(x) && x > 0 && x < 1);
addParameter(p, 'LateralityBasis', 'signed', @(s) ischar(s) || isstring(s));
addParameter(p, 'SaveResults',     true, @islogical);
addParameter(p, 'ResultsDir',      '',   @(s) ischar(s) || isstring(s));
addParameter(p, 'RunParameters',   [],   @(x) isempty(x) || isstruct(x));
addParameter(p, 'CasePrefix',      'Anderson', @(s) ischar(s) || isstring(s));
addParameter(p, 'ControlPrefix',   'Arnold',   @(s) ischar(s) || isstring(s));
parse(p, varargin{:});
opts = p.Results;
scheme = char(opts.Normalisation);
casePrefix = char(opts.CasePrefix);
controlPrefix = char(opts.ControlPrefix);

basis = lower(strtrim(char(opts.LateralityBasis)));
if ~any(strcmp(basis, {'norm', 'signed'}))
    error('compareLeftRightAsymmetry:BadBasis', ...
        'LateralityBasis must be ''signed'' or ''norm'' (got ''%s'').', basis);
end
opts.LateralityBasis = basis;
useSigned = strcmp(basis, 'signed');

resultsDir = resolveMouseResultsDir(opts.ResultsDir);
outDir = mouseResultsDir(resultsDir, 'lr_asymmetry');

%% Load per-mouse results
files = listPerMouseResultFiles(resultsDir, scheme);
if isempty(files)
    error('compareLeftRightAsymmetry:NoResults', ...
        'No stabilityCentralities_*_%s_results.mat found under %s.', ...
        scheme, resultsDir);
end

mice   = cell(0, 1);
loaded = cell(0, 1);
for k = 1:numel(files)
    s = load(fullfile(files(k).folder, files(k).name));
    if ~isfield(s, 'results'); continue; end
    r = s.results;
    if ~ischar(r.mouseId) && ~isstring(r.mouseId); continue; end
    mice{end+1, 1}   = char(r.mouseId); %#ok<AGROW>
    loaded{end+1, 1} = r;                %#ok<AGROW>
end
if isempty(mice)
    error('compareLeftRightAsymmetry:EmptyResults', ...
        'Loaded files contained no usable results structs.');
end

[isCase, isControl] = mouseCohortGroupMask(mice, casePrefix, controlPrefix);
isAnderson = isCase;
isArnold   = isControl;

if ~any(isCase)
    error('compareLeftRightAsymmetry:NoCaseMice', ...
        'No mice with case prefix "%s" found in %s (mice: %s).', ...
        casePrefix, resultsDir, strjoin(mice, ', '));
end
if ~any(isControl)
    error('compareLeftRightAsymmetry:NoControlMice', ...
        'No mice with control prefix "%s" found in %s (mice: %s).', ...
        controlPrefix, resultsDir, strjoin(mice, ', '));
end

labels = loaded{1}.labels;
N      = numel(labels);
for k = 2:numel(loaded)
    if numel(loaded{k}.labels) ~= N || ...
       ~all(string(loaded{k}.labels) == string(labels))
        warning('compareLeftRightAsymmetry:LabelMismatch', ...
            'Mouse %s labels differ from %s -- assuming index alignment.', ...
            mice{k}, mice{1});
    end
end

%% L-R pairing
[pairLeftIdx, pairRightIdx, pairLabels] = pairLRNodes(labels);
nPairs = numel(pairLeftIdx);
if nPairs == 0
    error('compareLeftRightAsymmetry:NoPairs', ...
        'No L-R region pairs found in labels.');
end

% Exclude pairs trivial in any Arnold mouse
trivialAcrossArnold = false(N, 1);
arnoldIdx = find(isArnold);
for ii = 1:numel(arnoldIdx)
    trivialAcrossArnold = trivialAcrossArnold | loaded{arnoldIdx(ii)}.isTrivial(:);
end
pairExcluded = trivialAcrossArnold(pairLeftIdx) | trivialAcrossArnold(pairRightIdx);
pairValid = ~pairExcluded;

pairLeftIdx  = pairLeftIdx(pairValid);
pairRightIdx = pairRightIdx(pairValid);
pairLabels   = pairLabels(pairValid);
nPairs = numel(pairLeftIdx);

%% Pack node-level metrics (N x nMice); mask indexes mice, not nodes
allMiceMask = true(numel(loaded), 1);
D_susc_all = packField(loaded, allMiceMask, 'D_susceptibility', N);
D_infl_all = packField(loaded, allMiceMask, 'D_influence',      N);
BC_all     = packCentrality(loaded, allMiceMask, 'betweenness', N);
hasBC = any(~isnan(BC_all(:)));
if ~hasBC
    warning('compareLeftRightAsymmetry:NoBC', ...
        'Betweenness not available; BC panels skipped.');
end

metricNames = {'D_susc', 'D_infl'};
nodeMats = {D_susc_all, D_infl_all};
ylabs = {'D(\rightarrow i)', 'D(k \rightarrow)'};
if hasBC
    metricNames{end+1} = 'BC';
    nodeMats{end+1} = BC_all;
    ylabs{end+1} = 'BC';
end
nMetrics = numel(metricNames);

%% Per-mouse laterality indices (nPairs x nMice)
LI_norm  = struct();
LI_signed = struct();
for m = 1:nMetrics
    mn = metricNames{m};
    M = nodeMats{m};
    LI_norm.(mn)  = computeLaterality(M, pairLeftIdx, pairRightIdx, 'norm');
    LI_signed.(mn) = computeLaterality(M, pairLeftIdx, pairRightIdx, 'signed');
end

%% Arnold cohort stats per pair x metric
LI_norm_arnold_mean  = struct();
LI_norm_arnold_std   = struct();
LI_signed_arnold_mean = struct();
LI_signed_arnold_std  = struct();
for m = 1:nMetrics
    mn = metricNames{m};
    Arn = LI_norm.(mn)(:, isArnold);
    LI_norm_arnold_mean.(mn)  = mean(Arn, 2, 'omitnan');
    LI_norm_arnold_std.(mn)   = std(Arn, 0, 2, 'omitnan');
    ArnS = LI_signed.(mn)(:, isArnold);
    LI_signed_arnold_mean.(mn) = mean(ArnS, 2, 'omitnan');
    LI_signed_arnold_std.(mn)  = std(ArnS, 0, 2, 'omitnan');
end

andersonNames = mice(isAnderson);
arnoldNames   = mice(isArnold);
nArnold       = numel(arnoldNames);

fprintf('compareLeftRightAsymmetry: %d Anderson, %d Arnold, %d L-R pairs.\n', ...
    numel(andersonNames), nArnold, nPairs);
fprintf('  Anderson: %s\n', strjoin(andersonNames, ', '));
fprintf('  Arnold  : %s\n', strjoin(arnoldNames, ', '));
fprintf('  Per-Anderson basis: LI_%s.\n', basis);

%% Per-Anderson z-score view (basis = LateralityBasis)
% LI_table / arnMean_table / arnStd_table point at the LI_norm or
% LI_signed structs depending on the chosen basis. Group-level Welch and
% the systematic sign-rank below always run on LI_norm, independent of
% this branch (see plan: docs/swap-LR-asymmetry-to-LI-signed_*.plan.md).
if useSigned
    LI_table       = LI_signed;
    arnMean_table  = LI_signed_arnold_mean;
    arnStd_table   = LI_signed_arnold_std;
    basisTex       = 'LI_{signed}';
    zeroLineLabel  = 'L - R = 0';
    clampToUnit    = false;
else
    LI_table       = LI_norm;
    arnMean_table  = LI_norm_arnold_mean;
    arnStd_table   = LI_norm_arnold_std;
    basisTex       = 'LI_{norm}';
    zeroLineLabel  = 'LI_{norm} = 0';
    clampToUnit    = true;
end

perAnderson = struct();
outlierTables = struct();
for a = 1:numel(andersonNames)
    aname = andersonNames{a};
    aIdx = find(strcmp(mice, aname), 1);
    devStruct = struct('mouseId', aname);

    nPanels = nMetrics;
    figHeight = 380 * nPanels + 80;
    fig = figure('Name', sprintf('L-R asymmetry %s vs Arnold', aname), ...
        'Position', [60 60 1700 figHeight]);

    outlierMaskAll = false(nPairs, nMetrics);
    zAll = cell(nMetrics, 1);
    pBonfAll = cell(nMetrics, 1);
    zAllVec = [];
    for m = 1:nMetrics
        mn = metricNames{m};
        aLI = LI_table.(mn)(:, aIdx);
        z = safeZ(aLI, arnMean_table.(mn), arnStd_table.(mn));
        zAll{m} = z;
        zAllVec = [zAllVec; z(:)]; %#ok<AGROW>
    end
    nTests = sum(~isnan(zAllVec));
    if nTests == 0; nTests = 1; end
    zThresh_bonf = norminv(1 - opts.Alpha / (2 * nTests));

    for m = 1:nMetrics
        mn = metricNames{m};
        aLI = LI_table.(mn)(:, aIdx);
        arnLI = LI_table.(mn)(:, isArnold);
        z = zAll{m};
        p = 2 * (1 - normcdf(abs(z)));
        pBonf = min(1, p * nTests);
        pBonfAll{m} = pBonf;
        outlierMask = ~isnan(pBonf) & pBonf < opts.Alpha;
        outlierMaskAll(:, m) = outlierMask;
        plotLRPanel(nPanels, m, aLI, arnLI, ...
            arnMean_table.(mn), arnStd_table.(mn), z, ...
            pairLabels, outlierMask, ...
            sprintf('%s (%s): %s vs Arnold (n=%d)  [Bonf |z|\\geq%.2f]', ...
                basisTex, ylabs{m}, strrep(aname, '_', '\_'), nArnold, zThresh_bonf), ...
            sprintf('%s (%s)', basisTex, ylabs{m}), ...
            clampToUnit, zeroLineLabel);
    end

    sgtitle(sprintf('Left-right laterality: %s vs Arnold  --  %s  (basis=%s)', ...
        strrep(aname, '_', '\_'), scheme, basis), 'Interpreter', 'tex');

    devStruct.basis = basis;
    devStruct.z = zAll;
    devStruct.p_bonf = pBonfAll;
    devStruct.outlierMask = outlierMaskAll;
    perAnderson.(matlab.lang.makeValidName(aname)) = devStruct;

    outMask = any(outlierMaskAll, 2);
    outIdx = find(outMask);
    if ~isempty(outIdx)
        outlierTables.(matlab.lang.makeValidName(aname)) = buildOutlierTable( ...
            outIdx, pairLabels, pairLeftIdx, pairRightIdx, ...
            LI_norm, LI_signed, LI_norm_arnold_mean, LI_norm_arnold_std, ...
            LI_signed_arnold_mean, LI_signed_arnold_std, ...
            metricNames, aIdx, zAll, pBonfAll, basis);
    else
        outlierTables.(matlab.lang.makeValidName(aname)) = table();
    end

    if opts.SaveResults
        baseName = sprintf('compare_LR_asymmetry_%s_%s', aname, scheme);
        savefig(fig, fullfile(outDir, [baseName '.fig']));
        try
            exportgraphics(fig, fullfile(outDir, [baseName '.png']), 'Resolution', 200);
        catch
            saveas(fig, fullfile(outDir, [baseName '.png']));
        end
        close(fig);

        if ~isempty(outIdx)
            T = outlierTables.(matlab.lang.makeValidName(aname));
            writetable(T, fullfile(outDir, [baseName '_outliers.csv']));
            fprintf('  %s: %d outlier pair(s) -> %s_outliers.csv\n', ...
                aname, numel(outIdx), baseName);
        else
            fprintf('  %s: no Bonferroni outliers at alpha=%.3f.\n', aname, opts.Alpha);
        end
    else
        close(fig);
    end
end

%% Group-level Welch t-tests (per pair x metric)
groupTest = struct();
groupRows = {};
canWelch = exist('ttest2', 'file') == 2;
pValAll = [];

for m = 1:nMetrics
    mn = metricNames{m};
    AndLI = LI_norm.(mn)(:, isAnderson);
    ArnLI = LI_norm.(mn)(:, isArnold);

    tStat = NaN(nPairs, 1);
    dfVal = NaN(nPairs, 1);
    pVal  = NaN(nPairs, 1);
    mAnd  = mean(AndLI, 2, 'omitnan');
    mArn  = mean(ArnLI, 2, 'omitnan');

    for p = 1:nPairs
        xA = AndLI(p, :)';
        xR = ArnLI(p, :)';
        xA = xA(~isnan(xA));
        xR = xR(~isnan(xR));
        if canWelch && numel(xA) >= 2 && numel(xR) >= 2
            [~, pVal(p), ~, st] = ttest2(xA, xR, 'Vartype', 'unequal');
            tStat(p) = st.tstat;
            dfVal(p) = st.df;
        end
    end

    groupTest.(mn) = struct( ...
        'tStat', tStat, 'df', dfVal, 'p', pVal, ...
        'LI_anderson_mean', mAnd, 'LI_arnold_mean', mArn);
    pValAll = [pValAll; pVal(:)]; %#ok<AGROW>
end

nGT = sum(~isnan(pValAll));
if nGT == 0; nGT = 1; end
for m = 1:nMetrics
    mn = metricNames{m};
    pVal = groupTest.(mn).p;
    pBonf = min(1, pVal * nGT);
    bonfSig = ~isnan(pBonf) & pBonf < opts.Alpha;
    groupTest.(mn).p_bonf = pBonf;
    groupTest.(mn).bonf_significant = bonfSig;
    mAnd = groupTest.(mn).LI_anderson_mean;
    mArn = groupTest.(mn).LI_arnold_mean;
    tStat = groupTest.(mn).tStat;
    dfVal = groupTest.(mn).df;
    for p = 1:nPairs
        groupRows{end+1, 1} = {pairLabels{p}, mn, mAnd(p), mArn(p), ...
            tStat(p), dfVal(p), pVal(p), pBonf(p), bonfSig(p)}; %#ok<AGROW>
    end
end

Tgroup = cell2table(vertcat(groupRows{:}), ...
    'VariableNames', {'Region', 'Metric', 'LI_anderson_mean', 'LI_arnold_mean', ...
    't_stat', 'df', 'p', 'p_bonf', 'bonf_significant'});

%% Systematic direction (sign-rank vs zero per cohort x metric)
systematic = struct();
sysRows = {};
for m = 1:nMetrics
    mn = metricNames{m};
    for cohortTag = {'Anderson', 'Arnold'}
        cname = cohortTag{1};
        if strcmp(cname, 'Anderson')
            mask = isAnderson;
        else
            mask = isArnold;
        end
        vals = LI_norm.(mn)(:, mask);
        v = vals(:);
        v = v(~isnan(v));
        nObs = numel(v);
        mLI = mean(v, 'omitnan');
        sLI = std(v, 0, 'omitnan');
        pSR = NaN;
        if exist('signrank', 'file') == 2 && nObs >= 1
            try
                pSR = signrank(v, 0, 'method', 'approximate');
            catch
                pSR = NaN;
            end
        end
        domSide = inferDominantSide(mLI, pSR, opts.Alpha);
        systematic.(mn).(cname) = struct( ...
            'n_obs', nObs, 'mean_LI', mLI, 'sd_LI', sLI, ...
            'signrank_p', pSR, 'dominant_side', domSide);
        sysRows{end+1, 1} = {cname, mn, nObs, mLI, sLI, pSR, domSide}; %#ok<AGROW>
    end
end
Tsys = cell2table(vertcat(sysRows{:}), ...
    'VariableNames', {'cohort', 'metric', 'n_obs', 'mean_LI', 'sd_LI', ...
    'signrank_p', 'inferred_dominant_side'});

%% Group-test and systematic figures
if opts.SaveResults
    plotGroupTestFigure(outDir, scheme, pairLabels, metricNames, ylabs, ...
        LI_norm, isAnderson, isArnold, groupTest, opts.Alpha);
    plotSystematicFigure(outDir, scheme, metricNames, ylabs, systematic, Tsys);
    writetable(Tgroup, fullfile(outDir, sprintf('LR_asymmetry_groupTest_%s.csv', scheme)));
    writetable(Tsys, fullfile(outDir, sprintf('LR_asymmetry_systematic_%s.csv', scheme)));
    fprintf('Wrote LR_asymmetry_groupTest_%s.csv, LR_asymmetry_systematic_%s.csv\n', ...
        scheme, scheme);
end

%% Summary struct
summary = struct();
summary.scheme = scheme;
summary.mice = mice;
summary.casePrefix = casePrefix;
summary.controlPrefix = controlPrefix;
summary.isCase = isCase;
summary.isControl = isControl;
summary.isAnderson = isAnderson;
summary.isArnold = isArnold;
summary.labels = labels;
summary.pairLabels = pairLabels;
summary.pairLeftIdx = pairLeftIdx;
summary.pairRightIdx = pairRightIdx;
summary.pairExcluded = pairExcluded;
summary.metricNames = metricNames;
summary.hasBC = hasBC;
summary.LI_norm = LI_norm;
summary.LI_signed = LI_signed;
summary.LI_norm_arnold_mean = LI_norm_arnold_mean;
summary.LI_norm_arnold_std = LI_norm_arnold_std;
summary.LI_signed_arnold_mean = LI_signed_arnold_mean;
summary.LI_signed_arnold_std = LI_signed_arnold_std;
summary.perAnderson = perAnderson;
summary.outlierTables = outlierTables;
summary.groupTest = groupTest;
summary.systematic = systematic;
summary.zThreshold = opts.ZThreshold;
summary.alpha = opts.Alpha;
summary.lateralityBasis = basis;

if isempty(opts.RunParameters)
    summary.runParameters = mouseExperimentRunParameters('buildFromCompareLR', ...
        opts, scheme, resultsDir);
else
    summary.runParameters = mouseExperimentRunParameters('merge', opts.RunParameters, ...
        mouseExperimentRunParameters('buildFromCompareLR', opts, scheme, resultsDir));
end

if opts.SaveResults
    save(fullfile(outDir, sprintf('LR_asymmetry_%s.mat', scheme)), ...
        '-struct', 'summary');
end

end

%% ==================================================================
function [leftIdx, rightIdx, pairLabels] = pairLRNodes(labels)
%PAIRLRNODES  Match L-<suffix> with R-<suffix> by shared suffix.
labels = string(labels);
n = numel(labels);
leftIdx = [];
rightIdx = [];
pairLabels = strings(0, 1);

for i = 1:n
    lab = labels(i);
    if ~startsWith(lab, "L-")
        continue;
    end
    suffix = extractAfter(lab, "L-");
    rLab = "R-" + suffix;
    j = find(labels == rLab, 1);
    if isempty(j)
        warning('pairLRNodes:NoMatch', 'No R pair for %s', lab);
        continue;
    end
    leftIdx(end+1, 1) = i;      %#ok<AGROW>
    rightIdx(end+1, 1) = j;     %#ok<AGROW>
    pairLabels(end+1, 1) = suffix; %#ok<AGROW>
end
leftIdx = leftIdx(:);
rightIdx = rightIdx(:);
pairLabels = cellstr(pairLabels);
end

%% ------------------------------------------------------------------
function LI = computeLaterality(M, leftIdx, rightIdx, mode)
%COMPUTELATERALITY  Per pair: LI from left/right columns of M (N x nMice).
nPairs = numel(leftIdx);
nMice = size(M, 2);
LI = NaN(nPairs, nMice);
for p = 1:nPairs
    L = M(leftIdx(p), :);
    R = M(rightIdx(p), :);
    switch mode
        case 'norm'
            denom = L + R;
            ok = denom > 0 & ~isnan(L) & ~isnan(R);
            LI(p, ok) = (L(ok) - R(ok)) ./ denom(ok);
        case 'signed'
            ok = ~isnan(L) & ~isnan(R);
            LI(p, ok) = L(ok) - R(ok);
        otherwise
            error('computeLaterality:BadMode', 'mode must be norm or signed.');
    end
end
end

%% ------------------------------------------------------------------
function T = buildOutlierTable(outIdx, pairLabels, pairLeftIdx, pairRightIdx, ...
    LI_norm, LI_signed, muN, sdN, muS, sdS, metricNames, aIdx, zAll, pBonfAll, ...
    basis)
% basis : 'signed' -> z/p_bonf columns named z_signed_<mn> / p_bonf_signed_<mn>
%         'norm'   -> z/p_bonf columns named z_<mn>        / p_bonf_<mn>
% LI_norm_* and LI_signed_* raw columns are emitted in both branches.
if nargin < 15 || isempty(basis)
    basis = 'norm';
end
if strcmp(basis, 'signed')
    zPrefix     = 'z_signed';
    pBonfPrefix = 'p_bonf_signed';
else
    zPrefix     = 'z';
    pBonfPrefix = 'p_bonf';
end

nOut = numel(outIdx);
nMet = numel(metricNames);
cols = {'PairIdx', 'Region', 'L_NodeIdx', 'R_NodeIdx'};
data = {outIdx, string(pairLabels(outIdx)), pairLeftIdx(outIdx), pairRightIdx(outIdx)};

for m = 1:nMet
    mn = metricNames{m};
    aLI = LI_norm.(mn)(outIdx, aIdx);
    aLIs = LI_signed.(mn)(outIdx, aIdx);
    data = [data, {aLI, muN.(mn)(outIdx), sdN.(mn)(outIdx), ...
        aLIs, muS.(mn)(outIdx), sdS.(mn)(outIdx), ...
        zAll{m}(outIdx), pBonfAll{m}(outIdx)}]; %#ok<AGROW>
    cols = [cols, { ...
        sprintf('LI_norm_anderson_%s', mn), sprintf('LI_norm_arnold_mean_%s', mn), ...
        sprintf('LI_norm_arnold_std_%s', mn), ...
        sprintf('LI_signed_anderson_%s', mn), sprintf('LI_signed_arnold_mean_%s', mn), ...
        sprintf('LI_signed_arnold_std_%s', mn), ...
        sprintf('%s_%s', zPrefix, mn), sprintf('%s_%s', pBonfPrefix, mn)}]; %#ok<AGROW>
end

T = table(data{:}, 'VariableNames', cols);
minP = inf(nOut, 1);
for m = 1:nMet
    pcol = T.(sprintf('%s_%s', pBonfPrefix, metricNames{m}));
    minP = min(minP, pcol, 'omitnan');
end
[~, ord] = sort(minP, 'ascend');
T = T(ord, :);
end

%% ------------------------------------------------------------------
function plotLRPanel(nPanels, panelIdx, anderVals, arnoldVals, ...
    arnMean, arnStd, zScore, labels, outlierMask, panelTitle, ylab, ...
    clampToUnit, zeroLineLabel)
% clampToUnit   : true to clamp the y-axis to [-1.05, 1.05] (only sensible
%                 for LI_norm, which is bounded). LI_signed is in the
%                 metric's native units and must not be clamped.
% zeroLineLabel : DisplayName for the y=0 reference line.

if nargin < 12 || isempty(clampToUnit)
    clampToUnit = true;
end
if nargin < 13 || isempty(zeroLineLabel)
    zeroLineLabel = 'LI = 0';
end

ax = subplot(nPanels, 1, panelIdx); hold(ax, 'on');
nP = numel(anderVals);
x = 1:nP;

yl = [-1.05 1.05];
if any(isfinite(anderVals)) || any(isfinite(arnMean))
    yl = [min([anderVals(:); arnMean - arnStd; arnoldVals(:)], [], 'omitnan'), ...
          max([anderVals(:); arnMean + arnStd; arnoldVals(:)], [], 'omitnan')];
    if ~all(isfinite(yl))
        yl = [-1 1];
    end
    pad = 0.05 * max(diff(yl), 0.1);
    yl = yl + [-pad pad];
    if clampToUnit
        yl(1) = max(yl(1), -1.05);
        yl(2) = min(yl(2), 1.05);
    end
end

yline(ax, 0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, ...
    'DisplayName', zeroLineLabel);

valid = ~isnan(arnMean) & ~isnan(arnStd);
if any(valid)
    xv = x(valid);
    yLo = arnMean(valid) - arnStd(valid);
    yHi = arnMean(valid) + arnStd(valid);
    fill(ax, [xv fliplr(xv)], [yLo' fliplr(yHi')], ...
        [0.80 0.84 0.95], 'EdgeColor', 'none', 'FaceAlpha', 0.65, ...
        'DisplayName', 'Arnold mean \pm SD');
end

plot(ax, x, arnMean, '-', 'Color', [0.20 0.30 0.70], ...
    'LineWidth', 1.2, 'DisplayName', 'Arnold mean');

for j = 1:size(arnoldVals, 2)
    plot(ax, x, arnoldVals(:, j), '.', 'Color', [0.55 0.60 0.85], ...
        'MarkerSize', 6, 'HandleVisibility', 'off');
end

plot(ax, x, anderVals, 'o', 'Color', [0.65 0.10 0.10], ...
    'MarkerFaceColor', [0.85 0.33 0.10], 'MarkerSize', 5, ...
    'DisplayName', 'Anderson');

outlierIdx = find(outlierMask & ~isnan(zScore));
if ~isempty(outlierIdx)
    plot(ax, x(outlierIdx), anderVals(outlierIdx), 'o', ...
        'Color', 'k', 'MarkerSize', 10, 'LineWidth', 1.5, ...
        'HandleVisibility', 'off');
    for k = 1:numel(outlierIdx)
        i = outlierIdx(k);
        text(ax, i, anderVals(i), ...
            sprintf('  %s (z=%+.1f)', labels{i}, zScore(i)), ...
            'Interpreter', 'none', 'FontSize', 7, 'Rotation', 25, ...
            'VerticalAlignment', 'bottom');
    end
end

xlim(ax, [0.5 nP + 0.5]);
ylim(ax, yl);
set(ax, 'XTick', 1:nP, 'XTickLabel', labels, ...
    'TickLabelInterpreter', 'none', 'FontSize', 6);
xtickangle(ax, 45);
grid(ax, 'on'); box(ax, 'on');
ylabel(ax, ylab, 'Interpreter', 'tex');
title(ax, panelTitle, 'Interpreter', 'tex');
legend(ax, 'Location', 'best');
end

%% ------------------------------------------------------------------
function plotGroupTestFigure(outDir, scheme, pairLabels, metricNames, ylabs, ...
    LI_norm, isAnderson, isArnold, groupTest, alpha)

nMet = numel(metricNames);
fig = figure('Name', sprintf('L-R group test (%s)', scheme), ...
    'Position', [80 80 520 * nMet 480]);
nPairs = numel(pairLabels);

for m = 1:nMet
    mn = metricNames{m};
    ax = subplot(1, nMet, m); hold(ax, 'on');
    gt = groupTest.(mn);

    AndLI = LI_norm.(mn)(:, isAnderson);
    ArnLI = LI_norm.(mn)(:, isArnold);
    xAll = [AndLI(:); ArnLI(:)];
    xAll = xAll(~isnan(xAll));
    if isempty(xAll)
        xlims = [-1 1];
    else
        pad = 0.1 * max(range(xAll), 0.1);
        xlims = [min(xAll) - pad, max(xAll) + pad];
    end

    for p = 1:nPairs
        xA = mean(AndLI(p, :), 'omitnan');
        xR = mean(ArnLI(p, :), 'omitnan');
        if isnan(xA) || isnan(xR)
            continue;
        end
        col = [0.4 0.4 0.4];
        mk = 4;
        if gt.bonf_significant(p)
            col = [0.85 0.2 0.1];
            mk = 8;
        end
        plot(ax, xR, xA, 'o', 'Color', col, 'MarkerFaceColor', col, ...
            'MarkerSize', mk, 'HandleVisibility', 'off');
    end

    plot(ax, xlims, xlims, 'k--', 'HandleVisibility', 'off');
    xlabel(ax, sprintf('Arnold mean LI (%s)', ylabs{m}), 'Interpreter', 'tex');
    ylabel(ax, sprintf('Anderson mean LI (%s)', ylabs{m}), 'Interpreter', 'tex');
    title(ax, sprintf('%s: pair means (red = Bonf p<%.2f)', ylabs{m}, alpha), ...
        'Interpreter', 'tex');
    axis(ax, 'equal');
    xlim(ax, xlims);
    ylim(ax, xlims);
    grid(ax, 'on');
end

sgtitle(sprintf('Group-level L-R laterality: Anderson vs Arnold  --  %s', scheme), ...
    'Interpreter', 'tex');

base = fullfile(outDir, sprintf('LR_asymmetry_groupTest_%s', scheme));
savefig(fig, [base '.fig']);
try
    exportgraphics(fig, [base '.png'], 'Resolution', 200);
catch
    saveas(fig, [base '.png']);
end
close(fig);
end

%% ------------------------------------------------------------------
function plotSystematicFigure(outDir, scheme, metricNames, ylabs, systematic, Tsys)

nMet = numel(metricNames);
fig = figure('Name', sprintf('L-R systematic (%s)', scheme), ...
    'Position', [100 100 720 420]);
ax = axes(fig); hold(ax, 'on');

cohorts = {'Anderson', 'Arnold'};
cohortColors = [0.85 0.33 0.10; 0.20 0.30 0.70];
xBase = 1:nMet;
barW = 0.35;

for c = 1:2
    cname = cohorts{c};
    means = zeros(nMet, 1);
    sds = zeros(nMet, 1);
    for m = 1:nMet
        mn = metricNames{m};
        if isfield(systematic, mn) && isfield(systematic.(mn), cname)
            means(m) = systematic.(mn).(cname).mean_LI;
            sds(m) = systematic.(mn).(cname).sd_LI;
        else
            means(m) = NaN;
            sds(m) = NaN;
        end
    end
    xPos = xBase + (c - 1.5) * barW;
    errorbar(ax, xPos, means, sds, 'o', 'Color', cohortColors(c, :), ...
        'MarkerFaceColor', cohortColors(c, :), 'LineWidth', 1.2, ...
        'CapSize', 8, 'DisplayName', cname);
end

yline(ax, 0, '--', 'Color', [0.5 0.5 0.5]);
set(ax, 'XTick', xBase, 'XTickLabel', ylabs, 'TickLabelInterpreter', 'tex');
ylabel(ax, 'Mean LI_{norm} (all pairs \times mice)');
title(ax, sprintf('Systematic L-R bias by cohort  --  %s', scheme), 'Interpreter', 'tex');
legend(ax, 'Location', 'best');
grid(ax, 'on');

% Annotate sign-rank p-values
for m = 1:nMet
    mn = metricNames{m};
    for c = 1:2
        cname = cohorts{c};
        row = Tsys(strcmp(Tsys.cohort, cname) & strcmp(Tsys.metric, mn), :);
        if ~isempty(row) && isfinite(row.signrank_p(1))
            xPos = xBase(m) + (c - 1.5) * barW;
            yPos = row.mean_LI(1) + row.sd_LI(1) + 0.02;
            text(ax, xPos, yPos, sprintf('p=%.3f', row.signrank_p(1)), ...
                'FontSize', 7, 'HorizontalAlignment', 'center');
        end
    end
end

base = fullfile(outDir, sprintf('LR_asymmetry_systematic_%s', scheme));
savefig(fig, [base '.fig']);
try
    exportgraphics(fig, [base '.png'], 'Resolution', 200);
catch
    saveas(fig, [base '.png']);
end
close(fig);
end

%% ------------------------------------------------------------------
function side = inferDominantSide(meanLI, pSR, alpha)
if isnan(meanLI) || isnan(pSR) || pSR >= alpha
    side = 'none';
elseif meanLI > 0
    side = 'L';
else
    side = 'R';
end
end

%% ------------------------------------------------------------------
function M = packField(loaded, mask, fieldName, N)
idx = find(mask);
M = NaN(N, numel(idx));
for ii = 1:numel(idx)
    v = loaded{idx(ii)}.(fieldName);
    M(:, ii) = v(:);
end
end

%% ------------------------------------------------------------------
function M = packCentrality(loaded, mask, fieldName, N)
idx = find(mask);
M = NaN(N, numel(idx));
for ii = 1:numel(idx)
    r = loaded{idx(ii)};
    if ~isfield(r, 'centralities') || ~isstruct(r.centralities)
        continue;
    end
    if ~isfield(r.centralities, fieldName)
        continue;
    end
    v = r.centralities.(fieldName);
    if numel(v) ~= N
        continue;
    end
    M(:, ii) = v(:);
end
end

%% ------------------------------------------------------------------
function z = safeZ(x, mu, sigma)
z = NaN(size(x));
ok = ~isnan(sigma) & sigma > 0 & ~isnan(mu) & ~isnan(x);
z(ok) = (x(ok) - mu(ok)) ./ sigma(ok);
end
