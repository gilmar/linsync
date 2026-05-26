function content = deriveRobustSuggestiveContent(andersonSummary, lrSummary, overrideCfg)
%DERIVEROBUSTSUGGESTIVECONTENT  ROBUST vs SUGGESTIVE buckets from comparison summaries.

if nargin < 3
    overrideCfg = struct();
end
alpha = andersonSummary.alpha;

%% ROBUST — p_bonf_susc < alpha in both Anderson mice (same base region)
suscRows = collectSuscOutlierRows(andersonSummary);
robustBases = cohortRobustBases(suscRows);
robustItems = {};
for b = 1:numel(robustBases)
    bn = robustBases{b};
    effs = [suscRows(strcmp({suscRows.baseName}, bn)).eff];
    effs = effs(isfinite(effs) & effs > 0);
    effs = sort(effs, 'descend');
    if numel(effs) >= 2
        effStr = sprintf('+%.1f%% / +%.1f%%', effs(1), effs(2));
    elseif numel(effs) == 1
        effStr = sprintf('+%.1f%%', effs(1));
    else
        effStr = '';
    end
    robustItems{end+1} = {prettyRegionName(bn), ...
        sprintf('Susceptibility D(\\rightarrow i)  %s', effStr)}; %#ok<AGROW>
end

%% SUGGESTIVE
suggestItems = {};

% (i) Laterality D(k->): Bonferroni-significant in exactly one Anderson mouse
lateralityRegs = lateralityOneMouseRegions(lrSummary, alpha);
if ~isempty(lateralityRegs)
    regList = strjoin(lateralityRegs, ', ');
    if numel(lateralityRegs) >= 2
        hdr = 'Laterality (motor, olfactory) — one mouse only';
    else
        hdr = sprintf('Laterality (%s) — one mouse only', lateralityRegs{1});
    end
    suggestItems{end+1} = {hdr, sprintf('Influence D(k \\rightarrow): %s', regList)}; %#ok<AGROW>
end

% (ii) Betweenness-only (BC Bonferroni, D-channels not)
bcRegs = betweennessOnlyRegions(andersonSummary, alpha);
if ~isempty(bcRegs)
    suggestItems{end+1} = {'Betweenness hubs (amygdala, cranial nerves) — atlas-sensitive', ...
        strjoin(bcRegs, ', ')}; %#ok<AGROW>
end

% (iii) Single-mouse susceptibility outliers
singleSusc = singleMouseSuscRegions(suscRows, robustBases);
if ~isempty(singleSusc)
    if numel(singleSusc) >= 2
        hdr = 'Single-mouse susceptibility (somatosensory, temporal-assoc)';
    else
        hdr = sprintf('Single-mouse susceptibility (%s)', singleSusc{1});
    end
    suggestItems{end+1} = {hdr, strjoin(singleSusc, ', ')}; %#ok<AGROW>
end

content = struct();
content.robust = struct('items', {robustItems});
content.suggestive = struct('items', {suggestItems});

if isfield(overrideCfg, 'robust')
    content.robust = overrideCfg.robust;
end
if isfield(overrideCfg, 'suggestive')
    content.suggestive = overrideCfg.suggestive;
end
end

%% ------------------------------------------------------------------
function bases = cohortRobustBases(suscRows)
baseNames = {suscRows.baseName};
allBases = unique(baseNames, 'stable');
bases = {};
for b = 1:numel(allBases)
    bn = allBases{b};
    miceHere = unique({suscRows(strcmp(baseNames, bn)).mouseId});
    if numel(miceHere) >= 2
        bases{end+1} = bn; %#ok<AGROW>
    end
end
end

%% ------------------------------------------------------------------
function regs = lateralityOneMouseRegions(lrSummary, alpha)
regs = {};
if nargin < 2 || isempty(lrSummary)
    return;
end
rows = collectLateralityDinflRows(lrSummary, alpha);
if isempty(rows)
    return;
end
pairSig = containers.Map('KeyType', 'char', 'ValueType', 'any');
for k = 1:numel(rows)
    reg = rows(k).region;
    mid = rows(k).mouseId;
    if pairSig.isKey(reg)
        pairSig(reg) = unique([pairSig(reg), {mid}], 'stable');
    else
        pairSig(reg) = {mid};
    end
end
keysLR = pairSig.keys;
for k = 1:numel(keysLR)
    if numel(pairSig(keysLR{k})) == 1
        regs{end+1} = prettyRegionName(keysLR{k}); %#ok<AGROW>
    end
end
regs = unique(regs, 'stable');
end

%% ------------------------------------------------------------------
function regs = betweennessOnlyRegions(andersonSummary, alpha)
regs = {};
andersonNames = andersonSummary.mice(andersonSummary.isAnderson);
for a = 1:numel(andersonNames)
    aname = andersonNames{a};
    fld = matlab.lang.makeValidName(aname);
    dev = andersonSummary.deviations.(fld);
    if ~dev.hasBC
        continue;
    end
    bcOnly = dev.p_bonf_BC < alpha & dev.p_bonf_susc >= alpha & dev.p_bonf_infl >= alpha;
    bcOnly = bcOnly & ~isnan(dev.p_bonf_BC);
    if isfield(andersonSummary, 'trivialAcrossArnold')
        bcOnly = bcOnly & ~andersonSummary.trivialAcrossArnold;
    end
    idx = find(bcOnly);
    for k = 1:numel(idx)
        [bn, ~] = stripHemispherePrefix(andersonSummary.labels{idx(k)});
        regs{end+1} = prettyRegionName(bn); %#ok<AGROW>
    end
end
regs = unique(regs, 'stable');
end

%% ------------------------------------------------------------------
function regs = singleMouseSuscRegions(suscRows, robustBases)
baseNames = {suscRows.baseName};
regs = {};
for k = 1:numel(suscRows)
    bn = suscRows(k).baseName;
    if any(strcmp(robustBases, bn))
        continue;
    end
    miceHere = unique({suscRows(strcmp(baseNames, bn)).mouseId});
    if numel(miceHere) == 1
        regs{end+1} = prettyRegionName(bn); %#ok<AGROW>
    end
end
regs = unique(regs, 'stable');
end

%% ------------------------------------------------------------------
function s = prettyRegionName(baseName)
s = strrep(char(baseName), '_', ' ');
s = lower(s);
end
