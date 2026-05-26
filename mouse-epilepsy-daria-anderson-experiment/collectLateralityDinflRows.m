function rows = collectLateralityDinflRows(lrSummary, alpha)
%COLLECTLATERALITYDINFLROWS  Bonferroni-significant signed D(k->) laterality pairs.

rows = struct('mouseId', {}, 'region', {}, 'raw', {}, 'arn_mean', {}, ...
    'arn_std', {}, 'z', {});

if isfield(lrSummary, 'outlierTables') && ~isempty(fieldnames(lrSummary.outlierTables))
    tblFields = fieldnames(lrSummary.outlierTables);
    info = cohortGroupInfo(lrSummary);
    for m = 1:numel(tblFields)
        fld = tblFields{m};
        T = lrSummary.outlierTables.(fld);
        if isempty(T)
            continue;
        end
        mouseId = resolveMouseIdFromField(fld, info.caseIds);
        pCol = 'p_bonf_signed_D_infl';
        if ~ismember(pCol, T.Properties.VariableNames)
            continue;
        end
        sig = T.(pCol) < alpha;
        for r = 1:height(T)
            if ~sig(r)
                continue;
            end
            rows(end+1) = rowFromTable(T, r, mouseId); %#ok<AGROW>
        end
    end
    if ~isempty(rows)
        return;
    end
end

% Legacy .mat: scan per-case z/p on D_infl signed LI
if ~isfield(lrSummary, 'perAnderson') || ~isfield(lrSummary, 'LI_signed')
    return;
end
mn = 'D_infl';
if ~isfield(lrSummary.LI_signed, mn)
    return;
end
info = cohortGroupInfo(lrSummary);
mice = lrSummary.mice;
caseIds = info.caseIds;
muS = lrSummary.LI_signed_arnold_mean.(mn);
sdS = lrSummary.LI_signed_arnold_std.(mn);
pairLabels = lrSummary.pairLabels;
metricNames = lrSummary.metricNames;
mInfl = find(strcmp(metricNames, mn), 1);

for a = 1:numel(caseIds)
    aname = caseIds{a};
    fld = matlab.lang.makeValidName(aname);
    aIdx = find(strcmp(mice, aname), 1);
    if ~isfield(lrSummary.perAnderson, fld)
        continue;
    end
    pa = lrSummary.perAnderson.(fld);
    if isempty(mInfl) || numel(pa.p_bonf) < mInfl
        continue;
    end
    pBonf = pa.p_bonf{mInfl};
    z = pa.z{mInfl};
    aLI = lrSummary.LI_signed.(mn)(:, aIdx);
    sig = pBonf < alpha & ~isnan(pBonf);
    idx = find(sig);
    for k = 1:numel(idx)
        p = idx(k);
        rows(end+1).mouseId = aname; %#ok<AGROW>
        rows(end).region = char(pairLabels{p});
        rows(end).raw = aLI(p);
        rows(end).arn_mean = muS(p);
        rows(end).arn_std = sdS(p);
        rows(end).z = z(p);
    end
end
end

%% ------------------------------------------------------------------
function mouseId = resolveMouseIdFromField(fld, caseIds)
mouseId = fld;
for k = 1:numel(caseIds)
    if strcmp(matlab.lang.makeValidName(caseIds{k}), fld)
        mouseId = caseIds{k};
        return;
    end
end
end

%% ------------------------------------------------------------------
function row = rowFromTable(T, r, mouseId)
row.mouseId = mouseId;
row.region = char(T.Region(r));
row.raw = T.LI_signed_anderson_D_infl(r);
row.arn_mean = T.LI_signed_arnold_mean_D_infl(r);
row.arn_std = T.LI_signed_arnold_std_D_infl(r);
row.z = T.z_signed_D_infl(r);
end
