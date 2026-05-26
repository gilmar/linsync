function rows = collectSuscOutlierRows(caseSummary)
%COLLECTSUSCOUTLIERROWS  Per-case-mouse Bonferroni-significant D(->i) outliers as struct rows.
%
%   rows(k) fields: mouseId, region, baseName, hemi, eff, z, D_case, D_control_mean

info = cohortGroupInfo(caseSummary);
alpha = info.alpha;
labels = caseSummary.labels;
caseIds = info.caseIds;
mu = caseSummary.D_susc_arnold_mean;
sd = caseSummary.D_susc_arnold_std; %#ok<ASGLU>

rows = struct('mouseId', {}, 'region', {}, 'baseName', {}, 'hemi', {}, ...
    'eff', {}, 'z', {}, 'D_case', {}, 'D_control_mean', {});

for a = 1:numel(caseIds)
    aname = caseIds{a};
    fld = matlab.lang.makeValidName(aname);
    dev = caseSummary.deviations.(fld);
    mask = dev.p_bonf_susc < alpha & ~isnan(dev.p_bonf_susc);
    if isfield(caseSummary, 'trivialAcrossArnold')
        mask = mask & ~caseSummary.trivialAcrossArnold;
    end
    idx = find(mask);
    for k = 1:numel(idx)
        i = idx(k);
        dCtrl = mu(i);
        dCase = dev.D_susc_anderson(i);
        if isnan(dCtrl) || dCtrl == 0
            eff = NaN;
        else
            eff = 100 * (dCase - dCtrl) / dCtrl;
        end
        [baseName, hemi] = stripHemispherePrefix(labels{i});
        rows(end+1).mouseId = aname; %#ok<AGROW>
        rows(end).region = char(labels{i});
        rows(end).baseName = baseName;
        rows(end).hemi = hemi;
        rows(end).eff = eff;
        rows(end).z = dev.z_susc(i);
        rows(end).D_case = dCase;
        rows(end).D_control_mean = dCtrl;
    end
end
end
