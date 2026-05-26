function rows = collectSuscOutlierRows(andersonSummary)
%COLLECTSUSCOUTLIERROWS  Per-Anderson Bonferroni-significant D(->i) outliers as struct rows.
%
%   rows(k) fields: mouseId, region, baseName, hemi, eff, z, D_and, D_arn_mean

P = struct(); %#ok<NASGU>
alpha = andersonSummary.alpha;
labels = andersonSummary.labels;
andersonNames = andersonSummary.mice(andersonSummary.isAnderson);
mu = andersonSummary.D_susc_arnold_mean;
sd = andersonSummary.D_susc_arnold_std; %#ok<ASGLU>

rows = struct('mouseId', {}, 'region', {}, 'baseName', {}, 'hemi', {}, ...
    'eff', {}, 'z', {}, 'D_and', {}, 'D_arn_mean', {});

for a = 1:numel(andersonNames)
    aname = andersonNames{a};
    fld = matlab.lang.makeValidName(aname);
    dev = andersonSummary.deviations.(fld);
    mask = dev.p_bonf_susc < alpha & ~isnan(dev.p_bonf_susc);
    if isfield(andersonSummary, 'trivialAcrossArnold')
        mask = mask & ~andersonSummary.trivialAcrossArnold;
    end
    idx = find(mask);
    for k = 1:numel(idx)
        i = idx(k);
        dArn = mu(i);
        dAnd = dev.D_susc_anderson(i);
        if isnan(dArn) || dArn == 0
            eff = NaN;
        else
            eff = 100 * (dAnd - dArn) / dArn;
        end
        [baseName, hemi] = stripHemispherePrefix(labels{i});
        rows(end+1).mouseId = aname; %#ok<AGROW>
        rows(end).region = char(labels{i});
        rows(end).baseName = baseName;
        rows(end).hemi = hemi;
        rows(end).eff = eff;
        rows(end).z = dev.z_susc(i);
        rows(end).D_and = dAnd;
        rows(end).D_arn_mean = dArn;
    end
end
end
