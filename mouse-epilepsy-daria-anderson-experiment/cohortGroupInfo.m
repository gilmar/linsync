function info = cohortGroupInfo(summary)
%COHORTGROUPINFO  Case/control cohort metadata from a comparison summary struct.
%
%   info = cohortGroupInfo(andersonSummary)
%
% Accepts summaries from compareAndersonVsArnold / compareLeftRightAsymmetry.
% Uses isCase/isControl when present; otherwise isAnderson/isArnold.

mice = summary.mice(:);
if isfield(summary, 'isCase')
    isCase = summary.isCase(:);
    isControl = summary.isControl(:);
elseif isfield(summary, 'isAnderson')
    isCase = summary.isAnderson(:);
    isControl = summary.isArnold(:);
else
    error('cohortGroupInfo:NoMask', 'Summary must include isCase or isAnderson masks.');
end

if isfield(summary, 'casePrefix') && ~isempty(summary.casePrefix)
    info.casePrefix = char(summary.casePrefix);
else
    info.casePrefix = inferPrefixFromMice(mice(isCase));
end
if isfield(summary, 'controlPrefix') && ~isempty(summary.controlPrefix)
    info.controlPrefix = char(summary.controlPrefix);
else
    info.controlPrefix = inferPrefixFromMice(mice(isControl));
end

info.caseIds = cellstr(mice(isCase));
info.controlIds = cellstr(mice(isControl));
info.nCase = numel(info.caseIds);
info.nControl = numel(info.controlIds);
info.caseLabel = prefixDisplayLabel(info.casePrefix);
info.controlLabel = prefixDisplayLabel(info.controlPrefix);
info.caseTex = texEscape(info.casePrefix);
info.controlTex = texEscape(info.controlPrefix);

if isfield(summary, 'scheme')
    info.scheme = char(summary.scheme);
else
    info.scheme = '';
end
if isfield(summary, 'alpha')
    info.alpha = summary.alpha;
else
    info.alpha = 0.05;
end
end

%% ------------------------------------------------------------------
function prefix = inferPrefixFromMice(mouseIds)
if isempty(mouseIds)
    prefix = 'case';
    return;
end
id = char(mouseIds{1});
tok = regexp(id, '^([A-Za-z]+)', 'tokens', 'once');
if isempty(tok)
    prefix = id;
else
    prefix = tok{1};
end
end

%% ------------------------------------------------------------------
function label = prefixDisplayLabel(prefix)
label = strrep(char(prefix), '_', ' ');
if isempty(label)
    label = 'case';
end
end

%% ------------------------------------------------------------------
function s = texEscape(str)
s = strrep(char(str), '_', '\_');
end
