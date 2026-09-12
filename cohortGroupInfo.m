function info = cohortGroupInfo(summary)
%COHORTGROUPINFO  Case/control metadata pulled out of a comparison summary.
%
%   info = cohortGroupInfo(summary)
%
%   Accepts the summary structs returned by compareCohortGroups and
%   compareHemisphericAsymmetry, and returns everything a figure or report
%   needs to describe the cohort without hard-coding group names:
%
%     casePrefix, controlPrefix : the ID prefixes in use
%     caseIds, controlIds       : cellstr of subject IDs per group
%     nCase, nControl           : group sizes
%     caseLabel, controlLabel   : display names (underscores -> spaces)
%     caseTex, controlTex       : TeX-safe names (underscores escaped)
%     scheme, alpha             : passed through when present
%
%   Prefixes are taken from the summary when recorded, and otherwise
%   inferred from the leading alphabetic run of the first member ID.
%
%   See also COHORTGROUPMASK, COMPARECOHORTGROUPS.

if ~isstruct(summary) || ~isfield(summary, 'subjects')
    error('cohortGroupInfo:BadSummary', ...
        'summary must be a struct with a "subjects" field.');
end

subjects = summary.subjects(:);
if ~isfield(summary, 'isCase') || ~isfield(summary, 'isControl')
    error('cohortGroupInfo:NoMask', ...
        'summary must include isCase and isControl masks.');
end
isCase = logical(summary.isCase(:));
isControl = logical(summary.isControl(:));

info = struct();
info.casePrefix    = resolvePrefix(summary, 'casePrefix', subjects(isCase), 'case');
info.controlPrefix = resolvePrefix(summary, 'controlPrefix', subjects(isControl), 'control');

info.caseIds    = cellstr(subjects(isCase));
info.controlIds = cellstr(subjects(isControl));
info.nCase      = numel(info.caseIds);
info.nControl   = numel(info.controlIds);

info.caseLabel    = displayLabel(info.casePrefix, 'case');
info.controlLabel = displayLabel(info.controlPrefix, 'control');
info.caseTex      = texEscape(info.casePrefix);
info.controlTex   = texEscape(info.controlPrefix);

info.scheme = fieldOr(summary, 'scheme', '');
info.alpha  = fieldOr(summary, 'alpha', 0.05);
end

%% ------------------------------------------------------------------
function prefix = resolvePrefix(summary, fieldName, memberIds, fallback)
if isfield(summary, fieldName) && ~isempty(summary.(fieldName))
    prefix = char(summary.(fieldName));
    return;
end
if isempty(memberIds)
    prefix = fallback;
    return;
end
id = char(memberIds{1});
tok = regexp(id, '^([A-Za-z]+)', 'tokens', 'once');
if isempty(tok)
    prefix = id;
else
    prefix = tok{1};
end
end

%% ------------------------------------------------------------------
function label = displayLabel(prefix, fallback)
label = strrep(char(prefix), '_', ' ');
if isempty(label)
    label = fallback;
end
end

%% ------------------------------------------------------------------
function s = texEscape(str)
s = strrep(char(str), '_', '\_');
end

%% ------------------------------------------------------------------
function v = fieldOr(s, name, defaultValue)
if isfield(s, name)
    v = s.(name);
else
    v = defaultValue;
end
end
