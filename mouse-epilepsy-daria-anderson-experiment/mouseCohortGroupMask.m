function [isCase, isControl] = mouseCohortGroupMask(mice, casePrefix, controlPrefix)
%MOUSECOHORTGROUPMASK  Logical masks for case vs control mice by ID prefix.
%
%   [isCase, isControl] = mouseCohortGroupMask(mice, 'Anderson', 'Arnold')
%
% A mouse matching both prefixes is excluded from both masks.

if nargin < 2 || isempty(casePrefix)
    casePrefix = 'Anderson';
end
if nargin < 3 || isempty(controlPrefix)
    controlPrefix = 'Arnold';
end
casePrefix = char(casePrefix);
controlPrefix = char(controlPrefix);

miceStr = string(mice(:));
isCase = startsWith(miceStr, casePrefix);
isControl = startsWith(miceStr, controlPrefix);
both = isCase & isControl;
isCase(both) = false;
isControl(both) = false;
end
