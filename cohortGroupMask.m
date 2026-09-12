function [isCase, isControl] = cohortGroupMask(subjectIds, casePrefix, controlPrefix)
%COHORTGROUPMASK  Logical case/control masks for subject IDs, by ID prefix.
%
%   [isCase, isControl] = cohortGroupMask(subjectIds, 'Case', 'Ctrl')
%
%   The toolkit identifies cohort membership from the subject ID itself so
%   that a study needs no separate group table: name the folders
%   <casePrefix>_1, <casePrefix>_2, <controlPrefix>_1, ... and the whole
%   pipeline picks the grouping up from the data directory.
%
%   A subject matching both prefixes (possible when one prefix is a prefix
%   of the other) is excluded from both masks rather than silently assigned.
%
%   See also COHORTGROUPINFO, LISTCONNECTOMESUBJECTS.

if nargin < 2 || isempty(casePrefix)
    error('cohortGroupMask:NoCasePrefix', 'casePrefix is required.');
end
if nargin < 3 || isempty(controlPrefix)
    error('cohortGroupMask:NoControlPrefix', 'controlPrefix is required.');
end

ids = string(subjectIds(:));
isCase    = startsWith(ids, char(casePrefix));
isControl = startsWith(ids, char(controlPrefix));

both = isCase & isControl;
isCase(both) = false;
isControl(both) = false;
end
