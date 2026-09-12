function colorMap = cohortSubjectColors(subjectIds, P)
%COHORTSUBJECTCOLORS  Stable colour per subject ID, for multi-subject figures.
%
%   colorMap = cohortSubjectColors({'Case_1', 'Case_2'})
%   c = colorMap('Case_1');
%
%   Assigns each subject a colour from the shared palette in ID order and
%   returns a containers.Map from ID to RGB. Because the assignment is
%   driven by the sorted ID list rather than by plotting order, a subject
%   keeps the same colour across every figure in a run -- which is what
%   lets a reader carry a colour from one panel to the next.
%
%   The palette wraps if there are more subjects than colours.
%
%   See also COHORTFIGSTYLE, PLOTCOHORTOUTLIERSUMMARY.

if nargin < 2 || isempty(P)
    P = cohortFigStyle();
end

palette = [P.a1; P.red; P.teal; P.navy; P.amber; P.grey];

colorMap = containers.Map('KeyType', 'char', 'ValueType', 'any');
if isempty(subjectIds)
    return;
end

subjectIds = cellstr(subjectIds);
for k = 1:numel(subjectIds)
    colorMap(subjectIds{k}) = palette(mod(k - 1, size(palette, 1)) + 1, :);
end
end
