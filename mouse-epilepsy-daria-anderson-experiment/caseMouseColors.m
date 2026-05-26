function colorMap = caseMouseColors(caseIds, P)
%CASEMOUSECOLORS  Stable face colours for each case-mouse ID (struct map).
%
%   colorMap = caseMouseColors({'Anderson_1','Anderson_2'}, cohortFigStyle())

if nargin < 2 || isempty(P)
    P = cohortFigStyle();
end
palette = [P.a1; P.red; P.teal; P.navy; P.amber; P.grey];
if isempty(caseIds)
    colorMap = containers.Map('KeyType', 'char', 'ValueType', 'any');
    return;
end
caseIds = cellstr(caseIds);
colorMap = containers.Map('KeyType', 'char', 'ValueType', 'any');
for k = 1:numel(caseIds)
    colorMap(caseIds{k}) = palette(mod(k - 1, size(palette, 1)) + 1, :);
end
end
