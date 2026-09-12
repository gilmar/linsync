function pairs = pairHemisphereNodes(labels, varargin)
%PAIRHEMISPHERENODES  Match left-hemisphere nodes to their right-hemisphere twins.
%
%   pairs = pairHemisphereNodes(labels)
%   pairs = pairHemisphereNodes(labels, 'LeftPrefix', 'lh.', 'RightPrefix', 'rh.')
%
%   Homotopic pairing by shared base name: every label carrying the left
%   prefix is matched with the label that carries the right prefix and the
%   same remainder. This is what makes a laterality analysis possible on a
%   bilateral atlas without a hand-maintained pairing table.
%
%   Name-value options
%     'LeftPrefix'  : 'L-'  -- prefix marking the left hemisphere
%     'RightPrefix' : 'R-'  -- prefix marking the right hemisphere
%     'WarnUnpaired': true  -- warn about left nodes with no right twin
%
%   Output struct
%     leftIdx    : nPairs x 1 indices into labels (left member)
%     rightIdx   : nPairs x 1 indices into labels (right member)
%     baseLabels : nPairs x 1 cellstr of the shared base names
%     nPairs     : number of pairs found
%     unpaired   : cellstr of hemisphere-tagged labels with no counterpart
%
%   See also STRIPHEMISPHEREPREFIX, COMPAREHEMISPHERICASYMMETRY.

p = inputParser;
addParameter(p, 'LeftPrefix',   'L-', @(s) ischar(s) || isstring(s));
addParameter(p, 'RightPrefix',  'R-', @(s) ischar(s) || isstring(s));
addParameter(p, 'WarnUnpaired', true, @islogical);
parse(p, varargin{:});
leftPrefix  = string(p.Results.LeftPrefix);
rightPrefix = string(p.Results.RightPrefix);

labelStr = string(labels(:));
n = numel(labelStr);

leftIdx = zeros(0, 1);
rightIdx = zeros(0, 1);
baseLabels = strings(0, 1);
unpaired = strings(0, 1);

for i = 1:n
    lab = labelStr(i);
    if ~startsWith(lab, leftPrefix)
        continue;
    end
    base = extractAfter(lab, leftPrefix);
    j = find(labelStr == rightPrefix + base, 1);
    if isempty(j)
        unpaired(end+1, 1) = lab; %#ok<AGROW>
        if p.Results.WarnUnpaired
            warning('pairHemisphereNodes:NoMatch', ...
                'No "%s%s" counterpart for "%s".', rightPrefix, base, lab);
        end
        continue;
    end
    leftIdx(end+1, 1)    = i;    %#ok<AGROW>
    rightIdx(end+1, 1)   = j;    %#ok<AGROW>
    baseLabels(end+1, 1) = base; %#ok<AGROW>
end

% Right-hemisphere labels that never got claimed.
for i = 1:n
    lab = labelStr(i);
    if startsWith(lab, rightPrefix) && ~ismember(i, rightIdx)
        unpaired(end+1, 1) = lab; %#ok<AGROW>
    end
end

pairs = struct();
pairs.leftIdx    = leftIdx;
pairs.rightIdx   = rightIdx;
pairs.baseLabels = cellstr(baseLabels);
pairs.nPairs     = numel(leftIdx);
pairs.unpaired   = cellstr(unpaired);
end
