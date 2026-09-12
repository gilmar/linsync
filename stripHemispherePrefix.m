function [baseName, hemisphere] = stripHemispherePrefix(regionLabel, leftPrefix, rightPrefix)
%STRIPHEMISPHEREPREFIX  Split a node label into hemisphere tag and base region.
%
%   [baseName, hemisphere] = stripHemispherePrefix('L-CORTEX_GUSTATORY')
%   -> 'CORTEX_GUSTATORY', 'L'
%
%   [baseName, hemisphere] = stripHemispherePrefix('lh.precentral', 'lh.', 'rh.')
%   -> 'precentral', 'L'
%
%   Defaults match the common 'L-'/'R-' atlas convention. A label matching
%   neither prefix is returned unchanged with an empty hemisphere tag.
%
%   See also PAIRHEMISPHERENODES.

if nargin < 2 || isempty(leftPrefix)
    leftPrefix = 'L-';
end
if nargin < 3 || isempty(rightPrefix)
    rightPrefix = 'R-';
end

lab = char(regionLabel);
leftPrefix = char(leftPrefix);
rightPrefix = char(rightPrefix);

if startsWith(lab, leftPrefix)
    hemisphere = 'L';
    baseName = lab(numel(leftPrefix) + 1:end);
elseif startsWith(lab, rightPrefix)
    hemisphere = 'R';
    baseName = lab(numel(rightPrefix) + 1:end);
else
    hemisphere = '';
    baseName = lab;
end
end
