function [baseName, hemi] = stripHemispherePrefix(regionLabel)
%STRIPHEMISPHEREPREFIX  Remove leading L-/R- from a node label; return base + hemisphere.
%   [baseName, hemi] = stripHemispherePrefix('L-CORTEX_GUSTATORY')
%   -> 'CORTEX_GUSTATORY', 'L'

lab = char(regionLabel);
hemi = '';
if startsWith(lab, 'L-')
    hemi = 'L';
    baseName = lab(3:end);
elseif startsWith(lab, 'R-')
    hemi = 'R';
    baseName = lab(3:end);
else
    baseName = lab;
end
end
