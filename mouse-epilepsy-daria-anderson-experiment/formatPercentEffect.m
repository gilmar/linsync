function s = formatPercentEffect(eff)
%FORMATPERCENTEFFECT  Signed percent string for bar annotations.
if isnan(eff)
    s = 'NaN';
elseif eff >= 0
    s = sprintf('+%.1f%%', eff);
else
    s = sprintf('%.1f%%', eff);
end
end
