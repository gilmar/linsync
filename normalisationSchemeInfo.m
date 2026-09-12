function info = normalisationSchemeInfo(scheme)
%NORMALISATIONSCHEMEINFO  Canonical name and properties of a normalisation scheme.
%
%   info = normalisationSchemeInfo('col')
%
%   Resolves the aliases accepted across the toolkit ('col'/'colnorm' ->
%   'column') and reports what the scheme implies for the rest of the
%   pipeline, so callers do not have to re-derive it with string tests.
%
%   Fields
%     scheme        : canonical name -- 'tvb' | 'none' | 'parkes' | 'column'
%     isLinear      : true when the normalised K is used directly as the
%                     coupling matrix C (no Epileptor fixed point, and
%                     therefore no critical-excitability sweep).
%     usesEpileptor : ~isLinear; the scheme solves a 1-D Epileptor fixed
%                     point and linearises around it.
%     isSymmetric   : true when the scheme preserves symmetry of a
%                     symmetric input, i.e. D(k->) == D(->i) by construction.
%     description   : one-line human-readable summary.
%
%   See also APPLYCONNECTOMENORMALISATION.

if nargin < 1 || isempty(scheme)
    error('normalisationSchemeInfo:NoScheme', 'scheme is required.');
end
raw = lower(strtrim(char(scheme)));

switch raw
    case {'tvb'}
        info = mk('tvb', false, true, ...
            'TVB 95th-percentile truncation and rescale to [0,1]; Epileptor fixed point');
    case {'none', 'raw'}
        info = mk('none', false, true, ...
            'raw weights, no rescaling; Epileptor fixed point');
    case {'parkes'}
        info = mk('parkes', true, true, ...
            'A / (|lambda|_max + c) after Parkes et al. 2024; linear, symmetric');
    case {'column', 'col', 'colnorm'}
        info = mk('column', true, false, ...
            'column sums rescaled to a fixed target; linear, asymmetric');
    otherwise
        error('normalisationSchemeInfo:UnknownScheme', ...
            'Unknown normalisation scheme "%s". Use tvb, none, parkes, or column.', raw);
end
end

%% ------------------------------------------------------------------
function info = mk(scheme, isLinear, isSymmetric, description)
info = struct( ...
    'scheme',        scheme, ...
    'isLinear',      isLinear, ...
    'usesEpileptor', ~isLinear, ...
    'isSymmetric',   isSymmetric, ...
    'description',   description);
end
