function tf = isAbsolutePath(p)
%ISABSOLUTEPATH  True for an absolute path on either Windows or POSIX.
%
%   isAbsolutePath('C:\studies\mouse')   % true
%   isAbsolutePath('/home/me/study')     % true
%   isAbsolutePath('data/Case_1')        % false
%
%   Accepts both separators on both platforms, so a config file written on
%   one operating system still resolves on the other.

tf = false;
p = char(p);
if isempty(p)
    return;
end
if p(1) == '/' || p(1) == filesep
    tf = true;
elseif numel(p) >= 2 && isletter(p(1)) && p(2) == ':'
    tf = true;
end
end
