function [folderName, stamp] = experimentFolderName(configName, when)
%EXPERIMENTFOLDERNAME  Results subfolder name: config name plus a run timestamp.
%
%   [folderName, stamp] = experimentFolderName('baseline_column')
%   -> 'baseline_column_2026-05-23_1430'
%
%   Stamping the folder means re-running the same config never overwrites
%   an earlier run, so a parameter sweep leaves an auditable trail instead
%   of a single mutable output directory.
%
%   Stamp format: yyyy-mm-dd_HHMM (minute resolution).
%
%   See also RUNCONNECTOMEEXPERIMENT, EXPERIMENTPROVENANCEFILE.

if nargin < 2 || isempty(when)
    when = now; %#ok<TNOW1>
end
stamp = datestr(when, 'yyyy-mm-dd_HHMM'); %#ok<DATST>
folderName = sprintf('%s_%s', char(configName), stamp);
end
