function [folderName, stamp] = mouseExperimentFolderName(configName)
%MOUSEEXPERIMENTFOLDERNAME  Results subfolder name = config name + run timestamp.
%
%   [folderName, stamp] = mouseExperimentFolderName('initial_column')
%   -> folderName = 'initial_column_2026-05-23_1430'
%
%   stamp format: yyyy-mm-dd_HHMM (date, hour, minute; no seconds).

configName = char(configName);
stamp = datestr(now, 'yyyy-mm-dd_HHMM');
folderName = sprintf('%s_%s', configName, stamp);
end
