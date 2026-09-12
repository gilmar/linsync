function info = setupConnectomePaths(varargin)
%SETUPCONNECTOMEPATHS  Put the linsync toolkit (and optionally BCT) on the MATLAB path.
%
%   info = setupConnectomePaths()
%   info = setupConnectomePaths('ProjectRoot', '/path/to/my-study')
%   info = setupConnectomePaths('BctRoot', '/path/to/2019_03_03_BCT')
%
%   Call once at the start of a session or driver script. Adds the folder
%   containing this file (the linsync toolkit root), so every helper the
%   connectome pipeline needs resolves regardless of the current folder.
%
%   The Brain Connectivity Toolbox is optional: it supplies betweenness_wei
%   and distance_wei for computeNetworkCentralities. If it is not found the
%   stability centralities still run and the BCT-only measures come back as
%   NaN. See docs/BCT.md for install instructions.
%
%   Name-value options
%     'ProjectRoot' : study folder holding configs/, data/ and results/.
%                     Default: the current folder. Added to the path too,
%                     so study-specific drivers stay callable.
%     'BctRoot'     : explicit BCT location. Default: look for a folder
%                     whose name contains 'BCT' beside the toolkit root and
%                     beside the project root.
%     'Quiet'       : true to suppress the missing-BCT warning.
%
%   Output struct
%     toolkitRoot, projectRoot, dataRoot, resultsRoot, configsRoot, bctRoot
%     hasBct : whether BCT functions are now resolvable
%
%   See also RUNCONNECTOMEEXPERIMENT, COMPUTENETWORKCENTRALITIES.

p = inputParser;
addParameter(p, 'ProjectRoot', '', @(s) ischar(s) || isstring(s));
addParameter(p, 'BctRoot',     '', @(s) ischar(s) || isstring(s));
addParameter(p, 'Quiet',       false, @islogical);
parse(p, varargin{:});

toolkitRoot = fileparts(mfilename('fullpath'));
addpath(toolkitRoot);

projectRoot = char(p.Results.ProjectRoot);
if isempty(projectRoot)
    projectRoot = pwd;
elseif ~isAbsolutePath(projectRoot)
    projectRoot = fullfile(pwd, projectRoot);
end
if ~strcmp(projectRoot, toolkitRoot) && isfolder(projectRoot)
    addpath(projectRoot);
end

bctRoot = char(p.Results.BctRoot);
if isempty(bctRoot)
    bctRoot = findBct({toolkitRoot, projectRoot});
end
if ~isempty(bctRoot) && isfolder(bctRoot)
    addpath(bctRoot);
end

hasBct = exist('betweenness_wei', 'file') == 2 && exist('distance_wei', 'file') == 2;
if ~hasBct && ~p.Results.Quiet
    warning('setupConnectomePaths:NoBCT', ...
        ['Brain Connectivity Toolbox not found. Betweenness and closeness ' ...
         'centralities will be NaN; stability centralities are unaffected. ' ...
         'See %s'], fullfile(toolkitRoot, 'docs', 'BCT.md'));
end

info = struct();
info.toolkitRoot = toolkitRoot;
info.projectRoot = projectRoot;
info.dataRoot    = fullfile(projectRoot, 'data');
info.resultsRoot = fullfile(projectRoot, 'results');
info.configsRoot = fullfile(projectRoot, 'configs');
info.bctRoot     = bctRoot;
info.hasBct      = hasBct;
end

%% ------------------------------------------------------------------
function bctRoot = findBct(searchRoots)
%FINDBCT  Look for a *BCT* folder in, or next to, each search root.
bctRoot = '';
candidates = {};
for k = 1:numel(searchRoots)
    if isempty(searchRoots{k})
        continue;
    end
    candidates{end+1} = searchRoots{k};               %#ok<AGROW>
    candidates{end+1} = fileparts(searchRoots{k});    %#ok<AGROW>
end
for k = 1:numel(candidates)
    base = candidates{k};
    if isempty(base) || ~isfolder(base)
        continue;
    end
    d = dir(fullfile(base, '*BCT*'));
    d = d([d.isdir]);
    if ~isempty(d)
        bctRoot = fullfile(base, d(1).name);
        return;
    end
end
end
