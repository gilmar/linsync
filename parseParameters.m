%% Parses varargin into the relevant parameters
% either based on a structure, filename, or
% each parameter supplied individually.
%
%% Linear Sync Toolkit (linsync)
% Copyright (C) 2023 Joseph T. Lizier
% Distributed under GNU General Public License v3

if (length(varargin) == 1)
    % User has provided a parameters object directly or a string specifying
    % the filename to load a parameters object in.
    parameters = varargin{1};
    if ischar(parameters)
        % Assume that this string contains a filename which when run will load
        % a properties object for this run
        eval(['run ', parameters]);
    end
    % Postcondition: parameters are in the parameters object
    
    % Assign all of the relevant variables:
    N = parameters.N;
    b = parameters.b;
    c = parameters.c;
    undirected = parameters.undirected;
    discretized = parameters.discretized;
    repeats = parameters.repeats;
    generateNetworkFunction = parameters.generateNetworkFunction;
    p = parameters.p;
    d = parameters.d;
    S = parameters.S;
    SRangeToPlot = parameters.SRangeToPlot;
    if (~isfield(parameters, 'motifLengthsToCheck'))
        parameters.motifLengthsToCheck = 1:parameters.maxMotifLength;
    end
    motifLengthsToCheck = parameters.motifLengthsToCheck;
    parameters.maxMotifLength = max(parameters.motifLengthsToCheck);
    maxMotifLength = parameters.maxMotifLength;
    folder = parameters.folder;
    MaxK = parameters.MaxK;
    dt = parameters.dt;
    randSeed = parameters.randSeed;
    if (~isfield(parameters, 'checkDiagonalizable'))
        parameters.checkDiagonalizable = true;
    end
    checkDiagonalizable = parameters.checkDiagonalizable;
    originalParameters = parameters; % Store for later
elseif (length(varargin) < 12)
    fprintf('Not enough arguments supplied, see code for details');
else
    % All parameters have been supplied individually
    %  We are *deprecating* this method so do not rely on it!
    N = varargin{1};
    b = varargin{2};
    c = varargin{3};
    undirected = varargin{4};
    discretized = varargin{5};
    repeats = varargin{6};
    generateNetworkFunction = varargin{7};
    p = varargin{8};
    d = varargin{9};
    S = varargin{10};
    SRangeToPlot = S; % only used by plot scripts in position of S
    maxMotifLength = varargin{11}; % This parameter is deprecated, see above.
    folder = varargin{12};
    if (length(varargin) > 12)
        MaxK = varargin{13};
    else
        MaxK = 2000;
    end
    if (length(varargin) > 13)
        dt = varargin{14};
    else
        dt = 1;
    end
    if (length(varargin) > 14)
        randSeed = varargin{15};
    else
        randSeed = 'shuffle';
    end
    if (length(varargin) > 15)
        checkDiagonalizable = varargin{16};
    else
        checkDiagonalizable = true;
    end
end

% Now pull out some useful strings (for filenames and debug prints)
% from these properties:
% Generate string for the boolean arguments ready for file names
if (undirected)
    undirString = 'un';
else
    undirString = 'dir';
end
if (discretized)
    discString = 'disc';
else
    discString = 'cont';
end
