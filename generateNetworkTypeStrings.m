function [netType, netTypeSuffix] = generateNetworkTypeStrings(parameters)
%%
% Generate strings for filenames etc for a given network type and
% parameters.
%
% Inputs:
% - parameters - an object containing the expected properties (as outlined for option 2)
%
% Outputs
% - netType - a string describing the network type
% - netTypeSuffix - a suffix capturing relevant parameters for the network
% type
%
%% Linear Sync Toolkit (linsync)
% Copyright (C) 2023 Joseph T. Lizier
% Distributed under GNU General Public License v3

if (strcmp(parameters.generateNetworkFunction, 'generateNewRandomMatrix'))
    netType = 'rand';
    netTypeSuffix = '';
elseif (strcmp(parameters.generateNetworkFunction, 'generateNewRandomFixedDMatrix'))
    netType = 'randFixedD';
    netTypeSuffix = '';
elseif (strcmp(parameters.generateNetworkFunction, 'generateNewRandomRingMatrix'))
    netType = 'randRing';
    netTypeSuffix = sprintf('-d%d', parameters.d);
else
    error('networkType \"%s\" not recognised\n', networkType);
end

end