% Compute <sigma^2> analytically and from numerical simulations
% for discrete-time VAR process through a small-world
% transition, with standard interaction delay.
%
%
%% Linear Sync Toolkit (linsync)
% Copyright (C) 2023 Joseph T. Lizier
% Distributed under GNU General Public License v3

%% Preliminaries:

% Load in the parameters object:
exp1Parameters;

% Add the toolkit to the path:
addpath(parameters.syncToolkitPath);

% Make sure the folder for the results to be stored in exists:
if (exist(parameters.folder, 'dir') == 0)
    % Need to create the results folder:
    [success, message] = mkdir(parameters.folder);
    if (success == 0)
        error('Creating the results folder failed: %s\n', message);
    end
end

% To speed up the experiment when not on a cluster you can set:
parameters.repeats = 10; % instead of 2000
parameters.S = [100, 1000, 10000, 100000]; % Removing the longest runs
parameters.SRangeToPlot = parameters.S;
% You can get a decent idea of the main trends for this number already,
% though there are significant fluctuations / std error differences
% compared to the longer runs.

%% Experiment 1 and plots:
% Here, the numerical simulations are the limiting factor
for S = parameters.SRangeToPlot
    parameters.S = S;
    computeSyncResults(parameters);
end
% Takes about 3 mins on my machine for 10 repeats for S = [100,1000,10000,100000]
% will take much longer for 2000 repeats and S up to 1000000 -- 
% Better to do on cluster for that many.

%% Now plot the results:
plotErrorInEmpiricalSyncResults(parameters);
figure(2); % Select the correct plot we're keeping
print('-depsc', [parameters.folder, '/figS1.eps'])
saveas(gca, [parameters.folder, '/figS1.fig'], 'fig');

fprintf('Figure S1 finished, using %d network samples (is this the same as the 2000 used for the paper?)\n', ...
    parameters.repeats);

