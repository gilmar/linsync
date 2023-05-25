function plotErrorInEmpiricalSyncResults(varargin)
%% function plotErrorInEmpiricalSyncResults(N, b, c, undirected, discretized, repeats, networkType, p, d, S, maxMotifLength, folder, MaxK, dt)
%function plotErrorInEmpiricalSyncResults(parameters)
%
% Plot the error between analytic and empirical synchronisation versus S and p
%  for a number of samples networks
%  (discrete time AR or continuous time Ornstein-Uhlenbeck specified in input).
% without running them, by retrieving results from files saved previously by computeSyncResults
%
% Inputs:
% See the specs of inputs in the parametersTemplate.m script.
% Parameters can be supplied in one of two ways:
% - Option 1 -- (a filename string or object)
%   - parameters - an object containing the expected properties (as outlined below), or a string
%     describing the filename to run load this object in
% - Option 2 -- DEPRECATED -- (all supplied individually)
%   - definitions of individual variables are as specified in the parametersTemplate.m script.
%     We require the following to be defined in order:
%     N, b, c, undirected, discretized, repeats, networkType, p, d, S, maxMotifLength, folder
%     The following are optional:
%     MaxK, dt
%  
% BUT NOTE: Parameter S (number of time samples for empirical results) if supplied 
%  via Option 2 should be an array for this call (else can use parameters
%  script to specify S and SRangeToPlot separately)
%
%% Linear Sync Toolkit (linsync)
% Copyright (C) 2023 Joseph T. Lizier
% Distributed under GNU General Public License v3

parseParameters;

tic;

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

% Pull out the array of parameters that we will sweep through here:
paramsToRunThrough = parameters.(parameters.tosweep);
% But replace the swept parameter with its final value (for generating the
% results file names)
parameters.(parameters.tosweep) = parameters.(parameters.tosweep)(end);

% Load the processed results here:
[networkType, netTypeSuffix] = generateNetworkTypeStrings(parameters);

% Now plot the errors of the empirical values as a function of S and the swept parameter:

avAbsoluteEmpiricalErrors = zeros(length(SRangeToPlot), length(paramsToRunThrough));
avRelativeEmpiricalErrors = zeros(length(SRangeToPlot), length(paramsToRunThrough));
rmsEmpiricalErrors = zeros(length(SRangeToPlot), length(paramsToRunThrough));
rmsLogEmpiricalErrors = zeros(length(SRangeToPlot), length(paramsToRunThrough));
stdAbsoluteEmpiricalErrors = zeros(length(SRangeToPlot), length(paramsToRunThrough));
stdRelativeEmpiricalErrors = zeros(length(SRangeToPlot), length(paramsToRunThrough));

for sIndex = 1:length(SRangeToPlot)
    s = SRangeToPlot(sIndex);

    % This will generate filename prefix with the *final* value of the swept
    %  parameter, but parameter being swept will be indicated
    fileNamePrefix = sprintf('%s/N%d-%s%s-b%.2f-c%.2f-p%.4f-sweep_%s-%s-k%d-%s-S%d-repeats%d', ...
                folder, parameters.N, networkType, netTypeSuffix, parameters.b, parameters.c, parameters.p, ...
                parameters.tosweep, undirString, maxMotifLength, discString, s, repeats);

    load([fileNamePrefix, '.mat'], '-mat', 'repeats', ...
        'syncWidths', 'syncWidthEmpirical');
    % Now grab the relative error against s for each value of the parameter being
    % swept
    avRelativeEmpiricalErrors(sIndex,:) = mean(abs(syncWidthEmpirical - syncWidths) ./ syncWidths, 2)';
    stdRelativeEmpiricalErrors(sIndex,:) = std(abs(syncWidthEmpirical - syncWidths) ./ syncWidths, 0, 2)';
    avAbsoluteEmpiricalErrors(sIndex,:) = mean(abs(syncWidthEmpirical - syncWidths), 2)';
    stdAbsoluteEmpiricalErrors(sIndex,:) = std(abs(syncWidthEmpirical - syncWidths), 0, 2)';
    rmsEmpiricalErrors(sIndex,:) = sqrt(mean((syncWidthEmpirical - syncWidths).^2, 2))';
    rmsLogEmpiricalErrors(sIndex,:) = sqrt(10.^mean(log10((syncWidthEmpirical - syncWidths).^2), 2))';
    fprintf('Gathered results from %d repeats for S=%d\n', repeats, s);
end

toc

% Save the plots like so:
% print('-depsc', [fileNamePrefix, '.eps'])
% saveas(gca, [fileNamePrefix, '.fig'], 'fig');

% Now plot Relative errors for each parameter and S: (this is the one for the
% paper)
figure(2)
hold off;
h = zeros(length(paramsToRunThrough), 1);
labels = {};
for pIndex = 1:length(paramsToRunThrough)
    % h(pIndex) = loglog(SRangeToPlot, avRelativeEmpiricalErrors(:,pIndex)', '-x', 'markersize', 10);
    % Correct symmetric error bars as described at ...
    % https://faculty.washington.edu/stuve/log_error.pdf -- ...
    logErrorBars = 1./log(10).*(stdRelativeEmpiricalErrors(:,pIndex)')./(avRelativeEmpiricalErrors(:,pIndex)');
    h(pIndex) = errorbar(SRangeToPlot, avRelativeEmpiricalErrors(:,pIndex)', ...
        avRelativeEmpiricalErrors(:,pIndex)'.*(1-10.^(-logErrorBars)), ...
        avRelativeEmpiricalErrors(:,pIndex)'.*(10.^(logErrorBars) - 1), ...
        'x', 'markersize', 10);
    labels{pIndex} = sprintf('%s=%.3f', parameters.tosweep_label, paramsToRunThrough(pIndex));
    hold on;
end
palette = jet (length(paramsToRunThrough));
% palette = hsv (length(paramsToRunThrough)); % Good for smaller numbers of parameters on the plot
for i =1:length(paramsToRunThrough)
    set(h(i),'color',palette(i,:))
end
% h
set(gca, 'XScale','log', 'YScale','log');
% labels
legend(h, labels, 'NumColumns', 2);
hold off;
title('Relative error of empirical \langle\sigma^2\rangle estimate versus sample length');
xlabel('Sample length L');
ylabel('Relative error');
% Put a bit of space on the side of the axes
a = axis;
a(1:2) = [a(1)/1.2, a(2)*1.2];
axis(a)

avRelativeEmpiricalErrors % debug print this
stdRelativeEmpiricalErrors % debug print this

% Now plot absolute errors for each parameter and S:
figure(3)
hold off;
h = zeros(length(paramsToRunThrough), 1);
labels = {};
for pIndex = 1:length(paramsToRunThrough)
    % h(pIndex) = loglog(SRangeToPlot, avAbsoluteEmpiricalErrors(:,pIndex)', '-x', 'markersize', 10);
    % h(pIndex) = errorbar(SRangeToPlot, avAbsoluteEmpiricalErrors(:,pIndex)', stdAbsoluteEmpiricalErrors(:,pIndex)', ...
    %    '-x', 'markersize', 10);
    logErrorBars = 1./log(10).*(stdAbsoluteEmpiricalErrors(:,pIndex)')./(avAbsoluteEmpiricalErrors(:,pIndex)');
    h(pIndex) = errorbar(SRangeToPlot, avAbsoluteEmpiricalErrors(:,pIndex)', ...
        avAbsoluteEmpiricalErrors(:,pIndex)'.*(1-10.^(-logErrorBars)), ...
        avAbsoluteEmpiricalErrors(:,pIndex)'.*(10.^(logErrorBars) - 1), ...
        '-x', 'markersize', 10);
    labels{pIndex} = sprintf('%s=%.3f', parameters.tosweep_label, paramsToRunThrough(pIndex));
    hold on;
end
palette = jet (length(paramsToRunThrough));
for i =1:length(paramsToRunThrough)
    set(h(i),'color',palette(i,:))
end
% h;
set(gca, 'XScale','log', 'YScale','log');
% labels;
legend(h, labels);
hold off;
title('Absolute error of empirical \langle\sigma^2\rangle estimate versus sample length');
xlabel('Sample length L');
ylabel('Absolute error');

avAbsoluteEmpiricalErrors

% Now plot RMS errors for each parameter and S:
figure(4)
hold off;
h = zeros(length(paramsToRunThrough), 1);
labels = {};
for pIndex = 1:length(paramsToRunThrough)
    h(pIndex) = loglog(SRangeToPlot, rmsEmpiricalErrors(:,pIndex)', '-x', 'markersize', 10);
    labels{pIndex} = sprintf('%s=%.3f', parameters.tosweep_label, paramsToRunThrough(pIndex));
    hold on;
end
palette = jet (length(paramsToRunThrough));
for i =1:length(paramsToRunThrough)
    set(h(i),'color',palette(i,:))
end
% h;
% labels;
legend(h, labels);
hold off;
title('RMS error of empirical \langle\sigma^2\rangle estimate versus sample length');
xlabel('Sample length L');
ylabel('RMS error');

rmsEmpiricalErrors

end

