function plotSyncResults(varargin)
%% function plotSyncResults(N, b, c, undirected, discretized, repeats, networkType, p, d, S, maxMotifLength, folder, MaxK, dt)
% function plotSyncResults(parameters)
%
% Plot the synchronisation versus one parameter for a number of network samples
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
%  One parameter should be an array, which determines which we plot
%  against.
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
    deltaT = 1; % The natural time step
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
% This will generate filename prefix with the *final* value of the swept
%  parameter, but parameter being swept will be indicated
fileNamePrefix = sprintf('%s/N%d-%s%s-b%.2f-c%.2f-p%.4f-sweep_%s-%s-k%d-%s-S%d-repeats%d', ...
            folder, parameters.N, networkType, netTypeSuffix, parameters.b, parameters.c, parameters.p, ...
            parameters.tosweep, undirString, maxMotifLength, discString, S, repeats);
try
    load([fileNamePrefix, '.mat'], '-mat', 'N', 'd', 'b', 'c', 'p', 'undirected', 'maxMotifLength', 'discretized', ...
        'networkType', 'paramsToRunThrough', 'S', 'repeats', ...
        'syncWidths', 'syncWidthApproxes', 'syncWidthEmpirical', ...
        'dominantEigenvalues', 'secondEigenvalues');
catch ME
    error(['Error loading file ', [fileNamePrefix, '.mat'], ' - if plotting a cluster run did you change parameters.repeats to the total number of repeats rather than for each cluster job?']);
end
toc

avSyncWidths = mean(syncWidths, 2)';
stdSyncWidths = std(syncWidths, 0, 2)';
avSyncWidthsEmpirical = mean(syncWidthEmpirical, 2)';
stdSyncWidthsEmpirical = std(syncWidthEmpirical, 0, 2)';
if (discretized)
    % We want to plot the magnitudes
    avDominantEigenvalues = mean(abs(dominantEigenvalues), 2)';
    stdDominantEigenvalues = std(abs(dominantEigenvalues), 0, 2)';
    avSecondEigenvalues = mean(abs(secondEigenvalues), 2)';
    stdSecondEigenvalues = std(abs(secondEigenvalues), 0, 2)';
    labelDominantEigenvalue = '|\lambda_1|';
    labelSecondEigenvalue = '|\lambda_2|';
    labelEigenvalueAxis = '|\lambda|';
else
    % We want to plot the real components
    avDominantEigenvalues = mean(real(dominantEigenvalues), 2)';
    stdDominantEigenvalues = std(real(dominantEigenvalues), 0, 2)';
    avSecondEigenvalues = mean(real(secondEigenvalues), 2)';
    stdSecondEigenvalues = std(real(secondEigenvalues), 0, 2)';
    labelDominantEigenvalue = 'Re(\lambda_1)';
    labelSecondEigenvalue = 'Re(\lambda_2)';
    labelEigenvalueAxis = 'Re(\lambda)';
end
avSyncWidthApproxes = mean(syncWidthApproxes, 3)';
stdSyncWidthApproxes = std(syncWidthApproxes, 0, 3)';

%     if (varyingP)
%         syncWidthsNormaliser = avSyncWidths(1);
%         eigenvaluesNormaliser = avDominantEigenvalues(1);
%     else
%         % We're not normalising when plotting against c
%         syncWidthsNormaliser = 1;
%         eigenvaluesNormaliser = 1;
%     end
% We're not normalising anymore
syncWidthsNormaliser = 1;
eigenvaluesNormaliser = 1;

%  and make the plots
figure(1);
clf;
hold off;
% We're using separate axes for sigma^2 and eigenvalues
yyaxis left;
%h1 = semilogx(paramsToRunThrough, avSyncWidths ./ syncWidthsNormaliser, 'rx', 'markersize', 10);
h1 = errorbar(paramsToRunThrough, avSyncWidths ./ syncWidthsNormaliser, ...
    stdSyncWidths ./ syncWidthsNormaliser, ...
    'rx', 'markersize', 10);
hold on;
yyaxis right;
% h6 = semilogx(paramsToRunThrough, avDominantEigenvalues ./ eigenvaluesNormaliser, 'ks', 'markersize', 10);
h6 = errorbar(paramsToRunThrough, avDominantEigenvalues ./ eigenvaluesNormaliser, ...
    stdDominantEigenvalues ./ eigenvaluesNormaliser, ...
    'ks', 'markersize', 10);
% h7 = semilogx(paramsToRunThrough, avSecondEigenvalues ./ eigenvaluesNormaliser, 'ys', 'markersize', 10); % Plot relative to dominant
h7 = errorbar(paramsToRunThrough, avSecondEigenvalues ./ eigenvaluesNormaliser, ...
    stdSecondEigenvalues ./ eigenvaluesNormaliser, ...
    'ms', 'markersize', 10);
yyaxis left;
if (discretized)
    % Default: print the cummulative first three terms in the sum
    ordersToPlot = [1,2,3];
else
    % Default: print the cummulative 2nd, 3rd and fourth terms in the sum
    % ordersToPlot = [2,3,4];
    ordersToPlot = [2,10,50];
end
% h2 = semilogx(paramsToRunThrough, avSyncWidthApproxes(ordersToPlot(1),:) ./ syncWidthsNormaliser, 'bo', 'markersize', 10); % Normalise to the total sync width
h2 = errorbar(paramsToRunThrough, avSyncWidthApproxes(ordersToPlot(1),:) ./ syncWidthsNormaliser, ...
    stdSyncWidthApproxes(ordersToPlot(1),:) ./ syncWidthsNormaliser, ...
    'bo', 'markersize', 10);
% h3 = semilogx(paramsToRunThrough, avSyncWidthApproxes(ordersToPlot(2),:) ./ syncWidthsNormaliser, 'go', 'markersize', 10); % Normalise to the total sync width
h3 = errorbar(paramsToRunThrough, avSyncWidthApproxes(ordersToPlot(2),:) ./ syncWidthsNormaliser, ...
    stdSyncWidthApproxes(ordersToPlot(2),:) ./ syncWidthsNormaliser, ...
    'go', 'markersize', 10);
% h4 = semilogx(paramsToRunThrough, avSyncWidthApproxes(ordersToPlot(3),:) ./ syncWidthsNormaliser, 'co', 'markersize', 10); % Normalise to the total sync width
h4 = errorbar(paramsToRunThrough, avSyncWidthApproxes(ordersToPlot(3),:) ./ syncWidthsNormaliser, ...
    stdSyncWidthApproxes(ordersToPlot(3),:) ./ syncWidthsNormaliser, ...
    'co', 'markersize', 10);
% syncWidthApproxLabels = {'\left\langle \sigma_2 \right\rangle', '\left\langle \sigma_3 \right\rangle', '\left\langle \sigma_4 \right\rangle'};
syncWidthApproxLabels = {sprintf('< \\sigma^2_{%d} >', ordersToPlot(3)), ...
    sprintf('< \\sigma^2_{%d} >', ordersToPlot(2)), sprintf('< \\sigma^2_{%d} >', ordersToPlot(1))};
if (S > 0)
    % Add the empirical results:
    % h5 = semilogx(paramsToRunThrough, avSyncWidthsEmpirical ./ syncWidthsNormaliser, 'm+', 'markersize', 10); % Normalise to the total sync width
    h5 = errorbar(paramsToRunThrough, avSyncWidthsEmpirical ./ syncWidthsNormaliser, ...
        stdSyncWidthsEmpirical ./ syncWidthsNormaliser, ...
        'y+', 'markersize', 10);
    % legend({'\left\langle \sigma \right\rangle', syncWidthApproxLabels{1}, syncWidthApproxLabels{2}, syncWidthApproxLabels{3}, ...
    %     '\left\langle \sigma_{emp} \right\rangle', '\lambda_1', '\lambda_2'});
    legend([h1; h4; h3; h2; h5; h6; h7], ...
        {'< \sigma^2 >', ...
        syncWidthApproxLabels{1}, syncWidthApproxLabels{2}, syncWidthApproxLabels{3}, ...
         '< \sigma^2_{emp} >', labelDominantEigenvalue, labelSecondEigenvalue});
else
    legend([h1; h4; h3; h2; h6; h7], ...
        {'< \sigma^2 >', ...
        syncWidthApproxLabels{1}, syncWidthApproxLabels{2}, syncWidthApproxLabels{3}, ...
        labelDominantEigenvalue, labelSecondEigenvalue});
end
hold off;
% Set properties of sigma^2 axis:
yyaxis left
set(gca, 'YColor', [0 0 0]); % otherwise it is blue
ylabel('< \sigma^2 >');
if (strcmp('p', parameters.tosweep))
    % Have a log axis if we're sweeping the p parameter for the network
    set(gca, 'XScale','log');
else
    set(gca, 'XScale','linear');
end
xlabel(parameters.tosweep_label);
a = axis;
% Hard code the sigma^2 limits based on the max deviation from sync
a(3:4) = [0, avSyncWidths(1) ./ syncWidthsNormaliser * 1.2];
maxForSigma = a(4); % Max y value on sigma^2 axis
if (strcmp('p', parameters.tosweep))
    % Put a bit of space on the side of the axes
    a(1:2) = [a(1)/1.2, a(2)*1.2];
else
    % Put a bit of space on the side of the axes
    a(1:2) = [paramsToRunThrough(1) - (paramsToRunThrough(2)-paramsToRunThrough(1)), ...
        paramsToRunThrough(end)+(paramsToRunThrough(end)-paramsToRunThrough(end-1))];
end
axis(a)
% Now adjust eigenvalue axis:
yyaxis right
a = axis;
% Work out the max value on second axis to make
%  the first points on each line up:
maxForEigenvalues = avDominantEigenvalues(1) ./ eigenvaluesNormaliser ...
    .* maxForSigma ./ (avSyncWidths(1) ./ syncWidthsNormaliser);
a(3:4) = [0, maxForEigenvalues];
axis(a)
ylabel(labelEigenvalueAxis);
set(gca, 'YColor', [0 0 0]); % otherwise it is red
% title('Width of sync landscape');

fprintf('Mean syncWidths(%.4f)=%.4f\n', paramsToRunThrough(1), avSyncWidths(1));
fprintf('Mean dominantEigenvalue(%.4f)=%.4f\n', paramsToRunThrough(1), avDominantEigenvalues(1));

% Save the plots like so:
% print('-depsc', [fileNamePrefix, '.eps'])
% saveas(gca, [fileNamePrefix, '.fig'], 'fig');
end

