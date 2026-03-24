function plotStabilityResults(varargin)
%% function plotStabilityResults(parameters)
% function plotStabilityResults(N, b, c, undirected, discretized, repeats, networkType, p, d, S, maxMotifLength, folder, MaxK, dt)
%
% Plot deviation from stability (D_st) versus one swept parameter for a number of network samples
%  (discrete time AR or continuous time Ornstein-Uhlenbeck specified in input),
%  by retrieving results from files saved previously by computeSyncResults
%  with parameters.computationMode = 'stability'.
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
%  One parameter should be an array, which determines which we plot against.
%
%% Linear Sync Toolkit (linsync)
% Copyright (C) 2023 Joseph T. Lizier
% Distributed under GNU General Public License v3

parseParameters;

tic;

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
    load([fileNamePrefix, '.mat'], '-mat', 'N', 'd', 'b', 'c', 'p', 'undirected', ...
        'motifLengthsToCheck', 'maxMotifLength', 'discretized', ...
        'networkType', 'paramsToRunThrough', 'S', 'repeats', ...
        'stabilityWidths', 'dominantEigenvalues', 'secondEigenvalues');
catch ME
    error(['Error loading file ', [fileNamePrefix, '.mat'], ' - if plotting a cluster run did you change parameters.repeats to the total number of repeats rather than for each cluster job?']);
end
toc

avStabilityWidths = mean(stabilityWidths, 2)';
stdStabilityWidths = std(stabilityWidths, 0, 2)';

% dominantEigenvalues and secondEigenvalues are already stored as the
% criterion-relevant scalar at save time: abs() for discrete, real() for
% continuous. No need to re-apply those transforms here.
avDominantEigenvalues = mean(dominantEigenvalues, 2)';
stdDominantEigenvalues = std(dominantEigenvalues, 0, 2)';
avSecondEigenvalues = mean(secondEigenvalues, 2)';
stdSecondEigenvalues = std(secondEigenvalues, 0, 2)';
if (discretized)
    labelDominantEigenvalue = '|\lambda_1|';
    labelSecondEigenvalue = '|\lambda_2|';
    labelEigenvalueAxis = '|\lambda|';
else
    labelDominantEigenvalue = 'Re(\lambda_1)';
    labelSecondEigenvalue = 'Re(\lambda_2)';
    labelEigenvalueAxis = 'Re(\lambda)';
end

figure(1);
clf;
hold off;

% Left: D_st (linear y). Log-Y on the left combined with yyaxis/errorbar can map
%  right-axis series incorrectly so all traces appear to share the left scale.
yyaxis left;
hold on;
hDst = errorbar(paramsToRunThrough, avStabilityWidths, stdStabilityWidths, ...
    'rx', 'markersize', 10);
set(gca, 'YColor', [0 0 0]); % otherwise it is blue
ylabel('D_{st}', 'Interpreter', 'tex');
if (strcmp('p', parameters.tosweep))
    set(gca, 'XScale', 'log');
else
    set(gca, 'XScale', 'linear');
end
xlabel(parameters.tosweep_label);

maxAvDst = max(avStabilityWidths);
minAvDst = min(avStabilityWidths);
if (maxAvDst > 0)
    dstPad = 0.05 * (maxAvDst - minAvDst);
    if (~(dstPad > 0))
        dstPad = 0.05 * maxAvDst;
    end
    yyaxis left;
    ylim([max(0, minAvDst - dstPad), maxAvDst + dstPad]);
end
if (strcmp('p', parameters.tosweep))
    xlim([paramsToRunThrough(1) / 1.2, paramsToRunThrough(end) * 1.2]);
else
    xlim([paramsToRunThrough(1) - (paramsToRunThrough(2) - paramsToRunThrough(1)), ...
        paramsToRunThrough(end) + (paramsToRunThrough(end) - paramsToRunThrough(end - 1))]);
end

hold off;

% Right: Re(lambda); keep linear and independent y-limits from D_st
yyaxis right;
set(gca, 'YScale', 'linear');
hold on;
hL1 = errorbar(paramsToRunThrough, avDominantEigenvalues, stdDominantEigenvalues, ...
    'ks', 'markersize', 10);
hL2 = errorbar(paramsToRunThrough, avSecondEigenvalues, stdSecondEigenvalues, ...
    'ms', 'markersize', 10);
eLo = min([avDominantEigenvalues - stdDominantEigenvalues, avSecondEigenvalues - stdSecondEigenvalues]);
eHi = max([avDominantEigenvalues + stdDominantEigenvalues, avSecondEigenvalues + stdSecondEigenvalues]);
pad = 0.05 * (eHi - eLo);
if (~(pad > 0))
    pad = 0.1;
end
if (discretized)
    ylim([eLo - pad, eHi + pad]);
else
    % Match paper Figure 2: Re(lambda) axis extends up to 1 (continuous-time).
    ylim([eLo - pad, 1]);
end
ylabel(labelEigenvalueAxis, 'Interpreter', 'tex');
set(gca, 'YColor', [0 0 0]); % otherwise it is red

hold off;

legend([hDst; hL1; hL2], {'D_{st}', labelDominantEigenvalue, labelSecondEigenvalue}, ...
    'Interpreter', 'tex', 'Location', 'southwest');

title('Deviation from Stability vs. p (Watts-Strogatz)');

fprintf('Mean stabilityWidths(%.4f)=%.4f\n', paramsToRunThrough(1), avStabilityWidths(1));
fprintf('Mean dominantEigenvalue(%.4f)=%.4f\n', paramsToRunThrough(1), avDominantEigenvalues(1));

end
