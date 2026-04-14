function [K, weightsFile] = loadPatientWeights(patient)
%LOADPATIENTWEIGHTS  Load structural weights for a patient ID (cross-platform paths).
%
%   Looks for (in order), relative to the epileptor-experiment directory:
%     data/connectivity_<patient>/weights.txt
%     weights.txt   (fallback in experiment folder)

experimentRoot = fileparts(mfilename('fullpath'));
primary = fullfile(experimentRoot, 'data', ['connectivity_' patient], 'weights.txt');
fallback = fullfile(experimentRoot, 'weights.txt');
if exist(primary, 'file')
    weightsFile = primary;
elseif exist(fallback, 'file')
    weightsFile = fallback;
else
    error('loadPatientWeights:MissingFile', ...
        ['Could not find weights for patient %s. Tried:\n  %s\n  %s'], ...
        patient, primary, fallback);
end
K = load(weightsFile);
end
