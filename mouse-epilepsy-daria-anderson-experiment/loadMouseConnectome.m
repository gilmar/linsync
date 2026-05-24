function [K, labels, info] = loadMouseConnectome(mouseId)
%LOADMOUSECONNECTOME  Load a coarse-grained weighted adjacency matrix
%for a Daria/Anderson mouse from data/<mouseId>/fine_family_labelled_coarse.csv.
%
%   Inputs
%     mouseId : string -- e.g. 'Anderson_1', 'Arnold_3'
%
%   Outputs
%     K      : NxN double, weighted (symmetric) adjacency matrix
%     labels : Nx1 cell array of region names (e.g. 'L-CORTEX_VISUAL')
%     info   : struct with mouseId, csvFile, refImage, isTrivial,
%              symmetryError fields. isTrivial(i) is true for rows that
%              are non-physical atlas placeholders and should be
%              excluded from per-node rankings downstream. A row counts
%              as trivial if either:
%                (a) its entire row + column (excluding diagonal) is
%                    exactly zero, or
%                (b) its label is a placeholder by convention --
%                    '*BACKGROUND' or '*_MASK' (e.g. ISOCORTEX_MASK,
%                    GLOBAL_MASKS). The latter rule catches the MASK
%                    rows whose connectivity is non-zero but tiny, and
%                    which the TVB 95th-percentile rescale otherwise
%                    pushes to the top of the D(->i) ranking.

experimentRoot = fileparts(mfilename('fullpath'));
csvFile = fullfile(experimentRoot, 'data', mouseId, 'fine_family_labelled_coarse.csv');
if ~exist(csvFile, 'file')
    error('loadMouseConnectome:MissingFile', ...
        'Could not find connectome CSV for %s at %s', mouseId, csvFile);
end

T = readtable(csvFile, 'VariableNamingRule', 'preserve');

% First column is the row label ("node"); remaining columns are region names
varNames  = T.Properties.VariableNames;
labels    = varNames(2:end)';
rowLabels = T{:, 1};

K = T{:, 2:end};
if ~isnumeric(K)
    K = table2array(T(:, 2:end));
end
K = double(K);

[N, M] = size(K);
if N ~= M
    error('loadMouseConnectome:NotSquare', ...
        'K is %dx%d for %s', N, M, mouseId);
end
if N ~= length(labels)
    error('loadMouseConnectome:LabelMismatch', ...
        '%d rows but %d column labels for %s', N, length(labels), mouseId);
end

% Sanity-check label ordering: row label name should equal column label
if iscell(rowLabels) || isstring(rowLabels)
    rowLabels = string(rowLabels);
    colLabels = string(labels);
    mismatch  = find(rowLabels(:) ~= colLabels(:), 1, 'first');
    if ~isempty(mismatch)
        warning('loadMouseConnectome:LabelOrder', ...
            'row/column label mismatch at index %d (row="%s", col="%s") for %s', ...
            mismatch, rowLabels(mismatch), colLabels(mismatch), mouseId);
    end
end

% Identify trivial nodes. Two rules combined (logical OR):
%   (a) connectivity-based: row + column (excluding diagonal) is exactly zero
%   (b) label-based: row name ends in 'BACKGROUND' or '_MASK' / 'MASKS'
% The label rule catches atlas placeholders such as ISOCORTEX_MASK and
% GLOBAL_MASKS that carry tiny non-zero weights but are not physical brain
% regions; without it, the TVB normalisation amplifies those rows into the
% top of the D(->i) ranking even though they have no biological meaning.
Knd = K - diag(diag(K));   % non-diagonal contributions only
isTrivialByWeight = (sum(abs(Knd), 1)' + sum(abs(Knd), 2)) == 0;
labelStr = string(labels);
isTrivialByLabel  = endsWith(labelStr, 'BACKGROUND') | ...
                    endsWith(labelStr, '_MASK')      | ...
                    endsWith(labelStr, '_MASKS');
isTrivial = isTrivialByWeight | isTrivialByLabel(:);

% Symmetry error (structural connectomes should be symmetric up to noise)
symmetryError = max(abs(K - K'), [], 'all');

info = struct();
info.mouseId       = mouseId;
info.csvFile       = csvFile;
info.refImage      = fullfile(experimentRoot, 'data', mouseId, ...
                               sprintf('coarse_connectome_%s.png', lower(mouseId)));
info.labels        = labels;
info.isTrivial     = isTrivial;
info.symmetryError = symmetryError;
info.N             = N;
end
