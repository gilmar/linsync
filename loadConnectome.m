function [K, labels, info] = loadConnectome(target, varargin)
%LOADCONNECTOME  Load a labelled weighted connectome for one subject.
%
%   [K, labels, info] = loadConnectome('data/Case_1/connectome.csv')
%   [K, labels, info] = loadConnectome('Case_1', 'DataRoot', dataRoot, ...
%                                      'FileName', 'connectome.csv')
%
%   Accepts either a direct path to a connectome file, or a subject ID plus
%   'DataRoot', in which case the file is read from
%   <DataRoot>/<subjectId>/<FileName>. The second form is what the cohort
%   pipeline uses, so a study is laid out as one folder per subject.
%
%   Two file shapes are understood:
%     * Labelled CSV -- header row of region names, first column holding
%       the matching row labels, remaining NxN block the weights. Row and
%       column labels are cross-checked and a mismatch warns.
%     * Plain numeric matrix -- whitespace- or comma-delimited NxN block
%       with no labels (e.g. a TVB weights.txt). Labels are generated as
%       node_001 ... node_NNN.
%
%   Name-value options
%     'DataRoot'             : '' -- when set, `target` is read as a subject ID
%     'FileName'             : 'connectome.csv' -- file inside the subject folder
%     'SubjectId'            : '' -- override the ID recorded in info
%     'TrivialLabelSuffixes' : {} -- label suffixes marking non-physical
%                              placeholder rows (e.g. {'BACKGROUND','_MASK'})
%     'WarnAsymmetry'        : true -- warn when the matrix is far from symmetric
%
%   Outputs
%     K      : NxN double weighted adjacency matrix, as stored
%     labels : Nx1 cellstr of node labels
%     info   : struct with fields
%                subjectId, sourceFile, N, labels
%                isTrivial     Nx1 logical, see below
%                symmetryError max(abs(K - K')), 0 for a symmetric matrix
%                isSymmetric   symmetryError below a relative tolerance
%
%   Trivial nodes
%     A node is flagged trivial when EITHER its off-diagonal row and column
%     are exactly zero (it is disconnected, so per-node measures are
%     meaningless), OR its label ends with one of 'TrivialLabelSuffixes'.
%     The label rule exists because atlas placeholder rows -- background,
%     masks -- often carry tiny non-zero weights that survive normalisation
%     and then dominate a ranking despite having no anatomical meaning.
%     Trivial nodes are excluded from rankings, sweeps and statistical
%     families downstream; they are never deleted, so indices stay aligned
%     with the source file.
%
%   See also LISTCONNECTOMESUBJECTS, RUNSTABILITYCENTRALITIES.

p = inputParser;
addParameter(p, 'DataRoot',             '', @(s) ischar(s) || isstring(s));
addParameter(p, 'FileName',             'connectome.csv', @(s) ischar(s) || isstring(s));
addParameter(p, 'SubjectId',            '', @(s) ischar(s) || isstring(s));
addParameter(p, 'TrivialLabelSuffixes', {}, @(c) iscell(c) || isstring(c) || ischar(c));
addParameter(p, 'WarnAsymmetry',        true, @islogical);
parse(p, varargin{:});
opts = p.Results;

target = char(target);
dataRoot = char(opts.DataRoot);
if isempty(dataRoot)
    sourceFile = target;
    subjectId = char(opts.SubjectId);
    if isempty(subjectId)
        subjectId = subjectIdFromPath(sourceFile);
    end
else
    subjectId = char(opts.SubjectId);
    if isempty(subjectId)
        subjectId = target;
    end
    sourceFile = fullfile(dataRoot, target, char(opts.FileName));
end

if ~isfile(sourceFile)
    error('loadConnectome:MissingFile', ...
        'Connectome file not found for "%s": %s', subjectId, sourceFile);
end

[K, labels] = readConnectomeFile(sourceFile, subjectId);
K = double(K);

[N, M] = size(K);
if N ~= M
    error('loadConnectome:NotSquare', ...
        'Connectome for %s is %dx%d; expected a square matrix (%s).', ...
        subjectId, N, M, sourceFile);
end
if numel(labels) ~= N
    error('loadConnectome:LabelMismatch', ...
        'Connectome for %s has %d rows but %d labels (%s).', ...
        subjectId, N, numel(labels), sourceFile);
end

trivialSuffixes = normaliseSuffixes(opts.TrivialLabelSuffixes);
isTrivial = detectTrivialNodes(K, labels, trivialSuffixes);

symmetryError = max(abs(K - K.'), [], 'all');
scale = max(abs(K), [], 'all');
isSymmetric = ~(scale > 0) || (symmetryError <= 1e-9 * scale);
if opts.WarnAsymmetry && ~isSymmetric
    warning('loadConnectome:Asymmetric', ...
        ['Connectome for %s is not symmetric (max |K - K''| = %.3g). ' ...
         'That is fine for a directed study; check the source if the data ' ...
         'was meant to be an undirected structural connectome.'], ...
        subjectId, symmetryError);
end

info = struct();
info.subjectId     = subjectId;
info.sourceFile    = sourceFile;
info.N             = N;
info.labels        = labels;
info.isTrivial     = isTrivial;
info.symmetryError = symmetryError;
info.isSymmetric   = isSymmetric;
info.trivialLabelSuffixes = trivialSuffixes;
end

%% ------------------------------------------------------------------
function [K, labels] = readConnectomeFile(sourceFile, subjectId)
%READCONNECTOMEFILE  Labelled CSV if we can, plain numeric matrix otherwise.
[~, ~, ext] = fileparts(sourceFile);

if strcmpi(ext, '.csv')
    T = readtable(sourceFile, 'VariableNamingRule', 'preserve');
    varNames = T.Properties.VariableNames;
    firstCol = T{:, 1};
    if iscell(firstCol) || isstring(firstCol) || ischar(firstCol)
        labels = varNames(2:end).';
        K = table2array(T(:, 2:end));
        checkLabelOrder(firstCol, labels, subjectId);
        return;
    end
    K = table2array(T);
    labels = varNames(:);
    if all(cellfun(@(s) ~isempty(regexp(s, '^(Var|Column)\d+$', 'once')), labels))
        labels = defaultLabels(size(K, 1));
    end
    return;
end

K = readmatrix(sourceFile);
labels = defaultLabels(size(K, 1));
end

%% ------------------------------------------------------------------
function checkLabelOrder(rowLabels, colLabels, subjectId)
rowLabels = string(rowLabels);
colLabels = string(colLabels);
if numel(rowLabels) ~= numel(colLabels)
    return;
end
mismatch = find(rowLabels(:) ~= colLabels(:), 1, 'first');
if ~isempty(mismatch)
    warning('loadConnectome:LabelOrder', ...
        ['Row/column label mismatch at index %d for %s (row="%s", col="%s"). ' ...
         'Row i and column i must name the same node.'], ...
        mismatch, subjectId, rowLabels(mismatch), colLabels(mismatch));
end
end

%% ------------------------------------------------------------------
function labels = defaultLabels(N)
labels = cell(N, 1);
for i = 1:N
    labels{i} = sprintf('node_%03d', i);
end
end

%% ------------------------------------------------------------------
function suffixes = normaliseSuffixes(suffixes)
if ischar(suffixes)
    suffixes = {suffixes};
end
suffixes = cellstr(suffixes);
suffixes = suffixes(:).';
suffixes = suffixes(~cellfun(@isempty, suffixes));
end

%% ------------------------------------------------------------------
function isTrivial = detectTrivialNodes(K, labels, trivialSuffixes)
Koff = K - diag(diag(K));
isTrivialByWeight = (sum(abs(Koff), 1).' + sum(abs(Koff), 2)) == 0;

labelStr = string(labels(:));
isTrivialByLabel = false(numel(labelStr), 1);
for s = 1:numel(trivialSuffixes)
    isTrivialByLabel = isTrivialByLabel | endsWith(labelStr, trivialSuffixes{s});
end

isTrivial = isTrivialByWeight | isTrivialByLabel;
end

%% ------------------------------------------------------------------
function subjectId = subjectIdFromPath(sourceFile)
[folder, base] = fileparts(sourceFile);
[~, leaf] = fileparts(folder);
if isempty(leaf)
    subjectId = base;
else
    subjectId = leaf;
end
end
