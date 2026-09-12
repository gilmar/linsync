function info = makeSyntheticCohort(dataRoot, varargin)
%MAKESYNTHETICCOHORT  Generate a synthetic bilateral connectome cohort.
%
%   info = makeSyntheticCohort('data')
%   info = makeSyntheticCohort(dataRoot, 'NCase', 2, 'NControl', 4, 'Seed', 7)
%
%   Writes one folder per subject under `dataRoot`, each containing a
%   labelled connectome CSV in the format LOADCONNECTOME expects. Use it to
%   exercise the whole pipeline end to end without any real data: for a
%   demo, for a regression test, or to sanity-check that an analysis can
%   recover an effect you planted yourself.
%
%   Structure of the generated networks
%     * Bilateral atlas: every region appears as <LeftPrefix><NAME> and
%       <RightPrefix><NAME>, so the laterality analysis has homotopic pairs
%       to work with.
%     * Modular: within-hemisphere edges are dense, cross-hemisphere edges
%       are sparse, and each region has a strong homotopic link to its twin
%       -- the coarse structure real connectomes show.
%     * Placeholder rows: a disconnected <prefix>BACKGROUND node per
%       hemisphere, so the trivial-node handling is exercised too.
%     * Symmetric and non-negative, as a structural connectome should be.
%
%   The planted effect
%     Case subjects get their LEFT-hemisphere copies of the regions named
%     by 'AffectedRegions' up-weighted by 'EffectSize'. That should show up
%     downstream as elevated D(->i) for those nodes and as a left-right
%     asymmetry -- which is exactly what the cohort comparison and the
%     laterality analysis are meant to detect.
%
%   Name-value options
%     'NCase'           : 2      case subjects
%     'NControl'        : 4      control subjects
%     'CasePrefix'      : 'Case'
%     'ControlPrefix'   : 'Ctrl'
%     'Regions'         : {}     region base names (default: 12 generic names)
%     'LeftPrefix'      : 'L-'
%     'RightPrefix'     : 'R-'
%     'AffectedRegions' : {'HIPPOCAMPUS', 'CORTEX_TEMPORAL'}
%     'EffectSize'      : 1.8    multiplier on the affected left-hemisphere edges
%     'SubjectNoise'    : 0.15   relative spread between subjects
%     'FileName'        : 'connectome.csv'
%     'Seed'            : 42     RNG seed, so the cohort is reproducible
%     'Overwrite'       : true
%
%   Output struct
%     dataRoot, subjects, caseIds, controlIds, labels, N, affectedNodes
%
%   See also LOADCONNECTOME, RUNCONNECTOMEEXPERIMENT, LISTCONNECTOMESUBJECTS.

if nargin < 1 || isempty(dataRoot)
    dataRoot = fullfile(pwd, 'data');
end

defaultRegions = { ...
    'CORTEX_VISUAL', 'CORTEX_MOTOR', 'CORTEX_SOMATOSENSORY', 'CORTEX_AUDITORY', ...
    'CORTEX_TEMPORAL', 'CORTEX_PREFRONTAL', 'HIPPOCAMPUS', 'AMYGDALA', ...
    'THALAMUS', 'STRIATUM', 'HYPOTHALAMUS', 'CEREBELLUM'};

p = inputParser;
addParameter(p, 'NCase',           2,  @(x) isscalar(x) && x >= 1);
addParameter(p, 'NControl',        4,  @(x) isscalar(x) && x >= 1);
addParameter(p, 'CasePrefix',      'Case', @(s) ischar(s) || isstring(s));
addParameter(p, 'ControlPrefix',   'Ctrl', @(s) ischar(s) || isstring(s));
addParameter(p, 'Regions',         {}, @(c) iscell(c) || isstring(c));
addParameter(p, 'LeftPrefix',      'L-', @(s) ischar(s) || isstring(s));
addParameter(p, 'RightPrefix',     'R-', @(s) ischar(s) || isstring(s));
addParameter(p, 'AffectedRegions', {'HIPPOCAMPUS', 'CORTEX_TEMPORAL'}, ...
    @(c) iscell(c) || isstring(c));
addParameter(p, 'EffectSize',      1.8,  @isscalar);
addParameter(p, 'SubjectNoise',    0.15, @isscalar);
addParameter(p, 'FileName',        'connectome.csv', @(s) ischar(s) || isstring(s));
addParameter(p, 'Seed',            42,   @isscalar);
addParameter(p, 'Overwrite',       true, @islogical);
parse(p, varargin{:});
opts = p.Results;

regions = opts.Regions;
if isempty(regions)
    regions = defaultRegions;
else
    regions = cellstr(regions);
end
regions = regions(:).';

leftPrefix = char(opts.LeftPrefix);
rightPrefix = char(opts.RightPrefix);

% Labels: all left regions, then all right regions, each hemisphere ending
% in a disconnected BACKGROUND placeholder row.
leftLabels  = [strcat(leftPrefix,  regions), {[leftPrefix  'BACKGROUND']}];
rightLabels = [strcat(rightPrefix, regions), {[rightPrefix 'BACKGROUND']}];
labels = [leftLabels, rightLabels].';
N = numel(labels);
nPerHemi = numel(regions) + 1;
backgroundIdx = [nPerHemi, 2 * nPerHemi];

rng(opts.Seed, 'twister');

baseK = buildBaseNetwork(nPerHemi, backgroundIdx);

affected = cellstr(opts.AffectedRegions);
affectedNodes = zeros(0, 1);
for k = 1:numel(affected)
    idx = find(strcmp(regions, affected{k}), 1);
    if isempty(idx)
        warning('makeSyntheticCohort:UnknownRegion', ...
            'Affected region "%s" is not in the region list; ignoring it.', affected{k});
        continue;
    end
    affectedNodes(end+1, 1) = idx; %#ok<AGROW>   index within the left hemisphere
end

subjects = cell(0, 1);
caseIds = cell(0, 1);
controlIds = cell(0, 1);

for s = 1:opts.NControl
    id = sprintf('%s_%d', char(opts.ControlPrefix), s);
    writeSubject(dataRoot, id, opts.FileName, ...
        perturb(baseK, opts.SubjectNoise), labels, opts.Overwrite);
    subjects{end+1, 1} = id;    %#ok<AGROW>
    controlIds{end+1, 1} = id;  %#ok<AGROW>
end

for s = 1:opts.NCase
    id = sprintf('%s_%d', char(opts.CasePrefix), s);
    K = applyCaseEffect(perturb(baseK, opts.SubjectNoise), affectedNodes, opts.EffectSize);
    writeSubject(dataRoot, id, opts.FileName, K, labels, opts.Overwrite);
    subjects{end+1, 1} = id; %#ok<AGROW>
    caseIds{end+1, 1} = id;  %#ok<AGROW>
end

info = struct();
info.dataRoot = char(dataRoot);
info.subjects = sort(subjects);
info.caseIds = caseIds;
info.controlIds = controlIds;
info.labels = labels;
info.N = N;
info.affectedNodes = affectedNodes;
info.affectedLabels = leftLabels(affectedNodes);

fprintf('Synthetic cohort written to %s\n', info.dataRoot);
fprintf('  %d subject(s): %s\n', numel(info.subjects), strjoin(info.subjects.', ', '));
fprintf('  %d nodes, planted effect on: %s\n', N, strjoin(info.affectedLabels, ', '));
end

%% ------------------------------------------------------------------
function K = buildBaseNetwork(nPerHemi, backgroundIdx)
%BUILDBASENETWORK  Symmetric two-hemisphere template with homotopic links.
N = 2 * nPerHemi;
K = zeros(N);

withinDensity = 0.55;
crossDensity  = 0.12;

for i = 1:N
    for j = i + 1:N
        sameHemisphere = (i <= nPerHemi) == (j <= nPerHemi);
        if sameHemisphere
            density = withinDensity;
            scale = 1.0;
        else
            density = crossDensity;
            scale = 0.5;
        end
        if rand() < density
            K(i, j) = scale * (0.2 + 0.8 * rand());
        end
    end
end

% Strong homotopic connections: each left region to its right twin.
for i = 1:nPerHemi
    j = i + nPerHemi;
    K(i, j) = 0.9 + 0.2 * rand();
end

K = K + K.';

% Placeholder rows carry no connectivity at all.
K(backgroundIdx, :) = 0;
K(:, backgroundIdx) = 0;
end

%% ------------------------------------------------------------------
function K = perturb(K, noise)
%PERTURB  Per-subject multiplicative jitter that preserves symmetry and sign.
N = size(K, 1);
J = 1 + noise * randn(N);
J = (J + J.') / 2;
K = max(K .* J, 0);
K = (K + K.') / 2;
end

%% ------------------------------------------------------------------
function K = applyCaseEffect(K, affectedNodes, effectSize)
%APPLYCASEEFFECT  Up-weight the left-hemisphere copies of the affected regions.
% Applied to the left side only, so the effect shows up both as elevated
% D(->i) and as a left-right asymmetry.
for k = 1:numel(affectedNodes)
    i = affectedNodes(k);
    K(i, :) = K(i, :) * effectSize;
    K(:, i) = K(:, i) * effectSize;
end
K = (K + K.') / 2;
end

%% ------------------------------------------------------------------
function writeSubject(dataRoot, subjectId, fileName, K, labels, overwrite)
%WRITESUBJECT  Write one labelled connectome CSV.
folder = fullfile(char(dataRoot), subjectId);
if ~isfolder(folder)
    mkdir(folder);
end
csvPath = fullfile(folder, char(fileName));
if isfile(csvPath) && ~overwrite
    return;
end

T = array2table(K, 'VariableNames', matlab.lang.makeValidName(labels));
T.Properties.VariableNames = labels;      % keep the original names verbatim
T = addvars(T, string(labels), 'Before', 1, 'NewVariableNames', 'node');
writetable(T, csvPath);
end
