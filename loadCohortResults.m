function cohort = loadCohortResults(resultsDir, scheme, varargin)
%LOADCOHORTRESULTS  Load every per-subject stability result of a run.
%
%   cohort = loadCohortResults(resultsDir, 'column')
%   cohort = loadCohortResults(resultsDir, scheme, 'CasePrefix', 'Case', ...
%                              'ControlPrefix', 'Ctrl', 'RequireGroups', true)
%
%   Shared entry point for every cohort-level comparison. Discovers the
%   per-subject result files of a run, loads them, checks that they share
%   one parcellation, and applies the case/control split by ID prefix.
%   Doing this in one place means every comparison sees exactly the same
%   cohort and the same node ordering.
%
%   Name-value options
%     'CasePrefix'    : 'Case'    -- subject-ID prefix for the case group
%     'ControlPrefix' : 'Control' -- subject-ID prefix for the control group
%     'RequireGroups' : false     -- error when either group is empty
%
%   Output struct
%     scheme            : the normalisation scheme loaded
%     resultsDir        : the run directory
%     subjects          : nSubjects x 1 cellstr of IDs
%     results           : nSubjects x 1 cell of per-subject result structs
%     labels            : N x 1 cellstr of node labels (from the first subject)
%     N                 : number of nodes
%     isCase, isControl : nSubjects x 1 logical masks
%     caseIds, controlIds
%     isTrivialAny      : N x 1 logical, trivial in any subject
%     isTrivialControl  : N x 1 logical, trivial in any control subject --
%                         the mask to exclude when comparing against a
%                         control mean, since a control-trivial node makes
%                         that mean meaningless
%
%   See also LISTSUBJECTRESULTFILES, PACKCOHORTMETRIC, COHORTGROUPMASK.

p = inputParser;
addParameter(p, 'CasePrefix',    'Case',    @(s) ischar(s) || isstring(s));
addParameter(p, 'ControlPrefix', 'Control', @(s) ischar(s) || isstring(s));
addParameter(p, 'RequireGroups', false, @islogical);
parse(p, varargin{:});
opts = p.Results;

scheme = normalisationSchemeInfo(scheme).scheme;

files = listSubjectResultFiles(resultsDir, scheme);
if isempty(files)
    error('loadCohortResults:NoResults', ...
        ['No %s_*_%s_results.mat found under %s. ' ...
         'Run runCohortStabilityCentralities for this scheme first.'], ...
        experimentResultPrefix(), scheme, resultsDir);
end

subjects = cell(0, 1);
loaded = cell(0, 1);
for k = 1:numel(files)
    s = load(fullfile(files(k).folder, files(k).name));
    if ~isfield(s, 'results')
        continue;
    end
    r = s.results;
    if ~isfield(r, 'subjectId')
        continue;
    end
    subjects{end+1, 1} = char(r.subjectId); %#ok<AGROW>
    loaded{end+1, 1} = r;                   %#ok<AGROW>
end
if isempty(subjects)
    error('loadCohortResults:EmptyResults', ...
        'The result files under %s contained no usable results structs.', resultsDir);
end

[subjects, order] = sort(subjects);
loaded = loaded(order);

labels = loaded{1}.labels(:);
N = numel(labels);
for k = 2:numel(loaded)
    if numel(loaded{k}.labels) ~= N || ...
            ~all(string(loaded{k}.labels(:)) == string(labels))
        warning('loadCohortResults:LabelMismatch', ...
            ['Node labels for %s differ from %s. Comparisons assume one ' ...
             'shared parcellation and align subjects by index.'], ...
            subjects{k}, subjects{1});
    end
end

[isCase, isControl] = cohortGroupMask(subjects, opts.CasePrefix, opts.ControlPrefix);
if opts.RequireGroups
    if ~any(isCase)
        error('loadCohortResults:NoCaseSubjects', ...
            'No subject ID starts with the case prefix "%s". Present: %s.', ...
            char(opts.CasePrefix), strjoin(subjects, ', '));
    end
    if ~any(isControl)
        error('loadCohortResults:NoControlSubjects', ...
            'No subject ID starts with the control prefix "%s". Present: %s.', ...
            char(opts.ControlPrefix), strjoin(subjects, ', '));
    end
end

isTrivialAny = false(N, 1);
isTrivialControl = false(N, 1);
for k = 1:numel(loaded)
    if ~isfield(loaded{k}, 'isTrivial') || isempty(loaded{k}.isTrivial)
        continue;
    end
    triv = logical(loaded{k}.isTrivial(:));
    isTrivialAny = isTrivialAny | triv;
    if isControl(k)
        isTrivialControl = isTrivialControl | triv;
    end
end

cohort = struct();
cohort.scheme           = scheme;
cohort.resultsDir       = char(resultsDir);
cohort.subjects         = subjects;
cohort.results          = loaded;
cohort.labels           = cellstr(labels);
cohort.N                = N;
cohort.isCase           = isCase;
cohort.isControl        = isControl;
cohort.caseIds          = subjects(isCase);
cohort.controlIds       = subjects(isControl);
cohort.casePrefix       = char(opts.CasePrefix);
cohort.controlPrefix    = char(opts.ControlPrefix);
cohort.isTrivialAny     = isTrivialAny;
cohort.isTrivialControl = isTrivialControl;
end
