function prefix = experimentResultPrefix()
%EXPERIMENTRESULTPREFIX  Filename prefix for per-subject result artefacts.
%
%   Every per-subject and cohort-summary file the stability pipeline writes
%   starts with this prefix, which is what lets listSubjectResultFiles
%   discover a run's outputs without being told the subject IDs.

prefix = 'stabilityCentralities';
end
