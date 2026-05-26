function P = cohortFigStyle()
%COHORTFIGSTYLE  Shared palette and default figure/axes styling for cohort comparison figures.
%
%   P = cohortFigStyle();
%   fig = figure(...); ax = axes(fig); applyCohortAxesStyle(ax, P);

P.red       = [0.902 0.275 0.149];
P.a1        = [0.949 0.627 0.494];
P.navy      = [0.102 0.204 0.369];
P.teal      = [0.247 0.561 0.482];
P.ink       = [0.169 0.169 0.169];
P.mute      = [0.541 0.541 0.541];
P.grid      = [0.886 0.855 0.827];
P.cream     = [0.984 0.929 0.910];
P.amber     = [0.690 0.416 0.118];
P.grey      = [0.420 0.420 0.420];
P.creamFill = [0.945 0.969 0.957];

set(groot, 'defaultAxesFontName', 'Arial', ...
    'defaultTextFontName', 'Arial', ...
    'defaultFigureColor', 'w');
end
