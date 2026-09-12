function paths = saveFigureBoth(fig, basePath, resolution)
%SAVEFIGUREBOTH  Save a figure as both .fig and .png next to each other.
%
%   paths = saveFigureBoth(fig, fullfile(outDir, 'my_figure'))
%   paths = saveFigureBoth(fig, basePath, 300)
%
%   basePath is the output path WITHOUT an extension. The .fig keeps the
%   figure editable for later relabelling; the .png is what goes into a
%   draft, a slide or a README. Every pipeline figure is written both ways
%   so no result is trapped in a format someone cannot open.
%
%   Falls back to saveas when exportgraphics is unavailable (pre-R2020a) or
%   refuses the figure.
%
%   Returns a struct with the `fig` and `png` paths written.

if nargin < 3 || isempty(resolution)
    resolution = 200;
end
if ~isgraphics(fig, 'figure')
    error('saveFigureBoth:InvalidFigure', 'First argument must be a figure handle.');
end

basePath = char(basePath);
outDir = fileparts(basePath);
if ~isempty(outDir) && ~exist(outDir, 'dir')
    mkdir(outDir);
end

paths = struct('fig', [basePath '.fig'], 'png', [basePath '.png']);

savefig(fig, paths.fig);
try
    exportgraphics(fig, paths.png, 'Resolution', resolution);
catch
    saveas(fig, paths.png);
end
end
