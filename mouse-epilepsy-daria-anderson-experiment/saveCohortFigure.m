function saveCohortFigure(fig, outDir, baseName)
%SAVECOHORTFIGURE  Export 1500x988 PNG (200 dpi) plus matching PDF and SVG.
%
%   saveCohortFigure(fig, outDir, 'region_susc_outliers')

if nargin < 3 || isempty(baseName)
    error('saveCohortFigure:NeedName', 'baseName is required.');
end
if nargin < 2 || isempty(outDir)
    outDir = 'comparison_figures';
end
if ~isfolder(outDir)
    mkdir(outDir);
end

if ~isgraphics(fig, 'figure')
    error('saveCohortFigure:InvalidFig', 'Figure handle is not valid.');
end

% 7.5 in x 4.94 in @ 200 dpi -> 1500 x 988 px
wIn = 7.5;
hIn = 4.94;
dpi = 200;
fig.Units = 'inches';
fig.Position = [1 1 wIn hIn];
fig.PaperUnits = 'inches';
fig.PaperSize = [wIn hIn];
fig.PaperPosition = [0 0 wIn hIn];
fig.PaperPositionMode = 'manual';

pngFile = fullfile(outDir, [baseName '.png']);
pdfFile = fullfile(outDir, [baseName '.pdf']);
svgFile = fullfile(outDir, [baseName '.svg']);

drawnow;
print(fig, pngFile, '-dpng', sprintf('-r%d', dpi), '-painters');
drawnow;
exportgraphics(fig, pdfFile, 'ContentType', 'vector');
drawnow;
try
    exportgraphics(fig, svgFile, 'ContentType', 'vector');
catch
    print(fig, svgFile, '-dsvg', '-painters');
end
end
