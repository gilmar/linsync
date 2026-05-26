function fig = plotRobustVsSuggestiveSummary(andersonSummary, lrSummary, varargin)
%PLOTROBUSTVSSUGGESTIVESUMMARY  Two-card ROBUST vs SUGGESTIVE summary from comparison data.

p = inputParser;
addParameter(p, 'OutputDir', '', @(s) ischar(s) || isstring(s));
addParameter(p, 'Save', false, @islogical);
addParameter(p, 'OverrideCfg', struct(), @isstruct);
parse(p, varargin{:});
opts = p.Results;

P = cohortFigStyle();
content = deriveRobustSuggestiveContent(andersonSummary, lrSummary, opts.OverrideCfg);

fig = figure('Units', 'inches', 'Position', [1 1 7.5 4.94], ...
    'Color', 'w', 'Visible', 'off', 'MenuBar', 'none', 'ToolBar', 'none');
ax = axes(fig, 'Position', [0 0 1 1]);
axis(ax, 'off');
xlim(ax, [0 10]);
ylim(ax, [0 10]);

rectangle(ax, 'Position', [0.35 1.0 4.2 7.8], 'Curvature', 0.12, ...
    'EdgeColor', P.teal, 'LineWidth', 2, 'FaceColor', P.creamFill);
text(ax, 0.65, 8.2, 'ROBUST', 'FontWeight', 'bold', 'Color', P.teal, 'FontSize', 16);
text(ax, 0.65, 7.55, 'large effect, flagged in both mice', 'Color', P.teal, 'FontSize', 9);
yL = 6.9;
for k = 1:numel(content.robust.items)
    item = content.robust.items{k};
    text(ax, 0.65, yL, item{1}, 'Color', P.ink, 'FontSize', 10, ...
        'FontWeight', 'bold', 'Interpreter', 'none');
    text(ax, 0.85, yL - 0.55, [char(10003) ' ' item{2}], 'Color', P.teal, ...
        'FontSize', 9, 'Interpreter', 'tex');
    yL = yL - 1.35;
end
text(ax, 0.65, 1.35, 'Large effect sizes relative to the control band — the dependable signals.', ...
    'Color', [0.35 0.35 0.35], 'FontSize', 8.5);

rectangle(ax, 'Position', [5.45 1.0 4.2 7.8], 'Curvature', 0.12, ...
    'EdgeColor', P.amber, 'LineWidth', 2, 'FaceColor', P.cream);
text(ax, 5.75, 8.2, 'SUGGESTIVE', 'FontWeight', 'bold', 'Color', P.amber, 'FontSize', 16);
text(ax, 5.75, 7.55, 'interpret with caution at this n', 'Color', P.amber, 'FontSize', 9);
yR = 6.9;
for k = 1:numel(content.suggestive.items)
    item = content.suggestive.items{k};
    text(ax, 5.75, yR, ['\bullet ' item{1}], 'Color', P.ink, 'FontSize', 9.5, 'Interpreter', 'none');
    text(ax, 5.95, yR - 0.5, item{2}, 'Color', P.mute, 'FontSize', 8.5, 'Interpreter', 'tex');
    yR = yR - 1.25;
end
text(ax, 5.75, 1.35, 'Small effects sitting near the Bonferroni line at n = 2 vs 4.', ...
    'Color', [0.35 0.35 0.35], 'FontSize', 8.5);

text(ax, 5.0, 0.35, ...
    'Bonferroni significance is fragile at this sample size — effect sizes, not z alone, are the dependable signal.', ...
    'FontAngle', 'italic', 'Color', P.mute, 'FontSize', 8, ...
    'HorizontalAlignment', 'center', 'Interpreter', 'none');

if opts.Save
    outDir = resolveFigureOutDir(opts.OutputDir, andersonSummary);
    saveCohortFigure(fig, outDir, 'robust_vs_suggestive_summary');
end
end
