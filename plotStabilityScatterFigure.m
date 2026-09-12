function fig = plotStabilityScatterFigure(x0_crit, D_susceptibility, D_influence, subjectLabel, x0_base, figName)
%PLOTSTABILITYSCATTERFIGURE  Two-panel scatter of stability centralities vs x0^c.

if nargin < 6 || isempty(figName)
    figName = sprintf('Stability centralities (%s, healthy)', subjectLabel);
end

fig = figure('Name', figName, 'Position', [100 100 1200 480]);

subplot(1, 2, 1);
scatter(x0_crit, D_susceptibility, 48, 'filled');
xlabel('critical excitability  x^c_{0,i}');
ylabel('D(\rightarrow i)');
title('Stability susceptibility centrality');
grid on;

subplot(1, 2, 2);
scatter(x0_crit, D_influence, 48, 'filled', 'MarkerFaceColor', [0.85 0.33 0.10]);
xlabel('critical excitability  x^c_{0,i}');
ylabel('D(k \rightarrow)');
title('Stability influence centrality');
grid on;

sgtitle(sprintf('%s, fully healthy network (x_{0,j\\neqi} = %.2f)', ...
    subjectLabel, x0_base), 'Interpreter', 'tex');
end
