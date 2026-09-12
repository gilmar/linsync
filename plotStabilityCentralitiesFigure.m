function fig = plotStabilityCentralitiesFigure(results, figName)
%PLOTSTABILITYCENTRALITIESFIGURE  Per-subject figure for a stability-centralities run.
%
%   fig = plotStabilityCentralitiesFigure(results)
%   fig = plotStabilityCentralitiesFigure(results, 'My title')
%
%   Picks the layout that actually carries information for the scheme that
%   produced `results`:
%
%     * Epileptor schemes -- two-panel scatter of D(->i) and D(k->) against
%       the critical excitability x0^c, which is the validation plot: if the
%       stability centralities predict epileptogenicity, these should trend.
%     * Linear asymmetric schemes (column) -- two bar panels, since D(->i)
%       and D(k->) differ.
%     * Linear symmetric schemes (parkes) -- one bar panel, because a
%       symmetric K makes D(k->) identical to D(->i) by construction and a
%       second panel would just repeat it.
%
%   Takes the struct returned by RUNSTABILITYCENTRALITIES.
%
%   See also RUNSTABILITYCENTRALITIES, PLOTSTABILITYSCATTERFIGURE.

if nargin < 2
    figName = '';
end

subjectId = char(results.subjectId);
subjectTex = strrep(subjectId, '_', '\_');
scheme = char(results.normalisation);
N = results.N;

if isfield(results, 'schemeInfo') && isstruct(results.schemeInfo)
    schemeInfo = results.schemeInfo;
else
    schemeInfo = normalisationSchemeInfo(scheme);
end

hasCriticalX0 = isfield(results, 'x0_crit') && any(isfinite(results.x0_crit));

if isempty(figName)
    figName = sprintf('Stability centralities - %s (%s)', subjectId, scheme);
end

if schemeInfo.usesEpileptor && hasCriticalX0
    fig = plotStabilityScatterFigure(results.x0_crit, ...
        results.D_susceptibility, results.D_influence, ...
        sprintf('%s, normalisation = %s', subjectTex, scheme), ...
        results.x0_base, figName);
    return;
end

if schemeInfo.isSymmetric
    fig = figure('Name', figName, 'Position', [100 100 860 520]);
    bar(results.D_susceptibility, 'FaceColor', [0.20 0.40 0.80]);
    xlabel('node index');
    ylabel('D(\rightarrow i)');
    xlim([0.5 N + 0.5]);
    grid on;
    sgtitle({sprintf('%s  --  %s normalisation, D_{st} = %.4f', ...
                     subjectTex, scheme, results.D_st_healthy), ...
             'D(k\rightarrow) \equiv D(\rightarrow i) for a symmetric coupling matrix'}, ...
            'Interpreter', 'tex');
    return;
end

fig = figure('Name', figName, 'Position', [100 100 1300 520]);

subplot(1, 2, 1);
bar(results.D_susceptibility, 'FaceColor', [0.20 0.40 0.80]);
xlabel('node index');
ylabel('D(\rightarrow i)');
xlim([0.5 N + 0.5]);
grid on;
title('Stability susceptibility centrality');

subplot(1, 2, 2);
bar(results.D_influence, 'FaceColor', [0.85 0.33 0.10]);
xlabel('node index');
ylabel('D(k \rightarrow)');
xlim([0.5 N + 0.5]);
grid on;
title('Stability influence centrality');

sgtitle(sprintf('%s  --  %s normalisation, D_{st} = %.4f', ...
                subjectTex, scheme, results.D_st_healthy), 'Interpreter', 'tex');
end
