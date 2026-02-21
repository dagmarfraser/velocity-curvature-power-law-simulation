%% CCC Demo: Constellation Comparison (Pre-Reg Section 7.5 Use Case)
% For ONE trial: 6 pipelines → 15 pairwise differences
% CCC tests: Does simulation predict the PATTERN of inter-pipeline disagreements?
%
% The key insight: β_true cancels in pairwise differences!
%   Simulation: β_rec,P = β_gen + bias_P  → Δβ_PQ = bias_P - bias_Q
%   Empirical:  β_rec,P = β_true + bias_P → Δβ_PQ = bias_P - bias_Q
% So we can validate the framework without knowing ground truth.
%
% Author: D. S. Fraser
% Date: January 2026

addpath('/Users/dsfraser/Dropbox/Brain2Bee/PowerLawSimulationPreReg/src/functions');

rng(42);

% Pipeline names (short versions for labels)
pipelines = {'BWFD-OLS', 'BWFD-LMLS', 'BWFD-IRLS', 'SG-OLS', 'SG-LMLS', 'SG-IRLS'};
pipeShort = {'BO', 'BL', 'BI', 'SO', 'SL', 'SI'};  % Abbreviated
nPipelines = 6;
nPairs = nchoosek(nPipelines, 2);  % 15 unique pairs

% Generate pair indices and labels
pairIdx = zeros(nPairs, 2);
pairLabels = cell(nPairs, 1);
pairType = zeros(nPairs, 1);  % 1=within-BWFD, 2=within-SG, 3=cross-method
k = 1;
for i = 1:nPipelines-1
    for j = i+1:nPipelines
        pairIdx(k,:) = [i, j];
        pairLabels{k} = sprintf('%s vs %s', pipeShort{i}, pipeShort{j});
        if i <= 3 && j <= 3
            pairType(k) = 1;  % Within BWFD
        elseif i >= 4 && j >= 4
            pairType(k) = 2;  % Within SG
        else
            pairType(k) = 3;  % Cross-method
        end
        k = k + 1;
    end
end

% Colors for each type
colors = [0.8 0.3 0.3;   % Within BWFD - red
          0.3 0.6 0.8;   % Within SG - blue  
          0.5 0.5 0.5];  % Cross-method - grey

figure('Position', [100 100 1600 550], 'Color', 'w');

%% Scenario A: Good framework validation
subplot(1,3,1);

beta_true = 0.35;
pipeline_bias_sim = [0.08, 0.04, 0.02, 0.01, 0.005, 0.002];
beta_pred = beta_true + pipeline_bias_sim + randn(1,6)*0.002;
pipeline_bias_emp = pipeline_bias_sim + randn(1,6)*0.01;
beta_obs = beta_true + pipeline_bias_emp + randn(1,6)*0.002;

delta_pred = zeros(nPairs, 1);
delta_obs = zeros(nPairs, 1);
for k = 1:nPairs
    i = pairIdx(k,1); j = pairIdx(k,2);
    delta_pred(k) = beta_pred(i) - beta_pred(j);
    delta_obs(k) = beta_obs(i) - beta_obs(j);
end

res = linCCC_v001(delta_pred, delta_obs);

hold on;
for k = 1:nPairs
    scatter(delta_pred(k), delta_obs(k), 60, colors(pairType(k),:), 'filled', 'MarkerFaceAlpha', 0.8);
    text(delta_pred(k)+0.002, delta_obs(k), pairLabels{k}, 'FontSize', 7, 'Color', colors(pairType(k),:)*0.7);
end
lims = [min([delta_pred; delta_obs])-0.02, max([delta_pred; delta_obs])+0.03];
plot(lims, lims, 'k--', 'LineWidth', 1.5);
xlim(lims); ylim(lims);
xlabel('Δβ predicted (simulation)');
ylabel('Δβ observed (empirical)');
title(sprintf('Good Validation\nr = %.3f, CCC = %.3f', res.r, res.ccc), 'FontSize', 12);
axis square;

% Add key in corner
text(0.02, 0.98, {'Key:', 'BO = BWFD-OLS', 'BL = BWFD-LMLS', 'BI = BWFD-IRLS', ...
    'SO = SG-OLS', 'SL = SG-LMLS', 'SI = SG-IRLS'}, ...
    'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 7, ...
    'BackgroundColor', [1 1 1 0.8], 'EdgeColor', [0.7 0.7 0.7]);

%% Scenario B: Scale mismatch
subplot(1,3,2);

rng(43);
pipeline_bias_sim = [0.08, 0.04, 0.02, 0.01, 0.005, 0.002];
beta_pred = beta_true + pipeline_bias_sim + randn(1,6)*0.002;
pipeline_bias_emp = pipeline_bias_sim * 1.5;
beta_obs = beta_true + pipeline_bias_emp + randn(1,6)*0.002;

for k = 1:nPairs
    i = pairIdx(k,1); j = pairIdx(k,2);
    delta_pred(k) = beta_pred(i) - beta_pred(j);
    delta_obs(k) = beta_obs(i) - beta_obs(j);
end

res = linCCC_v001(delta_pred, delta_obs);

hold on;
for k = 1:nPairs
    scatter(delta_pred(k), delta_obs(k), 60, colors(pairType(k),:), 'filled', 'MarkerFaceAlpha', 0.8);
    text(delta_pred(k)+0.002, delta_obs(k), pairLabels{k}, 'FontSize', 7, 'Color', colors(pairType(k),:)*0.7);
end
lims = [min([delta_pred; delta_obs])-0.02, max([delta_pred; delta_obs])+0.03];
plot(lims, lims, 'k--', 'LineWidth', 1.5);
xlim(lims); ylim(lims);
xlabel('Δβ predicted (simulation)');
ylabel('Δβ observed (empirical)');
title(sprintf('Scale Mismatch\nr = %.3f, CCC = %.3f', res.r, res.ccc), 'FontSize', 12);
axis square;
text(0.05, 0.95, 'Simulation under-predicts', 'Units', 'normalized', 'FontSize', 9, ...
    'Color', [0.7 0.2 0.2], 'VerticalAlignment', 'top');

%% Scenario C: Pattern mismatch
subplot(1,3,3);

rng(44);
pipeline_bias_sim = [0.08, 0.04, 0.02, 0.01, 0.005, 0.002];
beta_pred = beta_true + pipeline_bias_sim + randn(1,6)*0.002;
pipeline_bias_emp = [0.03, 0.02, 0.01, 0.06, 0.02, 0.01];  % SG-OLS surprisingly bad
beta_obs = beta_true + pipeline_bias_emp + randn(1,6)*0.002;

for k = 1:nPairs
    i = pairIdx(k,1); j = pairIdx(k,2);
    delta_pred(k) = beta_pred(i) - beta_pred(j);
    delta_obs(k) = beta_obs(i) - beta_obs(j);
end

res = linCCC_v001(delta_pred, delta_obs);

hold on;
for k = 1:nPairs
    scatter(delta_pred(k), delta_obs(k), 60, colors(pairType(k),:), 'filled', 'MarkerFaceAlpha', 0.8);
    text(delta_pred(k)+0.002, delta_obs(k), pairLabels{k}, 'FontSize', 7, 'Color', colors(pairType(k),:)*0.7);
end
lims = [min([delta_pred; delta_obs])-0.03, max([delta_pred; delta_obs])+0.03];
plot(lims, lims, 'k--', 'LineWidth', 1.5);
xlim(lims); ylim(lims);
xlabel('Δβ predicted (simulation)');
ylabel('Δβ observed (empirical)');
title(sprintf('Pattern Mismatch\nr = %.3f, CCC = %.3f', res.r, res.ccc), 'FontSize', 12);
axis square;
text(0.05, 0.95, 'Wrong pipeline ranking', 'Units', 'normalized', 'FontSize', 9, ...
    'Color', [0.7 0.2 0.2], 'VerticalAlignment', 'top');

sgtitle({'CCC for Framework Validation: 15 Pipeline-Pair Differences per Trial', ...
    '\color[rgb]{0.8,0.3,0.3}Within BWFD  \color[rgb]{0.3,0.6,0.8}Within SG  \color[rgb]{0.5,0.5,0.5}Cross-method (BWFD vs SG)'}, ...
    'FontSize', 13, 'FontWeight', 'bold');

%% Save figure
saveas(gcf, '/Users/dsfraser/Dropbox/Brain2Bee/PowerLawSimulationPreReg/docs/CCC_demo_constellation.png');
fprintf('Figure saved to docs/CCC_demo_constellation.png\n');

%% Print pair breakdown
fprintf('\n=== PAIR CLASSIFICATION ===\n');
fprintf('\nWithin BWFD (red, 3 pairs):\n');
for k = find(pairType==1)'
    fprintf('  %s vs %s\n', pipelines{pairIdx(k,1)}, pipelines{pairIdx(k,2)});
end
fprintf('\nWithin SG (blue, 3 pairs):\n');
for k = find(pairType==2)'
    fprintf('  %s vs %s\n', pipelines{pairIdx(k,1)}, pipelines{pairIdx(k,2)});
end
fprintf('\nCross-method BWFD vs SG (grey, 9 pairs):\n');
for k = find(pairType==3)'
    fprintf('  %s vs %s\n', pipelines{pairIdx(k,1)}, pipelines{pairIdx(k,2)});
end

fprintf('\n=== VALIDATION LOGIC ===\n');
fprintf('For each empirical trial:\n');
fprintf('  1. Characterise noise profile (α, σ, fs)\n');
fprintf('  2. Query simulation → 15 predicted Δβ values\n');
fprintf('  3. Apply 6 pipelines → 15 observed Δβ values\n');
fprintf('  4. CCC compares these 15-element vectors\n');
fprintf('\nβ_true cancels in pairwise differences:\n');
fprintf('  Δβ_PQ = (β_true + bias_P) - (β_true + bias_Q) = bias_P - bias_Q\n');
fprintf('\nSo we validate the framework without knowing ground truth!\n');
