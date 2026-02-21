%% CCC vs Pearson's r Demo: Why Correlation ≠ Agreement
% For power law simulation pre-registration context
%
% Demonstrates why Lin's Concordance Correlation Coefficient (1989) is
% essential for method comparison studies. Pearson's r can be perfect
% even when methods systematically disagree.
%
% USAGE:
%   CCC_demo_correlation_vs_agreement
%
% OUTPUT:
%   Figure saved to docs/CCC_demo_correlation_vs_agreement.png
%
% Dagmar Fraser, University of Birmingham
% 2026-01-21

function CCC_demo_correlation_vs_agreement()

% Add functions path (adjust if needed)
scriptDir = fileparts(mfilename('fullpath'));
addpath(fullfile(scriptDir, '..', 'src', 'functions'));

figure('Position', [100, 100, 1400, 900], 'Color', 'w');

%% Scenario 1: Perfect correlation, location shift (constant bias)
% e.g., Two methods estimating β, one systematically reads higher
rng(42);
n = 50;
beta_true = 0.33 + randn(n,1)*0.05;  % "True" β estimates
beta_method2 = beta_true + 0.05;      % Method 2 reads 0.05 higher (systematic bias)

subplot(2,3,1);
scatter(beta_true, beta_method2, 50, [0.2 0.4 0.8], 'filled', 'MarkerFaceAlpha', 0.7);
hold on;
plot([0.15 0.55], [0.15 0.55], 'k--', 'LineWidth', 1.5);  % Identity line
plot([0.15 0.55], [0.20 0.60], 'r-', 'LineWidth', 1.5);   % Actual relationship
hold off;
xlabel('Method 1: β estimate');
ylabel('Method 2: β estimate');
title('Location Shift (+0.05)');
axis equal; xlim([0.15 0.55]); ylim([0.15 0.60]);
legend('Data', '45° line (perfect agreement)', 'Actual relationship', 'Location', 'southeast');

r1 = corr(beta_true, beta_method2);
ccc1 = linCCC_v001(beta_true, beta_method2);
text(0.17, 0.57, sprintf('r = %.3f\nCCC = %.3f', r1, ccc1.ccc), 'FontSize', 11, 'FontWeight', 'bold');

%% Scenario 2: Perfect correlation, scale shift
% e.g., One protocol consistently estimates β 20% higher
beta_method2_scale = beta_true * 1.20;  % 20% scale inflation

subplot(2,3,2);
scatter(beta_true, beta_method2_scale, 50, [0.2 0.4 0.8], 'filled', 'MarkerFaceAlpha', 0.7);
hold on;
plot([0.15 0.55], [0.15 0.55], 'k--', 'LineWidth', 1.5);
plot([0.15 0.55], [0.15 0.55]*1.2, 'r-', 'LineWidth', 1.5);
hold off;
xlabel('Method 1: β estimate');
ylabel('Method 2: β estimate');
title('Scale Shift (×1.2)');
axis equal; xlim([0.15 0.55]); ylim([0.15 0.70]);
legend('Data', '45° line', 'Actual relationship', 'Location', 'southeast');

r2 = corr(beta_true, beta_method2_scale);
ccc2 = linCCC_v001(beta_true, beta_method2_scale);
text(0.17, 0.65, sprintf('r = %.3f\nCCC = %.3f', r2, ccc2.ccc), 'FontSize', 11, 'FontWeight', 'bold');

%% Scenario 3: Both location and scale shift
% Realistic: different protocols give different β ranges
beta_method2_both = beta_true * 1.15 + 0.03;

subplot(2,3,3);
scatter(beta_true, beta_method2_both, 50, [0.2 0.4 0.8], 'filled', 'MarkerFaceAlpha', 0.7);
hold on;
plot([0.15 0.55], [0.15 0.55], 'k--', 'LineWidth', 1.5);
xfit = [0.15 0.55];
plot(xfit, xfit*1.15 + 0.03, 'r-', 'LineWidth', 1.5);
hold off;
xlabel('Method 1: β estimate');
ylabel('Method 2: β estimate');
title('Location + Scale Shift');
axis equal; xlim([0.15 0.55]); ylim([0.15 0.70]);
legend('Data', '45° line', 'Actual relationship', 'Location', 'southeast');

r3 = corr(beta_true, beta_method2_both);
ccc3 = linCCC_v001(beta_true, beta_method2_both);
text(0.17, 0.65, sprintf('r = %.3f\nCCC = %.3f', r3, ccc3.ccc), 'FontSize', 11, 'FontWeight', 'bold');

%% Scenario 4: Good agreement (what we want!)
beta_method2_good = beta_true + randn(n,1)*0.01;  % Small random error, no bias

subplot(2,3,4);
scatter(beta_true, beta_method2_good, 50, [0.2 0.7 0.3], 'filled', 'MarkerFaceAlpha', 0.7);
hold on;
plot([0.15 0.55], [0.15 0.55], 'k--', 'LineWidth', 1.5);
hold off;
xlabel('Method 1: β estimate');
ylabel('Method 2: β estimate');
title('Good Agreement (no systematic bias)');
axis equal; xlim([0.15 0.55]); ylim([0.15 0.55]);

r4 = corr(beta_true, beta_method2_good);
ccc4 = linCCC_v001(beta_true, beta_method2_good);
text(0.17, 0.52, sprintf('r = %.3f\nCCC = %.3f', r4, ccc4.ccc), 'FontSize', 11, 'FontWeight', 'bold', 'Color', [0 0.5 0]);

%% Scenario 5: Power law context - Legacy vs Vetted protocols
% Simulating what Fraser et al. 2025 found: Butterworth compresses toward 1/3
rng(123);
beta_vetted = 0.25 + randn(n,1)*0.08;  % Vetted protocol: wider spread, true variation
beta_legacy = beta_vetted * 0.4 + 0.20; % Legacy: compressed toward ~0.33

subplot(2,3,5);
scatter(beta_vetted, beta_legacy, 50, [0.8 0.3 0.2], 'filled', 'MarkerFaceAlpha', 0.7);
hold on;
plot([0.0 0.5], [0.0 0.5], 'k--', 'LineWidth', 1.5);
yline(1/3, 'b:', 'LineWidth', 1.5);
xfit = [0.0 0.5];
plot(xfit, xfit*0.4 + 0.20, 'r-', 'LineWidth', 1.5);
hold off;
xlabel('Vetted protocol: β estimate');
ylabel('Legacy protocol: β estimate');
title('Legacy vs Vetted (compression toward β=1/3)');
xlim([0.0 0.5]); ylim([0.2 0.45]);
legend('Data', '45° line', 'β = 1/3', 'Compression', 'Location', 'southeast');

r5 = corr(beta_vetted, beta_legacy);
ccc5 = linCCC_v001(beta_vetted, beta_legacy);
text(0.02, 0.43, sprintf('r = %.3f\nCCC = %.3f', r5, ccc5.ccc), 'FontSize', 11, 'FontWeight', 'bold', 'Color', [0.7 0 0]);

%% Summary bar chart
subplot(2,3,6);
scenarios = {'Location', 'Scale', 'Both', 'Good', 'Legacy/Vetted'};
r_vals = [r1, r2, r3, r4, r5];
ccc_vals = [ccc1.ccc, ccc2.ccc, ccc3.ccc, ccc4.ccc, ccc5.ccc];

x = 1:5;
b = bar(x, [r_vals; ccc_vals]', 'grouped');
b(1).FaceColor = [0.7 0.7 0.7];
b(2).FaceColor = [0.2 0.5 0.8];
set(gca, 'XTickLabel', scenarios);
ylabel('Coefficient value');
title('Pearson r vs CCC: Summary');
legend({'Pearson r', 'CCC'}, 'Location', 'southwest');
ylim([0 1.15]);
yline(1, 'k:', 'LineWidth', 1);

% Add value labels
for i = 1:5
    text(i-0.15, r_vals(i)+0.03, sprintf('%.2f', r_vals(i)), 'FontSize', 9, 'HorizontalAlignment', 'center');
    text(i+0.15, ccc_vals(i)+0.03, sprintf('%.2f', ccc_vals(i)), 'FontSize', 9, 'HorizontalAlignment', 'center', 'Color', [0 0 0.7]);
end

sgtitle('Why Correlation ≠ Agreement: The Case for Lin''s CCC', 'FontSize', 14, 'FontWeight', 'bold');

% Save figure
outputPath = fullfile(scriptDir, 'CCC_demo_correlation_vs_agreement.png');
saveas(gcf, outputPath);
fprintf('Figure saved to %s\n', outputPath);

end
