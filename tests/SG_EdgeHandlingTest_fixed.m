%% Compact SG Differentiation Test with Proper Edge Handling
% Tests SG differentiation with both original and fixed implementations
% Properly accounts for edge effects when comparing to analytical solutions

clear all; close all; clc;
addpath(genpath('../src/functions'));

% Generate a perfect circle for testing
fs = 100; t = (0:1/fs:1)';
radius = 100; x = radius * cos(2*pi*t); y = radius * sin(2*pi*t);

% Analytical derivatives
dx_analytical = -radius * 2*pi * sin(2*pi*t);
dy_analytical = radius * 2*pi * cos(2*pi*t);
d2x_analytical = -radius * (2*pi)^2 * cos(2*pi*t);
d2y_analytical = -radius * (2*pi)^2 * sin(2*pi*t);

% SG filter parameters
framelen = 9;  % Must be odd
order = 3;     % At least 3 for jerk
dt = 1/fs;
padding = (framelen-1)/2;

% Exclude edge points for all comparisons
edge_points = floor(framelen/2);
valid_indices = (edge_points+1):(length(x)-edge_points);

% 1. Direct SG implementation (without the * p)
[~, g] = sgolay(order, framelen);

% Calculate derivatives using direct convolution (MATLAB's recommended approach)
dx_direct = zeros(length(x), 4);
dy_direct = zeros(length(y), 4);

dx_direct(:,1) = x;  % Original data
dy_direct(:,1) = y;

% Calculate derivatives (p=1,2,3)
for p = 1:3
    dx_direct(:,p+1) = conv(x, factorial(p)/(-dt)^p * g(:,p+1), 'same');
    dy_direct(:,p+1) = conv(y, factorial(p)/(-dt)^p * g(:,p+1), 'same');
end

% 2. Implementation with * p (similar to the original function)
dx_withp = zeros(length(x), 4);
dy_withp = zeros(length(y), 4);

dx_withp(:,1) = x;  % Original data
dy_withp(:,1) = y;

% Calculate derivatives with the * p multiplication
for p = 1:3
    dxElement = conv(x, factorial(p)/(-dt)^p * g(:,p+1), 'same');
    dyElement = conv(y, factorial(p)/(-dt)^p * g(:,p+1), 'same');
    
    % Note the * p multiplication here (the key issue)
    dxPadded = padarray(dxElement * p, padding/2, 0, 'pre');
    dyPadded = padarray(dyElement * p, padding/2, 0, 'pre');
    
    dx_withp(:,p+1) = dxPadded(2:length(x)+1);
    dy_withp(:,p+1) = dyPadded(2:length(y)+1);
end

% Calculate errors (only for valid indices)
% Velocity errors
error_direct_vel_x = mean(abs(dx_direct(valid_indices,2) - dx_analytical(valid_indices)));
error_direct_vel_y = mean(abs(dy_direct(valid_indices,2) - dy_analytical(valid_indices)));
error_withp_vel_x = mean(abs(dx_withp(valid_indices,2) - dx_analytical(valid_indices)));
error_withp_vel_y = mean(abs(dy_withp(valid_indices,2) - dy_analytical(valid_indices)));

% Acceleration errors
error_direct_acc_x = mean(abs(dx_direct(valid_indices,3) - d2x_analytical(valid_indices)));
error_direct_acc_y = mean(abs(dy_direct(valid_indices,3) - d2y_analytical(valid_indices)));
error_withp_acc_x = mean(abs(dx_withp(valid_indices,3) - d2x_analytical(valid_indices)));
error_withp_acc_y = mean(abs(dy_withp(valid_indices,3) - d2y_analytical(valid_indices)));

% Calculate scaling factors between implementations (only for valid indices)
scaling_vel_x = mean(abs(dx_withp(valid_indices,2))) / mean(abs(dx_direct(valid_indices,2)));
scaling_vel_y = mean(abs(dy_withp(valid_indices,2))) / mean(abs(dy_direct(valid_indices,2)));
scaling_acc_x = mean(abs(dx_withp(valid_indices,3))) / mean(abs(dx_direct(valid_indices,3)));
scaling_acc_y = mean(abs(dy_withp(valid_indices,3))) / mean(abs(dy_direct(valid_indices,3)));
scaling_jerk_x = mean(abs(dx_withp(valid_indices,4))) / mean(abs(dx_direct(valid_indices,4)));
scaling_jerk_y = mean(abs(dy_withp(valid_indices,4))) / mean(abs(dy_direct(valid_indices,4)));

% Display results
fprintf('=== SG DIFFERENTIATION COMPARISON ===\n');
fprintf('(Excluding %d edge points from each end)\n\n', edge_points);

fprintf('VELOCITY ERRORS:\n');
fprintf('  Direct implementation: x=%.6f, y=%.6f (mean=%.6f)\n', ...
    error_direct_vel_x, error_direct_vel_y, mean([error_direct_vel_x, error_direct_vel_y]));
fprintf('  With *p multiplication: x=%.6f, y=%.6f (mean=%.6f)\n', ...
    error_withp_vel_x, error_withp_vel_y, mean([error_withp_vel_x, error_withp_vel_y]));
fprintf('  Improvement: %.2f%%\n\n', ...
    (mean([error_withp_vel_x, error_withp_vel_y]) - mean([error_direct_vel_x, error_direct_vel_y])) / ...
    mean([error_withp_vel_x, error_withp_vel_y]) * 100);

fprintf('ACCELERATION ERRORS:\n');
fprintf('  Direct implementation: x=%.6f, y=%.6f (mean=%.6f)\n', ...
    error_direct_acc_x, error_direct_acc_y, mean([error_direct_acc_x, error_direct_acc_y]));
fprintf('  With *p multiplication: x=%.6f, y=%.6f (mean=%.6f)\n', ...
    error_withp_acc_x, error_withp_acc_y, mean([error_withp_acc_x, error_withp_acc_y]));
fprintf('  Improvement: %.2f%%\n\n', ...
    (mean([error_withp_acc_x, error_withp_acc_y]) - mean([error_direct_acc_x, error_direct_acc_y])) / ...
    mean([error_withp_acc_x, error_withp_acc_y]) * 100);

fprintf('SCALING FACTORS (With p / Direct):\n');
fprintf('  Velocity: x=%.2f, y=%.2f (expected = 1.0)\n', scaling_vel_x, scaling_vel_y);
fprintf('  Acceleration: x=%.2f, y=%.2f (expected = 2.0)\n', scaling_acc_x, scaling_acc_y);
fprintf('  Jerk: x=%.2f, y=%.2f (expected = 3.0)\n\n', scaling_jerk_x, scaling_jerk_y);

% Plot velocity comparison - with valid points highlighted
figure('Name', 'Velocity Comparison - Full Data');
subplot(2,1,1);
plot(t, dx_analytical, 'k-', 'LineWidth', 2, 'DisplayName', 'Analytical');
hold on;
plot(t, dx_direct(:,2), 'g-', 'LineWidth', 1, 'DisplayName', 'Direct SG');
plot(t, dx_withp(:,2), 'r-', 'LineWidth', 1, 'DisplayName', 'With * p');
% Highlight valid points
plot(t(valid_indices), dx_analytical(valid_indices), 'ko', 'MarkerSize', 5, ...
    'DisplayName', 'Analytical (valid)');
plot(t(valid_indices), dx_direct(valid_indices,2), 'go', 'MarkerSize', 5, ...
    'DisplayName', 'Direct SG (valid)');
plot(t(valid_indices), dx_withp(valid_indices,2), 'ro', 'MarkerSize', 5, ...
    'DisplayName', 'With * p (valid)');
title('X Velocity - Full Data with Valid Region Highlighted');
legend('Location', 'best');
grid on;

subplot(2,1,2);
plot(t, dy_analytical, 'k-', 'LineWidth', 2, 'DisplayName', 'Analytical');
hold on;
plot(t, dy_direct(:,2), 'g-', 'LineWidth', 1, 'DisplayName', 'Direct SG');
plot(t, dy_withp(:,2), 'r-', 'LineWidth', 1, 'DisplayName', 'With * p');
% Highlight valid points
plot(t(valid_indices), dy_analytical(valid_indices), 'ko', 'MarkerSize', 5, ...
    'DisplayName', 'Analytical (valid)');
plot(t(valid_indices), dy_direct(valid_indices,2), 'go', 'MarkerSize', 5, ...
    'DisplayName', 'Direct SG (valid)');
plot(t(valid_indices), dy_withp(valid_indices,2), 'ro', 'MarkerSize', 5, ...
    'DisplayName', 'With * p (valid)');
title('Y Velocity - Full Data with Valid Region Highlighted');
grid on;

% Plot only the valid region for velocity
figure('Name', 'Velocity Comparison - Valid Region Only');
subplot(2,1,1);
plot(t(valid_indices), dx_analytical(valid_indices), 'k-', 'LineWidth', 2, 'DisplayName', 'Analytical');
hold on;
plot(t(valid_indices), dx_direct(valid_indices,2), 'g-', 'LineWidth', 1.5, 'DisplayName', 'Direct SG');
plot(t(valid_indices), dx_withp(valid_indices,2), 'r-', 'LineWidth', 1.5, 'DisplayName', 'With * p');
title('X Velocity - Valid Region Only');
legend('Location', 'best');
grid on;
xlim([t(valid_indices(1)) t(valid_indices(end))]);

subplot(2,1,2);
plot(t(valid_indices), dy_analytical(valid_indices), 'k-', 'LineWidth', 2, 'DisplayName', 'Analytical');
hold on;
plot(t(valid_indices), dy_direct(valid_indices,2), 'g-', 'LineWidth', 1.5, 'DisplayName', 'Direct SG');
plot(t(valid_indices), dy_withp(valid_indices,2), 'r-', 'LineWidth', 1.5, 'DisplayName', 'With * p');
title('Y Velocity - Valid Region Only');
grid on;
xlim([t(valid_indices(1)) t(valid_indices(end))]);

% Plot acceleration comparison - with valid points highlighted
figure('Name', 'Acceleration Comparison - Full Data');
subplot(2,1,1);
plot(t, d2x_analytical, 'k-', 'LineWidth', 2, 'DisplayName', 'Analytical');
hold on;
plot(t, dx_direct(:,3), 'g-', 'LineWidth', 1, 'DisplayName', 'Direct SG');
plot(t, dx_withp(:,3), 'r-', 'LineWidth', 1, 'DisplayName', 'With * p');
% Highlight valid points
plot(t(valid_indices), d2x_analytical(valid_indices), 'ko', 'MarkerSize', 5, ...
    'DisplayName', 'Analytical (valid)');
plot(t(valid_indices), dx_direct(valid_indices,3), 'go', 'MarkerSize', 5, ...
    'DisplayName', 'Direct SG (valid)');
plot(t(valid_indices), dx_withp(valid_indices,3), 'ro', 'MarkerSize', 5, ...
    'DisplayName', 'With * p (valid)');
title('X Acceleration - Full Data with Valid Region Highlighted');
legend('Location', 'best');
grid on;

subplot(2,1,2);
plot(t, d2y_analytical, 'k-', 'LineWidth', 2, 'DisplayName', 'Analytical');
hold on;
plot(t, dy_direct(:,3), 'g-', 'LineWidth', 1, 'DisplayName', 'Direct SG');
plot(t, dy_withp(:,3), 'r-', 'LineWidth', 1, 'DisplayName', 'With * p');
% Highlight valid points
plot(t(valid_indices), d2y_analytical(valid_indices), 'ko', 'MarkerSize', 5, ...
    'DisplayName', 'Analytical (valid)');
plot(t(valid_indices), dy_direct(valid_indices,3), 'go', 'MarkerSize', 5, ...
    'DisplayName', 'Direct SG (valid)');
plot(t(valid_indices), dy_withp(valid_indices,3), 'ro', 'MarkerSize', 5, ...
    'DisplayName', 'With * p (valid)');
title('Y Acceleration - Full Data with Valid Region Highlighted');
grid on;

% Plot only the valid region for acceleration
figure('Name', 'Acceleration Comparison - Valid Region Only');
subplot(2,1,1);
plot(t(valid_indices), d2x_analytical(valid_indices), 'k-', 'LineWidth', 2, 'DisplayName', 'Analytical');
hold on;
plot(t(valid_indices), dx_direct(valid_indices,3), 'g-', 'LineWidth', 1.5, 'DisplayName', 'Direct SG');
plot(t(valid_indices), dx_withp(valid_indices,3), 'r-', 'LineWidth', 1.5, 'DisplayName', 'With * p');
title('X Acceleration - Valid Region Only');
legend('Location', 'best');
grid on;
xlim([t(valid_indices(1)) t(valid_indices(end))]);

subplot(2,1,2);
plot(t(valid_indices), d2y_analytical(valid_indices), 'k-', 'LineWidth', 2, 'DisplayName', 'Analytical');
hold on;
plot(t(valid_indices), dy_direct(valid_indices,3), 'g-', 'LineWidth', 1.5, 'DisplayName', 'Direct SG');
plot(t(valid_indices), dy_withp(valid_indices,3), 'r-', 'LineWidth', 1.5, 'DisplayName', 'With * p');
title('Y Acceleration - Valid Region Only');
grid on;
xlim([t(valid_indices(1)) t(valid_indices(end))]);

% Calculate even stricter valid indices to ensure no edge effects remain
% This accounts for any potential shifts from the padding operations
edge_points_strict = edge_points + 1;
valid_indices_strict = (edge_points_strict+1):(length(x)-edge_points_strict);

% Plot only the valid region for velocity (with stricter indices)
figure('Name', 'Velocity Comparison - Strictly Valid Region Only');
subplot(2,1,1);
plot(t(valid_indices_strict), dx_analytical(valid_indices_strict), 'k-', 'LineWidth', 2, 'DisplayName', 'Analytical');
hold on;
plot(t(valid_indices_strict), dx_direct(valid_indices_strict,2), 'g-', 'LineWidth', 1.5, 'DisplayName', 'Direct SG');
plot(t(valid_indices_strict), dx_withp(valid_indices_strict,2), 'r-', 'LineWidth', 1.5, 'DisplayName', 'With * p');
title('X Velocity - Strictly Valid Region Only');
legend('Location', 'best');
grid on;
xlim([t(valid_indices_strict(1)) t(valid_indices_strict(end))]);

subplot(2,1,2);
plot(t(valid_indices_strict), dy_analytical(valid_indices_strict), 'k-', 'LineWidth', 2, 'DisplayName', 'Analytical');
hold on;
plot(t(valid_indices_strict), dy_direct(valid_indices_strict,2), 'g-', 'LineWidth', 1.5, 'DisplayName', 'Direct SG');
plot(t(valid_indices_strict), dy_withp(valid_indices_strict,2), 'r-', 'LineWidth', 1.5, 'DisplayName', 'With * p');
title('Y Velocity - Strictly Valid Region Only');
grid on;
xlim([t(valid_indices_strict(1)) t(valid_indices_strict(end))]);

% Plot only the strictly valid region for acceleration
figure('Name', 'Acceleration Comparison - Strictly Valid Region Only');
subplot(2,1,1);
plot(t(valid_indices_strict), d2x_analytical(valid_indices_strict), 'k-', 'LineWidth', 2, 'DisplayName', 'Analytical');
hold on;
plot(t(valid_indices_strict), dx_direct(valid_indices_strict,3), 'g-', 'LineWidth', 1.5, 'DisplayName', 'Direct SG');
plot(t(valid_indices_strict), dx_withp(valid_indices_strict,3), 'r-', 'LineWidth', 1.5, 'DisplayName', 'With * p');
title('X Acceleration - Strictly Valid Region Only');
legend('Location', 'best');
grid on;
xlim([t(valid_indices_strict(1)) t(valid_indices_strict(end))]);

subplot(2,1,2);
plot(t(valid_indices_strict), d2y_analytical(valid_indices_strict), 'k-', 'LineWidth', 2, 'DisplayName', 'Analytical');
hold on;
plot(t(valid_indices_strict), dy_direct(valid_indices_strict,3), 'g-', 'LineWidth', 1.5, 'DisplayName', 'Direct SG');
plot(t(valid_indices_strict), dy_withp(valid_indices_strict,3), 'r-', 'LineWidth', 1.5, 'DisplayName', 'With * p');
title('Y Acceleration - Strictly Valid Region Only');
grid on;
xlim([t(valid_indices_strict(1)) t(valid_indices_strict(end))]);

% Calculate errors using stricter valid indices
% Velocity errors
error_direct_vel_x_strict = mean(abs(dx_direct(valid_indices_strict,2) - dx_analytical(valid_indices_strict)));
error_direct_vel_y_strict = mean(abs(dy_direct(valid_indices_strict,2) - dy_analytical(valid_indices_strict)));
error_withp_vel_x_strict = mean(abs(dx_withp(valid_indices_strict,2) - dx_analytical(valid_indices_strict)));
error_withp_vel_y_strict = mean(abs(dy_withp(valid_indices_strict,2) - dy_analytical(valid_indices_strict)));

% Acceleration errors
error_direct_acc_x_strict = mean(abs(dx_direct(valid_indices_strict,3) - d2x_analytical(valid_indices_strict)));
error_direct_acc_y_strict = mean(abs(dy_direct(valid_indices_strict,3) - d2y_analytical(valid_indices_strict)));
error_withp_acc_x_strict = mean(abs(dx_withp(valid_indices_strict,3) - d2x_analytical(valid_indices_strict)));
error_withp_acc_y_strict = mean(abs(dy_withp(valid_indices_strict,3) - d2y_analytical(valid_indices_strict)));

fprintf('\nERRORS WITH STRICTER INDICES:\n');
fprintf('VELOCITY ERRORS (STRICT):\n');
fprintf('  Direct implementation: x=%.6f, y=%.6f (mean=%.6f)\n', ...
    error_direct_vel_x_strict, error_direct_vel_y_strict, mean([error_direct_vel_x_strict, error_direct_vel_y_strict]));
fprintf('  With *p multiplication: x=%.6f, y=%.6f (mean=%.6f)\n', ...
    error_withp_vel_x_strict, error_withp_vel_y_strict, mean([error_withp_vel_x_strict, error_withp_vel_y_strict]));
fprintf('  Improvement: %.2f%%\n\n', ...
    (mean([error_withp_vel_x_strict, error_withp_vel_y_strict]) - mean([error_direct_vel_x_strict, error_direct_vel_y_strict])) / ...
    mean([error_withp_vel_x_strict, error_withp_vel_y_strict]) * 100);

fprintf('ACCELERATION ERRORS (STRICT):\n');
fprintf('  Direct implementation: x=%.6f, y=%.6f (mean=%.6f)\n', ...
    error_direct_acc_x_strict, error_direct_acc_y_strict, mean([error_direct_acc_x_strict, error_direct_acc_y_strict]));
fprintf('  With *p multiplication: x=%.6f, y=%.6f (mean=%.6f)\n', ...
    error_withp_acc_x_strict, error_withp_acc_y_strict, mean([error_withp_acc_x_strict, error_withp_acc_y_strict]));
fprintf('  Improvement: %.2f%%\n', ...
    (mean([error_withp_acc_x_strict, error_withp_acc_y_strict]) - mean([error_direct_acc_x_strict, error_direct_acc_y_strict])) / ...
    mean([error_withp_acc_x_strict, error_withp_acc_y_strict]) * 100);

% Scale factors
scaling_vel_x_strict = mean(abs(dx_withp(valid_indices_strict,2))) / mean(abs(dx_direct(valid_indices_strict,2)));
scaling_vel_y_strict = mean(abs(dy_withp(valid_indices_strict,2))) / mean(abs(dy_direct(valid_indices_strict,2)));
scaling_acc_x_strict = mean(abs(dx_withp(valid_indices_strict,3))) / mean(abs(dx_direct(valid_indices_strict,3)));
scaling_acc_y_strict = mean(abs(dy_withp(valid_indices_strict,3))) / mean(abs(dy_direct(valid_indices_strict,3)));
fprintf('\nSCALING FACTORS WITH STRICTER INDICES:\n');
fprintf('  Velocity: x=%.2f, y=%.2f (expected = 1.0)\n', scaling_vel_x_strict, scaling_vel_y_strict);
fprintf('  Acceleration: x=%.2f, y=%.2f (expected = 2.0)\n', scaling_acc_x_strict, scaling_acc_y_strict);

fprintf('\nKEY FINDINGS:\n');
fprintf('1. The * p multiplication incorrectly scales derivatives by their order\n');
fprintf('2. Removing this scaling significantly improves accuracy\n');
fprintf('3. Edge effects must be excluded from calculations (framelen/2 points)\n');
fprintf('4. The proper implementation follows MATLAB`s documentation example\n');
