%% Simplified Test for Savitzky-Golay Implementation
% This script focuses only on testing the SG filter implementation

clear all;
close all;
clc;

% Generate test signal - simple sine wave
t = linspace(0,2*pi,101)';
x = sin(t);
y = cos(t);

% Analytical derivatives
dx_analytical = cos(t);
dy_analytical = -sin(t);
d2x_analytical = -sin(t);
d2y_analytical = -cos(t);
d3x_analytical = -cos(t);
d3y_analytical = sin(t);

% Setup
fs = 100;  % Hz
dt = 1/fs;

% Different implementations of SG filtering
disp('Testing different SG implementations...');

% Parameters
order = 3;
framelen = 9;
padding = (framelen-1)/2;

% 1. Direct implementation (based on MATLAB documentation)
[b, g] = sgolay(order, framelen);
dx_direct = zeros(length(x), 4);
dy_direct = zeros(length(y), 4);

% Position (smoothed signal)
dx_direct(:,1) = conv(x, g(:,1), 'same');
dy_direct(:,1) = conv(y, g(:,1), 'same');

% Velocity (1st derivative)
dx_direct(:,2) = conv(x, factorial(1)/(-dt)^1 * g(:,1+1), 'same');
dy_direct(:,2) = conv(y, factorial(1)/(-dt)^1 * g(:,1+1), 'same');

% Acceleration (2nd derivative)
dx_direct(:,3) = conv(x, factorial(2)/(-dt)^2 * g(:,2+1), 'same');
dy_direct(:,3) = conv(y, factorial(2)/(-dt)^2 * g(:,2+1), 'same');

% Jerk (3rd derivative)
dx_direct(:,4) = conv(x, factorial(3)/(-dt)^3 * g(:,3+1), 'same');
dy_direct(:,4) = conv(y, factorial(3)/(-dt)^3 * g(:,3+1), 'same');

% 2. Implementation with padding but no * p
dx_pad = zeros(length(x), 4);
dy_pad = zeros(length(y), 4);

% Position (original signal)
dx_pad(:,1) = x;
dy_pad(:,1) = y;

for p = 1:3
    % Calculate derivative
    dxElement = conv(x, factorial(p)/(-dt)^p * g(:,p+1), 'same');
    dyElement = conv(y, factorial(p)/(-dt)^p * g(:,p+1), 'same');
    
    % Pad array (without * p)
    dxPadded = padarray(dxElement, padding/2, 0, 'pre');
    dyPadded = padarray(dyElement, padding/2, 0, 'pre');
    
    % Extract subset
    dx_pad(:,p+1) = dxPadded(2:length(x)+1);
    dy_pad(:,p+1) = dyPadded(2:length(y)+1);
end

% 3. Original implementation with * p (for comparison)
dx_orig = zeros(length(x), 4);
dy_orig = zeros(length(y), 4);

% Position (original signal)
dx_orig(:,1) = x;
dy_orig(:,1) = y;

for p = 1:3
    % Calculate derivative
    dxElement = conv(x, factorial(p)/(-dt)^p * g(:,p+1), 'same');
    dyElement = conv(y, factorial(p)/(-dt)^p * g(:,p+1), 'same');
    
    % Pad array (with * p)
    dxPadded = padarray(dxElement * p, padding/2, 0, 'pre');
    dyPadded = padarray(dyElement * p, padding/2, 0, 'pre');
    
    % Extract subset
    dx_orig(:,p+1) = dxPadded(2:length(x)+1);
    dy_orig(:,p+1) = dyPadded(2:length(y)+1);
end

% Calculate errors
% Direct SG
error_dx_direct_vel = mean(abs(dx_direct(:,2) - dx_analytical));
error_dy_direct_vel = mean(abs(dy_direct(:,2) - dy_analytical));
error_dx_direct_acc = mean(abs(dx_direct(:,3) - d2x_analytical));
error_dy_direct_acc = mean(abs(dy_direct(:,3) - d2y_analytical));
error_dx_direct_jerk = mean(abs(dx_direct(:,4) - d3x_analytical));
error_dy_direct_jerk = mean(abs(dy_direct(:,4) - d3y_analytical));

% Padded without * p
error_dx_pad_vel = mean(abs(dx_pad(:,2) - dx_analytical));
error_dy_pad_vel = mean(abs(dy_pad(:,2) - dy_analytical));
error_dx_pad_acc = mean(abs(dx_pad(:,3) - d2x_analytical));
error_dy_pad_acc = mean(abs(dy_pad(:,3) - d2y_analytical));
error_dx_pad_jerk = mean(abs(dx_pad(:,4) - d3x_analytical));
error_dy_pad_jerk = mean(abs(dy_pad(:,4) - d3y_analytical));

% Original with * p
error_dx_orig_vel = mean(abs(dx_orig(:,2) - dx_analytical));
error_dy_orig_vel = mean(abs(dy_orig(:,2) - dy_analytical));
error_dx_orig_acc = mean(abs(dx_orig(:,3) - d2x_analytical));
error_dy_orig_acc = mean(abs(dy_orig(:,3) - d2y_analytical));
error_dx_orig_jerk = mean(abs(dx_orig(:,4) - d3x_analytical));
error_dy_orig_jerk = mean(abs(dy_orig(:,4) - d3y_analytical));

% Display results
disp('-------------------------------------------------');
disp('COMPARISON OF SG IMPLEMENTATIONS');
disp('-------------------------------------------------');

disp('1. Direct SG (MATLAB documentation style):');
fprintf('  Velocity error: %.6f\n', mean([error_dx_direct_vel, error_dy_direct_vel]));
fprintf('  Acceleration error: %.6f\n', mean([error_dx_direct_acc, error_dy_direct_acc]));
fprintf('  Jerk error: %.6f\n', mean([error_dx_direct_jerk, error_dy_direct_jerk]));

disp('2. SG with padding (without * p):');
fprintf('  Velocity error: %.6f\n', mean([error_dx_pad_vel, error_dy_pad_vel]));
fprintf('  Acceleration error: %.6f\n', mean([error_dx_pad_acc, error_dy_pad_acc]));
fprintf('  Jerk error: %.6f\n', mean([error_dx_pad_jerk, error_dy_pad_jerk]));

disp('3. SG with padding (with * p):');
fprintf('  Velocity error: %.6f\n', mean([error_dx_orig_vel, error_dy_orig_vel]));
fprintf('  Acceleration error: %.6f\n', mean([error_dx_orig_acc, error_dy_orig_acc]));
fprintf('  Jerk error: %.6f\n', mean([error_dx_orig_jerk, error_dy_orig_jerk]));

% Plot results
figure('Name', 'Velocity Comparison');
plot(t, dx_analytical, 'k-', 'LineWidth', 2);
hold on;
plot(t, dx_direct(:,2), 'b-', 'LineWidth', 1);
plot(t, dx_pad(:,2), 'g--', 'LineWidth', 1);
plot(t, dx_orig(:,2), 'r-.', 'LineWidth', 1);
title('Velocity (1st Derivative)');
legend('Analytical', 'Direct SG', 'Padded (no *p)', 'Padded (with *p)', 'Location', 'best');
grid on;

figure('Name', 'Acceleration Comparison');
plot(t, d2x_analytical, 'k-', 'LineWidth', 2);
hold on;
plot(t, dx_direct(:,3), 'b-', 'LineWidth', 1);
plot(t, dx_pad(:,3), 'g--', 'LineWidth', 1);
plot(t, dx_orig(:,3), 'r-.', 'LineWidth', 1);
title('Acceleration (2nd Derivative)');
legend('Analytical', 'Direct SG', 'Padded (no *p)', 'Padded (with *p)', 'Location', 'best');
grid on;

% Check padding behavior in detail
disp('-------------------------------------------------');
disp('DETAILED EXAMINATION OF PADDING');
disp('-------------------------------------------------');

% Create a very small sine wave
ts = linspace(0, 2*pi, 11)';
xs = sin(ts);
dxs_analytical = cos(ts);

% Parameters
order_s = 3;
framelen_s = 9;  % 9 is large compared to 11 data points
padding_s = (framelen_s-1)/2;

% Get SG coefficients
[b_s, g_s] = sgolay(order_s, framelen_s);

% 1. Direct SG
dxs_direct = conv(xs, factorial(1)/(-dt)^1 * g_s(:,1+1), 'same');

% 2. Pad without * p
dxsElement = conv(xs, factorial(1)/(-dt)^1 * g_s(:,1+1), 'same');
dxsPadded = padarray(dxsElement, padding_s/2, 0, 'pre');
dxs_pad = dxsPadded(2:length(xs)+1);

% 3. Both methods with detailed output
disp('Original signal:');
disp(xs);
disp('Analytical derivative:');
disp(dxs_analytical);
disp('Direct SG derivative:');
disp(dxs_direct);
disp('Padded (no *p) derivative:');
disp(dxs_pad);

% Check if padding/2 is actually an integer
disp(['Is padding/2 an integer? ' num2str(padding_s/2) ' - ' num2str(mod(padding_s/2, 1) == 0)]);

% Check the size of the data vs frame length
disp(['Data size: ' num2str(length(xs)) ', Frame length: ' num2str(framelen_s)]);
disp(['Is frame length > data size? ' num2str(framelen_s > length(xs))]);

% Examine the end effects more carefully
disp('-------------------------------------------------');
disp('EXAMINATION OF EDGE EFFECTS');
disp('-------------------------------------------------');

% Generate a new signal with 101 points for better visualization
te = linspace(0, 2*pi, 101)';
xe = sin(te);
dxe_analytical = cos(te);

% Calculate direct SG without padding
dxe_direct = conv(xe, factorial(1)/(-dt)^1 * g(:,1+1), 'same');

% Calculate errors at edges vs center
edge_width = round(framelen/2);
center_indices = (edge_width+1):(length(xe)-edge_width);
edge_indices = [1:edge_width, (length(xe)-edge_width+1):length(xe)];

dxe_direct_center_error = mean(abs(dxe_direct(center_indices) - dxe_analytical(center_indices)));
dxe_direct_edge_error = mean(abs(dxe_direct(edge_indices) - dxe_analytical(edge_indices)));

% Display edge vs center errors
fprintf('Direct SG center error: %.6f\n', dxe_direct_center_error);
fprintf('Direct SG edge error: %.6f\n', dxe_direct_edge_error);
fprintf('Edge/center error ratio: %.2f\n', dxe_direct_edge_error/dxe_direct_center_error);

% Plot to visualize
figure('Name', 'Edge Effects');
plot(te, dxe_analytical, 'k-', 'LineWidth', 2);
hold on;
plot(te, dxe_direct, 'b-', 'LineWidth', 1);
plot(te(edge_indices), dxe_direct(edge_indices), 'ro', 'LineWidth', 1);
title('Edge Effects in SG Filtering');
legend('Analytical', 'SG Derivative', 'Edge Points', 'Location', 'best');
grid on;

disp('-------------------------------------------------');
disp('CONCLUSION');
disp('-------------------------------------------------');
disp('Based on this analysis, we can determine:');
disp('1. Which implementation performs best for velocity, acceleration, and jerk');
disp('2. The impact of edge effects on SG filtering accuracy');
disp('3. The correct way to handle padding in this context');
disp('-------------------------------------------------');
