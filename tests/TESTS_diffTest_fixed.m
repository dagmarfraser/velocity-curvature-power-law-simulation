%% Test Script for differentiateKinematicsEBR function
% This script tests all methods of the differentiateKinematicsEBR function
% to verify they produce expected results for known trajectories
%
% Created: May 2025
% Author: Dagmar Scott Fraser
% d.s.fraser@bham.ac.uk

% Clear workspace
clear all;
close all;
clc;

% Add necessary paths
addpath(genpath('../src/functions'));
addpath(genpath('../src/req'));
addpath(genpath('../src/utils'));

disp('Running tests for differentiateKinematicsEBR function...');


%% Test Setup
% We'll test with a circular trajectory where the analytical derivatives are known

% Set up parameters
radius = 100; % mm
frequency = 1; % Hz, one full circle per second
fs = 100; % Hz, sampling frequency
duration = 1; % seconds
phi = 0; % initial phase

% Generate time vector
t = linspace(0, duration, fs*duration)';

% Generate circle trajectory
x = radius * cos(2*pi*frequency*t + phi);
y = radius * sin(2*pi*frequency*t + phi);

% Analytical derivatives (known values for this case)
dx_analytical = -radius * 2*pi*frequency * sin(2*pi*frequency*t + phi);
dy_analytical = radius * 2*pi*frequency * cos(2*pi*frequency*t + phi);
d2x_analytical = -radius * (2*pi*frequency)^2 * cos(2*pi*frequency*t + phi);
d2y_analytical = -radius * (2*pi*frequency)^2 * sin(2*pi*frequency*t + phi);
d3x_analytical = radius * (2*pi*frequency)^3 * sin(2*pi*frequency*t + phi);
d3y_analytical = -radius * (2*pi*frequency)^3 * cos(2*pi*frequency*t + phi);

%% Test Type 1: Finite Differences Differentiation
disp('Testing Method 1: Finite Differences Differentiation');
filterType = 1;
filterParams = [0, 0, 0]; % Not used for type 1
[dx1, dy1] = differentiateKinematicsEBR(x, y, filterType, filterParams, fs);

% Calculate errors
error_dx1_vel = mean(abs(dx1(2:end,2) - dx_analytical(1:end-1)));
error_dy1_vel = mean(abs(dy1(2:end,2) - dy_analytical(1:end-1)));
error_dx1_acc = mean(abs(dx1(2:end-1,3) - d2x_analytical(1:end-2)));
error_dy1_acc = mean(abs(dy1(2:end-1,3) - d2y_analytical(1:end-2)));
error_dx1_jerk = mean(abs(dx1(2:end-2,4) - d3x_analytical(1:end-3)));
error_dy1_jerk = mean(abs(dy1(2:end-2,4) - d3y_analytical(1:end-3)));

% Report
fprintf('Method 1 Mean Absolute Errors:\n');
fprintf('  Velocity - x: %.6f, y: %.6f\n', error_dx1_vel, error_dy1_vel);
fprintf('  Acceleration - x: %.6f, y: %.6f\n', error_dx1_acc, error_dy1_acc);
fprintf('  Jerk - x: %.6f, y: %.6f\n', error_dx1_jerk, error_dy1_jerk);

%% Test Type 2: Low-pass filter then Finite Differences
disp('Testing Method 2: Low-pass filter then Finite Differences');
filterType = 2;
filterParams = [2, 10, 1]; % 2nd order, 10 Hz cutoff, zero lag
[dx2, dy2] = differentiateKinematicsEBR(x, y, filterType, filterParams, fs);

% Calculate errors
error_dx2_vel = mean(abs(dx2(2:end,2) - dx_analytical(1:end-1)));
error_dy2_vel = mean(abs(dy2(2:end,2) - dy_analytical(1:end-1)));
error_dx2_acc = mean(abs(dx2(2:end-1,3) - d2x_analytical(1:end-2)));
error_dy2_acc = mean(abs(dy2(2:end-1,3) - d2y_analytical(1:end-2)));
error_dx2_jerk = mean(abs(dx2(2:end-2,4) - d3x_analytical(1:end-3)));
error_dy2_jerk = mean(abs(dy2(2:end-2,4) - d3y_analytical(1:end-3)));

% Report
fprintf('Method 2 Mean Absolute Errors:\n');
fprintf('  Velocity - x: %.6f, y: %.6f\n', error_dx2_vel, error_dy2_vel);
fprintf('  Acceleration - x: %.6f, y: %.6f\n', error_dx2_acc, error_dy2_acc);
fprintf('  Jerk - x: %.6f, y: %.6f\n', error_dx2_jerk, error_dy2_jerk);

%% Test Type 3: Finite Differences then Low-pass filter
disp('Testing Method 3: Finite Differences then Low-pass filter');
filterType = 3;
filterParams = [2, 10, 1]; % 2nd order, 10 Hz cutoff, zero lag
[dx3, dy3] = differentiateKinematicsEBR(x, y, filterType, filterParams, fs);

% Calculate errors
error_dx3_vel = mean(abs(dx3(2:end,2) - dx_analytical(1:end-1)));
error_dy3_vel = mean(abs(dy3(2:end,2) - dy_analytical(1:end-1)));
error_dx3_acc = mean(abs(dx3(2:end-1,3) - d2x_analytical(1:end-2)));
error_dy3_acc = mean(abs(dy3(2:end-1,3) - d2y_analytical(1:end-2)));
error_dx3_jerk = mean(abs(dx3(2:end-2,4) - d3x_analytical(1:end-3)));
error_dy3_jerk = mean(abs(dy3(2:end-2,4) - d3y_analytical(1:end-3)));

% Report
fprintf('Method 3 Mean Absolute Errors:\n');
fprintf('  Velocity - x: %.6f, y: %.6f\n', error_dx3_vel, error_dy3_vel);
fprintf('  Acceleration - x: %.6f, y: %.6f\n', error_dx3_acc, error_dy3_acc);
fprintf('  Jerk - x: %.6f, y: %.6f\n', error_dx3_jerk, error_dy3_jerk);

%% Test Type 4: Savitzky-Golay smoothing differential filter
disp('Testing Method 4: Savitzky-Golay smoothing differential filter');
filterType = 4;
order = 3;  % Use 3rd order (can calculate up to 3rd derivative)
framelen = 9;
filterParams = [order, framelen];
[dx4, dy4] = differentiateKinematicsEBR(x, y, filterType, filterParams, fs);

% Calculate errors
error_dx4_vel = mean(abs(dx4(:,2) - dx_analytical));
error_dy4_vel = mean(abs(dy4(:,2) - dy_analytical));
error_dx4_acc = mean(abs(dx4(:,3) - d2x_analytical));
error_dy4_acc = mean(abs(dy4(:,3) - d2y_analytical));
error_dx4_jerk = mean(abs(dx4(:,4) - d3x_analytical));
error_dy4_jerk = mean(abs(dy4(:,4) - d3y_analytical));

% Report
fprintf('Method 4 Mean Absolute Errors:\n');
fprintf('  Velocity - x: %.6f, y: %.6f\n', error_dx4_vel, error_dy4_vel);
fprintf('  Acceleration - x: %.6f, y: %.6f\n', error_dx4_acc, error_dy4_acc);
fprintf('  Jerk - x: %.6f, y: %.6f\n', error_dx4_jerk, error_dy4_jerk);

%% Testing with added noise
disp('Testing all methods with added noise...');
noise_level = 0.5; % mm, standard deviation of noise
rng(42); % For reproducibility
x_noisy = x + noise_level * randn(size(x));
y_noisy = y + noise_level * randn(size(y));

%% Test Type 1 with noise
fprintf('\nMethod 1 with noise (%.2f mm std):\n', noise_level);
[dx1n, dy1n] = differentiateKinematicsEBR(x_noisy, y_noisy, 1, filterParams, fs);
error_dx1n_vel = mean(abs(dx1n(2:end,2) - dx_analytical(1:end-1)));
error_dy1n_vel = mean(abs(dy1n(2:end,2) - dy_analytical(1:end-1)));
fprintf('  Velocity - x: %.6f, y: %.6f\n', error_dx1n_vel, error_dy1n_vel);

%% Test Type 2 with noise
fprintf('\nMethod 2 with noise (%.2f mm std):\n', noise_level);
[dx2n, dy2n] = differentiateKinematicsEBR(x_noisy, y_noisy, 2, [2, 10, 1], fs);
error_dx2n_vel = mean(abs(dx2n(2:end,2) - dx_analytical(1:end-1)));
error_dy2n_vel = mean(abs(dy2n(2:end,2) - dy_analytical(1:end-1)));
fprintf('  Velocity - x: %.6f, y: %.6f\n', error_dx2n_vel, error_dy2n_vel);

%% Test Type 3 with noise
fprintf('\nMethod 3 with noise (%.2f mm std):\n', noise_level);
[dx3n, dy3n] = differentiateKinematicsEBR(x_noisy, y_noisy, 3, [2, 10, 1], fs);
error_dx3n_vel = mean(abs(dx3n(2:end,2) - dx_analytical(1:end-1)));
error_dy3n_vel = mean(abs(dy3n(2:end,2) - dy_analytical(1:end-1)));
fprintf('  Velocity - x: %.6f, y: %.6f\n', error_dx3n_vel, error_dy3n_vel);

%% Test Type 4 with noise
fprintf('\nMethod 4 with noise (%.2f mm std):\n', noise_level);
[dx4n, dy4n] = differentiateKinematicsEBR(x_noisy, y_noisy, 4, [order, framelen], fs);
error_dx4n_vel = mean(abs(dx4n(:,2) - dx_analytical));
error_dy4n_vel = mean(abs(dy4n(:,2) - dy_analytical));
fprintf('  Velocity - x: %.6f, y: %.6f\n', error_dx4n_vel, error_dy4n_vel);

%% Plot results for visualization
figure('Name', 'Differentiation Methods Comparison - Velocity (Clean Data)');

subplot(2,1,1);
plot(t, dx_analytical, 'b-', 'LineWidth', 2);
hold on;
plot(t(2:end), dx1(2:end,2), 'r--', 'LineWidth', 1);
plot(t(2:end), dx2(2:end,2), 'g-.', 'LineWidth', 1);
plot(t(2:end), dx3(2:end,2), 'm:', 'LineWidth', 1);
plot(t, dx4(:,2), 'c-', 'LineWidth', 1);
title('X Velocity - Clean Data');
legend('Analytical', 'Method 1 (FD)', 'Method 2 (Filter→FD)', 'Method 3 (FD→Filter)', 'Method 4 (SG)', 'Location', 'best');
grid on;

subplot(2,1,2);
plot(t, dy_analytical, 'b-', 'LineWidth', 2);
hold on;
plot(t(2:end), dy1(2:end,2), 'r--', 'LineWidth', 1);
plot(t(2:end), dy2(2:end,2), 'g-.', 'LineWidth', 1);
plot(t(2:end), dy3(2:end,2), 'm:', 'LineWidth', 1);
plot(t, dy4(:,2), 'c-', 'LineWidth', 1);
title('Y Velocity - Clean Data');
grid on;

% Create a second figure for acceleration
figure('Name', 'Differentiation Methods Comparison - Acceleration (Clean Data)');

subplot(2,1,1);
plot(t, d2x_analytical, 'b-', 'LineWidth', 2);
hold on;
plot(t(2:end-1), dx1(2:end-1,3), 'r--', 'LineWidth', 1);
plot(t(2:end-1), dx2(2:end-1,3), 'g-.', 'LineWidth', 1);
plot(t(2:end-1), dx3(2:end-1,3), 'm:', 'LineWidth', 1);
plot(t, dx4(:,3), 'c-', 'LineWidth', 1);
title('X Acceleration - Clean Data');
legend('Analytical', 'Method 1 (FD)', 'Method 2 (Filter→FD)', 'Method 3 (FD→Filter)', 'Method 4 (SG)', 'Location', 'best');
grid on;

subplot(2,1,2);
plot(t, d2y_analytical, 'b-', 'LineWidth', 2);
hold on;
plot(t(2:end-1), dy1(2:end-1,3), 'r--', 'LineWidth', 1);
plot(t(2:end-1), dy2(2:end-1,3), 'g-.', 'LineWidth', 1);
plot(t(2:end-1), dy3(2:end-1,3), 'm:', 'LineWidth', 1);
plot(t, dy4(:,3), 'c-', 'LineWidth', 1);
title('Y Acceleration - Clean Data');
grid on;

% Create a third figure for jerk
figure('Name', 'Differentiation Methods Comparison - Jerk (Clean Data)');

subplot(2,1,1);
plot(t, d3x_analytical, 'b-', 'LineWidth', 2);
hold on;
plot(t(2:end-2), dx1(2:end-2,4), 'r--', 'LineWidth', 1);
plot(t(2:end-2), dx2(2:end-2,4), 'g-.', 'LineWidth', 1);
plot(t(2:end-2), dx3(2:end-2,4), 'm:', 'LineWidth', 1);
plot(t, dx4(:,4), 'c-', 'LineWidth', 1);
title('X Jerk - Clean Data');
legend('Analytical', 'Method 1 (FD)', 'Method 2 (Filter→FD)', 'Method 3 (FD→Filter)', 'Method 4 (SG)', 'Location', 'best');
grid on;

subplot(2,1,2);
plot(t, d3y_analytical, 'b-', 'LineWidth', 2);
hold on;
plot(t(2:end-2), dy1(2:end-2,4), 'r--', 'LineWidth', 1);
plot(t(2:end-2), dy2(2:end-2,4), 'g-.', 'LineWidth', 1);
plot(t(2:end-2), dy3(2:end-2,4), 'm:', 'LineWidth', 1);
plot(t, dy4(:,4), 'c-', 'LineWidth', 1);
title('Y Jerk - Clean Data');
grid on;

% Create a fourth figure for noisy data
figure('Name', 'Differentiation Methods Comparison - Velocity (Noisy Data)');

subplot(2,1,1);
plot(t, dx_analytical, 'b-', 'LineWidth', 2);
hold on;
plot(t(2:end), dx1n(2:end,2), 'r--', 'LineWidth', 1);
plot(t(2:end), dx2n(2:end,2), 'g-.', 'LineWidth', 1);
plot(t(2:end), dx3n(2:end,2), 'm:', 'LineWidth', 1);
plot(t, dx4n(:,2), 'c-', 'LineWidth', 1);
title(sprintf('X Velocity - Noisy Data (%.2f mm std)', noise_level));
legend('Analytical', 'Method 1 (FD)', 'Method 2 (Filter→FD)', 'Method 3 (FD→Filter)', 'Method 4 (SG)', 'Location', 'best');
grid on;

subplot(2,1,2);
plot(t, dy_analytical, 'b-', 'LineWidth', 2);
hold on;
plot(t(2:end), dy1n(2:end,2), 'r--', 'LineWidth', 1);
plot(t(2:end), dy2n(2:end,2), 'g-.', 'LineWidth', 1);
plot(t(2:end), dy3n(2:end,2), 'm:', 'LineWidth', 1);
plot(t, dy4n(:,2), 'c-', 'LineWidth', 1);
title(sprintf('Y Velocity - Noisy Data (%.2f mm std)', noise_level));
grid on;

%% Summary of results
disp('-------------------------------------------------');
disp('SUMMARY OF RESULTS');
disp('-------------------------------------------------');
fprintf('Error Summary (Clean Data):\n');
disp('Method 1 (Finite Differences):');
fprintf('  Velocity: %.6f\n', mean([error_dx1_vel, error_dy1_vel]));
fprintf('  Acceleration: %.6f\n', mean([error_dx1_acc, error_dy1_acc]));
fprintf('  Jerk: %.6f\n', mean([error_dx1_jerk, error_dy1_jerk]));

disp('Method 2 (Filter then Differentiate):');
fprintf('  Velocity: %.6f\n', mean([error_dx2_vel, error_dy2_vel]));
fprintf('  Acceleration: %.6f\n', mean([error_dx2_acc, error_dy2_acc]));
fprintf('  Jerk: %.6f\n', mean([error_dx2_jerk, error_dy2_jerk]));

disp('Method 3 (Differentiate then Filter):');
fprintf('  Velocity: %.6f\n', mean([error_dx3_vel, error_dy3_vel]));
fprintf('  Acceleration: %.6f\n', mean([error_dx3_acc, error_dy3_acc]));
fprintf('  Jerk: %.6f\n', mean([error_dx3_jerk, error_dy3_jerk]));

disp('Method 4 (Savitzky-Golay):');
fprintf('  Velocity: %.6f\n', mean([error_dx4_vel, error_dy4_vel]));
fprintf('  Acceleration: %.6f\n', mean([error_dx4_acc, error_dy4_acc]));
fprintf('  Jerk: %.6f\n', mean([error_dx4_jerk, error_dy4_jerk]));

fprintf('\nError Summary (Noisy Data, %.2f mm std):\n', noise_level);
fprintf('  Method 1 Velocity: %.6f\n', mean([error_dx1n_vel, error_dy1n_vel]));
fprintf('  Method 2 Velocity: %.6f\n', mean([error_dx2n_vel, error_dy2n_vel]));
fprintf('  Method 3 Velocity: %.6f\n', mean([error_dx3n_vel, error_dy3n_vel]));
fprintf('  Method 4 Velocity: %.6f\n', mean([error_dx4n_vel, error_dy4n_vel]));

% Calculate improvement percentages
m1_to_m4_clean = (mean([error_dx1_vel, error_dy1_vel]) - mean([error_dx4_vel, error_dy4_vel])) / mean([error_dx1_vel, error_dy1_vel]) * 100;
m1_to_m4_noisy = (mean([error_dx1n_vel, error_dy1n_vel]) - mean([error_dx4n_vel, error_dy4n_vel])) / mean([error_dx1n_vel, error_dy1n_vel]) * 100;

fprintf('\nRelative Performance:\n');
fprintf('  SG method improves over basic FD by %.2f%% (clean data)\n', m1_to_m4_clean);
fprintf('  SG method improves over basic FD by %.2f%% (noisy data)\n', m1_to_m4_noisy);

disp('-------------------------------------------------');
disp('TEST COMPLETED');
disp('-------------------------------------------------');


