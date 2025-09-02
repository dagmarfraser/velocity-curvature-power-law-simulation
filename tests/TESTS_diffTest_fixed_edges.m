%% Test Script for differentiateKinematicsEBR function with proper edge handling
% This script tests all methods and properly handles edge effects for SG filter
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

disp('Running tests for differentiateKinematicsEBR function with proper edge handling...');

%% Create a fixed version of differentiateKinematicsEBR function for testing
% Define the fixed function without the * p multiplication
differentiateKinematicsEBR_fixed = @(x, y, filterType, filterParams, fs) fixed_differentiateKinematics(x, y, filterType, filterParams, fs);

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
[dx1, dy1] = differentiateKinematicsEBR_fixed(x, y, filterType, filterParams, fs);

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
[dx2, dy2] = differentiateKinematicsEBR_fixed(x, y, filterType, filterParams, fs);

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
[dx3, dy3] = differentiateKinematicsEBR_fixed(x, y, filterType, filterParams, fs);

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
[dx4, dy4] = differentiateKinematicsEBR_fixed(x, y, filterType, filterParams, fs);

% Exclude edge points for SG method
edge_points = floor(framelen/2);
valid_indices = (edge_points+1):(length(x)-edge_points);

% Calculate errors for SG method (excluding edge points)
error_dx4_vel = mean(abs(dx4(valid_indices,2) - dx_analytical(valid_indices)));
error_dy4_vel = mean(abs(dy4(valid_indices,2) - dy_analytical(valid_indices)));
error_dx4_acc = mean(abs(dx4(valid_indices,3) - d2x_analytical(valid_indices)));
error_dy4_acc = mean(abs(dy4(valid_indices,3) - d2y_analytical(valid_indices)));
error_dx4_jerk = mean(abs(dx4(valid_indices,4) - d3x_analytical(valid_indices)));
error_dy4_jerk = mean(abs(dy4(valid_indices,4) - d3y_analytical(valid_indices)));

% Report
fprintf('Method 4 Mean Absolute Errors (excluding %d edge points):\n', edge_points);
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
[dx1n, dy1n] = differentiateKinematicsEBR_fixed(x_noisy, y_noisy, 1, filterParams, fs);
error_dx1n_vel = mean(abs(dx1n(2:end,2) - dx_analytical(1:end-1)));
error_dy1n_vel = mean(abs(dy1n(2:end,2) - dy_analytical(1:end-1)));
fprintf('  Velocity - x: %.6f, y: %.6f\n', error_dx1n_vel, error_dy1n_vel);

%% Test Type 2 with noise
fprintf('\nMethod 2 with noise (%.2f mm std):\n', noise_level);
[dx2n, dy2n] = differentiateKinematicsEBR_fixed(x_noisy, y_noisy, 2, [2, 10, 1], fs);
error_dx2n_vel = mean(abs(dx2n(2:end,2) - dx_analytical(1:end-1)));
error_dy2n_vel = mean(abs(dy2n(2:end,2) - dy_analytical(1:end-1)));
fprintf('  Velocity - x: %.6f, y: %.6f\n', error_dx2n_vel, error_dy2n_vel);

%% Test Type 3 with noise
fprintf('\nMethod 3 with noise (%.2f mm std):\n', noise_level);
[dx3n, dy3n] = differentiateKinematicsEBR_fixed(x_noisy, y_noisy, 3, [2, 10, 1], fs);
error_dx3n_vel = mean(abs(dx3n(2:end,2) - dx_analytical(1:end-1)));
error_dy3n_vel = mean(abs(dy3n(2:end,2) - dy_analytical(1:end-1)));
fprintf('  Velocity - x: %.6f, y: %.6f\n', error_dx3n_vel, error_dy3n_vel);

%% Test Type 4 with noise
fprintf('\nMethod 4 with noise (%.2f mm std):\n', noise_level);
[dx4n, dy4n] = differentiateKinematicsEBR_fixed(x_noisy, y_noisy, 4, [order, framelen], fs);
error_dx4n_vel = mean(abs(dx4n(valid_indices,2) - dx_analytical(valid_indices)));
error_dy4n_vel = mean(abs(dy4n(valid_indices,2) - dy_analytical(valid_indices)));
fprintf('  Velocity - x: %.6f, y: %.6f (excluding %d edge points)\n', error_dx4n_vel, error_dy4n_vel, edge_points);

%% Plot results for visualization
figure('Name', 'Differentiation Methods Comparison - Velocity (Clean Data)');

subplot(2,1,1);
plot(t, dx_analytical, 'b-', 'LineWidth', 2);
hold on;
plot(t(2:end), dx1(2:end,2), 'r--', 'LineWidth', 1);
plot(t(2:end), dx2(2:end,2), 'g-.', 'LineWidth', 1);
plot(t(2:end), dx3(2:end,2), 'm:', 'LineWidth', 1);
plot(t, dx4(:,2), 'c-', 'LineWidth', 1);
plot(t(valid_indices), dx4(valid_indices,2), 'c-', 'LineWidth', 2);  % Highlight valid points
title('X Velocity - Clean Data');
legend('Analytical', 'Method 1 (FD)', 'Method 2 (Filter→FD)', 'Method 3 (FD→Filter)', 'Method 4 (SG)', 'SG Valid Points', 'Location', 'best');
grid on;

subplot(2,1,2);
plot(t, dy_analytical, 'b-', 'LineWidth', 2);
hold on;
plot(t(2:end), dy1(2:end,2), 'r--', 'LineWidth', 1);
plot(t(2:end), dy2(2:end,2), 'g-.', 'LineWidth', 1);
plot(t(2:end), dy3(2:end,2), 'm:', 'LineWidth', 1);
plot(t, dy4(:,2), 'c-', 'LineWidth', 1);
plot(t(valid_indices), dy4(valid_indices,2), 'c-', 'LineWidth', 2);  % Highlight valid points
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
plot(t(valid_indices), dx4(valid_indices,3), 'c-', 'LineWidth', 2);  % Highlight valid points
title('X Acceleration - Clean Data');
legend('Analytical', 'Method 1 (FD)', 'Method 2 (Filter→FD)', 'Method 3 (FD→Filter)', 'Method 4 (SG)', 'SG Valid Points', 'Location', 'best');
grid on;

subplot(2,1,2);
plot(t, d2y_analytical, 'b-', 'LineWidth', 2);
hold on;
plot(t(2:end-1), dy1(2:end-1,3), 'r--', 'LineWidth', 1);
plot(t(2:end-1), dy2(2:end-1,3), 'g-.', 'LineWidth', 1);
plot(t(2:end-1), dy3(2:end-1,3), 'm:', 'LineWidth', 1);
plot(t, dy4(:,3), 'c-', 'LineWidth', 1);
plot(t(valid_indices), dy4(valid_indices,3), 'c-', 'LineWidth', 2);  % Highlight valid points
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
plot(t(valid_indices), dx4(valid_indices,4), 'c-', 'LineWidth', 2);  % Highlight valid points
title('X Jerk - Clean Data');
legend('Analytical', 'Method 1 (FD)', 'Method 2 (Filter→FD)', 'Method 3 (FD→Filter)', 'Method 4 (SG)', 'SG Valid Points', 'Location', 'best');
grid on;

subplot(2,1,2);
plot(t, d3y_analytical, 'b-', 'LineWidth', 2);
hold on;
plot(t(2:end-2), dy1(2:end-2,4), 'r--', 'LineWidth', 1);
plot(t(2:end-2), dy2(2:end-2,4), 'g-.', 'LineWidth', 1);
plot(t(2:end-2), dy3(2:end-2,4), 'm:', 'LineWidth', 1);
plot(t, dy4(:,4), 'c-', 'LineWidth', 1);
plot(t(valid_indices), dy4(valid_indices,4), 'c-', 'LineWidth', 2);  % Highlight valid points
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
plot(t(valid_indices), dx4n(valid_indices,2), 'c-', 'LineWidth', 2);  % Highlight valid points
title(sprintf('X Velocity - Noisy Data (%.2f mm std)', noise_level));
legend('Analytical', 'Method 1 (FD)', 'Method 2 (Filter→FD)', 'Method 3 (FD→Filter)', 'Method 4 (SG)', 'SG Valid Points', 'Location', 'best');
grid on;

subplot(2,1,2);
plot(t, dy_analytical, 'b-', 'LineWidth', 2);
hold on;
plot(t(2:end), dy1n(2:end,2), 'r--', 'LineWidth', 1);
plot(t(2:end), dy2n(2:end,2), 'g-.', 'LineWidth', 1);
plot(t(2:end), dy3n(2:end,2), 'm:', 'LineWidth', 1);
plot(t, dy4n(:,2), 'c-', 'LineWidth', 1);
plot(t(valid_indices), dy4n(valid_indices,2), 'c-', 'LineWidth', 2);  % Highlight valid points
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

disp('Method 4 (Savitzky-Golay, excluding edge points):');
fprintf('  Velocity: %.6f\n', mean([error_dx4_vel, error_dy4_vel]));
fprintf('  Acceleration: %.6f\n', mean([error_dx4_acc, error_dy4_acc]));
fprintf('  Jerk: %.6f\n', mean([error_dx4_jerk, error_dy4_jerk]));

fprintf('\nError Summary (Noisy Data, %.2f mm std):\n', noise_level);
fprintf('  Method 1 Velocity: %.6f\n', mean([error_dx1n_vel, error_dy1n_vel]));
fprintf('  Method 2 Velocity: %.6f\n', mean([error_dx2n_vel, error_dy2n_vel]));
fprintf('  Method 3 Velocity: %.6f\n', mean([error_dx3n_vel, error_dy3n_vel]));
fprintf('  Method 4 Velocity: %.6f (excluding edge points)\n', mean([error_dx4n_vel, error_dy4n_vel]));

% Calculate improvement percentages
m1_to_m4_clean = (mean([error_dx1_vel, error_dy1_vel]) - mean([error_dx4_vel, error_dy4_vel])) / mean([error_dx1_vel, error_dy1_vel]) * 100;
m1_to_m4_noisy = (mean([error_dx1n_vel, error_dy1n_vel]) - mean([error_dx4n_vel, error_dy4n_vel])) / mean([error_dx1n_vel, error_dy1n_vel]) * 100;

fprintf('\nRelative Performance:\n');
fprintf('  SG method improves over basic FD by %.2f%% (clean data)\n', m1_to_m4_clean);
fprintf('  SG method improves over basic FD by %.2f%% (noisy data)\n', m1_to_m4_noisy);

disp('-------------------------------------------------');
disp('RECOMMENDATIONS:');
disp('-------------------------------------------------');
disp('1. Fix the Savitzky-Golay implementation by:');
disp('   a. Removing the "* p" multiplication');
disp('   b. Using appropriate polynomial order (order ≥ 3 for computing jerk)');
disp('   c. Excluding edge points from analysis (framelen/2 points from each end)');
disp('2. Use the fixed implementation for power law analysis');
disp('-------------------------------------------------');

disp('TEST COMPLETED');
disp('-------------------------------------------------');

%% Define the fixed version of differentiateKinematicsEBR
function [dx, dy] = fixed_differentiateKinematics(x, y, filterType, filterParams, fs)
    % This is a copy of differentiateKinematicsEBR with the * p multiplication removed
    % and support for polynomial order >= 3 for computing jerk
    
    dx = zeros(length(x),4);
    dy = zeros(length(y),4);
    dt = 1/fs;

    switch filterType

        case 1 % Finite Differences Differentation
            dx(:,1) = x; % not smoothed!
            dy(:,1) = y;

            dx(2:end,2) = diff(x,1) * fs; % raw diff -  scaled by the fs
            dy(2:end,2) = diff(y,1) * fs; % velocity

            dx(2:end-1,3) = diff(x,2) * fs * fs; % raw diff -  scaled by the fs^2
            dy(2:end-1,3) = diff(y,2) * fs * fs; % acceleration

            dx(2:end-2,4) = diff(x,3) * fs * fs * fs; % raw diff -  scaled by the fs^3
            dy(2:end-2,4) = diff(y,3) * fs * fs * fs; % jerk

        case 2 % Nth Order Fp Hz Low pass filter, filtfilt for zero lag followed by Finite Differences 
            dx(:,1) = x; % not smoothed!
            dy(:,1) = y;
            
            Fp = filterParams(2);
            N = filterParams(1);
            zeroLag = filterParams(3);
            %% make an Nth order Butterworth zero lag low pass filter with corner frequency Fp
            fc = Fp;
            bOrder = N; % this is filter order, will be 2*bOrder when used with filtfilt below
            [b,a] = butter(bOrder,fc/(fs/2));

            % use zerophase digital filtering i.e. filtfilt
            % diff the output of the low pass filter
            if ~zeroLag
                dx(2:end,2) = diff(filter(b,a, (x))) * fs; % raw diff -  scaled by the fs, and now diminshed by the filter
                dy(2:end,2) = diff(filter(b,a, (y))) * fs; % velocity
            else
                dx(2:end,2) = diff(filtfilt(b,a, (x))) * fs; % raw diff -  scaled by the fs, and now diminshed by the filter
                dy(2:end,2) = diff(filtfilt(b,a, (y))) * fs; % velocity
            end

            dx(2:end-1,3) = diff(dx(2:end,2),1) * fs; % raw diff -  scaled by the fs
            dy(2:end-1,3) = diff(dy(2:end,2),1) * fs; % acceleration

            dx(2:end-2,4) = diff(dx(2:end,2),2) * fs * fs; % raw diff - scaled by the fs
            dy(2:end-2,4) = diff(dy(2:end,2),2) * fs * fs; % jerk

        case 3 % Finite Differences with Nth Order Fp Hz Low pass filter
            dx(:,1) = x; % not smoothed!
            dy(:,1) = y;

            Fp = filterParams(2);
            N = filterParams(1);
            zeroLag = filterParams(3);

            %% make an Nth order Butterworth zero lag low pass filter with corner frequency Fp
            fc = Fp;
            bOrder = N; % this is filter order, will be 2*bOrder when used with filtfilt below
            [b,a] = butter(bOrder,fc/(fs/2)); % we pass the order and the cut off freq / Nyquist Frequency (i.e. half sampling rate)

            % use zerophase digital filtering i.e. filtfilt
            % low pass filter the output of the diff
            if ~zeroLag
                dx(2:end,2) = filter(b,a, (diff(x,1)*fs)); % raw diff -  scaled by the fs, and now diminshed by the filter
                dy(2:end,2) = filter(b,a, (diff(y,1)*fs));
            else
                dx(2:end,2) = filtfilt(b,a, (diff(x,1)*fs)); % raw diff -  scaled by the fs, and now diminshed by the filter
                dy(2:end,2) = filtfilt(b,a, (diff(y,1)*fs));
            end
            dx(2:end-1,3) = diff(dx(2:end,2),1) * fs; % raw diff -  scaled by the fs
            dy(2:end-1,3) = diff(dy(2:end,2),1) * fs;

            dx(2:end-2,4) = diff(dx(2:end,2),2) * fs * fs; % raw diff - scaled by the fs
            dy(2:end-2,4) = diff(dy(2:end,2),2) * fs * fs;

        case 4 % Savitzky-Golay smooth differentiation
            % adapted from https://uk.mathworks.com/help/signal/ref/sgolay.html
            % example

            order = filterParams(1);
            framelen = filterParams(2);
            padding = (framelen-1)/2;

            % Check if order is sufficient for computing up to jerk
            maxDerivative = min(order, 3);
            
            % Initialize with original data
            dx(:,1) = x;
            dy(:,1) = y;
            
            % Get Savitzky-Golay coefficients
            [b,g] = sgolay(order, framelen);

            % Direct implementation from MATLAB documentation
            % Calculate derivatives directly using convolution
            for p = 1:maxDerivative
                dx(:,p+1) = conv(x, factorial(p)/(-dt)^p * g(:,p+1), 'same');
                dy(:,p+1) = conv(y, factorial(p)/(-dt)^p * g(:,p+1), 'same');
            end

        case 5
            disp('NOT IMPLEMNTED')
    end
end
