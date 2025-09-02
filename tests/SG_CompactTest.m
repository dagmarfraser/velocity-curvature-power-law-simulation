%% Compact Test Script for differentiateKinematicsEBR
% Tests the function with a minimal, correct implementation of SG differentiation 
% that removes the * p multiplication and properly handles edge effects

clear all; close all; clc;
addpath(genpath('../src/functions'));

% Generate a perfect circle for testing
fs = 100; t = (0:1/fs:1)';
radius = 100; x = radius * cos(2*pi*t); y = radius * sin(2*pi*t);

% Analytical derivatives
dx_analytical = -radius * 2*pi * sin(2*pi*t);
dy_analytical = radius * 2*pi * cos(2*pi*t);

% Fixed implementation
function [dx, dy] = fixedSG(x, y, order, framelen, fs)
    dx = zeros(length(x), 4);
    dy = zeros(length(y), 4);
    dt = 1/fs;
    
    % Store original data
    dx(:,1) = x; dy(:,1) = y;
    
    % Get SG coefficients
    [~, g] = sgolay(order, framelen);
    
    % Direct SG implementation without the * p
    for p = 1:min(order, 3)
        dx(:,p+1) = conv(x, factorial(p)/(-dt)^p * g(:,p+1), 'same');
        dy(:,p+1) = conv(y, factorial(p)/(-dt)^p * g(:,p+1), 'same');
    end
end

% Original implementation (simplified for clarity)
function [dx, dy] = originalSG(x, y, order, framelen, fs)
    dx = zeros(length(x), 4);
    dy = zeros(length(y), 4);
    dt = 1/fs;
    padding = (framelen-1)/2;
    
    % Store original data
    dx(:,1) = x; dy(:,1) = y;
    
    % Get SG coefficients
    [~, g] = sgolay(order, framelen);
    
    % Implementation with the problematic * p
    for p = 1:min(order, 3)
        dxElement = conv(x, factorial(p)/(-dt)^p * g(:,p+1), 'same');
        dyElement = conv(y, factorial(p)/(-dt)^p * g(:,p+1), 'same');
        
        % Note the * p multiplication here!
        dxPadded = padarray(dxElement * p, padding/2, 0, 'pre');
        dyPadded = padarray(dyElement * p, padding/2, 0, 'pre');
        
        dx(:,p+1) = dxPadded(2:length(x)+1);
        dy(:,p+1) = dyPadded(2:length(y)+1);
    end
end

% Run tests
disp('TESTING SG DIFFERENTIATION IMPLEMENTATIONS');
framelen = 9;
order = 3;  % 3rd order needed for jerk

% Get derivatives
[dx_fixed, dy_fixed] = fixedSG(x, y, order, framelen, fs);
[dx_orig, dy_orig] = originalSG(x, y, order, framelen, fs);

% Exclude edge points to handle boundary effects
edge_points = floor(framelen/2);
valid_indices = (edge_points+1):(length(x)-edge_points);

% Calculate errors for fixed implementation
error_fixed_vel = mean(abs([dx_fixed(valid_indices,2) - dx_analytical(valid_indices); 
                          dy_fixed(valid_indices,2) - dy_analytical(valid_indices)]));
                          
% Calculate errors for original implementation
error_orig_vel = mean(abs([dx_orig(valid_indices,2) - dx_analytical(valid_indices); 
                         dy_orig(valid_indices,2) - dy_analytical(valid_indices)]));

% Calculate scaling factors between implementations
scaling_vel = mean(abs(dx_orig(valid_indices,2))) / mean(abs(dx_fixed(valid_indices,2)));
scaling_acc = mean(abs(dx_orig(valid_indices,3))) / mean(abs(dx_fixed(valid_indices,3)));
scaling_jerk = mean(abs(dx_orig(valid_indices,4))) / mean(abs(dx_fixed(valid_indices,4)));

% Display results
fprintf('RESULTS SUMMARY (excluding %d edge points)\n', edge_points);
fprintf('Original implementation error: %.6f\n', error_orig_vel);
fprintf('Fixed implementation error: %.6f\n', error_fixed_vel);
fprintf('Improvement: %.2f%%\n', (error_orig_vel - error_fixed_vel)/error_orig_vel * 100);

fprintf('\nSCALING FACTORS (Original/Fixed)\n');
fprintf('Velocity scaling: %.2f (should be ~1.0 if only scaling is p=1)\n', scaling_vel);
fprintf('Acceleration scaling: %.2f (should be ~2.0 if only scaling is p=2)\n', scaling_acc);
fprintf('Jerk scaling: %.2f (should be ~3.0 if only scaling is p=3)\n', scaling_jerk);

% Plot results - velocity
figure('Name', 'Velocity Comparison');
plot(t, dx_analytical, 'k-', 'LineWidth', 2);
hold on;
plot(t, dx_fixed(:,2), 'g-', 'LineWidth', 1);
plot(t, dx_orig(:,2), 'r-', 'LineWidth', 1);
plot(t(valid_indices), dx_fixed(valid_indices,2), 'g-', 'LineWidth', 2);
plot(t(valid_indices), dx_orig(valid_indices,2), 'r-', 'LineWidth', 2);
legend('Analytical', 'Fixed SG', 'Original SG', 'Fixed Valid', 'Original Valid');
title('Velocity Comparison (X-component)');
grid on;

disp('-------------------------------------------------');
disp('RECOMMENDATIONS:');
disp('-------------------------------------------------');
disp('1. Fix the Savitzky-Golay implementation by:');
disp('   a. Removing the "* p" multiplication');
disp('   b. Using the direct convolution approach from MATLAB documentation');
disp('   c. Always exclude edge points from analysis');
disp('2. Use at least 3rd order polynomial for computing up to jerk');
disp('-------------------------------------------------');
