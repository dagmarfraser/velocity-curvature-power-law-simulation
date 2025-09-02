function [method, a, b, a_confint, b_confint] = power_analysis(x, y, varargin)
% Dagmar Fraser 2024 d.s.fraser@bham.ac.uk
% MATLAB Reimplementation of Xiao et al.'s 2024 
% POWER_ANALYSIS Perform power law analysis on input data
%   [METHOD, A, B, A_CONFINT, B_CONFINT] = POWER_ANALYSIS(X, Y) performs
%   power law analysis on the input vectors X and Y.
%
%   [...] = POWER_ANALYSIS(..., 'CI_boot', true/false) specifies whether to
%   compute confidence intervals using bootstrapping for model averaging.
%   Default is true.
%
%   [...] = POWER_ANALYSIS(..., 'diagno', true/false) specifies whether to
%   create residual plots as visual diagnostics. Default is true.
%
%   [...] = POWER_ANALYSIS(..., 'output_plot', true/false) specifies whether
%   to create a scatter plot with fitted relationship. Default is false.

% Input parsing and validation
p = inputParser;
addRequired(p, 'x', @(x) isnumeric(x) && isvector(x) && ~any(isnan(x) | isinf(x)));
addRequired(p, 'y', @(y) isnumeric(y) && isvector(y) && ~any(isnan(y) | isinf(y)));
addParameter(p, 'CI_boot', true, @islogical);
addParameter(p, 'diagno', true, @islogical);
addParameter(p, 'output_plot', false, @islogical);
parse(p, x, y, varargin{:});

CI_boot = p.Results.CI_boot;
diagno = p.Results.diagno;
output_plot = p.Results.output_plot;

% Ensure x and y are column vectors
x = x(:);
y = y(:);

if length(x) ~= length(y)
    error('Input vectors x and y must have the same length.');
end

% Step 1: Likelihood analysis
model_lr = fitlm(log(x), log(y));
a_lr = exp(model_lr.Coefficients.Estimate(1));
b_lr = model_lr.Coefficients.Estimate(2);
sd_lr = std(log(y) - (log(a_lr) + b_lr * log(x)));

% Define the nonlinear model function
powerlaw = @(b,x) b(1) * x.^b(2);

% Fit the nonlinear model
model_nlr = fitnlm(x, y, powerlaw, [a_lr, b_lr]);
a_nlr = model_nlr.Coefficients.Estimate(1);
b_nlr = model_nlr.Coefficients.Estimate(2);
sd_nlr = std(y - a_nlr * x.^b_nlr);

l_logn = sum(log(lognpdf(y, log(a_lr * x.^b_lr), sd_lr)));
l_norm = sum(log(normpdf(y, a_nlr * x.^b_nlr, sd_nlr)));

n = length(x);
k = 3;
AICc_logn = 2 * k - 2 * l_logn + 2 * k * (k + 1) / (n - k - 1);
AICc_norm = 2 * k - 2 * l_norm + 2 * k * (k + 1) / (n - k - 1);
delta_AICc = AICc_norm - AICc_logn;
fprintf('AICc_logn: %f\nAICc_norm: %f\n', AICc_logn, AICc_norm);

w_logn = exp(-(AICc_logn - min(AICc_logn, AICc_norm)) / 2);
w_norm = exp(-(AICc_norm - min(AICc_logn, AICc_norm)) / 2);
weight_logn = w_logn / (w_logn + w_norm);
weight_norm = w_norm / (w_logn + w_norm);

% Step 2: Determine analysis method
if delta_AICc < -2
    fprintf('The assumption of additive normal error is better supported.\nProceed with NLR.\n');
    method = 'NLR';
    a = a_nlr;
    b = b_nlr;
    a_confint = coefCI(model_nlr);
    a_confint = a_confint(1, :);
    b_confint = coefCI(model_nlr);
    b_confint = b_confint(2, :);
elseif delta_AICc > 2
    fprintf('The assumption of multiplicative log-normal error is better supported.\nProceed with LR.\n');
    method = 'LR';
    a = a_lr;
    b = b_lr;
    a_confint = coefCI(model_lr);
    a_confint = exp(a_confint(1, :));
    b_confint = coefCI(model_lr);
    b_confint = b_confint(2, :);
else
    fprintf('The two error distributions have similar support.\nProceed with model averaging.\n');
    method = 'Model Averaging';
    a = a_lr * weight_logn + a_nlr * weight_norm;
    b = b_lr * weight_logn + b_nlr * weight_norm;
    
    if ~CI_boot
        a_confint = [NaN, NaN];
        b_confint = [NaN, NaN];
    else
        % Implement bootstrapping for confidence intervals
        % Note: This is a simplified version and may need further refinement
        nboot = 1000;
        boot_estimates = zeros(nboot, 2);
        for i = 1:nboot
            boot_indices = randi(n, n, 1);
            x_boot = x(boot_indices);
            y_boot = y(boot_indices);
            
            model_lr_boot = fitlm(log(x_boot), log(y_boot));
            a_lr_boot = exp(model_lr_boot.Coefficients.Estimate(1));
            b_lr_boot = model_lr_boot.Coefficients.Estimate(2);
            
            model_nlr_boot = fitnlm(x_boot, y_boot, powerlaw, [a_lr_boot, b_lr_boot]);
            a_nlr_boot = model_nlr_boot.Coefficients.Estimate(1);
            b_nlr_boot = model_nlr_boot.Coefficients.Estimate(2);
            
            a_boot = a_lr_boot * weight_logn + a_nlr_boot * weight_norm;
            b_boot = b_lr_boot * weight_logn + b_nlr_boot * weight_norm;
            
            boot_estimates(i, :) = [a_boot, b_boot];
        end
        
        a_confint = prctile(boot_estimates(:, 1), [2.5, 97.5]);
        b_confint = prctile(boot_estimates(:, 2), [2.5, 97.5]);
    end
end

fprintf('a: %f\nb: %f\n', a, b);

% Step 3: Residual plots as visual diagnostics
if diagno
    figure(1001);
    clf
    subplot(2, 2, 1);
    histogram(model_lr.Residuals.Raw, 'Normalization', 'pdf');
    hold on;
    x_range = linspace(min(model_lr.Residuals.Raw), max(model_lr.Residuals.Raw), 100);
    plot(x_range, normpdf(x_range, mean(model_lr.Residuals.Raw), std(model_lr.Residuals.Raw)));
    title('LR Residuals');
    
    subplot(2, 2, 2);
    scatter(model_lr.Fitted, model_lr.Residuals.Raw);
    xlabel('Predicted y');
    ylabel('Residuals');
    title('LR Homoscedasticity');
    
    subplot(2, 2, 3);
    histogram(model_nlr.Residuals.Raw, 'Normalization', 'pdf');
    hold on;
    x_range = linspace(min(model_nlr.Residuals.Raw), max(model_nlr.Residuals.Raw), 100);
    plot(x_range, normpdf(x_range, mean(model_nlr.Residuals.Raw), std(model_nlr.Residuals.Raw)));
    title('NLR Residuals');
    
    subplot(2, 2, 4);
    scatter(model_nlr.Fitted, model_nlr.Residuals.Raw);
    xlabel('Predicted y');
    ylabel('Residuals');
    title('NLR Homoscedasticity');
end

% Output plot
if output_plot
    figure(1002);
    clf
    subplot(1, 2, 1);
    loglog(x, y, 'o');
    hold on;
    x_range = logspace(log10(min(x)), log10(max(x)), 100);
    loglog(x_range, a * x_range.^b, '--');
    xlabel('x');
    ylabel('y');
    title('Logarithmic Scale');
    
    subplot(1, 2, 2);
    scatter(x, y);
    hold on;
    x_range = linspace(min(x), max(x), 100);
    plot(x_range, a * x_range.^b, '--');
    xlabel('x');
    ylabel('y');
    title('Arithmetic Scale');
    
    suptitle(['Fitting Power-Law Data with ', method]);
end

end