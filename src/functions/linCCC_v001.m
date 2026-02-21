function cccResults = linCCC_v001(y1, y2, alpha)
% Lin's Concordance Correlation Coefficient (CCC)
%
% Measures agreement between two continuous variables by evaluating
% deviation from the 45° line (perfect concordance). Combines precision
% (Pearson's r) and accuracy (bias correction factor Cb).
%
% USAGE:
%   cccResults = linCCC_v001(y1, y2)
%   cccResults = linCCC_v001(y1, y2, alpha)
%
% INPUTS:
%   y1    : First measurement vector (n x 1)
%   y2    : Second measurement vector (n x 1) - same length as y1
%   alpha : Significance level for CI (default 0.05)
%
% OUTPUTS:
%   cccResults : struct with fields
%     .ccc      - Concordance correlation coefficient (rho_c)
%     .r        - Pearson correlation coefficient (precision)
%     .Cb       - Bias correction factor (accuracy)
%     .LowerCI  - Lower 95% confidence bound (Z-transformed)
%     .UpperCI  - Upper 95% confidence bound (Z-transformed)
%     .n        - Sample size
%     .u        - Location shift (standardised mean difference)
%     .v        - Scale shift (SD ratio)
%     .Z        - Z-transformed CCC (for asymptotic inference)
%     .seZ      - Standard error of Z
%
% INTERPRETATION (McBride, 2005):
%   CCC > 0.99  : Almost perfect agreement
%   CCC 0.95-0.99: Substantial agreement
%   CCC 0.90-0.95: Moderate agreement
%   CCC < 0.90  : Poor agreement
%
% FORMULA:
%   rho_c = 2*sigma_12 / (sigma_1^2 + sigma_2^2 + (mu_1 - mu_2)^2)
%         = rho * Cb
%
%   where Cb = 2 / (v + 1/v + u^2)
%         v  = sigma_1 / sigma_2  (scale shift)
%         u  = (mu_1 - mu_2) / sqrt(sigma_1 * sigma_2)  (location shift)
%
% REFERENCES:
%   Lin LI-K (1989). A concordance correlation coefficient to evaluate
%   reproducibility. Biometrics 45(1):255-268.
%
%   Lin LI-K (2000). A note on the concordance correlation coefficient.
%   Biometrics 56(1):324-325. [Correction to variance formula]
%
%   McBride GB (2005). A proposal for strength-of-agreement criteria
%   for Lin's concordance correlation coefficient. NIWA Client Report.
%
% SEE ALSO: corrcoef, computeCCC_R
%
% VERSION: v001 (2026-01-15)
% AUTHOR: D.S. Fraser, University of Birmingham

%% Input validation
if nargin < 2
    error('linCCC:TooFewInputs', 'Requires two input vectors y1 and y2');
end
if nargin < 3
    alpha = 0.05;
end

% Ensure column vectors
y1 = y1(:);
y2 = y2(:);

% Check equal length
if length(y1) ~= length(y2)
    error('linCCC:UnequalLength', 'y1 and y2 must have the same length');
end

% Remove NaN pairs (pairwise deletion)
validIdx = ~isnan(y1) & ~isnan(y2);
y1 = y1(validIdx);
y2 = y2(validIdx);
n = length(y1);

if n < 3
    error('linCCC:InsufficientData', 'Need at least 3 valid pairs');
end

%% Compute sample statistics
% Means
mu1 = mean(y1);
mu2 = mean(y2);

% Population variances (1/n normalisation per Lin 1989)
% Note: Lin's original uses N, not N-1 in denominator
var1 = var(y1, 1);  % 1 = normalise by N (population variance)
var2 = var(y2, 1);
sigma1 = sqrt(var1);
sigma2 = sqrt(var2);

% Covariance (also population, 1/n)
sigma12 = sum((y1 - mu1) .* (y2 - mu2)) / n;

%% Concordance Correlation Coefficient
% Direct formula from Lin (1989) Eq. 1
numerator = 2 * sigma12;
denominator = var1 + var2 + (mu1 - mu2)^2;
ccc = numerator / denominator;

%% Decomposition into precision (r) and accuracy (Cb)
% Pearson correlation coefficient
r = sigma12 / (sigma1 * sigma2);

% Location shift (u) and scale shift (v)
if sigma1 > 0 && sigma2 > 0
    u = (mu1 - mu2) / sqrt(sigma1 * sigma2);
    v = sigma1 / sigma2;
else
    u = 0;
    v = 1;
end

% Bias correction factor
% Cb = 2 / (v + 1/v + u^2)
Cb = 2 / (v + 1/v + u^2);

%% Confidence Interval via Z-transformation (Lin 1989, Appendix)
% Z = tanh^(-1)(ccc) for better asymptotic normality
if abs(ccc) < 1
    Z = atanh(ccc);  % = 0.5 * log((1+ccc)/(1-ccc))
else
    Z = sign(ccc) * Inf;
end

% Asymptotic variance of Z (Lin 1989, corrected in Lin 2000)
% sigma_Z^2 = (1/(n-2)) * [...complex expression...]
% Simplified approximation for large n:
%   seZ^2 ≈ (1-r^2)*ccc^2 / (n*r^2) + 2*ccc^3*(1-ccc)*u^2 / (r^3*n)
%         + ccc^4*u^4 / (2*r^4*n)
%
% More robust approximation (used in DescTools/epiR):
if abs(r) > eps && n > 2
    r2 = r^2;
    ccc2 = ccc^2;
    
    % Variance components (Lin 1989 Appendix, simplified)
    term1 = (1 - r2) * ccc2 / (1 - ccc2);
    term2 = 4 * ccc^3 * (1 - ccc) * u^2 / r;
    term3 = ccc^4 * u^4 / (2 * r2);
    
    varZ = (1 / (n - 2)) * (term1 + term2 + term3 + 2);
    seZ = sqrt(varZ);
else
    seZ = NaN;
end

% CI on Z scale
z_crit = norminv(1 - alpha/2);
Z_lower = Z - z_crit * seZ;
Z_upper = Z + z_crit * seZ;

% Back-transform to CCC scale
LowerCI = tanh(Z_lower);
UpperCI = tanh(Z_upper);

%% Package results
cccResults.ccc = ccc;
cccResults.r = r;
cccResults.Cb = Cb;
cccResults.LowerCI = LowerCI;
cccResults.UpperCI = UpperCI;
cccResults.n = n;
cccResults.u = u;
cccResults.v = v;
cccResults.Z = Z;
cccResults.seZ = seZ;
cccResults.alpha = alpha;

%% Display if no output requested
if nargout == 0
    fprintf('\nLin''s Concordance Correlation Coefficient\n');
    fprintf('=========================================\n');
    fprintf('CCC (rho_c)    = %.4f [%.4f, %.4f]\n', ccc, LowerCI, UpperCI);
    fprintf('Precision (r)  = %.4f\n', r);
    fprintf('Accuracy (Cb)  = %.4f\n', Cb);
    fprintf('Location shift = %.4f\n', u);
    fprintf('Scale shift    = %.4f\n', v);
    fprintf('N              = %d\n', n);
    fprintf('\nInterpretation (McBride 2005):\n');
    if ccc > 0.99
        fprintf('  Almost perfect agreement\n');
    elseif ccc > 0.95
        fprintf('  Substantial agreement\n');
    elseif ccc > 0.90
        fprintf('  Moderate agreement\n');
    else
        fprintf('  Poor agreement\n');
    end
end

end
