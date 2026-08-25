function cccResults = linCCC(y1, y2, alpha, options)
%LINCCC Lin's Concordance Correlation Coefficient (CCC).
%
% Measures agreement between two continuous variables by evaluating
% deviation from the 45-degree line (perfect concordance). Combines
% precision (Pearson's r) and accuracy (bias correction factor Cb).
%
% SYNTAX:
%   cccResults = linCCC(y1, y2)
%   cccResults = linCCC(y1, y2, alpha)
%   cccResults = linCCC(y1, y2, alpha, CompareR=true)
%
% INPUTS:
%   y1    - First measurement vector (n x 1)
%   y2    - Second measurement vector (n x 1), same length as y1
%   alpha - Significance level for the confidence interval, in (0, 1)
%           exclusive. Default 0.05. NOTE: this is the standard
%           statistical significance-level "alpha" (Lin 1989's own
%           notation) -- it has no connection whatsoever to the
%           "alpha" (spectral exponent) terminology used in the
%           unrelated fractalnoise toolbox, if you have both installed.
%   CompareR - Name-value, default false. If true, also computes CCC via
%           R's DescTools::CCC() (Signorell, 2024) on the same data for
%           independent cross-validation, returned in rInfo. See
%           R DEPENDENCY POLICY below.
%
% OUTPUTS:
%   cccResults : struct with fields
%     .ccc      - Concordance correlation coefficient (rho_c)
%     .r        - Pearson correlation coefficient (precision)
%     .Cb       - Bias correction factor (accuracy)
%     .LowerCI  - Lower confidence bound (Z-transformed)
%     .UpperCI  - Upper confidence bound (Z-transformed)
%     .n        - Sample size (after NaN-pair removal)
%     .u        - Location shift (standardised mean difference)
%     .v        - Scale shift (SD ratio)
%     .Z        - Z-transformed CCC (for asymptotic inference)
%     .seZ      - Standard error of Z
%     .alpha    - Significance level used
%     .rInfo    - Struct, populated only when CompareR=true:
%         .available - logical, whether R + DescTools were found
%         .ccc       - R's computed CCC (NaN if unavailable)
%         .gap       - abs(cccResults.ccc - rInfo.ccc) (NaN if unavailable)
%         .message   - string, empty on success; explains why
%                      available=false otherwise
%         When CompareR=false (default), rInfo has available=false and
%         message="CompareR not requested".
%
% INTERPRETATION (McBride, 2005):
%   CCC > 0.99   : Almost perfect agreement
%   CCC 0.95-0.99: Substantial agreement
%   CCC 0.90-0.95: Moderate agreement
%   CCC < 0.90   : Poor agreement
%
% FORMULA:
%   rho_c = 2*sigma_12 / (sigma_1^2 + sigma_2^2 + (mu_1 - mu_2)^2)
%         = rho * Cb
%   where Cb = 2 / (v + 1/v + u^2)
%         v  = sigma_1 / sigma_2  (scale shift)
%         u  = (mu_1 - mu_2) / sqrt(sigma_1 * sigma_2)  (location shift)
%
% R DEPENDENCY POLICY:
%   CompareR defaults to false, so linCCC has no hard R dependency for
%   its core (homebrew) computation. Requesting CompareR=true without R
%   and the DescTools package installed does NOT error -- it warns
%   visibly and returns rInfo.available=false with an explanatory
%   message (Fail Loud, Never Fake: the degraded mode is disclosed, not
%   silently skipped). Two distinct failure modes are checked and
%   reported separately: R itself absent (no Rscript on the system
%   path), versus R present but the DescTools package not installed.
%
% VALIDATION:
%   This implementation was cross-validated against R's DescTools::CCC()
%   across 7 test cases (|MATLAB - R| < 0.001 for CCC in every case),
%   including the canonical Bland-Altman (1986) Peak Expiratory Flow
%   Rate dataset -- the same worked example used in DescTools' and
%   epiR's own documentation. See doc/WhyValidatedAgainstR.mlx for the
%   full comparison.
%
% REFERENCES:
%   Lin LI-K (1989). A concordance correlation coefficient to evaluate
%   reproducibility. Biometrics 45(1):255-268.
%   Lin LI-K (2000). A note on the concordance correlation coefficient.
%   Biometrics 56(1):324-325. [Correction to variance formula]
%   McBride GB (2005). A proposal for strength-of-agreement criteria
%   for Lin's concordance correlation coefficient. NIWA Client Report.
%   Bland JM, Altman DG (1986). Statistical methods for assessing
%   agreement between two methods of clinical measurement. Lancet
%   327:307-310.
%   Signorell A (2024). DescTools: Tools for Descriptive Statistics.
%   R package version 0.99.54.
%
% Toolbox: concordance
% Extracted from: src/functions/linCCC_v001.m (core computation
% unchanged; CompareR gated R cross-validation inlined as a local
% function below rather than depending on src/functions/computeCCC_R.m
% existing outside this toolbox's own tree)
% (PowerLawSimulationPreReg, session 52)

arguments
    y1 (:, 1) double
    y2 (:, 1) double
    alpha (1, 1) double {mustBeInRange(alpha, 0, 1, "exclusive")} = 0.05
    options.CompareR (1, 1) logical = false
end

if length(y1) ~= length(y2)
    error("linCCC:UnequalLength", "y1 and y2 must have the same length (got %d and %d).", ...
        length(y1), length(y2));
end

% Remove NaN pairs (pairwise deletion)
validIdx = ~isnan(y1) & ~isnan(y2);
y1 = y1(validIdx);
y2 = y2(validIdx);
n = length(y1);

if n < 3
    error("linCCC:InsufficientData", "Need at least 3 valid pairs, got %d.", n);
end

%% Compute sample statistics
mu1 = mean(y1);
mu2 = mean(y2);

% Population variances (1/n normalisation per Lin 1989)
var1 = var(y1, 1);
var2 = var(y2, 1);
sigma1 = sqrt(var1);
sigma2 = sqrt(var2);

sigma12 = sum((y1 - mu1) .* (y2 - mu2)) / n;

%% Concordance Correlation Coefficient (Lin 1989, Eq. 1)
numerator = 2 * sigma12;
denominator = var1 + var2 + (mu1 - mu2)^2;
ccc = numerator / denominator;

%% Decomposition into precision (r) and accuracy (Cb)
r = sigma12 / (sigma1 * sigma2);

if sigma1 > 0 && sigma2 > 0
    u = (mu1 - mu2) / sqrt(sigma1 * sigma2);
    v = sigma1 / sigma2;
else
    u = 0;
    v = 1;
end

Cb = 2 / (v + 1 / v + u^2);

%% Confidence interval via Z-transformation (Lin 1989 Appendix, Lin 2000 correction)
if abs(ccc) < 1
    Z = atanh(ccc);
else
    Z = sign(ccc) * Inf;
end

if abs(r) > eps && n > 2
    r2 = r^2;
    ccc2 = ccc^2;
    term1 = (1 - r2) * ccc2 / (1 - ccc2);
    term2 = 4 * ccc^3 * (1 - ccc) * u^2 / r;
    term3 = ccc^4 * u^4 / (2 * r2);
    varZ = (1 / (n - 2)) * (term1 + term2 + term3 + 2);
    seZ = sqrt(varZ);
else
    seZ = NaN;
end

zCrit = norminv(1 - alpha / 2);
zLower = Z - zCrit * seZ;
zUpper = Z + zCrit * seZ;
LowerCI = tanh(zLower);
UpperCI = tanh(zUpper);

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

if options.CompareR
    cccResults.rInfo = linCCCCompareR(y1, y2, alpha, ccc);
else
    cccResults.rInfo = struct("available", false, "ccc", NaN, "gap", NaN, ...
        "message", "CompareR not requested");
end

%% Display if no output requested
if nargout == 0
    fprintf("\nLin's Concordance Correlation Coefficient\n");
    fprintf("=========================================\n");
    fprintf("CCC (rho_c)    = %.4f [%.4f, %.4f]\n", ccc, LowerCI, UpperCI);
    fprintf("Precision (r)  = %.4f\n", r);
    fprintf("Accuracy (Cb)  = %.4f\n", Cb);
    fprintf("Location shift = %.4f\n", u);
    fprintf("Scale shift    = %.4f\n", v);
    fprintf("N              = %d\n", n);
    fprintf("\nInterpretation (McBride 2005):\n");
    if ccc > 0.99
        fprintf("  Almost perfect agreement\n");
    elseif ccc > 0.95
        fprintf("  Substantial agreement\n");
    elseif ccc > 0.90
        fprintf("  Moderate agreement\n");
    else
        fprintf("  Poor agreement\n");
    end
    clear cccResults
end

end

function rInfo = linCCCCompareR(y1, y2, alpha, cccMatlab)
%LINCCCCOMPARER Optional R (DescTools::CCC) cross-validation panel.
%   Two distinct failure modes are checked and reported separately:
%   R itself absent, versus R present but DescTools not installed.
%   Inlined here (not calling an external computeCCC_R.m) so this
%   toolbox has no dependency on a file outside its own tree -- see
%   README_MakingToolboxes_v001.md's lesson on this exact class of bug.

if isempty(which_rscript())
    warning("linCCC:NoR", "%s", ...
        "CompareR=true requested but Rscript not found on the system " + ...
        "path. Returning rInfo.available=false. Install R from " + ...
        "https://www.r-project.org to enable this comparison.");
    rInfo = struct("available", false, "ccc", NaN, "gap", NaN, ...
        "message", "Rscript not found on system path");
    return
end

[descToolsStatus, ~] = system('Rscript -e "library(DescTools)" 2>&1');
if descToolsStatus ~= 0
    warning("linCCC:NoDescTools", "%s", ...
        "CompareR=true requested and R was found, but the DescTools " + ...
        "package is not installed. Returning rInfo.available=false. " + ...
        "Install it with: Rscript -e ""install.packages('DescTools')""");
    rInfo = struct("available", false, "ccc", NaN, "gap", NaN, ...
        "message", "R found but DescTools package not installed");
    return
end

tempData = fullfile(tempdir, "matlab_linccc_data.csv");
tempResults = fullfile(tempdir, "r_linccc_results.csv");
tempScript = fullfile(tempdir, "compute_linccc.R");

writematrix([y1, y2], tempData);

nl = newline;
rScript = ['library(DescTools)' nl ...
    sprintf('data <- read.csv("%s", header=FALSE)', strrep(tempData, '\', '/')) nl ...
    'y1 <- data$V1' nl ...
    'y2 <- data$V2' nl ...
    sprintf('ccc_result <- CCC(y1, y2, ci="z-transform", conf.level=1-%f)', alpha) nl ...
    'ccc_val <- as.numeric(ccc_result$rho.c$est)' nl ...
    'result_df <- data.frame(ccc = ccc_val)' nl ...
    sprintf('write.csv(result_df, "%s", row.names=FALSE)', strrep(tempResults, '\', '/')) nl];

fid = fopen(tempScript, 'w');
fprintf(fid, "%s", rScript);
fclose(fid);

message = "";
cccR = NaN;
try
    [status, cmdout] = system(sprintf('Rscript "%s" 2>&1', tempScript));
    if status ~= 0
        message = "R execution failed: " + string(cmdout);
    else
        resultsTable = readtable(tempResults);
        cccR = resultsTable.ccc(1);
    end
catch ME
    message = "R call failed: " + string(ME.message);
end

if isfile(tempData); delete(tempData); end
if isfile(tempResults); delete(tempResults); end
if isfile(tempScript); delete(tempScript); end

if isfinite(cccR)
    gap = abs(cccMatlab - cccR);
else
    gap = NaN;
end

rInfo = struct("available", true, "ccc", cccR, "gap", gap, "message", message);

end

function p = which_rscript()
%WHICH_RSCRIPT Locate Rscript on the system path.
%   which() alone does not search the system PATH for external binaries
%   (only MATLAB's own path for .m/.mex files), so this uses the shell
%   directly, mirroring estimateIRASA's which('ft_defaults') gate
%   pattern but adapted for a non-MATLAB external dependency.

[status, out] = system('which Rscript 2>&1');
if status == 0 && ~isempty(strtrim(out))
    p = strtrim(out);
else
    p = '';
end

end
