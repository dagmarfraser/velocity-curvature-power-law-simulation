function cccResults = computeCCC_R(y1, y2, alpha)
% Compute Lin's CCC using R's DescTools package via system call
%
% Uses R's DescTools::CCC() function for validation against MATLAB
% implementation. Requires R with DescTools package installed.
%
% INPUTS:
%   y1    : First measurement vector
%   y2    : Second measurement vector
%   alpha : Confidence level (default 0.05)
%
% OUTPUTS:
%   cccResults: struct with fields
%     .ccc      - Concordance correlation coefficient
%     .r        - Pearson correlation (precision)
%     .Cb       - Bias correction factor (accuracy)
%     .LowerCI  - Lower confidence bound
%     .UpperCI  - Upper confidence bound
%     .n        - Sample size
%     .u        - Location shift
%     .v        - Scale shift
%
% REFERENCE:
%   Lin LI-K (1989). A concordance correlation coefficient to evaluate
%   reproducibility. Biometrics 45(1):255-268.
%
%   Signorell A (2024). DescTools: Tools for Descriptive Statistics.
%   R package version 0.99.54.
%
% SEE ALSO: linCCC_v001, test_R_CCC
%
% VERSION: v001 (2026-01-15)

if nargin < 3
    alpha = 0.05;
end

% Ensure column vectors
y1 = y1(:);
y2 = y2(:);

% Remove NaN pairs
validIdx = ~isnan(y1) & ~isnan(y2);
y1 = y1(validIdx);
y2 = y2(validIdx);

% Write data to temporary CSV
tempData = fullfile(tempdir, 'matlab_ccc_data.csv');
tempResults = fullfile(tempdir, 'r_ccc_results.csv');

% Two columns: y1, y2
dataMatrix = [y1, y2];
writematrix(dataMatrix, tempData);

% R script using DescTools package
% Note: s.shift, l.shift, C.b are atomic vectors; rho.c is a data frame
nl = char(10);
rScript = ['library(DescTools)' nl ...
    sprintf('data <- read.csv("%s", header=FALSE)', strrep(tempData, '\', '/')) nl ...
    'y1 <- data$V1' nl ...
    'y2 <- data$V2' nl ...
    nl ...
    '# Compute CCC using DescTools' nl ...
    sprintf('ccc_result <- CCC(y1, y2, ci="z-transform", conf.level=1-%f)', alpha) nl ...
    nl ...
    '# Extract components - rho.c is data frame, others are atomic' nl ...
    'ccc_val <- as.numeric(ccc_result$rho.c$est)' nl ...
    'lb_val <- as.numeric(ccc_result$rho.c$lwr.ci)' nl ...
    'ub_val <- as.numeric(ccc_result$rho.c$upr.ci)' nl ...
    nl ...
    '# Pearson r calculated separately' nl ...
    'r_pearson <- cor(y1, y2)' nl ...
    nl ...
    '# s.shift, l.shift, C.b are atomic vectors' nl ...
    'cb_val <- as.numeric(ccc_result$C.b)' nl ...
    'u_val <- as.numeric(ccc_result$l.shift)' nl ...
    'v_val <- as.numeric(ccc_result$s.shift)' nl ...
    'n_val <- length(y1)' nl ...
    nl ...
    'result_df <- data.frame(' nl ...
    '  ccc = ccc_val,' nl ...
    '  r = r_pearson,' nl ...
    '  Cb = cb_val,' nl ...
    '  LowerCI = lb_val,' nl ...
    '  UpperCI = ub_val,' nl ...
    '  n = n_val,' nl ...
    '  u = u_val,' nl ...
    '  v = v_val' nl ...
    ')' nl ...
    sprintf('write.csv(result_df, "%s", row.names=FALSE)', strrep(tempResults, '\', '/')) nl];

tempScript = fullfile(tempdir, 'compute_ccc.R');
fid = fopen(tempScript, 'w');
fprintf(fid, '%s', rScript);
fclose(fid);

% Execute R
[status, cmdout] = system(sprintf('Rscript "%s" 2>&1', tempScript));

if status ~= 0
    error('computeCCC_R:RFailed', ...
        'R execution failed:\n%s\n\nInstall DescTools: Rscript -e "install.packages(''DescTools'')"', ...
        cmdout);
end

% Read results
opts = detectImportOptions(tempResults);
resultsTable = readtable(tempResults, opts);

% Extract values
cccResults.ccc = resultsTable.ccc(1);
cccResults.r = resultsTable.r(1);
cccResults.Cb = resultsTable.Cb(1);
cccResults.LowerCI = resultsTable.LowerCI(1);
cccResults.UpperCI = resultsTable.UpperCI(1);
cccResults.n = resultsTable.n(1);
cccResults.u = resultsTable.u(1);
cccResults.v = resultsTable.v(1);
cccResults.alpha = alpha;

% Cleanup
delete(tempData, tempResults, tempScript);

% Display if no output requested
if nargout == 0
    fprintf('R DescTools::CCC() = %.4f [%.4f, %.4f] (n=%d)\n', ...
        cccResults.ccc, cccResults.LowerCI, cccResults.UpperCI, cccResults.n);
end

end
