%[text] # Why `linCCC` Is Validated Against R's `DescTools::CCC()`
%[text] Lin's Concordance Correlation Coefficient has a well-known formula (Lin, 1989), but the confidence-interval construction and the precision/accuracy decomposition (Cb, u, v) have enough room for subtle implementation differences -- normalisation conventions (N vs N-1), the exact variance-of-Z approximation used for the CI, sign conventions on the location shift -- that independent cross-validation against a trusted reference implementation is worth doing, not just assuming the formula transcription is correct.
%[text] `linCCC` was checked against R's `DescTools::CCC()` (Signorell, 2024), a widely-used, actively-maintained R package specifically for descriptive and agreement statistics, across 7 test cases spanning near-perfect agreement through independence, plus one canonical published dataset. This document shows that comparison and lets you reproduce it yourself if you have R installed.
%[text] ## Validation criterion
%[text] `abs(MATLAB_ccc - R_ccc) < 0.001` for every test case. This is a tight tolerance -- tight enough that agreement within it rules out anything but genuine numerical noise as the source of any residual difference.
%[text] ## The seven test cases
%[text] Six synthetic cases spanning the agreement spectrum, plus the canonical Bland-Altman (1986) Peak Expiratory Flow Rate dataset.
caseNames = ["Near-perfect agreement"; "Location shift (+5)"; ...
    "Scale shift (x1.2)"; "Moderate agreement"; "No agreement (independent)"; ...
    "Larger sample (n=200)"; "Bland-Altman (1986) PEFR"];
matlabCCC = zeros(7, 1);
rCCC = zeros(7, 1);
gap = zeros(7, 1);

rng(42); y1 = randn(50, 1) * 10 + 100; y2 = y1 + randn(50, 1) * 0.01;
res = linCCC(y1, y2, 0.05, CompareR=true);
matlabCCC(1) = res.ccc; rCCC(1) = res.rInfo.ccc; gap(1) = res.rInfo.gap;

rng(43); y1 = randn(50, 1) * 10 + 100; y2 = y1 + 5;
res = linCCC(y1, y2, 0.05, CompareR=true);
matlabCCC(2) = res.ccc; rCCC(2) = res.rInfo.ccc; gap(2) = res.rInfo.gap;

rng(44); y1 = randn(50, 1) * 10 + 100; y2 = y1 * 1.2;
res = linCCC(y1, y2, 0.05, CompareR=true);
matlabCCC(3) = res.ccc; rCCC(3) = res.rInfo.ccc; gap(3) = res.rInfo.gap;

rng(45); y1 = randn(50, 1) * 10 + 100; y2 = y1 + randn(50, 1) * 5;
res = linCCC(y1, y2, 0.05, CompareR=true);
matlabCCC(4) = res.ccc; rCCC(4) = res.rInfo.ccc; gap(4) = res.rInfo.gap;

rng(46); y1 = randn(50, 1); y2 = randn(50, 1);
res = linCCC(y1, y2, 0.05, CompareR=true);
matlabCCC(5) = res.ccc; rCCC(5) = res.rInfo.ccc; gap(5) = res.rInfo.gap;

rng(47); y1 = randn(200, 1) * 10 + 100; y2 = y1 + randn(200, 1) * 2;
res = linCCC(y1, y2, 0.05, CompareR=true);
matlabCCC(6) = res.ccc; rCCC(6) = res.rInfo.ccc; gap(6) = res.rInfo.gap;

% Bland-Altman (1986) PEFR data -- Wright peak flow meter measurements.
% Source: Bland JM, Altman DG (1986). Statistical methods for assessing
% agreement between two methods of clinical measurement. Lancet 327:307-310.
% This exact dataset also appears in DescTools' and epiR's own
% documentation, so it's independently checkable against those packages
% too, not just this toolbox's word for it.
y1 = [494,395,516,434,476,557,413,442,650,433,417,656,267,478,178,423,427]';
y2 = [512,430,520,428,500,600,364,380,658,445,432,626,260,477,259,350,451]';
res = linCCC(y1, y2, 0.05, CompareR=true);
matlabCCC(7) = res.ccc; rCCC(7) = res.rInfo.ccc; gap(7) = res.rInfo.gap;

comparisonTable = table(caseNames, matlabCCC, rCCC, gap, ...
    VariableNames=["Case", "MATLAB_ccc", "R_ccc", "gap"]);
disp(comparisonTable)

fprintf("Max gap across all 7 cases: %.6f (criterion: < 0.001)\n", max(gap));
fprintf("All cases pass: %d\n", all(gap < 0.001));
%%
%[text] If R and the DescTools package are installed on your machine, the table above was computed live, right here, via `linCCC`'s own `CompareR` option -- not from a static historical record. If R isn't available, `CompareR` degrades visibly (a warning, not a silent skip or an error) and the table above would show NaN in the `R_ccc` and `gap` columns instead.
%[text] ## What this means for using this toolbox
%[text] - Trust `linCCC`'s core computation on its own -- it's validated, and  `CompareR` doesn't need to be true every time you use it.
%[text] - `CompareR=true` gives you a live second opinion for any specific dataset you're worried about, at the cost of an R + DescTools dependency and the shell-out overhead of calling R.
%[text] - The Bland-Altman PEFR case is worth knowing by name: if you ever want to sanity-check a **different** CCC implementation (not this one), that dataset and its expected CCC value (~0.9427) is a portable, citable check anyone can run against any implementation. \
%[text] ## Cross-references
%[text] - `linCCC.m` docstring -- states this same validation criterion
%[text] - Lin LI-K (1989). A concordance correlation coefficient to evaluate reproducibility. Biometrics 45(1):255-268.
%[text] - Signorell A (2024). DescTools: Tools for Descriptive Statistics. R package version 0.99.54.
%[text] - Bland JM, Altman DG (1986). Statistical methods for assessing agreement between two methods of clinical measurement. Lancet 327:307-310. \

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright"}
%---
