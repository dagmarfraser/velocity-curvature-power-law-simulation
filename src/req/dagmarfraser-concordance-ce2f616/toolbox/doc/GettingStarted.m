%[text] # Getting Started with `concordance`
%[text] A single-function toolbox: Lin's Concordance Correlation Coefficient (CCC), a measure of agreement between two continuous variables.
%[text] ## Installation
%[text] Add the toolbox folder to your path:
%[text] ```matlabCodeExample
%[text] addpath('/path/to/concordance/toolbox');
%[text] ```
%[text] ## The basic call
rng(1);
methodA = randn(30, 1) * 10 + 100;
methodB = methodA + randn(30, 1) * 3;  % a second measurement method, close but not identical

result = linCCC(methodA, methodB);
fprintf("CCC = %.4f\n", result.ccc);
%%
%[text] ## Interpreting the output
%[text] `linCCC` returns a struct with everything you need to report and interpret the result:
disp(result)
%%
%[text] - `ccc` -- the concordance correlation coefficient itself, in \[-1, 1\]. This is the number you usually want.
%[text] - `r` -- Pearson correlation (precision: are the two variables  **linearly related**, regardless of whether they agree in absolute terms).
%[text] - `Cb` -- bias correction factor (accuracy: do the two variables fall on the 45-degree identity line, or just some other straight line).  `ccc = r * Cb`.
%[text] - `LowerCI`/|UpperCI| -- 95% confidence interval on `ccc` (or whatever level you specify via the third input argument).
%[text] - `u`, `v` -- location shift and scale shift, the two components that make up `Cb`. Useful for diagnosing **why** Cb is low, not just that it is.
%[text] - `n` -- sample size actually used (after NaN-pair removal). \
%[text] ## Interpreting the value (McBride, 2005)
%[text] - CCC \> 0.99: almost perfect agreement
%[text] - CCC 0.95-0.99: substantial agreement
%[text] - CCC 0.90-0.95: moderate agreement
%[text] - CCC \< 0.90: poor agreement \
%[text] Calling with no output argument prints this interpretation directly:
linCCC(methodA, methodB);
%%
%[text] ## Why CCC, not just Pearson r
%[text] Pearson r only measures whether two variables move together linearly -- it's blind to systematic offsets or scale differences. Two measurement methods can correlate at r=0.99 while one consistently reads 20% higher than the other; Pearson won't flag this, but CCC will, because CCC specifically penalises deviation from the 45-degree identity line, not just linear association.
biasedMethod = methodA * 1.3 + 15;  % same shape, offset and scaled
resultBiased = linCCC(methodA, biasedMethod);
fprintf("Pearson r = %.4f (looks great)\n", resultBiased.r);
fprintf("CCC       = %.4f (correctly flags the disagreement)\n", resultBiased.ccc);
fprintf("Cb        = %.4f (the accuracy penalty driving that gap)\n", resultBiased.Cb);
%%
%[text] ## Optional: cross-validate against R
%[text] If you have R and the DescTools package installed, `CompareR=true` computes the same CCC via R's `DescTools::CCC()` for an independent second opinion on any specific dataset you're checking:
resultChecked = linCCC(methodA, methodB, 0.05, CompareR=true);
if resultChecked.rInfo.available
    fprintf("MATLAB: %.4f, R: %.4f, gap: %.6f\n", ...
        resultChecked.ccc, resultChecked.rInfo.ccc, resultChecked.rInfo.gap);
else
    fprintf("R comparison unavailable: %s\n", resultChecked.rInfo.message);
end
%%
%[text] See `doc/WhyValidatedAgainstR.m` for the full validation story behind this feature -- 7 test cases, including a canonical published dataset, all matching R to within the stated tolerance.
%[text] ## Where to go next
%[text] - `doc/WhyValidatedAgainstR.m` -- the R cross-validation story
%[text] - `examples/AgreementAnalysisWorkedExample.m` -- a full method- comparison worked example
%[text] - `tests/` -- the test suite doubles as executable documentation of expected behaviour \

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright"}
%---
