%[text] # Why `mantelTest` Is Validated Against R's `vegan::mantel()`
%[text] The Mantel test's correlation statistic is simple to implement correctly, but its permutation-based significance test has a subtlety worth checking against a trusted reference: at small sample sizes, the **way** you generate the null distribution (random sampling vs exhaustive enumeration) materially affects the p-value, even when the underlying statistic is computed identically. This document shows both findings from validating `mantelTest` against R's `vegan::mantel()` (Oksanen et al., 2024) -- the standard, most-cited authoritative implementation.
%[text] ## Finding 1: the statistic matches exactly
%[text] `mantelTest`'s originating implementation (before generalisation) was checked against `vegan::mantel()` on a synthetic 6-item example.
bPred = [0.28, 0.31, 0.30, 0.33, 0.29, 0.34];
bObs  = [0.27, 0.33, 0.29, 0.35, 0.28, 0.36];

[mantelR, pValue, info] = mantelTest(bObs(:), bPred(:), InputType="observations", CompareR=true);

fprintf("MATLAB statistic: %.7f\n", mantelR);
if info.rInfo.available
    fprintf("R (vegan) statistic: %.7f\n", info.rInfo.statistic);
    fprintf("Gap: %.2e\n", info.rInfo.statisticGap);
else
    fprintf("R comparison unavailable: %s\n", info.rInfo.message);
end
%%
%[text] The gap is at machine-precision level -- the correlation statistic itself was never in question.
%[text] ## Finding 2: the p-value needs exact enumeration to match exactly
%[text] With only 6 items, there are only 6! = 720 distinct label permutations -- a small enough space to enumerate completely rather than approximate via random sampling. `vegan` detects this automatically and switches to exhaustive enumeration; `mantelTest`'s `Exact="auto"` (the default) replicates that behaviour.
fprintf("\nMATLAB p-value (Exact=auto, triggered since n=6<=8): %.9f\n", pValue);
if info.rInfo.available
    fprintf("R (vegan) p-value: %.9f\n", info.rInfo.pValue);
end
fprintf("Both equal 1/720 = %.9f\n", 1/720);
%%
%[text] Compare this against what happens if exact enumeration is deliberately turned off, forcing random-sampling approximation even at this small n:
rng(1);
[~, pValueRandom, infoRandom] = mantelTest(bObs(:), bPred(:), InputType="observations", ...
    Exact="off", NPerm=999);
fprintf("\nMATLAB p-value (Exact=off, 999 random permutations): %.9f\n", pValueRandom);
fprintf("Exact used: %d, nPermUsed: %d\n", infoRandom.exact, infoRandom.nPermUsed);
%%
%[text] This is NOT a bug -- both p-values are valid permutation-test estimates. The random-sampling version approximates the same null distribution the exact version characterises completely; at small n the two can differ by a non-trivial amount simply because 999 random draws (with replacement, from a 720-permutation space) won't necessarily hit exactly the same extreme cases exhaustive enumeration guarantees to find. If you ever compare `mantelTest` output against `vegan` (or any other exact-capable implementation) and see a p-value mismatch at small n, check whether `Exact` was actually triggered on both sides before suspecting a bug.
%[text] ## Why this matters practically
%[text] - At small n (\<=8 by default), always prefer `Exact="auto"` (the default) -- there's no reason to accept random-sampling approximation error when the exact answer is computable in milliseconds.
%[text] - At larger n, exhaustive enumeration becomes intractable (n! grows extremely fast: 10! = 3,628,800, 12! = 479,001,600) -- random sampling is the only practical option, and `NPerm=999` (the default, matching common practice) gives a p-value with a minimum resolution of 1/1000.
%[text] - If precise p-values near a significance threshold matter at larger n, increase `NPerm` rather than forcing `Exact="on"` beyond its hard safety limit (n\<=10). \
%[text] ## Cross-references
%[text] - `mantelTest.m` docstring -- states this same distinction
%[text] - Mantel, N. (1967). The detection of disease clustering and a generalized regression approach. Cancer Research, 27(2), 209-220.
%[text] - Oksanen, J. et al. (2024). vegan: Community Ecology Package. R package version 2.6-8.
%[text] - Legendre, P., Fortin, M.-J., & Borcard, D. (2015). Should the Mantel test be used in spatial analysis? Methods in Ecology and Evolution, 6(11), 1239-1247. \

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright"}
%---
