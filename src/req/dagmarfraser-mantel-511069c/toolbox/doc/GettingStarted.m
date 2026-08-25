%[text] # Getting Started with `mantel`
%[text] A single-function toolbox: the Mantel test, a permutation-based correlation test for comparing two distance/dissimilarity matrices on the same set of items.
%[text] ## Installation
%[text] ```matlabCodeExample
%[text] addpath('/path/to/mantel/toolbox');
%[text] ```
%[text] ## Why you can't just correlate two distance matrices normally
%[text] The entries within a single distance matrix are not independent -- each item participates in n-1 pairwise distances, so ordinary correlation significance tests (which assume independent observations) give invalid p-values if applied directly to distance-matrix entries. The Mantel test solves this with a permutation approach: shuffle which item is which, recompute, and see how often a random relabelling produces as strong a correlation as the one actually observed.
%[text] ## Basic usage: two distance matrices
%[text] Six sites, distances derived from two different measurements (e.g. geographic distance and community dissimilarity in an ecology context).
rng(1);
sites = 6;
geoDist = [0 12 8 25 30 18; 12 0 6 20 22 15; 8 6 0 18 24 10; ...
    25 20 18 0 15 22; 30 22 24 15 0 20; 18 15 10 22 20 0];
communityDiss = geoDist * 0.03 + 0.05 * abs(randn(sites)); % correlated with distance, plus noise
communityDiss = (communityDiss + communityDiss') / 2;
communityDiss(logical(eye(sites))) = 0;

[mantelR, pValue, info] = mantelTest(geoDist, communityDiss);
fprintf("Mantel r = %.4f, p = %.4f (n=%d, exact=%d)\n", mantelR, pValue, info.n, info.exact);
%%
%[text] ## Basic usage: raw observations instead of a pre-built distance matrix
%[text] If you have raw feature data rather than a pre-computed distance matrix, use `InputType="observations"` -- Euclidean distance is computed internally.
featuresA = randn(6, 3);
featuresB = featuresA + randn(6, 3) * 0.3;

[mantelR2, pValue2] = mantelTest(featuresA, featuresB, InputType="observations");
fprintf("Mantel r = %.4f, p = %.4f\n", mantelR2, pValue2);
%%
%[text] For any distance metric other than Euclidean, compute your own distance matrix (with your metric of choice) and pass it directly instead -- `InputType="observations"` deliberately only covers the Euclidean case, to keep this toolbox free of a hard Statistics and Machine Learning Toolbox dependency.
%[text] ## The `Exact` option
%[text] At small sample sizes, the exact p-value (from exhaustively enumerating every possible relabelling) is affordable to compute and strictly more precise than a random-sampling approximation. `mantelTest` does this automatically when n\<=8:
disp(info)
%%
%[text] `info.exact` tells you whether exact enumeration was used; `info.nPermUsed` tells you how many permutations were actually tested (n!-1 if exact, `NPerm` if not). See `doc/WhyValidatedAgainstR.m` for why this matters and how it was validated against R's `vegan::mantel()`.
%[text] ## Optional: cross-validate against R
%[text] If you have R and the vegan package installed, `CompareR=true` computes the same test via `vegan::mantel()` for an independent second opinion:
[~, ~, infoChecked] = mantelTest(geoDist, communityDiss, CompareR=true);
if infoChecked.rInfo.available
    fprintf("R agrees: statistic gap = %.2e\n", infoChecked.rInfo.statisticGap);
else
    fprintf("R comparison unavailable: %s\n", infoChecked.rInfo.message);
end
%%
%[text] ## Where to go next
%[text] - `doc/WhyValidatedAgainstR.m` -- the R cross-validation story and the exact-vs-random p-value distinction
%[text] - `examples/` -- a worked example
%[text] - `tests/` -- the test suite doubles as executable documentation \

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright"}
%---
