%[text] # Worked Example: Isolation by Distance
%[text] A classic application of the Mantel test, from the field it originated in: testing whether genetic differentiation between populations correlates with geographic distance between them ("isolation by distance" -- populations that are geographically closer tend to be more genetically similar, because gene flow is more likely between nearby populations).
%[text] This example also demonstrates the toolbox at a larger sample size (12 populations) than `doc/GettingStarted.mlx`'s 6-site example, where exact enumeration (12! ~ 479 million) is no longer tractable -- `Exact="auto"` correctly falls back to random-sampling permutation.
%[text] ## Simulated data: 12 populations along a landscape
%[text] Coordinates spread across a landscape; genetic distance constructed to genuinely correlate with geographic distance (simulating real gene flow limited by distance), plus realistic noise.
rng(7);
nPop = 12;
coords = rand(nPop, 2) * 100; % arbitrary spatial units

geoDist = zeros(nPop);
for i = 1:nPop
    geoDist(:, i) = sqrt(sum((coords - coords(i, :)) .^ 2, 2));
end

% Genetic distance: correlates with geographic distance (isolation by
% distance signal) plus substantial noise (genetic drift, sampling).
geneticDist = 0.01 * geoDist + 0.15 * abs(randn(nPop));
geneticDist = (geneticDist + geneticDist') / 2;
geneticDist(logical(eye(nPop))) = 0;
%%
%[text] ## Testing the isolation-by-distance hypothesis
[mantelR, pValue, info] = mantelTest(geoDist, geneticDist);
fprintf("n = %d populations\n", info.n);
fprintf("Exact enumeration used: %d (12! = %d permutations -- intractable, correctly falls back)\n", ...
    info.exact, factorial(12));
fprintf("Permutations actually used: %d\n", info.nPermUsed);
fprintf("Mantel r = %.4f, p = %.4f\n", mantelR, pValue);
%%
%[text] At n=12, exact enumeration (479,001,600 permutations) is well beyond the toolbox's tractable range -- `Exact="auto"` correctly falls back to the default 999 random-sampled permutations, giving a p-value with a minimum resolution of 1/1000. If a p-value close to your significance threshold matters here, increase `NPerm` rather than attempting `Exact="on"` (which would error immediately with a clear message, rather than hang attempting an intractable computation -- see `mantelTest.m`'s docstring for the hard n\<=10 safety limit).
%[text] ## Contrast: a null case (no real spatial signal)
%[text] For comparison, genetic distance with NO relationship to geography -- the pattern you'd see if gene flow were unrestricted by distance (e.g. a highly mobile species, or populations connected by frequent long-distance dispersal).
geneticDistNull = 0.15 * abs(randn(nPop));
geneticDistNull = (geneticDistNull + geneticDistNull') / 2;
geneticDistNull(logical(eye(nPop))) = 0;

[mantelRNull, pValueNull] = mantelTest(geoDist, geneticDistNull);
fprintf("\nNull case: Mantel r = %.4f, p = %.4f\n", mantelRNull, pValueNull);
%%
%[text] The isolation-by-distance case should show a clear positive correlation with a low p-value; the null case should show a weak correlation with a p-value that doesn't support rejecting the null hypothesis of no spatial structure. Comparing the two side by side is the point: a Mantel r alone, without its permutation p-value, doesn't tell you whether an observed correlation is distinguishable from chance restructuring of the same distance values.
%[text] ## Visualising the relationship
figure;
tiledlayout(1, 2);

nexttile;
lowerMask = tril(true(nPop), -1);
scatter(geoDist(lowerMask), geneticDist(lowerMask), 25, "filled");
xlabel("Geographic distance");
ylabel("Genetic distance");
title(sprintf("Isolation by distance\nr=%.3f, p=%.3f", mantelR, pValue));
lsline;

nexttile;
scatter(geoDist(lowerMask), geneticDistNull(lowerMask), 25, "filled");
xlabel("Geographic distance");
ylabel("Genetic distance");
title(sprintf("Null case\nr=%.3f, p=%.3f", mantelRNull, pValueNull));
lsline;
%%
%[text] ## Reporting the result
fprintf("\nExample reportable sentence:\n");
if pValue < 0.05
    verdict = "significantly";
else
    verdict = "not significantly";
end
fprintf("""Genetic distance was %s correlated with geographic distance " + ...
    "(Mantel r = %.3f, p = %.3f, %d permutations, n = %d populations).""\n", ...
    verdict, mantelR, pValue, info.nPermUsed, info.n);
%%
%[text] ## Cross-references
%[text] - `mantelTest.m` docstring
%[text] - `doc/GettingStarted.m` -- basic usage at a smaller, exact-mode n
%[text] - `doc/WhyValidatedAgainstR.m` -- validation against R's vegan::mantel()
%[text] - Classic references for this application: Wright, S. (1943). Isolation by distance. Genetics, 28(2), 114-138. \

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright"}
%---
