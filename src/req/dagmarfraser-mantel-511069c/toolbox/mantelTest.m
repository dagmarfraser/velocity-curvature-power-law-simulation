function [mantelR, pValue, info] = mantelTest(D1, D2, options)
%MANTELTEST Mantel test: permutation-based correlation of two distance matrices.
%
% Tests whether two distance/dissimilarity structures on the same n
% items are correlated, using a permutation test for significance
% (since entries within a distance matrix are not independent, standard
% parametric correlation p-values are invalid -- Mantel, 1967).
%
% SYNTAX:
%   [mantelR, pValue] = mantelTest(D1, D2)
%   [mantelR, pValue] = mantelTest(D1, D2, Method="spearman")
%   [mantelR, pValue] = mantelTest(X1, X2, InputType="observations")
%   [mantelR, pValue, info] = mantelTest(D1, D2, CompareR=true)
%
% INPUTS:
%   D1, D2 - Two n-by-n distance/dissimilarity matrices (symmetric,
%            zero diagonal) by default (InputType="distance"), OR two
%            sets of raw observations (InputType="observations") that
%            get converted to distance matrices internally.
%   InputType - Name-value, "distance" (default) | "observations".
%            "observations": D1/D2 are treated as n-by-p matrices (n
%            items, p features) or n-by-1 vectors, converted via
%            pairwise Euclidean distance. For any distance metric other
%            than Euclidean, compute your own distance matrix and use
%            InputType="distance" instead -- this function does not
%            bundle a general multi-metric distance library (that would
%            need Statistics and Machine Learning Toolbox's pdist; this
%            function has no hard dependency on that toolbox).
%   Method - Name-value, "pearson" (default) | "spearman". Correlation
%            type applied to the two matrices' lower-triangle vectors.
%   NPerm  - Name-value, default 999. Number of random permutations used
%            when exact enumeration is not used (see Exact below).
%   Exact  - Name-value, "auto" (default) | "on" | "off". When the
%            number of items n is small enough that all n! label
%            permutations can be enumerated exhaustively in reasonable
%            time, "auto" does so automatically, giving an EXACT p-value
%            instead of a randomly-sampled approximation. Default
%            threshold: n <= 8 (8! = 40320). "on" forces exact
%            enumeration (hard safety limit n <= 10, 10! = 3628800,
%            beyond which it errors rather than attempting an
%            intractable allocation). "off" always uses random sampling
%            (NPerm draws), matching this function's originating
%            implementation's behaviour at every n.
%   CompareR - Name-value logical, default false. If true, also computes
%            the Mantel test via R's vegan::mantel() on the same data
%            for independent cross-validation. See R DEPENDENCY POLICY.
%
% OUTPUTS:
%   mantelR - The Mantel correlation statistic (Pearson or Spearman r
%             between the two lower-triangle vectors). NaN if either
%             distance matrix has zero variance in its lower triangle
%             (a degenerate/constant structure -- not an error, this is
%             the same "correctly undefined" case as a standard
%             correlation between two constants).
%   pValue  - Two-sided permutation p-value, (nGE+1)/(nUsed+1), where
%             nGE counts permutations achieving |statistic| at least as
%             extreme as observed, and nUsed is the number of
%             permutations tested (n!-1 if exact, NPerm if random). NaN
%             when mantelR is NaN.
%   info    - Struct with fields:
%       .n         - Number of items (matrix dimension)
%       .exact     - Logical, whether exact enumeration was used
%       .nPermUsed - Number of permutations actually tested
%       .method    - Echoed Method
%       .rInfo     - Struct, populated only when CompareR=true (see
%                    R DEPENDENCY POLICY); default when not requested:
%                    available=false, message="CompareR not requested"
%
% DEPENDENCY POLICY:
%   No hard dependency on Statistics and Machine Learning Toolbox.
%   Pearson correlation uses corrcoef (base MATLAB); Spearman is
%   implemented via manual average-rank tie handling followed by Pearson
%   on the ranks (avoiding tiedrank, which needs the Statistics
%   Toolbox); the "observations" input path computes Euclidean distance
%   directly (avoiding pdist/squareform, likewise Statistics Toolbox
%   functions). perms (used for exact enumeration) is base MATLAB.
%
% R DEPENDENCY POLICY:
%   CompareR defaults to false, so mantelTest has no hard R dependency.
%   Requesting CompareR=true without R and the vegan package installed
%   does NOT error -- it warns visibly and returns rInfo.available=false
%   with an explanatory message (Fail Loud, Never Fake). Two distinct
%   failure modes are checked and reported separately: R itself absent,
%   versus R present but the vegan package not installed.
%
% VALIDATION:
%   The correlation statistic was cross-validated against R's
%   vegan::mantel() (statistic matched to ~1e-8). A precision
%   improvement over this function's originating implementation was
%   identified during that validation: the source always used random-
%   sampling permutation regardless of n, while vegan automatically
%   switches to exact enumeration when n! is small -- this function's
%   Exact="auto" (default) replicates that behaviour. See
%   doc/WhyValidatedAgainstR.mlx for the full comparison and an
%   explanation of why a p-value from Exact=false will differ slightly
%   from vegan's exact p-value at small n, which is expected, not a bug.
%
% REFERENCES:
%   Mantel, N. (1967). The detection of disease clustering and a
%   generalized regression approach. Cancer Research, 27(2), 209-220.
%   Legendre, P., Fortin, M.-J., & Borcard, D. (2015). Should the Mantel
%   test be used in spatial analysis? Methods in Ecology and Evolution,
%   6(11), 1239-1247.
%   Somers, K. M., & Jackson, D. A. (2022). Putting the Mantel test back
%   together again. Ecology, 103(10), e3780.
%   Oksanen, J. et al. (2024). vegan: Community Ecology Package. R
%   package version 2.6-8. (Reference implementation for CompareR.)
%
% Toolbox: mantel
% Extracted and generalised from: src/constellationMetrics_v002.m
% (mantel_local, a private local function; generalised from a narrow
% 6-scalar-pipeline-value construction to accept arbitrary distance
% matrices or raw observations, matching the standard textbook Mantel
% test interface. R-calling logic for CompareR is original to this
% toolbox, not present in the source.)
% (PowerLawSimulationPreReg, session 52)

arguments
    D1 (:, :) double
    D2 (:, :) double
    options.InputType (1, 1) string {mustBeMember(options.InputType, ["distance", "observations"])} = "distance"
    options.Distance (1, 1) string {mustBeMember(options.Distance, "euclidean")} = "euclidean"
    options.Method (1, 1) string {mustBeMember(options.Method, ["pearson", "spearman"])} = "pearson"
    options.NPerm (1, 1) double {mustBeInteger, mustBePositive} = 999
    options.Exact (1, 1) string {mustBeMember(options.Exact, ["auto", "on", "off"])} = "auto"
    options.CompareR (1, 1) logical = false
end

AUTO_EXACT_CAP_N = 8;
HARD_EXACT_LIMIT_N = 10;

if options.InputType == "observations"
    D1 = mantelPairwiseEuclidean(D1);
    D2 = mantelPairwiseEuclidean(D2);
end

mantelValidateDistanceMatrix(D1, "D1");
mantelValidateDistanceMatrix(D2, "D2");

if ~isequal(size(D1), size(D2))
    error("mantelTest:SizeMismatch", "D1 (%dx%d) and D2 (%dx%d) must be the same size.", ...
        size(D1, 1), size(D1, 2), size(D2, 1), size(D2, 2));
end

n = size(D1, 1);
if n < 3
    error("mantelTest:TooFewItems", "Need at least 3 items (got n=%d).", n);
end

lowerMask = tril(true(n), -1);
v1 = D1(lowerMask);
v2 = D2(lowerMask);

if std(v1) == 0 || std(v2) == 0
    mantelR = NaN;
    pValue = NaN;
    info = struct("n", n, "exact", false, "nPermUsed", 0, "method", options.Method, ...
        "rInfo", struct("available", false, "statistic", NaN, "pValue", NaN, ...
        "statisticGap", NaN, "message", "CompareR not evaluated: degenerate input"));
    return
end

mantelR = mantelCorr(v1, v2, options.Method);

switch options.Exact
    case "on"
        if n > HARD_EXACT_LIMIT_N
            error("mantelTest:ExactTooLarge", "%s", ...
                "Exact=""on"" requested for n=" + n + " (" + n + "! = " + ...
                factorial(n) + " permutations), exceeding the hard safety " + ...
                "limit of n<=" + HARD_EXACT_LIMIT_N + " (" + HARD_EXACT_LIMIT_N + ...
                "! = " + factorial(HARD_EXACT_LIMIT_N) + "). Use Exact=""off"" " + ...
                "for larger n.");
        end
        useExact = true;
    case "off"
        useExact = false;
    otherwise % "auto"
        useExact = n <= AUTO_EXACT_CAP_N;
end

if useExact
    allPerms = perms(1:n);
    isIdentity = all(allPerms == (1:n), 2);
    testPerms = allPerms(~isIdentity, :);
else
    testPerms = zeros(options.NPerm, n);
    for k = 1:options.NPerm
        testPerms(k, :) = randperm(n);
    end
end
nUsed = size(testPerms, 1);

% PERFORMANCE NOTE (added after real usage at n in the hundreds/thousands,
% not just the original n=6 validation case): v1 (the D1 side) never
% changes across permutations -- only D2 is permuted -- so for Spearman
% it is ranked ONCE here rather than being re-ranked on every one of
% potentially thousands of permutations inside the loop. This changes
% nothing about the statistic or its distribution (mantelCorr's rank
% transform is idempotent and depends only on v1 itself), it only avoids
% redundant work. Left as a fully serial optimisation -- no parallelism
% added to this toolbox; project-specific callers needing HPC-scale
% parfor should implement that in their own calling code, not here.
if options.Method == "spearman"
    v1ForPerm = mantelRank(v1);
else
    v1ForPerm = v1;
end

nGE = 0;
for k = 1:nUsed
    idx = testPerms(k, :);
    D2perm = D2(idx, idx);
    v2perm = D2perm(lowerMask);
    if std(v2perm) == 0
        continue
    end
    if options.Method == "spearman"
        v2ForPerm = mantelRank(v2perm);
    else
        v2ForPerm = v2perm;
    end
    C = corrcoef(v1ForPerm, v2ForPerm);
    statPerm = C(1, 2);
    if abs(statPerm) >= abs(mantelR)
        nGE = nGE + 1;
    end
end
pValue = (nGE + 1) / (nUsed + 1);

info.n = n;
info.exact = useExact;
info.nPermUsed = nUsed;
info.method = options.Method;

if options.CompareR
    info.rInfo = mantelCompareR(D1, D2, options.Method, options.NPerm, mantelR, pValue);
else
    info.rInfo = struct("available", false, "statistic", NaN, "pValue", NaN, ...
        "statisticGap", NaN, "message", "CompareR not requested");
end

end

function D = mantelPairwiseEuclidean(X)
%MANTELPAIRWISEEUCLIDEAN Pairwise Euclidean distance, no Statistics Toolbox needed.
if isvector(X)
    X = X(:);
end
n = size(X, 1);
D = zeros(n, n);
for i = 1:n
    diffs = X - X(i, :);
    D(:, i) = sqrt(sum(diffs .^ 2, 2));
end
end

function mantelValidateDistanceMatrix(D, name)
%MANTELVALIDATEDISTANCEMATRIX Fail Loud checks: square, symmetric, zero diagonal.
[nRows, nCols] = size(D);
if nRows ~= nCols
    error("mantelTest:NotSquare", "%s must be square (got %dx%d).", name, nRows, nCols);
end
if ~isequal(D, D') && max(abs(D - D'), [], "all") > 1e-9
    error("mantelTest:NotSymmetric", "%s must be symmetric (within 1e-9 tolerance).", name);
end
if max(abs(diag(D))) > 1e-9
    error("mantelTest:NonzeroDiagonal", "%s must have a zero diagonal (within 1e-9 tolerance).", name);
end
end

function rho = mantelCorr(x, y, method)
%MANTELCORR Pearson or Spearman correlation without Statistics Toolbox.
if method == "spearman"
    x = mantelRank(x);
    y = mantelRank(y);
end
C = corrcoef(x, y);
rho = C(1, 2);
end

function r = mantelRank(x)
%MANTELRANK Average-rank transform, handling ties, no Statistics Toolbox needed.
x = x(:);
n = numel(x);
[sortedX, sortIdx] = sort(x);
% Vectorized average-rank tie correction (patched 2026-07-27, Fraser D.S.):
% the previous per-unique-value loop (tiedMask = ic==k, looped once per
% unique value) is O(n^2) -- fine at small n, but at RDM lower-triangle
% scale (n ~ 9e6 for Pilot's full n=4237 trials) it never completes in
% practical time. This groups by consecutive equal sorted values and
% averages the sorted-order rank position (1:n) within each group via
% accumarray -- mathematically identical output, O(n log n) via sort.
isNewGroup = [true; diff(sortedX) ~= 0];
groupID = cumsum(isNewGroup);
positions = (1:n)';
meanRankPerGroup = accumarray(groupID, positions) ./ accumarray(groupID, 1);
r = zeros(n, 1);
r(sortIdx) = meanRankPerGroup(groupID);
end

function rInfo = mantelCompareR(D1, D2, method, nPerm, statMatlab, pMatlab)
%MANTELCOMPARER Optional R (vegan::mantel) cross-validation panel.
%   Inlined here (not calling an external file) for self-containment,
%   same discipline as linCCC's CompareR.

if isempty(mantelWhichRscript())
    warning("mantelTest:NoR", "%s", ...
        "CompareR=true requested but Rscript not found on the system " + ...
        "path. Returning rInfo.available=false. Install R from " + ...
        "https://www.r-project.org to enable this comparison.");
    rInfo = struct("available", false, "statistic", NaN, "pValue", NaN, ...
        "statisticGap", NaN, "message", "Rscript not found on system path");
    return
end

[veganStatus, ~] = system('Rscript -e "library(vegan)" 2>&1');
if veganStatus ~= 0
    warning("mantelTest:NoVegan", "%s", ...
        "CompareR=true requested and R was found, but the vegan " + ...
        "package is not installed. Returning rInfo.available=false. " + ...
        "Install it with: Rscript -e ""install.packages('vegan')""");
    rInfo = struct("available", false, "statistic", NaN, "pValue", NaN, ...
        "statisticGap", NaN, "message", "R found but vegan package not installed");
    return
end

tempD1 = fullfile(tempdir, "matlab_mantel_d1.csv");
tempD2 = fullfile(tempdir, "matlab_mantel_d2.csv");
tempResults = fullfile(tempdir, "r_mantel_results.csv");
tempScript = fullfile(tempdir, "compute_mantel.R");

writematrix(D1, tempD1);
writematrix(D2, tempD2);

nl = newline;
rScript = ['library(vegan)' nl ...
    sprintf('D1 <- as.dist(as.matrix(read.csv("%s", header=FALSE)))', strrep(tempD1, '\', '/')) nl ...
    sprintf('D2 <- as.dist(as.matrix(read.csv("%s", header=FALSE)))', strrep(tempD2, '\', '/')) nl ...
    sprintf('result <- mantel(D1, D2, method="%s", permutations=%d)', method, nPerm) nl ...
    'result_df <- data.frame(statistic = result$statistic, pvalue = result$signif)' nl ...
    sprintf('write.csv(result_df, "%s", row.names=FALSE)', strrep(tempResults, '\', '/')) nl];

fid = fopen(tempScript, 'w');
fprintf(fid, "%s", rScript);
fclose(fid);

message = "";
statR = NaN;
pR = NaN;
try
    [status, cmdout] = system(sprintf('Rscript "%s" 2>&1', tempScript));
    if status ~= 0
        message = "R execution failed: " + string(cmdout);
    else
        resultsTable = readtable(tempResults);
        statR = resultsTable.statistic(1);
        pR = resultsTable.pvalue(1);
    end
catch ME
    message = "R call failed: " + string(ME.message);
end

if isfile(tempD1); delete(tempD1); end
if isfile(tempD2); delete(tempD2); end
if isfile(tempResults); delete(tempResults); end
if isfile(tempScript); delete(tempScript); end

if isfinite(statR)
    statisticGap = abs(statMatlab - statR);
    if strlength(message) == 0
        message = sprintf( ...
            "Statistic gap %.2e. p-value comparison (MATLAB=%.4f, R=%.4f) may " + ...
            "differ even when both are correct if Exact was not used on the " + ...
            "MATLAB side and n is small enough that vegan used exact " + ...
            "enumeration -- see doc/WhyValidatedAgainstR.mlx.", statisticGap, pMatlab, pR);
    end
else
    statisticGap = NaN;
end

rInfo = struct("available", true, "statistic", statR, "pValue", pR, ...
    "statisticGap", statisticGap, "message", message);

end

function p = mantelWhichRscript()
%MANTELWHICHRSCRIPT Locate Rscript on the system path.
[status, out] = system('which Rscript 2>&1');
if status == 0 && ~isempty(strtrim(out))
    p = strtrim(out);
else
    p = '';
end
end
