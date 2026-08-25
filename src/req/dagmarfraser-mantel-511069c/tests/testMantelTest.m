classdef testMantelTest < matlab.unittest.TestCase
%TESTMANTELTEST Unit tests for mantelTest.m.
%
% Toolbox: mantel
% Covers: source-equivalence (mantel_local, reconstructed inline since
% it's a private local function), both input forms (distance matrix,
% observations), both correlation methods, the Exact auto/on/off modes,
% error paths, and both CompareR branches.

    methods (TestClassSetup)
        function addToolboxPath(testCase)
            toolboxDir = fullfile(fileparts(mfilename("fullpath")), "..", "toolbox");
            addpath(toolboxDir);
            testCase.addTeardown(@() rmpath(toolboxDir));
        end
    end

    methods (Test)
        function testMatchesSourceStatistic(testCase)
            % Reconstruct source mantel_local's construction inline
            % (private local function in constellationMetrics_v002.m,
            % not directly callable).
            bPred = [0.28, 0.31, 0.30, 0.33, 0.29, 0.34];
            bObs  = [0.27, 0.33, 0.29, 0.35, 0.28, 0.36];
            n = 6;
            pairP = zeros(1, n * (n - 1) / 2);
            pairQ = zeros(1, n * (n - 1) / 2);
            k = 0;
            for pp = 2:n
                for qq = 1:pp - 1
                    k = k + 1;
                    pairP(k) = pp;
                    pairQ(k) = qq;
                end
            end
            dPred = abs(bPred(pairP) - bPred(pairQ));
            dObs = abs(bObs(pairP) - bObs(pairQ));
            oC = dObs - mean(dObs);
            oN = sqrt(sum(oC .^ 2));
            pC = dPred - mean(dPred);
            rSource = sum(pC .* oC) / (sqrt(sum(pC .^ 2)) * oN);

            rNew = mantelTest(bObs(:), bPred(:), InputType="observations");
            testCase.verifyEqual(rNew, rSource, "AbsTol", 1e-9);
        end

        function testExactPValueMatchesKnownVeganResult(testCase)
            % Known result from vegan::mantel() on this exact data,
            % validated during development (session 52): statistic
            % 0.9073134, p-value 0.001388889 (=1/720, exact enumeration
            % over 6!=720 label permutations).
            bPred = [0.28, 0.31, 0.30, 0.33, 0.29, 0.34];
            bObs  = [0.27, 0.33, 0.29, 0.35, 0.28, 0.36];
            [r, p, info] = mantelTest(bObs(:), bPred(:), InputType="observations");
            testCase.verifyEqual(r, 0.9073134, "AbsTol", 1e-6);
            testCase.verifyEqual(p, 1 / 720, "AbsTol", 1e-9);
            testCase.verifyTrue(info.exact);
            testCase.verifyEqual(info.nPermUsed, 719);
        end

        function testDistanceMatrixInputPath(testCase)
            D1 = [0 1 2 3 4 5; 1 0 1 2 3 4; 2 1 0 1 2 3; ...
                3 2 1 0 1 2; 4 3 2 1 0 1; 5 4 3 2 1 0];
            [r, p, info] = mantelTest(D1, D1);
            testCase.verifyEqual(r, 1.0, "AbsTol", 1e-9, ...
                "A distance matrix correlated with itself should give r=1");
            testCase.verifyLessThan(p, 0.05);
            testCase.verifyEqual(info.n, 6);
        end

        function testSpearmanMethod(testCase)
            D1 = [0 1 2 3 4 5; 1 0 1 2 3 4; 2 1 0 1 2 3; ...
                3 2 1 0 1 2; 4 3 2 1 0 1; 5 4 3 2 1 0];
            [r, ~, info] = mantelTest(D1, D1, Method="spearman");
            testCase.verifyEqual(r, 1.0, "AbsTol", 1e-9);
            testCase.verifyEqual(info.method, "spearman");
        end

        function testExactOffForcesRandomSampling(testCase)
            D1 = [0 1 2 3 4 5; 1 0 1 2 3 4; 2 1 0 1 2 3; ...
                3 2 1 0 1 2; 4 3 2 1 0 1; 5 4 3 2 1 0];
            rng(1);
            [r, ~, info] = mantelTest(D1, D1, Exact="off", NPerm=200);
            testCase.verifyFalse(info.exact);
            testCase.verifyEqual(info.nPermUsed, 200);
            testCase.verifyEqual(r, 1.0, "AbsTol", 1e-9, ...
                "Statistic should be identical regardless of Exact mode -- only the p-value method changes");
        end

        function testExactAutoThresholdBehaviour(testCase)
            % n<=8 (AUTO_EXACT_CAP_N) triggers exact; n>8 falls back to random
            X8 = rand(8, 2);
            [~, ~, info8] = mantelTest(X8, X8, InputType="observations");
            testCase.verifyTrue(info8.exact);
            testCase.verifyEqual(info8.nPermUsed, factorial(8) - 1);

            X9 = rand(9, 2);
            [~, ~, info9] = mantelTest(X9, X9, InputType="observations", NPerm=100);
            testCase.verifyFalse(info9.exact);
            testCase.verifyEqual(info9.nPermUsed, 100);
        end

        function testExactOnHardLimitRejectsLargeN(testCase)
            X11 = rand(11, 2);
            testCase.verifyError(@() mantelTest(X11, X11, InputType="observations", Exact="on"), ...
                "mantelTest:ExactTooLarge");

            X10 = rand(10, 2);
            [~, ~, info10] = mantelTest(X10, X10, InputType="observations", Exact="on");
            testCase.verifyTrue(info10.exact);
            testCase.verifyEqual(info10.nPermUsed, factorial(10) - 1);
        end

        function testValidationErrorPaths(testCase)
            testCase.verifyError(@() mantelTest([1 2; 2 1], [1 2; 2 1]), "mantelTest:NonzeroDiagonal");
            testCase.verifyError(@() mantelTest(ones(2, 3), ones(2, 3)), "mantelTest:NotSquare");
            testCase.verifyError(@() mantelTest([0 1; 2 0], [0 1; 1 0]), "mantelTest:NotSymmetric");
            D1 = zeros(6); D1(1, 2) = 1; D1(2, 1) = 1;
            testCase.verifyError(@() mantelTest(D1, D1(1:5, 1:5)), "mantelTest:SizeMismatch");
            testCase.verifyError(@() mantelTest(zeros(2), zeros(2)), "mantelTest:TooFewItems");
        end

        function testDegenerateInputReturnsNaNNotError(testCase)
            D1 = zeros(4); % all-zero distance matrix, degenerate (zero variance)
            [r, p, info] = mantelTest(D1, D1);
            testCase.verifyTrue(isnan(r));
            testCase.verifyTrue(isnan(p));
            testCase.verifyFalse(info.rInfo.available);
        end

        function testDefaultRInfoWhenNotRequested(testCase)
            D1 = [0 1 2; 1 0 1; 2 1 0];
            [~, ~, info] = mantelTest(D1, D1);
            testCase.verifyFalse(info.rInfo.available);
            testCase.verifyEqual(info.rInfo.message, "CompareR not requested");
        end

        function testCompareRSuccessPathWhenAvailable(testCase)
            testCase.assumeTrue(~isempty(mantelTestWhichRscript()), ...
                "Rscript not on path on this machine; success-path test not applicable.");
            bPred = [0.28, 0.31, 0.30, 0.33, 0.29, 0.34];
            bObs  = [0.27, 0.33, 0.29, 0.35, 0.28, 0.36];
            [r, ~, info] = testCase.verifyWarningFree(...
                @() mantelTest(bObs(:), bPred(:), InputType="observations", CompareR=true));
            testCase.verifyTrue(info.rInfo.available);
            testCase.verifyEqual(r, info.rInfo.statistic, "AbsTol", 1e-6);
            testCase.verifyLessThan(info.rInfo.statisticGap, 1e-6);
        end

        function testCompareRDegradesGracefullyWhenRscriptUnavailable(testCase)
            originalPath = getenv("PATH");
            setenv("PATH", "/usr/bin:/bin");
            testCase.addTeardown(@() setenv("PATH", originalPath));

            testCase.assumeTrue(isempty(mantelTestWhichRscript()), ...
                "Rscript still resolvable with a stripped PATH on this machine; " + ...
                "degraded-mode simulation not applicable here.");

            D1 = [0 1 2; 1 0 1; 2 1 0];
            testCase.verifyWarning(@() mantelTest(D1, D1, CompareR=true), "mantelTest:NoR");
            [r, ~, info] = mantelTest(D1, D1, CompareR=true);
            testCase.verifyFalse(info.rInfo.available);
            testCase.verifyTrue(isfinite(r), ...
                "Core statistic must succeed regardless of CompareR availability");
        end

        function testSpearmanTiesMatchesRAverageRank(testCase)
            % Validates the fixed (O(n^2)->O(n log n)) mantelRank tie-averaging
            % against R's vegan::mantel(method="spearman"), which uses R's own
            % rank(ties.method="average") internally -- the exact behaviour
            % mantelRank is designed to replicate. Values are deliberately
            % coarse so the pairwise-distance lower triangle is tie-rich
            % (repeated |xi-xj| values), not just theoretically tie-capable.
            testCase.assumeTrue(~isempty(mantelTestWhichRscript()), ...
                "Rscript not on path on this machine; R cross-validation not applicable.");

            rng(42);
            n = 30;
            x = round(rand(n, 1) * 4);
            y = round(rand(n, 1) * 4);
            D1 = abs(x - x.');
            D2 = abs(y - y.');

            [rMatlab, ~, info] = mantelTest(D1, D2, Method="spearman", NPerm=5, CompareR=true);
            testCase.assumeTrue(info.rInfo.available, "CompareR reported unavailable despite Rscript check passing.");
            testCase.verifyEqual(rMatlab, info.rInfo.statistic, "AbsTol", 1e-9, ...
                "mantelRank's average-rank tie handling must match R's rank(ties.method=""average"") exactly.");
        end

        function testSpearmanTiesMatchesRAtModerateScale(testCase)
            % Same validation at n=300 (~44,850 pairs) -- large enough that the
            % old per-unique-value loop would be visibly slow, though far short
            % of the ~9e6-element RDM scale that motivated the fix (kept moderate
            % here to bound test-suite runtime; the full-scale case is exercised
            % operationally in PowerLawSimulationPreReg's RSA pipeline).
            testCase.assumeTrue(~isempty(mantelTestWhichRscript()), ...
                "Rscript not on path on this machine; R cross-validation not applicable.");

            rng(7);
            n = 300;
            x = round(rand(n, 1) * 20);
            y = round(rand(n, 1) * 20);
            D1 = abs(x - x.');
            D2 = abs(y - y.');

            [rMatlab, ~, info] = mantelTest(D1, D2, Method="spearman", NPerm=5, CompareR=true);
            testCase.assumeTrue(info.rInfo.available, "CompareR reported unavailable despite Rscript check passing.");
            testCase.verifyEqual(rMatlab, info.rInfo.statistic, "AbsTol", 1e-9);
        end
    end
end

function p = mantelTestWhichRscript()
%MANTELTESTWHICHRSCRIPT Local copy of the Rscript-detection logic, used
%   only to decide whether to run/skip R-dependent tests. Deliberately
%   duplicated (self-containment discipline), not calling into
%   mantelTest's own private local function.
[status, out] = system('which Rscript 2>&1');
if status == 0 && ~isempty(strtrim(out))
    p = strtrim(out);
else
    p = '';
end
end
