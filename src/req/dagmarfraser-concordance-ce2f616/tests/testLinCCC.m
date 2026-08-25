classdef testLinCCC < matlab.unittest.TestCase
%TESTLINCCC Unit tests for linCCC.m.
%
% Toolbox: concordance
% Covers: source-equivalence (linCCC_v001), the canonical Bland-Altman
% (1986) PEFR dataset, arguments-block validation, error paths, and
% both CompareR branches (available and degraded-mode).

    properties (Constant)
        % No hardcoded personal path here (see CHANGELOG.md 2026-07-20 fix) --
        % set CONCORDANCE_LEGACY_DIR as an environment variable if you want
        % to opt in to the source-equivalence check against the original
        % linCCC_v001.m on your own machine. Unset (the case for every
        % File Exchange/GitHub user, and previously silently true on every
        % machine here too, since the old hardcoded path was already
        % stale) -- the test below filters cleanly, by design.
    end

    methods (TestClassSetup)
        function addToolboxAndLegacyPaths(testCase)
            toolboxDir = fullfile(fileparts(mfilename("fullpath")), "..", "toolbox");
            addpath(toolboxDir);
            testCase.addTeardown(@() rmpath(toolboxDir));

            legacyDir = string(getenv("CONCORDANCE_LEGACY_DIR"));
            if strlength(legacyDir) > 0 && isfolder(legacyDir)
                addpath(legacyDir);
                testCase.addTeardown(@() rmpath(legacyDir));
            end
        end
    end

    methods (Test)
        function testMatchesSourceForSameData(testCase)
            testCase.assumeTrue(exist("linCCC_v001", "file") == 2, ...
                "linCCC_v001 not on path; skipping source-equivalence test.");
            rng(1);
            y1 = randn(50, 1) * 10 + 100;
            y2 = y1 + randn(50, 1) * 3;

            resOld = linCCC_v001(y1, y2);
            resNew = linCCC(y1, y2);

            fields = ["ccc", "r", "Cb", "LowerCI", "UpperCI", "n", "u", "v", "Z", "seZ", "alpha"];
            for f = fields
                testCase.verifyEqual(resNew.(f), resOld.(f), ...
                    sprintf("Field %s did not match source linCCC_v001", f));
            end
        end

        function testPEFRCanonicalDataset(testCase)
            % Bland-Altman (1986) PEFR data, also used in DescTools'/epiR's
            % own documentation. Known validated value from test_R_CCC.m.
            y1 = [494,395,516,434,476,557,413,442,650,433,417,656,267,478,178,423,427]';
            y2 = [512,430,520,428,500,600,364,380,658,445,432,626,260,477,259,350,451]';
            result = linCCC(y1, y2);
            testCase.verifyEqual(result.ccc, 0.942742, "AbsTol", 1e-6);
        end

        function testAlphaOutOfRangeErrors(testCase)
            y1 = randn(20, 1); y2 = y1 + randn(20, 1) * 0.1;
            testCase.verifyError(@() linCCC(y1, y2, 0), "MATLAB:validators:mustBeInRange");
            testCase.verifyError(@() linCCC(y1, y2, 1), "MATLAB:validators:mustBeInRange");
            testCase.verifyError(@() linCCC(y1, y2, 1.5), "MATLAB:validators:mustBeInRange");
        end

        function testUnequalLengthErrorsWithClearMessage(testCase)
            testCase.verifyError(@() linCCC([1; 2; 3], [1; 2]), "linCCC:UnequalLength");
        end

        function testInsufficientDataErrors(testCase)
            testCase.verifyError(@() linCCC([1; 2], [1; 2]), "linCCC:InsufficientData");
        end

        function testNaNPairwiseDeletion(testCase)
            rng(5);
            y1 = randn(20, 1);
            y2 = y1 + randn(20, 1) * 0.1;
            y1(3) = NaN;
            y2(7) = NaN;
            result = linCCC(y1, y2);
            testCase.verifyEqual(result.n, 18);
        end

        function testDefaultRInfoWhenNotRequested(testCase)
            y1 = randn(20, 1); y2 = y1 + randn(20, 1) * 0.1;
            result = linCCC(y1, y2);
            testCase.verifyFalse(result.rInfo.available);
            testCase.verifyEqual(result.rInfo.message, "CompareR not requested");
        end

        function testCompareRSuccessPathWhenAvailable(testCase)
            testCase.assumeTrue(~isempty(which_test_rscript()), ...
                "Rscript not on path on this machine; success-path test not applicable.");
            rng(1);
            y1 = randn(50, 1) * 10 + 100;
            y2 = y1 + randn(50, 1) * 3;
            result = testCase.verifyWarningFree(@() linCCC(y1, y2, 0.05, CompareR=true));
            testCase.verifyTrue(result.rInfo.available);
            testCase.verifyTrue(isfinite(result.rInfo.ccc));
            testCase.verifyEqual(result.ccc, result.rInfo.ccc, "AbsTol", 0.001, ...
                "MATLAB and R CCC should agree within the validated 0.001 threshold");
            testCase.verifyEqual(result.rInfo.message, "");
        end

        function testCompareRPEFRAgreesWithinValidatedThreshold(testCase)
            testCase.assumeTrue(~isempty(which_test_rscript()), ...
                "Rscript not on path on this machine; success-path test not applicable.");
            y1 = [494,395,516,434,476,557,413,442,650,433,417,656,267,478,178,423,427]';
            y2 = [512,430,520,428,500,600,364,380,658,445,432,626,260,477,259,350,451]';
            result = linCCC(y1, y2, 0.05, CompareR=true);
            testCase.verifyLessThan(result.rInfo.gap, 0.001);
        end

        function testCompareRDegradesGracefullyWhenRscriptUnavailable(testCase)
            % Simulates R absence by stripping PATH for the duration of
            % this test, rather than requiring a machine without R.
            % Restores PATH afterward via addTeardown regardless of
            % test outcome.
            originalPath = getenv("PATH");
            setenv("PATH", "/usr/bin:/bin");
            testCase.addTeardown(@() setenv("PATH", originalPath));

            testCase.assumeTrue(isempty(which_test_rscript()), ...
                "Rscript still resolvable even with a stripped PATH on this machine; " + ...
                "degraded-mode simulation not applicable here.");

            rng(2);
            y1 = randn(30, 1) * 5 + 50;
            y2 = y1 + randn(30, 1) * 2;

            testCase.verifyWarning(@() linCCC(y1, y2, 0.05, CompareR=true), "linCCC:NoR");
            result = linCCC(y1, y2, 0.05, CompareR=true);
            testCase.verifyFalse(result.rInfo.available);
            testCase.verifyTrue(isfinite(result.ccc), ...
                "Core CCC computation must succeed regardless of CompareR availability");
        end
    end
end

function p = which_test_rscript()
%WHICH_TEST_RSCRIPT Local copy of linCCC's Rscript-detection logic, used
%   only to decide whether to run/skip R-dependent tests. Deliberately
%   duplicated rather than calling into linCCC's private local function
%   (not accessible from outside the file) -- keeps this test file
%   self-contained per the same self-containment discipline the toolbox
%   itself follows.
[status, out] = system('which Rscript 2>&1');
if status == 0 && ~isempty(strtrim(out))
    p = strtrim(out);
else
    p = '';
end
end
