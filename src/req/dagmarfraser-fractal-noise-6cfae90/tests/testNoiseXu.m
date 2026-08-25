classdef testNoiseXu < matlab.unittest.TestCase
%TESTNOISEXU Unit tests for noiseXu.m.
%
% Toolbox: fractalnoise
% Covers: source-equivalence (XuNoise_v002, generateCustomNoise_v004),
% edge cases (zero amplitude, arguments-block validation), and a
% known-alpha recovery check (regression net per TODO Item 7).

    properties (Constant)
        LegacyFunctionsDir = "/Users/d.s.fraser/Library/CloudStorage/Dropbox/Brain2Bee/PowerLawSimulationPreReg/src/functions"
    end

    methods (TestClassSetup)
        function addToolboxAndLegacyPaths(testCase)
            toolboxDir = fullfile(fileparts(mfilename("fullpath")), "..", "toolbox");
            addpath(toolboxDir);
            testCase.addTeardown(@() rmpath(toolboxDir));
            if isfolder(testCase.LegacyFunctionsDir)
                addpath(testCase.LegacyFunctionsDir);
                testCase.addTeardown(@() rmpath(testCase.LegacyFunctionsDir));
            end
        end
    end

    methods (Test)
        function testMatchesXuNoiseV002ForSameSeed(testCase)
            testCase.assumeTrue(exist("XuNoise_v002", "file") == 2, ...
                "XuNoise_v002 not on path; skipping source-equivalence test.");
            rng(42);
            n1 = noiseXu(2000, 2.0, 1.5, 60);
            rng(42);
            n2 = XuNoise_v002(2000, 2.0, 1.5, 60, 0.99);
            testCase.verifyEqual(n1, n2);
        end

        function testMatchesGenerateCustomNoiseV004ForSameSeed(testCase)
            testCase.assumeTrue(exist("generateCustomNoise_v004", "file") == 2, ...
                "generateCustomNoise_v004 not on path; skipping wrapper-equivalence test.");
            rng(7);
            nOld = generateCustomNoise_v004(1000, 3.0, 1.0, 100);
            rng(7);
            nNew = noiseXu(1000, 3.0, 1.0, 100);
            testCase.verifyEqual(nOld, nNew);
        end

        function testZeroAmplitudeReturnsExactZeros(testCase)
            z = noiseXu(50, 1.0, 0);
            testCase.verifyTrue(all(z == 0));
            testCase.verifyEqual(numel(z), 50);
        end

        function testPhiOverrideRuns(testCase)
            n = noiseXu(500, 1.0, 1.0, 100, Phi=0.995);
            testCase.verifyEqual(numel(n), 500);
            testCase.verifyEqual(std(n), 1.0, "AbsTol", 1e-9);
        end

        function testDefaultTwoArgCall(testCase)
            n = noiseXu(200, 0.5);
            testCase.verifyEqual(numel(n), 200);
            testCase.verifyEqual(std(n), 1.0, "AbsTol", 1e-9);
        end

        function testRejectsNonPositiveN(testCase)
            testCase.verifyError(@() noiseXu(-5, 1.0), "MATLAB:validators:mustBePositive");
        end

        function testRejectsPhiOutOfRange(testCase)
            testCase.verifyError(@() noiseXu(100, 1.0, 1.0, 100, Phi=1.0), ...
                "MATLAB:validators:mustBeLessThan");
            testCase.verifyError(@() noiseXu(100, 1.0, 1.0, 100, Phi=0), ...
                "MATLAB:validators:mustBeGreaterThan");
        end

        function testKnownAlphaRecoveryWithinTolerance(testCase)
            % Regression net: a direct pmtm log-log slope fit on a long
            % realisation should recover the requested alpha within a
            % generous tolerance. Not a precision estimator test (see
            % estimateIRASA for that); this only guards against a gross
            % regression in the generator itself.
            rng(1);
            N = 4000;
            fs = 100;
            targetAlpha = 2.0;
            n = noiseXu(N, targetAlpha, 5.0, fs);
            [pxx, f] = pmtm(n, 4, [], fs);
            fitMask = f >= 1 & f <= fs / 4 & pxx > 0;
            fit = polyfit(log10(f(fitMask)), log10(pxx(fitMask)), 1);
            recoveredAlpha = -fit(1);
            testCase.verifyEqual(recoveredAlpha, targetAlpha, "AbsTol", 0.6);
        end
    end
end
