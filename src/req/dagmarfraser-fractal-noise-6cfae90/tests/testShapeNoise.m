classdef testShapeNoise < matlab.unittest.TestCase
%TESTSHAPENOISE Unit tests for shapeNoise.m.
%
% Toolbox: fractalnoise
% Covers: source-equivalence (shapeXu_local, reconstructed inline since
% it is a private local function not directly callable), all PhaseSource
% carriers, edge cases, and error paths.

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
        function testMatchesSourceShapeXuLocalForSameSeed(testCase)
            testCase.assumeTrue(exist("generateCustomNoise_v004", "file") == 2, ...
                "generateCustomNoise_v004 not on path; skipping source-equivalence test.");
            rng(11);
            resid = cumsum(randn(500, 1));
            resid = resid - mean(resid);
            fs = 100;
            alpha = 3.0;

            % Reconstruct source shapeXu_local inline (private local
            % function, not independently callable).
            rng(99);
            sigOld = std(resid, 0, 1);
            eEmpOld = abs(fft(resid));
            baseOld = generateCustomNoise_v004(numel(resid), alpha, sigOld, fs);
            fbOld = fft(baseOld);
            foutOld = eEmpOld .* exp(1i * angle(fbOld));
            nOld = real(ifft(foutOld));
            sNOld = std(nOld, 0, 1);
            if sNOld > 0
                nOld = nOld * (sigOld / sNOld);
            end

            rng(99);
            nNew = shapeNoise(resid, fs, alpha);

            testCase.verifyEqual(nOld, nNew);
        end

        function testWhiteCarrierPreservesStd(testCase)
            rng(21);
            resid = cumsum(randn(300, 1));
            n = shapeNoise(resid, 100, 2.0, PhaseSource="white");
            testCase.verifyEqual(std(n), std(resid), "RelTol", 1e-9);
        end

        function testCustomFunctionHandleCarrier(testCase)
            rng(31);
            resid = cumsum(randn(300, 1));
            customCarrier = @(M, fs, a) randn(M, 1) * 2;
            n = shapeNoise(resid, 100, 2.0, PhaseSource=customCarrier);
            testCase.verifyEqual(numel(n), numel(resid));
            testCase.verifyEqual(std(n), std(resid), "RelTol", 1e-9);
        end

        function testZeroVarianceResidualReturnsZeros(testCase)
            n = shapeNoise(zeros(50, 1), 100, 2.0);
            testCase.verifyTrue(all(n == 0));
        end

        function testUnknownPhaseSourceErrors(testCase)
            resid = randn(100, 1);
            testCase.verifyError(@() shapeNoise(resid, 100, 2.0, PhaseSource="bogus"), ...
                "shapeNoise:UnknownPhaseSource");
        end

        function testInvalidCarrierLengthErrors(testCase)
            resid = randn(100, 1);
            badCarrier = @(M, fs, a) randn(M + 1, 1);
            testCase.verifyError(@() shapeNoise(resid, 100, 2.0, PhaseSource=badCarrier), ...
                "shapeNoise:InvalidCarrier");
        end

        function testPinknoiseWarnsOnAlphaMismatchWhenAvailable(testCase)
            testCase.assumeTrue(exist("pinknoise", "file") == 2 || exist("pinknoise", "builtin") ~= 0, ...
                "pinknoise not available on this MATLAB installation; skipping.");
            resid = randn(200, 1);
            testCase.verifyWarning(@() shapeNoise(resid, 100, 3.0, PhaseSource="pinknoise"), ...
                "shapeNoise:PinknoiseAlphaMismatch");
        end

        function testDspInRangeAlphaRunsCleanly(testCase)
            testCase.assumeTrue(exist("dsp.ColoredNoise", "class") == 8, ...
                "dsp.ColoredNoise (DSP System Toolbox) not installed; skipping.");
            rng(41);
            resid = cumsum(randn(400, 1));
            resid = resid - mean(resid);
            n = testCase.verifyWarningFree(@() shapeNoise(resid, 100, 1.5, PhaseSource="dsp"));
            testCase.verifyEqual(numel(n), numel(resid));
            testCase.verifyEqual(std(n), std(resid), "RelTol", 1e-9);
        end

        function testDspOutOfRangeAlphaClampsWithWarning(testCase)
            testCase.assumeTrue(exist("dsp.ColoredNoise", "class") == 8, ...
                "dsp.ColoredNoise (DSP System Toolbox) not installed; skipping.");
            rng(42);
            resid = cumsum(randn(400, 1));
            resid = resid - mean(resid);
            testCase.verifyWarning(@() shapeNoise(resid, 100, 4.0, PhaseSource="dsp"), ...
                "shapeNoise:DspAlphaClamped");
            n = shapeNoise(resid, 100, 4.0, PhaseSource="dsp");
            testCase.verifyEqual(numel(n), numel(resid));
            testCase.verifyEqual(std(n), std(resid), "RelTol", 1e-9);
        end

        function testDspUnavailableErrorsCleanly(testCase)
            testCase.assumeFalse(exist("dsp.ColoredNoise", "class") == 8, ...
                "dsp.ColoredNoise is installed on this machine; unavailable-path test not applicable.");
            resid = randn(200, 1);
            testCase.verifyError(@() shapeNoise(resid, 100, 2.0, PhaseSource="dsp"), ...
                "shapeNoise:DspUnavailable");
        end
    end
end
