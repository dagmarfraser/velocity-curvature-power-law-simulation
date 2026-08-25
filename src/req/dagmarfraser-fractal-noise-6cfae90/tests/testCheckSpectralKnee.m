classdef testCheckSpectralKnee < matlab.unittest.TestCase
%TESTCHECKSPECTRALKNEE Unit tests for checkSpectralKnee.m.
%
% Toolbox: fractalnoise
% Covers: source-equivalence (checkLorentzianKnee_v001's per-trial core,
% reconstructed inline for comparison since the source is a script, not
% a callable function), default-starts behaviour, and error paths.

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
        function testMatchesSourcePerTrialCoreGivenSameStarts(testCase)
            testCase.assumeTrue(exist("generateCustomNoise_v004", "file") == 2 ...
                && exist("iraAlphaSigma_v003", "file") == 2, ...
                "Legacy source functions not on path; skipping source-equivalence test.");

            fs = 133;
            fitBandLow = 1.0;
            wideFitLow = 0.15;
            fHigh = fs / 2.2;
            hset = 1.1:0.05:1.9;
            phiXu = 0.99;
            fKneeXu = (1 - phiXu) * fs / (2 * pi);

            rng(5);
            x = generateCustomNoise_v004(900, 2.5, 3.0, fs);
            sig = x - mean(x);

            lorentzLog10 = @(th, f) th(1) - (th(3) / 2) .* log10(10 ^ (2 * th(2)) + f .^ 2);
            fminOpts = optimset("Display", "off", "TolFun", 1e-10, "TolX", 1e-10, "MaxFunEvals", 10000);

            [~, ~, pFrac, fVec] = iraAlphaSigma_v003(sig, fs, wideFitLow, fHigh, hset);
            fVec = fVec(:); pFrac = pFrac(:);

            maskS = fVec >= fitBandLow & fVec <= fHigh & ~isnan(pFrac) & pFrac > 0;
            pCS = polyfit(log10(fVec(maskS)), log10(pFrac(maskS)), 1);
            alphaStdOld = -pCS(1);

            maskW = fVec >= wideFitLow & fVec <= fHigh & ~isnan(pFrac) & pFrac > 0;
            fW = fVec(maskW); logPW = log10(pFrac(maskW)); nW = numel(fW);
            pCW = polyfit(log10(fW), logPW, 1);
            bW = pCW(2); alphaW = -pCW(1);
            rssF = sum((logPW - (bW - alphaW * log10(fW))) .^ 2);
            aicF = nW * log(rssF / nW) + 2 * 2;

            sseFn = @(th) sum((logPW - lorentzLog10(th, fW)) .^ 2);
            starts = [bW, log10(fKneeXu), alphaW; bW, log10(wideFitLow), alphaW; bW, log10(fitBandLow), alphaW];
            bestRss = Inf; thBest = starts(1, :);
            for si = 1:size(starts, 1)
                [thTry, rssTry] = fminsearch(sseFn, starts(si, :), fminOpts);
                if rssTry < bestRss
                    bestRss = rssTry; thBest = thTry;
                end
            end
            kneeHzOld = 10 ^ thBest(2);
            daicOld = nW * log(bestRss / nW) + 2 * 3 - aicF;

            result = checkSpectralKnee(x, fs, FitBandLow=fitBandLow, WideBandLow=wideFitLow, ...
                Hset=hset, ExtraKneeStarts=fKneeXu);

            testCase.verifyEqual(result.alphaStd, alphaStdOld, "AbsTol", 1e-9);
            testCase.verifyEqual(result.alphaWide, alphaW, "AbsTol", 1e-9);
            testCase.verifyEqual(result.kneeHz, kneeHzOld, "AbsTol", 1e-6);
            testCase.verifyEqual(result.deltaAIC, daicOld, "AbsTol", 1e-6);
        end

        function testDefaultArgsCallRuns(testCase)
            testCase.assumeTrue(exist("generateCustomNoise_v004", "file") == 2, ...
                "generateCustomNoise_v004 not on path; skipping.");
            rng(5);
            x = generateCustomNoise_v004(900, 2.5, 3.0, 133);
            result = checkSpectralKnee(x, 133);
            testCase.verifyEqual(result.nSamples, 900);
            testCase.verifyTrue(ismember(result.verdict, ["strong", "marginal", "none"]));
            testCase.verifyTrue(isfinite(result.alphaStd));
            testCase.verifyTrue(isfinite(result.kneeHz));
        end

        function testTooShortResidualErrors(testCase)
            testCase.verifyError(@() checkSpectralKnee(randn(50, 1), 133), ...
                "checkSpectralKnee:TooShort");
        end

        function testKneeInBandLogicConsistentWithBounds(testCase)
            testCase.assumeTrue(exist("generateCustomNoise_v004", "file") == 2, ...
                "generateCustomNoise_v004 not on path; skipping.");
            rng(8);
            x = generateCustomNoise_v004(1200, 2.0, 2.0, 100);
            result = checkSpectralKnee(x, 100);
            expected = result.kneeHz >= result.fBandLow && result.kneeHz <= result.fHigh;
            testCase.verifyEqual(result.kneeInBand, expected);
        end
    end
end
