classdef testEstimateIRASA < matlab.unittest.TestCase
%TESTESTIMATEIRASA Unit tests for estimateIRASA.m.
%
% Toolbox: fractalnoise
% Covers: output arities, Hset override, the CompareFieldTrip
% degraded-mode path, and the ComparePMTM resampling-free comparison
% (including pRaw plumbing source-equivalence against a direct pmtm
% call). All of the above use the toolbox's own noiseXu as the test
% signal generator, so they run standalone on any machine with only the
% toolbox directory on path -- no dependency on PowerLawSimulationPreReg
% legacy code or its machine-specific paths.
%
% testMatchesSourceForSameSeed is the one deliberate exception: its job
% is specifically to confirm extraction fidelity against
% iraAlphaSigma_v003/generateCustomNoise_v004, so it stays coupled to
% the legacy source and self-filters via assumeTrue when that source
% isn't on this machine's path.

    properties (Constant)
        LegacyFunctionsDir = "/Users/d.s.fraser/Library/CloudStorage/Dropbox/Brain2Bee/PowerLawSimulationPreReg/src/functions"
        FieldTripDir = "/Users/d.s.fraser/Library/CloudStorage/Dropbox/MATLAB/fieldtrip-20251201"
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
            if isempty(which("ft_defaults")) && isfolder(testCase.FieldTripDir)
                addpath(testCase.FieldTripDir);
                evalc("ft_defaults");
            end
        end
    end

    methods (Test)
        function testMatchesSourceForSameSeed(testCase)
            % Genuine source-equivalence test: this one must stay coupled
            % to the legacy functions, since its job is confirming
            % estimateIRASA matches iraAlphaSigma_v003 exactly. Every
            % other test below only needs *a* noise signal with known
            % alpha, which the toolbox's own noiseXu already provides --
            % see the vexed-about-hardcoded-paths discussion this
            % session, 2026-07-15.
            testCase.assumeTrue(exist("generateCustomNoise_v004", "file") == 2 ...
                && exist("iraAlphaSigma_v003", "file") == 2, ...
                "Legacy source functions not on path; skipping source-equivalence test.");
            rng(23);
            x = generateCustomNoise_v004(3000, 3.5, 5.0, 120);
            fLow = 1.0; fHigh = 30.0; fs = 120;

            [a1, s1, p1, f1] = iraAlphaSigma_v003(x, fs, fLow, fHigh);
            [a2, s2, p2, f2] = estimateIRASA(x, fs, fLow, fHigh);

            testCase.verifyEqual(a1, a2);
            testCase.verifyEqual(s1, s2);
            testCase.verifyTrue(isequaln(p1, p2), ...
                "pFractal must match source including NaN positions (isequaln, not isequal).");
            testCase.verifyEqual(f1, f2);
        end

        function testTwoOutputArity(testCase)
            rng(1);
            x = noiseXu(2000, 2.0, 3.0, 100);
            [a, s] = estimateIRASA(x, 100, 1.0, 25.0);
            testCase.verifyTrue(isscalar(a));
            testCase.verifyTrue(isscalar(s));
            testCase.verifyTrue(isfinite(a));
            testCase.verifyTrue(isfinite(s));
        end

        function testCustomHsetRuns(testCase)
            rng(2);
            x = noiseXu(2000, 2.0, 3.0, 100);
            [a, s] = estimateIRASA(x, 100, 1.0, 25.0, Hset=1.2:0.1:1.8);
            testCase.verifyTrue(isfinite(a));
            testCase.verifyTrue(isfinite(s));
        end

        function testDefaultFtInfoWhenNotRequested(testCase)
            rng(4);
            x = noiseXu(1500, 2.0, 3.0, 100);
            [~, ~, ~, ~, ftInfo] = estimateIRASA(x, 100, 1.0, 25.0);
            testCase.verifyFalse(ftInfo.available);
            testCase.verifyEqual(ftInfo.message, "CompareFieldTrip not requested");
        end

        function testCompareFieldTripDegradesGracefullyWhenUnavailable(testCase)
            testCase.assumeTrue(isempty(which("ft_defaults")), ...
                "FieldTrip is on the path on this machine; unavailable-path test not applicable.");
            rng(6);
            x = noiseXu(1500, 2.0, 3.0, 100);
            testCase.verifyWarning(@() estimateIRASA(x, 100, 1.0, 25.0, CompareFieldTrip=true), ...
                "estimateIRASA:NoFieldTrip");
            [~, ~, ~, ~, ftInfo] = estimateIRASA(x, 100, 1.0, 25.0, CompareFieldTrip=true); %#ok<ASGLU>
        end

        function testCompareFieldTripSuccessPathWhenAvailable(testCase)
            testCase.assumeTrue(~isempty(which("ft_defaults")), ...
                "FieldTrip not on path on this machine; success-path test not applicable.");
            rng(9);
            x = noiseXu(3000, 2.5, 4.0, 120);
            [alphaHB, ~, ~, ~, ftInfo] = testCase.verifyWarningFree(...
                @() estimateIRASA(x, 120, 1.0, 30.0, CompareFieldTrip=true));
            testCase.verifyTrue(ftInfo.available);
            testCase.verifyTrue(isfinite(ftInfo.alphaFT));
            testCase.verifyTrue(isfinite(ftInfo.gapHBFT));
            testCase.verifyEqual(ftInfo.gapHBFT, abs(alphaHB - ftInfo.alphaFT), "AbsTol", 1e-9);
            testCase.verifyEqual(ftInfo.message, "");
            % Sanity, not precision: two independent-ish IRASA implementations
            % on the same signal should not be wildly discrepant.
            testCase.verifyLessThan(ftInfo.gapHBFT, 1.0);
        end

        function testComparePMTMSuccessPath(testCase)
            % No external dependency (unlike CompareFieldTrip): should
            % always succeed when requested, regardless of machine.
            rng(9);
            x = noiseXu(3000, 2.5, 4.0, 120);
            [alphaHB, ~, ~, ~, ~, pmInfo] = testCase.verifyWarningFree(...
                @() estimateIRASA(x, 120, 1.0, 30.0, ComparePMTM=true));
            testCase.verifyTrue(pmInfo.available);
            testCase.verifyTrue(isfinite(pmInfo.alphaPM));
            testCase.verifyTrue(isfinite(pmInfo.gapHBPM));
            testCase.verifyEqual(pmInfo.gapHBPM, abs(alphaHB - pmInfo.alphaPM), "AbsTol", 1e-9);
            testCase.verifyEqual(pmInfo.message, "");
            % Sanity, not precision: HB and PM share the same PSD estimator
            % and differ only in the resampling-median step, so at this
            % alpha (comfortably inside the [-2,6] verified range) they
            % should not be wildly discrepant -- see
            % README_IRASAThreeWayAgreement_v001.md Section 3, alpha=2-3 rows.
            testCase.verifyLessThan(pmInfo.gapHBPM, 1.0);
        end

        function testDefaultPmInfoWhenNotRequested(testCase)
            rng(4);
            x = noiseXu(1500, 2.0, 3.0, 100);
            [~, ~, ~, ~, ~, pmInfo] = estimateIRASA(x, 100, 1.0, 25.0);
            testCase.verifyFalse(pmInfo.available);
            testCase.verifyEqual(pmInfo.message, "ComparePMTM not requested");
        end

        function testPRawMatchesDirectPmtmCall(testCase)
            % Source-equivalence for the pRaw plumbing added alongside
            % ComparePMTM: the h=1 PSD returned internally must be bit-
            % identical to calling pmtm(x,4,[],fs) directly -- this is
            % what ComparePMTM's alphaPM is fit on, so if this drifts,
            % ComparePMTM silently stops being the resampling-free
            % estimator it claims to be.
            rng(11);
            x = randn(2048, 1);
            fs = 100;
            [pDirect, fDirect] = pmtm(x, 4, [], fs);

            [~, ~, ~, fVec, ~, pmInfo] = estimateIRASA(x, fs, 1.0, 25.0, ComparePMTM=true); %#ok<ASGLU>
            testCase.verifyEqual(fVec, fDirect);

            fitMask = fDirect >= 1.0 & fDirect <= 25.0 & pDirect > 0;
            pFit = polyfit(log10(fDirect(fitMask)), log10(pDirect(fitMask)), 1);
            alphaDirect = -pFit(1);
            testCase.verifyEqual(pmInfo.alphaPM, alphaDirect, "AbsTol", 1e-12);
        end
    end
end
