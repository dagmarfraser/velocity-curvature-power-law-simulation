classdef testToolboxSelfContained < matlab.unittest.TestCase
%TESTTOOLBOXSELFCONTAINED Verifies the toolbox has no hidden path dependencies.
%
% Toolbox: fractalnoise
%
% This is deliberately separate from the per-function test files, which
% add PowerLawSimulationPreReg's legacy functions dir and FieldTrip to
% path for their own source-equivalence and comparison tests. Those
% convenience paths mask exactly the bug this file exists to catch: a
% toolbox function silently depending on something outside its own tree.
%
% This class caught a real bug during development: estimateIRASA.m
% originally called iraAlphaSigma_v003 (a PowerLawSimulationPreReg
% function never copied into the toolbox) rather than a bundled copy of
% the algorithm. Every per-function test still passed, because
% TestClassSetup in those files always added the legacy functions dir
% alongside the toolbox dir. Only a test that adds the toolbox path
% ALONE exposed it.

    properties (Constant)
        ToolboxDir = fullfile(fileparts(mfilename("fullpath")), "..", "toolbox")
    end

    methods (Test)
        function testCoreFunctionsWorkWithOnlyToolboxOnPath(testCase)
            originalPath = path();
            testCase.addTeardown(@() path(originalPath));
            restoredefaultpath();
            rehash();
            addpath(testCase.ToolboxDir);

            testCase.verifyEqual(exist("noiseXu", "file"), 2);
            testCase.verifyEqual(exist("shapeNoise", "file"), 2);
            testCase.verifyEqual(exist("estimateIRASA", "file"), 2);
            testCase.verifyEqual(exist("checkSpectralKnee", "file"), 2);

            rng(1);
            x = noiseXu(1500, 2.0, 5.0, 100);
            testCase.verifyEqual(numel(x), 1500);

            [a, s, p, f] = estimateIRASA(x, 100, 1.0, 25.0);
            testCase.verifyTrue(isfinite(a));
            testCase.verifyTrue(isfinite(s));
            testCase.verifyEqual(numel(p), numel(f));

            n = shapeNoise(x, 100, 2.0);
            testCase.verifyEqual(numel(n), numel(x));

            result = checkSpectralKnee(x, 100);
            testCase.verifyTrue(ismember(result.verdict, ["strong", "marginal", "none"]));
        end
    end
end
