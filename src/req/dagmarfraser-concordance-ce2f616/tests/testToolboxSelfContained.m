classdef testToolboxSelfContained < matlab.unittest.TestCase
%TESTTOOLBOXSELFCONTAINED Verifies the toolbox has no hidden path dependencies.
%
% Toolbox: concordance
%
% Deliberately separate from testLinCCC.m, which may add other project
% paths for its own CompareR / source-comparison purposes. This class
% adds ONLY the toolbox directory to a freshly restored default path and
% confirms the public function still resolves and runs correctly -- see
% docs/README_MakingToolboxes_v001.md ("Self-containment is not proven
% by tests that add convenience paths").

    properties (Constant)
        ToolboxDir = fullfile(fileparts(mfilename("fullpath")), "..", "toolbox")
    end

    methods (Test)
        function testCoreFunctionWorksWithOnlyToolboxOnPath(testCase)
            originalPath = path();
            testCase.addTeardown(@() path(originalPath));
            restoredefaultpath();
            rehash();
            addpath(testCase.ToolboxDir);

            testCase.verifyEqual(exist("linCCC", "file"), 2);

            rng(1);
            y1 = randn(50, 1);
            y2 = y1 + 0.1 * randn(50, 1);

            cccResults = linCCC(y1, y2);
            testCase.verifyClass(cccResults, "struct");
            testCase.verifyTrue(cccResults.ccc > 0.9);
            testCase.verifyTrue(cccResults.ccc <= 1);
            testCase.verifyEqual(cccResults.n, 50);
            testCase.verifyTrue(cccResults.LowerCI <= cccResults.ccc);
            testCase.verifyTrue(cccResults.UpperCI >= cccResults.ccc);
        end
    end
end
