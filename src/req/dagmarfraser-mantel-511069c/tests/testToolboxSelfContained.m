classdef testToolboxSelfContained < matlab.unittest.TestCase
%TESTTOOLBOXSELFCONTAINED Verifies the toolbox has no hidden path dependencies.
%
% Toolbox: mantel
%
% Deliberately separate from testMantelTest.m, which may add other
% project paths for its own CompareR / source-comparison purposes. This
% class adds ONLY the toolbox directory to a freshly restored default
% path and confirms the public function still resolves and runs
% correctly -- see docs/README_MakingToolboxes_v001.md ("Self-containment
% is not proven by tests that add convenience paths").

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

            testCase.verifyEqual(exist("mantelTest", "file"), 2);

            rng(1);
            n = 6;
            X = rand(n, 2);
            Y = rand(n, 2);
            D1 = squareform(pdist_local(X));
            D2 = squareform(pdist_local(Y));

            [mantelR, pValue, info] = mantelTest(D1, D2);
            testCase.verifyTrue(isfinite(mantelR));
            testCase.verifyTrue(pValue >= 0 && pValue <= 1);
            testCase.verifyClass(info, "struct");

            % Also exercise the "observations" input path directly (no
            % Statistics and Machine Learning Toolbox pdist involved).
            [mantelR2, pValue2] = mantelTest(X, Y, InputType="observations");
            testCase.verifyTrue(isfinite(mantelR2));
            testCase.verifyTrue(pValue2 >= 0 && pValue2 <= 1);
        end
    end
end

function D = pdist_local(X)
% Hand-rolled pairwise Euclidean distance vector (upper triangle order),
% avoiding a Statistics and Machine Learning Toolbox dependency in this
% self-containment test itself.
    n = size(X, 1);
    D = zeros(1, n * (n - 1) / 2);
    idx = 1;
    for i = 1:n-1
        for j = i+1:n
            D(idx) = sqrt(sum((X(i, :) - X(j, :)) .^ 2));
            idx = idx + 1;
        end
    end
end
