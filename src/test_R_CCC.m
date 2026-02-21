%% Test Lin's CCC: MATLAB vs R Implementation
% Verify linCCC_v001 matches R's DescTools::CCC()
%
% Follows the ICC validation pattern from MSc Analysis.
% CSV bridge between MATLAB and R for reproducibility.
%
% USAGE:
%   test_R_CCC          % Run all tests
%
% VERSION: v002 (2026-01-20) - Added Bland-Altman PEFR canonical test

function test_R_CCC()

fprintf('=== Testing Lin''s CCC: MATLAB vs R (DescTools) ===\n\n');

% Add functions path
addpath(fullfile(fileparts(mfilename('fullpath')), 'functions'));

%% Verify DescTools package
fprintf('Checking R DescTools package...\n');
[status, ~] = system('Rscript -e "library(DescTools); cat(''DescTools loaded'')" 2>&1');
if status ~= 0
    fprintf('Installing DescTools package...\n');
    [status2, out] = system('Rscript -e "install.packages(''DescTools'', repos=''https://cran.r-project.org'')" 2>&1');
    if status2 ~= 0
        error('Failed to install DescTools:\n%s', out);
    end
    fprintf('  ✓ DescTools installed\n\n');
else
    fprintf('  ✓ DescTools ready\n\n');
end

%% Run tests
results = struct('name', {}, 'pass', {}, 'diff_ccc', {});

% Test 1: Near-perfect agreement
rng(42);
y1 = randn(50, 1) * 10 + 100;
y2 = y1 + randn(50, 1) * 0.01;
results(end+1) = runTest(y1, y2, 'Near-perfect agreement (tiny noise)');

% Test 2: Good agreement with location shift
rng(43);
y1 = randn(50, 1) * 10 + 100;
y2 = y1 + 5;  % Constant shift
results(end+1) = runTest(y1, y2, 'Location shift (+5)');

% Test 3: Good agreement with scale shift
rng(44);
y1 = randn(50, 1) * 10 + 100;
y2 = y1 * 1.2;  % Scale shift
results(end+1) = runTest(y1, y2, 'Scale shift (×1.2)');

% Test 4: Moderate agreement
rng(45);
y1 = randn(50, 1) * 10 + 100;
y2 = y1 + randn(50, 1) * 5;
results(end+1) = runTest(y1, y2, 'Moderate agreement');

% Test 5: No agreement (independent)
rng(46);
y1 = randn(50, 1);
y2 = randn(50, 1);
results(end+1) = runTest(y1, y2, 'No agreement (independent)');

% Test 6: Larger sample
rng(47);
y1 = randn(200, 1) * 10 + 100;
y2 = y1 + randn(200, 1) * 2;
results(end+1) = runTest(y1, y2, 'Larger sample (n=200)');

% Test 7: Bland-Altman (1986) PEFR data — canonical test case
% Peak Expiratory Flow Rate measured by two methods (Wright peak flow meter)
% Source: Bland JM, Altman DG (1986). Statistical methods for assessing
%         agreement between two methods of clinical measurement. Lancet 327:307-310.
% This data appears in DescTools::CCC and epiR::epi.ccc documentation.
y1 = [494,395,516,434,476,557,413,442,650,433,417,656,267,478,178,423,427]';
y2 = [512,430,520,428,500,600,364,380,658,445,432,626,260,477,259,350,451]';
results(end+1) = runTest(y1, y2, 'Bland-Altman (1986) PEFR data');

%% Summary
fprintf('\n=== Summary ===\n');
n_pass = sum([results.pass]);
fprintf('Passed: %d/%d tests\n', n_pass, length(results));
fprintf('Validation criterion: |MATLAB - R| < 0.001 for CCC\n');

if n_pass == length(results)
    fprintf('\n✓ ALL TESTS PASSED - linCCC_v001 validated against R DescTools::CCC\n');
else
    fprintf('\n✗ SOME TESTS FAILED\n');
    for i = 1:length(results)
        if ~results(i).pass
            fprintf('  - %s (diff = %.6f)\n', results(i).name, results(i).diff_ccc);
        end
    end
end

end

%% Helper function
function result = runTest(y1, y2, testName)
    try
        % MATLAB implementation
        res_M = linCCC_v001(y1, y2);
        
        % R implementation
        res_R = computeCCC_R(y1, y2);
        
        % Calculate differences
        diff_ccc = abs(res_M.ccc - res_R.ccc);
        diff_r = abs(res_M.r - res_R.r);
        diff_Cb = abs(res_M.Cb - res_R.Cb);
        
        % Display
        fprintf('%s:\n', testName);
        fprintf('  MATLAB: CCC=%.4f, r=%.4f, Cb=%.4f [%.4f, %.4f]\n', ...
            res_M.ccc, res_M.r, res_M.Cb, res_M.LowerCI, res_M.UpperCI);
        fprintf('  R:      CCC=%.4f, r=%.4f, Cb=%.4f [%.4f, %.4f]\n', ...
            res_R.ccc, res_R.r, res_R.Cb, res_R.LowerCI, res_R.UpperCI);
        fprintf('  Diff:   CCC=%.6f, r=%.6f, Cb=%.6f\n', diff_ccc, diff_r, diff_Cb);
        
        % Pass/fail
        pass = diff_ccc < 0.001;
        if pass
            fprintf('  ✓ PASS\n\n');
        else
            fprintf('  ✗ FAIL\n\n');
        end
        
        result.name = testName;
        result.pass = pass;
        result.diff_ccc = diff_ccc;
        
    catch ME
        fprintf('%s:\n', testName);
        fprintf('  ✗ ERROR: %s\n\n', ME.message);
        result.name = testName;
        result.pass = false;
        result.diff_ccc = NaN;
    end
end
