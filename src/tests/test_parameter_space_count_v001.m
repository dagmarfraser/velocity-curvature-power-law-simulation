function test_parameter_space_count_v001()
% test_parameter_space_count_v001  defineParameterSpace returns the expected count.
%
% Calls defineParameterSpace(0) (full production run) and verifies:
%   (a) paramSpace.totalConfigCount == EXPECTED_TOTAL (14,764,680)
%   (b) Cartesian product of individual dimension vectors equals totalConfigCount
%   (c) Each dimension has the expected number of values
%
% Derivation of expected total:
%   3 fs x 1 shape x 21 beta x 14 VGF x 31 noiseType x 18 noiseMag
%   x 2 filter x 3 regress x 5 trials = 14,764,680
%
% Run from src/ (tests/ must be on the path).
%
% Dagmar Scott Fraser | PowerLawSimulationPreReg | 2026

addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..', 'functions')));

EXPECTED_TOTAL  = 14764680;
EXPECTED_TRIALS = 5;

%  { fieldName, expectedCount, description }
dims = {
    'samplingRates',   3,  'fs: [60, 120, 240] Hz';
    'shapes',          1,  'shape: ellipse only (index 6)';
    'generatedBetas',  21, 'beta: 0:(2/3)/20:2/3';
    'vgfValues',       14, 'VGF: exp(4.5:0.1:5.8)';
    'noiseTypes',      31, 'alpha: 0:0.1:3.0';
    'noiseMagnitudes', 18, 'sigma: 5+9+4 levels';
    'filterTypes',     2,  'filter: [BWFD=2, SG=6]';
    'regressTypes',    3,  'regress: [OLS=3, LMLS=4, IRLS=5]';
};

fprintf('=== test_parameter_space_count_v001 ===\n');
fprintf('Expected total configs: %d\n\n', EXPECTED_TOTAL);

%% ---- Call defineParameterSpace ----
try
    paramSpace = defineParameterSpace(0);
catch ME
    fprintf('FAIL: defineParameterSpace(0) threw: %s\n', ME.message);
    return;
end

if ~isfield(paramSpace, 'totalConfigCount')
    fprintf('FAIL: missing field totalConfigCount\n');
    return;
end

%% ---- (a) totalConfigCount ----
actual     = paramSpace.totalConfigCount;
countMatch = (actual == EXPECTED_TOTAL);
fprintf('  (a) totalConfigCount = %d  (expected %d)  %s\n\n', actual, EXPECTED_TOTAL, verdict(countMatch));

%% ---- (b) individual dimensions ----
fprintf('  %-22s  %-10s  %-10s  %s\n', 'field', 'expected', 'actual', 'verdict');
fprintf('  %s\n', repmat('-', 1, 68));

allDimPass = true;
nProduct   = 1;

for i = 1:size(dims, 1)
    field    = dims{i, 1};
    expected = dims{i, 2};
    n        = -1;
    if isfield(paramSpace, field), n = numel(paramSpace.(field)); end
    pass = (n == expected);
    if ~pass, allDimPass = false; end
    nProduct = nProduct * max(n, 0);
    fprintf('  %-22s  %-10d  %-10d  %s   %s\n', field, expected, n, verdict(pass), dims{i,3});
end

actualTrials = paramSpace.repeatTrial;
trialPass    = (actualTrials == EXPECTED_TRIALS);
if ~trialPass, allDimPass = false; end
nProduct = nProduct * actualTrials;
fprintf('  %-22s  %-10d  %-10d  %s   trials per combo\n', ...
    'repeatTrial', EXPECTED_TRIALS, actualTrials, verdict(trialPass));

%% ---- (c) independent Cartesian product ----
cartMatch = (nProduct == EXPECTED_TOTAL);
fprintf('\n  (c) product of all dims x trials = %d  %s\n', nProduct, verdict(cartMatch));

%% ---- Verdict ----
allPass = countMatch && allDimPass && cartMatch;
fprintf('\n');
if allPass
    fprintf('PASS: defineParameterSpace returns %d configs across all dimension checks\n', EXPECTED_TOTAL);
else
    fprintf('FAIL: count or dimension mismatch (see table above)\n');
end

end

function s = verdict(pass)
if pass, s = 'PASS'; else, s = 'FAIL'; end
end
