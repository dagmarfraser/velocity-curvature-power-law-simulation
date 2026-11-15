function T = checkV020Acceptance_v001(opts)
%CHECKV020ACCEPTANCE_V001  Acceptance test for runLoopClosureFftnoise_v020 (2026-11-03).
%   On the same machine, same seed, serial (no parfor), for a few 'Fraser' trials at a small grid:
%     A. v020 Template='endAnchored'  must equal  v019   (betaObs, betaRecCurveMed, betaGenStar, invertStatus, a_mm)
%     B. v020 Template='production'   must equal  v015   (the same fields)
%   isequaln, no tolerance: v020 changes no numerics at these settings, so anything but exact equality fails.
%   Each runner writes a per-trial file; the check refuses to start if any of those names already exists,
%   and deletes only the files it created.   Fraser, D.S. (2026)  v001
arguments
    opts.TrialIDs (1,:) string = ["F211_V2_S1_T049_sh3" "F101_V1_S1_T005_sh3"]
    opts.N_REPS   (1,1) double = 4
    opts.N_BETA   (1,1) double = 5
end
srcDir = fileparts(mfilename('fullpath'));  addpath(srcDir);
fields = ["betaObs" "betaRecCurveMed" "betaGenStar" "invertStatus" "a_mm" "nEbrErr"];
common = {'N_REPS', opts.N_REPS, 'N_BETA', opts.N_BETA, 'UseParfor', false, 'NoiseModel', 'shaped_xu'};
runs = { ...
    "v019",             @(tid) runLoopClosureFftnoise_v019('Fraser', 'TrialID', char(tid), common{:}), "_v019.mat"; ...
    "v020 endAnchored", @(tid) runLoopClosureFftnoise_v020('Fraser', 'TrialID', char(tid), common{:}, 'Template', 'endAnchored'), "_v020.mat"; ...
    "v015",             @(tid) runLoopClosureFftnoise_v015('Fraser', 'TrialID', char(tid), common{:}), "_v015.mat"; ...
    "v020 production",  @(tid) runLoopClosureFftnoise_v020('Fraser', 'TrialID', char(tid), common{:}, 'Template', 'production'), "_tplProd_v020.mat"};

% Output names each runner will write (fail loud if any exists: never overwrite)
outOf = @(tid, suffix) fullfile(srcDir, "loopClosureResults_Fraser_" + tid + "_shaped_xu" + suffix);
for tid = opts.TrialIDs
    for k = 1:size(runs, 1)
        f = outOf(tid, runs{k, 3});
        if isfile(f), error('checkV020Acceptance:exists', '%s', "Refusing to run: " + f + " exists."); end
    end
end

R = cell(numel(opts.TrialIDs), size(runs, 1));  made = strings(0);
for i = 1:numel(opts.TrialIDs)
    for k = 1:size(runs, 1)
        r = runs{k, 2}(opts.TrialIDs(i));
        f = outOf(opts.TrialIDs(i), runs{k, 3});
        if ~isfile(f), error('checkV020Acceptance:noOutput', '%s', "Runner did not write " + f); end
        made(end + 1) = f; %#ok<AGROW>
        R{i, k} = r;
    end
end

deleteMade_local(made);   % only the files this check created (none existed before it ran)

% Compare
pairs = [1 2; 3 4];  pairName = ["A: v020 endAnchored = v019" "B: v020 production = v015"];
rows = {};
for i = 1:numel(opts.TrialIDs)
    for q = 1:2
        a = R{i, pairs(q, 1)};  b = R{i, pairs(q, 2)};
        for fld = fields(1:end-1)
            rows(end + 1, :) = {opts.TrialIDs(i), pairName(q), fld, isequaln(a.(fld), b.(fld))}; %#ok<AGROW>
        end
    end
    rows(end + 1, :) = {opts.TrialIDs(i), "v020 EBR errors (endAnchored, production)", "nEbrErr", ...
        sprintf("%d, %d", R{i, 2}.nEbrErr, R{i, 4}.nEbrErr)}; %#ok<AGROW>
end
T = cell2table(rows, 'VariableNames', ["trialID" "check" "field" "equal"]);
disp(T)
eqCol = T.equal(T.field ~= "nEbrErr");
if ~all([eqCol{:}]), error('checkV020Acceptance:FAIL', '%s', 'v020 does not reproduce v019/v015 exactly (see table).'); end
% The two templates must differ somewhere, or the switch is not wired
anyDiff = false;
for i = 1:numel(opts.TrialIDs), anyDiff = anyDiff || ~isequaln(R{i, 2}.betaRecCurveMed, R{i, 4}.betaRecCurveMed); end
fprintf('Templates differ on at least one test trial: %d (0 means neither trial has an uncovered tail; pick another trial)\n', anyDiff);
fprintf('=== checkV020Acceptance_v001: PASS ===\n');
end

function deleteMade_local(made)
    for f = made, if isfile(f), delete(f); end, end
end
