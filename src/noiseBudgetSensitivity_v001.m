function S = noiseBudgetSensitivity_v001(opts)
%NOISEBUDGETSENSITIVITY_V001 How far does beta_gen* move if the end-matched noise carries x% more energy? (2026-11-15, M5)
%   Grounds the segment budget (Dagmar, 10:44): the budget bounds the extra derivative-noise energy the wrap can add; the rule
%   stated BEFORE running is: choose the largest budget whose beta_gen* shift stays below 0.003 (a tenth of the MDC, 0.03),
%   judged on the median and the 90th percentile of |shift| over trials, per arm.
%   0. Acceptance: v027 at NoiseScale 1 equals v026 exactly (Cook CTRL per_subject, small grid, endMatched, end-anchored,
%      budget 1%); and at least 80% of trials ran.
%   1. v027, NoiseModel 'endMatched', Template 'endAnchored', SegmentBudget 0.01, TrialSelection 'per_subject', N_BETA 13,
%      N_REPS 50, seed 1729, parfor, five core arms, NoiseScale 1, 1.01, 1.02, 1.05 (the same noise draws, scaled).
%   2. Per arm and scale: per-trial beta_gen* shift against scale 1 (betaGenStarMed, window basis; and SG-IRLS alone):
%      median, 90th percentile of |shift|, share of trials beyond 0.003 and beyond the MDC; trials that ran.
%   Writes results/templateComparison/noiseBudgetSensitivity_v001.csv (per trial) and _summary_v001.csv.
arguments
    opts.Arms (1,:) string = ["FraserN48", "Hickman PLAC", "Hickman HALO", "Cook CTRL", "Cook ASD"]
    opts.Scales (1,:) double = [1 1.01 1.02 1.05]
end
here = fileparts(mfilename('fullpath'));  addpath(here, fullfile(here, 'functions'));
od = fullfile(here, '..', 'results', 'templateComparison');  if ~isfolder(od), mkdir(od); end
a = {'N_BETA', 7, 'N_REPS', 3, 'N_R7', 2, 'ForceRerun', true, 'OutTag', '_acc27', 'NoiseModel', 'endMatched', 'Template', 'endAnchored', 'SegmentBudget', 0.01};
r26 = runLoopClosureFftnoise_v026('Cook CTRL', a{:});  r27 = runLoopClosureFftnoise_v027('Cook CTRL', a{:});
if nnz([r27.emStatus] == "ok") < 0.8 * numel(r27), error('nbs:ran', 'acceptance: too few trials ran under v027'); end
for k = 1:numel(r26)
    for f = ["betaObs", "betaGenStar", "betaRecCurveMed", "emSeg"]
        if ~isequaln(r26(k).(f), r27(k).(f)), error('nbs:acceptance', 'v027 at NoiseScale 1 differs from v026: %s, %s', r27(k).trialID, f); end
    end
end
fprintf('0. acceptance: v027 (NoiseScale 1) = v026 exactly on %d trials, %d ran\n', numel(r27), nnz([r27.emStatus] == "ok"));
ensurePool_v001();
P = table();
for arm = opts.Arms
    for sc = opts.Scales
        t0 = tic;
        R = runLoopClosureFftnoise_v027(char(arm), 'TrialSelection', 'per_subject', 'NoiseModel', 'endMatched', 'Template', 'endAnchored', ...
            'SegmentBudget', 0.01, 'NoiseScale', sc, 'N_BETA', 13, 'N_REPS', 50, 'N_R7', 2, 'RngSeed', 1729, 'UseParfor', true, 'ForceRerun', true, 'OutTag', '_nbs1');
        sg = arrayfun(@(r) r.betaGenStar(6), R);
        P = [P; table(repmat(arm, numel(R), 1), repmat(sc, numel(R), 1), string({R.trialID})', [R.betaGenStarMed]', sg(:), [R.emStatus]', ...
            'VariableNames', {'arm', 'scale', 'trialID', 'bgMed', 'bgSGIRLS', 'emStatus'})]; %#ok<AGROW>
        fprintf('%s scale %.2f: %d trials, %d ran, %.1f min\n', arm, sc, numel(R), nnz([R.emStatus] == "ok"), toc(t0) / 60);
    end
end
writetable(P, fullfile(od, 'noiseBudgetSensitivity_v001.csv'));
S = table();
for arm = opts.Arms
    B = P(P.arm == arm & P.scale == 1, :);
    for sc = opts.Scales(opts.Scales ~= 1)
        Q = P(P.arm == arm & P.scale == sc, :);  [~, ib, iq] = intersect(B.trialID, Q.trialID);
        d = Q.bgMed(iq) - B.bgMed(ib);  d = d(isfinite(d));  ds = Q.bgSGIRLS(iq) - B.bgSGIRLS(ib);  ds = ds(isfinite(ds));
        S = [S; table(arm, sc, numel(d), median(d), prctile(abs(d), 90), mean(abs(d) > 0.003), mean(abs(d) > 0.03), median(ds), prctile(abs(ds), 90), ...
            'VariableNames', {'arm', 'scale', 'n', 'medShift', 'p90AbsShift', 'shareOver0003', 'shareOverMDC', 'medShiftSGIRLS', 'p90AbsShiftSGIRLS'})]; %#ok<AGROW>
    end
end
writetable(S, fullfile(od, 'noiseBudgetSensitivity_summary_v001.csv'));  disp(S);
end
