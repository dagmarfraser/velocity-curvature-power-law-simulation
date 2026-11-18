function S = templateComparisonBeta_v001(opts)
%TEMPLATECOMPARISONBETA_V001 Template comparison, stage 2: how much does the template move beta_gen*? (2026-11-15, M5)
%   Stage 1 (templateComparison_v002) showed a trade-off and a production edge artefact; this measures the consequence.
%   0. Acceptance: runLoopClosureFftnoise_v025 at its defaults equals v024 exactly (Cook CTRL per_subject, endMatched,
%      production template, N_BETA 7, N_REPS 3); errors otherwise.
%   1. Runs v025, NoiseModel 'endMatched', TrialSelection 'per_subject' (one trial per participant), N_BETA 13 over [0 0.75],
%      N_REPS 50, seed 1729, parfor, on the five core arms, with three templates:
%        prod4  Template 'production',  NHarmonics 4     (= the running v024 runs)
%        end4   Template 'endAnchored', NHarmonics 4
%        endFc  Template 'endAnchored', NHarmonics 'fc'  (every harmonic up to the 10 Hz BWFD cutoff)
%      Outputs tagged _tc1 in src/ (loopClosureResults_<ds>_per_subject_endMatched[_tplProd][_nHfc]_tc1_v025.mat).
%   2. Per arm and template: trials, end-matched failures (emStatus not ok, reported as the JOSS thread asked), noise
%      replicates rejected, median emL, median betaGenStarMed on the window and the whole trial, and per-trial paired
%      differences (end4 - prod4, endFc - end4; window basis): median, IQR and share beyond the MDC (0.03).
%   Writes results/templateComparison/templateComparisonBeta_v001.csv (per trial) and _summary_v001.csv. Screening grid,
%   one seed: differences below about 0.01 are within the grid's resolution.
arguments
    opts.Arms (1,:) string = ["FraserN48", "Hickman PLAC", "Hickman HALO", "Cook CTRL", "Cook ASD"]
end
here = fileparts(mfilename('fullpath'));  addpath(here, fullfile(here, 'functions'));
od = fullfile(here, '..', 'results', 'templateComparison');  if ~isfolder(od), mkdir(od); end
a = {'N_BETA', 7, 'N_REPS', 3, 'N_R7', 2, 'ForceRerun', true, 'OutTag', '_acc25', 'NoiseModel', 'endMatched', 'Template', 'production'};
r24 = runLoopClosureFftnoise_v024('Cook CTRL', a{:});  r25 = runLoopClosureFftnoise_v025('Cook CTRL', a{:});
for k = 1:numel(r24)
    for f = ["betaObs", "betaObsWhole", "betaGenStar", "betaRecCurveMed", "emSeg"]
        if ~isequaln(r24(k).(f), r25(k).(f)), error('tcb:acceptance', 'v025 defaults differ from v024: trial %s field %s', r24(k).trialID, f); end
    end
end
fprintf('0. acceptance: v025 defaults = v024 exactly on %d Cook CTRL trials\n', numel(r25));
ensurePool_v001();
tpl = struct('name', {"prod4", "end4", "endFc"}, 'template', {'production', 'endAnchored', 'endAnchored'}, 'nh', {4, 4, 'fc'});
P = table();
for arm = opts.Arms
    for t = 1:numel(tpl)
        t0 = tic;
        R = runLoopClosureFftnoise_v025(char(arm), 'TrialSelection', 'per_subject', 'NoiseModel', 'endMatched', 'Template', tpl(t).template, ...
            'NHarmonics', tpl(t).nh, 'N_BETA', 13, 'N_REPS', 50, 'N_R7', 2, 'RngSeed', 1729, 'UseParfor', true, 'ForceRerun', true, 'OutTag', '_tc1');
        P = [P; table(repmat(arm, numel(R), 1), repmat(tpl(t).name, numel(R), 1), string({R.trialID})', [R.betaGenStarMed]', [R.betaGenStarMedWhole]', ...
            [R.emL]', [R.nNoiseRejected]', [R.emStatus]', [R.nHarm]', 'VariableNames', {'arm', 'template', 'trialID', 'bgWin', 'bgWhole', 'emL', 'nRej', 'emStatus', 'nHarm'})]; %#ok<AGROW>
        fprintf('%s %s: %d trials, %.1f min\n', arm, tpl(t).name, numel(R), toc(t0) / 60);
    end
end
writetable(P, fullfile(od, 'templateComparisonBeta_v001.csv'));
S = table();
for arm = opts.Arms
    Pa = P(P.arm == arm, :);
    for t = 1:numel(tpl)
        q = Pa(Pa.template == tpl(t).name, :);
        S = [S; table(arm, tpl(t).name, height(q), nnz(q.emStatus ~= "ok"), sum(q.nRej, 'omitnan'), median(q.nHarm, 'omitnan'), median(q.emL, 'omitnan'), ...
            median(q.bgWin, 'omitnan'), median(q.bgWhole, 'omitnan'), nnz(isfinite(q.bgWin)), ...
            'VariableNames', {'arm', 'template', 'nTrials', 'nEmFail', 'nNoiseRejected', 'medNHarm', 'medEmL', 'medBetaGenWin', 'medBetaGenWhole', 'nInverted'})]; %#ok<AGROW>
    end
    for pr = [["end4", "prod4"]; ["endFc", "end4"]]'
        A = Pa(Pa.template == pr(1), :);  B = Pa(Pa.template == pr(2), :);  [~, ia, ib] = intersect(A.trialID, B.trialID);
        d = A.bgWin(ia) - B.bgWin(ib);  d = d(isfinite(d));
        fprintf('%-13s %s - %s: n %d, median d %+.4f, IQR [%+.4f %+.4f], |d| > 0.03 in %.0f%%\n', arm, pr(1), pr(2), numel(d), median(d), prctile(d, 25), prctile(d, 75), 100 * mean(abs(d) > 0.03));
    end
end
writetable(S, fullfile(od, 'templateComparisonBeta_summary_v001.csv'));
disp(S);
end
