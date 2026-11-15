function checkV022Acceptance_v001(opts)
%CHECKV022ACCEPTANCE_V001 Acceptance of runLoopClosureFftnoise_v022 (makeRunLoopClosureV022_v001). Errors on any failure.
%   1. Defaults (shaped_xu): v022 equals v021 exactly (isequaln of betaObs, vgfObs, betaGenStar, betaGenStarMed, ciLo, ciHi,
%      invertStatus, betaRecCurveMed) on the same machine and seed, serial, small grid.
%   2. endMatched smoke: runs, every result carries the v022 fields, the template length equals emL where it ran, window
%      beta_obs differs from whole-trial beta_obs, and a second run with the trials in parallel (parfor) gives identical
%      betaGenStar (the per-cell substreams make results independent of execution order).
%   Small grid: N_BETA 7, N_REPS 3, N_R7 2, per_subject trials of one dataset. Outputs are tagged _acc20261115 in src/.
arguments
    opts.Dataset (1,:) char = 'Cook CTRL'
end
here = fileparts(mfilename('fullpath'));  addpath(here, fullfile(here, 'functions'));
common = {'N_BETA', 7, 'N_REPS', 3, 'N_R7', 2, 'ForceRerun', true};
r21 = runLoopClosureFftnoise_v021(opts.Dataset, common{:}, 'OutTag', '_acc20261115');
r22 = runLoopClosureFftnoise_v022(opts.Dataset, common{:}, 'OutTag', '_acc20261115');
f = ["betaObs", "vgfObs", "betaGenStar", "betaGenStarMed", "ciLo", "ciHi", "invertStatus", "betaRecCurveMed"];
if numel(r21) ~= numel(r22), error('acc:n', 'trial counts differ (%d vs %d)', numel(r21), numel(r22)); end
for k = 1:numel(r21)
    for q = f
        if ~isequaln(r21(k).(q), r22(k).(q)), error('acc:default', 'trial %s field %s differs between v021 and v022 at the defaults', r21(k).trialID, q); end
    end
end
fprintf('1. defaults: v022 = v021 exactly on %d trials (%s)\n', numel(r21), strjoin(f, ", "));
e1 = runLoopClosureFftnoise_v022(opts.Dataset, common{:}, 'NoiseModel', 'endMatched', 'OutTag', '_acc20261115');
need = ["betaObsWhole", "vgfObsWhole", "betaGenStarWhole", "betaGenStarMedWhole", "ciLoWhole", "ciHiWhole", "invertStatusWhole", "emL", "emLfrac", "emSeg", "emAlphaSeg", "emWindow", "emStatus", "nNoiseRejected", "betaObsBasis"];
miss = need(~isfield(e1, need));  if ~isempty(miss), error('acc:fields', 'missing fields: %s', strjoin(miss, ", ")); end
ran = [e1.emStatus] == "ok";
fprintf('2. endMatched: %d of %d trials ran (statuses: %s)\n', nnz(ran), numel(e1), strjoin(unique([e1.emStatus]), " | "));
if ~any(ran), error('acc:none', 'no trial ran under endMatched'); end
dW = arrayfun(@(r) max(abs(r.betaObs - r.betaObsWhole)), e1(ran));
fprintf('   window vs whole beta_obs, max |difference| per trial: median %.4f (range %.4f-%.4f)\n', median(dW), min(dW), max(dW));
if all(dW == 0), error('acc:window', 'window beta_obs equals whole-trial beta_obs everywhere: the window is not being used'); end
okW = arrayfun(@(r) r.emWindow(2) - r.emWindow(1) + 1 == r.emL, e1(ran));  if ~all(okW), error('acc:wlen', 'window length differs from emL'); end
fprintf('   noise replicates rejected (no accepted draw): %d of %d\n', sum([e1(ran).nNoiseRejected]), nnz(ran) * 7 * 3);
fprintf('   betaGenStarMed: window %s | whole %s\n', mat2str(round([e1(ran).betaGenStarMed], 3)), mat2str(round([e1(ran).betaGenStarMedWhole], 3)));
if isempty(gcp('nocreate')), ensurePool_v001(); end
e2 = runLoopClosureFftnoise_v022(opts.Dataset, common{:}, 'NoiseModel', 'endMatched', 'UseParfor', true, 'OutTag', '_acc20261115par');
for k = 1:numel(e1)
    if ~isequaln(e1(k).betaGenStar, e2(k).betaGenStar) || ~isequaln(e1(k).betaRecCurveMed, e2(k).betaRecCurveMed)
        error('acc:order', 'trial %s: serial and parfor endMatched runs differ', e1(k).trialID);
    end
end
fprintf('   serial and parfor endMatched runs identical on %d trials\nALL CHECKS PASSED\n', numel(e1));
end
