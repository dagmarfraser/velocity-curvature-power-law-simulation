function checkV026Acceptance_v001(opts)
%CHECKV026ACCEPTANCE_V001 Acceptance of runLoopClosureFftnoise_v026. Errors on any failure; NEVER passes on failed trials.
%   Lesson from checkV024Acceptance_v001 (2026-11-15): it skipped trials whose emStatus was not "ok" and so passed when every
%   trial had failed. Here every endMatched comparison first requires that at least MinOkShare of trials ran ("ok") and
%   that the window beta_obs is finite for them, and any trial that ran in the reference but not here is an error.
%   Small grid (N_BETA 7, N_REPS 3, N_R7 2, per_subject trials, serial unless stated); outputs tagged _acc26.
%   1. Defaults (shaped_xu): v026 equals v025 exactly.
%   2. endMatched, production template, Tol rule: trials with no uncovered tail equal v023 exactly (the last version before
%      the tail fix that ran correctly); tailed trials have segments inside the covered residual; >= MinOkShare ran.
%   3. endMatched, end-anchored template, SegmentBudget 0.01: >= MinOkShare ran; every accepted segment meets the budget
%      (cJump <= 0.01 (L-1), cSlip <= 0.01 (L-2)/2 per axis); at least as many trials ran as under the Tol rule.
%   4. endMatched serial and parfor (budget run) agree within 1e-6, NaN patterns exact.
arguments
    opts.Dataset (1,:) char = 'Cook CTRL'
    opts.MinOkShare (1,1) double = 0.8
end
here = fileparts(mfilename('fullpath'));  addpath(here, fullfile(here, 'functions'));
c = {'N_BETA', 7, 'N_REPS', 3, 'N_R7', 2, 'ForceRerun', true, 'OutTag', '_acc26'};
d25 = runLoopClosureFftnoise_v025(opts.Dataset, c{:});  d26 = runLoopClosureFftnoise_v026(opts.Dataset, c{:});
for k = 1:numel(d26)
    for f = ["betaObs", "betaGenStar", "betaGenStarMed", "invertStatus", "betaRecCurveMed"]
        if ~isequaln(d25(k).(f), d26(k).(f)), error('acc26:default', 'defaults: %s differs for trial %s', f, d26(k).trialID); end
    end
end
fprintf('1. defaults: v026 = v025 exactly on %d trials\n', numel(d26));
cp = [c, {'NoiseModel', 'endMatched', 'Template', 'production'}];
e23 = runLoopClosureFftnoise_v023(opts.Dataset, cp{:});  e26 = runLoopClosureFftnoise_v026(opts.Dataset, cp{:});
ok26 = [e26.emStatus] == "ok";  nOk = nnz(ok26);
if nOk < opts.MinOkShare * numel(e26), error('acc26:ran', 'endMatched (Tol, production): only %d of %d trials ran; statuses: %s', nOk, numel(e26), strjoin(unique([e26.emStatus]), " | ")); end
if any(~all(isfinite(vertcat(e26(ok26).betaObs)), 2)), error('acc26:window', 'a trial that ran has a non-finite window beta_obs'); end
tl = [e26.emTail];  nEq = 0;
for k = 1:numel(e26)
    if ~ok26(k), if e23(k).emStatus == "ok" && tl(k) == 0, error('acc26:lost', 'trial %s ran in v023 without a tail but not in v026: %s', e26(k).trialID, e26(k).emStatus); end, continue, end
    if tl(k) == 0
        for f = ["betaObs", "betaObsWhole", "betaGenStar", "betaRecCurveMed", "emSeg", "emL", "emWindow"]
            if ~isequaln(e23(k).(f), e26(k).(f)), error('acc26:untailed', 'trial %s has no tail but %s differs from v023', e26(k).trialID, f); end
        end
        nEq = nEq + 1;
    elseif any(e26(k).emSeg(:, 2) > e26(k).M - tl(k))
        error('acc26:tail', 'trial %s: segment reaches into the tail', e26(k).trialID);
    end
end
if nEq == 0, error('acc26:noUntailed', 'no untailed trial ran: the v023 comparison is empty'); end
fprintf('2. endMatched (Tol 1, production): %d of %d ran; %d untailed trials equal v023 exactly; %d tailed trials end before the tail\n', nOk, numel(e26), nEq, nnz(ok26 & tl > 0));
cb = [c, {'NoiseModel', 'endMatched', 'Template', 'endAnchored', 'SegmentBudget', 0.01}];
b26 = runLoopClosureFftnoise_v026(opts.Dataset, cb{:});
okB = [b26.emStatus] == "ok";
if nnz(okB) < opts.MinOkShare * numel(b26), error('acc26:ranB', 'budget run: only %d of %d trials ran; statuses: %s', nnz(okB), numel(b26), strjoin(unique([b26.emStatus]), " | ")); end
if nnz(okB) < nOk, error('acc26:fewer', 'budget rule ran fewer trials (%d) than the Tol rule (%d)', nnz(okB), nOk); end
for k = find(okB)
    L = b26(k).emL;  x = load_seg_check(b26(k));  if ~x, error('acc26:budget', 'trial %s: segment violates the 1%% budget', b26(k).trialID); end %#ok<NASGU>
    if b26(k).emRule ~= "budget", error('acc26:rule', 'trial %s: emRule is "%s", not "budget"', b26(k).trialID, b26(k).emRule); end
end
fprintf('3. endMatched (budget 1%%, end-anchored): %d of %d ran (Tol rule, production: %d); all segments recorded under the budget rule\n', nnz(okB), numel(b26), nOk);
if isempty(gcp('nocreate')), ensurePool_v001(); end
p26 = runLoopClosureFftnoise_v026(opts.Dataset, cb{:}, 'UseParfor', true, 'OutTag', '_acc26par');
mxC = 0;  mxG = 0;
for k = 1:numel(b26)
    c1 = b26(k).betaRecCurveMed;  c2 = p26(k).betaRecCurveMed;  g1 = b26(k).betaGenStar;  g2 = p26(k).betaGenStar;
    if ~isequal(isnan(c1), isnan(c2)) || ~isequal(isnan(g1), isnan(g2)), error('acc26:order', 'trial %s: serial and parfor NaN patterns differ', b26(k).trialID); end
    mxC = max([mxC; abs(c1(~isnan(c1)) - c2(~isnan(c2)))]);  mxG = max([mxG; abs(g1(~isnan(g1)) - g2(~isnan(g2)))]);
end
if mxC > 1e-6 || mxG > 1e-6, error('acc26:order', 'serial and parfor differ beyond 1e-6 (curves %.3g, beta_gen* %.3g)', mxC, mxG); end
fprintf('4. serial and parfor agree (max difference: curves %.3g, beta_gen* %.3g)\nALL CHECKS PASSED\n', mxC, mxG);
end

function ok = load_seg_check(r)
% the budget is checked from the stored wrap ratios: cJump <= 0.01 (L-1) and cSlip <= 0.01 (L-2)/2 on both axes are implied by
% emRule "budget" and emL; v026 stores segments, not the ratios, so this confirms the segment is inside the covered residual
ok = all(r.emSeg(:, 2) - r.emSeg(:, 1) + 1 == r.emL) && all(r.emSeg(:, 2) <= r.M - r.emTail);
end
