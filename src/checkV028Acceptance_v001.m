function checkV028Acceptance_v001(opts)
%CHECKV028ACCEPTANCE_V001 Acceptance of runLoopClosureFftnoise_v028. Errors on any failure; never passes on failed trials.
%   Small grid (N_BETA 7, N_REPS 3, N_R7 2, per_subject trials); outputs tagged _acc28.
%   1. Defaults (shaped_xu): v028 equals v027 exactly.
%   2. endMatched, end-anchored, SegmentBudget 0.01 (rule auto): v028 equals v027 exactly, and >= MinOkShare of trials ran.
%   3. endMatched, end-anchored, SegmentRule "detect": >= MinOkShare ran, every ran trial has emRule "detect", a finite
%      window beta_obs, a segment of length emL inside the covered residual, and at least as many trials ran as under 1%.
%   4. shapedXuSeg with "detect" shares segments, windows and emRule with endMatched "detect".
%   5. endMatched "detect": serial and parfor agree within 1e-6, NaN patterns exact.
arguments
    opts.Dataset (1,:) char = 'Cook CTRL'
    opts.MinOkShare (1,1) double = 0.8
end
here = fileparts(mfilename('fullpath'));  addpath(here, fullfile(here, 'functions'));
c = {'N_BETA', 7, 'N_REPS', 3, 'N_R7', 2, 'ForceRerun', true, 'OutTag', '_acc28'};
d27 = runLoopClosureFftnoise_v027(opts.Dataset, c{:});  d28 = runLoopClosureFftnoise_v028(opts.Dataset, c{:});
eqAll(d27, d28, ["betaObs", "betaGenStar", "betaGenStarMed", "invertStatus", "betaRecCurveMed"], 'defaults');
fprintf('1. defaults: v028 = v027 exactly on %d trials\n', numel(d28));
cb = [c, {'NoiseModel', 'endMatched', 'Template', 'endAnchored', 'SegmentBudget', 0.01}];
b27 = runLoopClosureFftnoise_v027(opts.Dataset, cb{:});  b28 = runLoopClosureFftnoise_v028(opts.Dataset, cb{:});
okB = [b28.emStatus] == "ok";  needOk(okB, opts.MinOkShare, '1% budget');
eqAll(b27, b28, ["betaObs", "betaObsWhole", "betaGenStar", "betaRecCurveMed", "emSeg", "emL", "emWindow"], '1% budget');
fprintf('2. 1%% budget (auto): v028 = v027 exactly on %d trials, %d ran\n', numel(b28), nnz(okB));
cdt = [c, {'NoiseModel', 'endMatched', 'Template', 'endAnchored', 'SegmentRule', 'detect'}];
e28 = runLoopClosureFftnoise_v028(opts.Dataset, cdt{:});
okD = [e28.emStatus] == "ok";  needOk(okD, opts.MinOkShare, 'detect');
if nnz(okD) < nnz(okB), error('acc28:fewer', 'detect ran fewer trials (%d) than the 1%% budget (%d)', nnz(okD), nnz(okB)); end
for k = find(okD)
    r = e28(k);
    if r.emRule ~= "detect", error('acc28:rule', 'trial %s: emRule "%s"', r.trialID, r.emRule); end
    if ~all(isfinite(r.betaObs)), error('acc28:window', 'trial %s: non-finite window beta_obs', r.trialID); end
    if any(r.emSeg(:, 2) - r.emSeg(:, 1) + 1 ~= r.emL) || any(r.emSeg(:, 2) > r.M - r.emTail), error('acc28:seg', 'trial %s: segment length or position wrong', r.trialID); end
end
fprintf('3. detect: %d of %d ran (1%% budget: %d); all with emRule "detect", finite window beta_obs and valid segments\n', nnz(okD), numel(e28), nnz(okB));
cs = [cdt, {'NoiseModel', 'shapedXuSeg'}];  s28 = runLoopClosureFftnoise_v028(opts.Dataset, cs{:});
eqAll(e28, s28, ["betaObs", "betaObsWhole", "emSeg", "emL", "emWindow", "emRule"], 'shapedXuSeg vs endMatched (detect)');
fprintf('4. shapedXuSeg (detect): same segments, windows and rule as endMatched on %d trials\n', numel(s28));
if isempty(gcp('nocreate')), ensurePool_v001(); end
cp = [cdt, {'UseParfor', true, 'OutTag', '_acc28par'}];  p28 = runLoopClosureFftnoise_v028(opts.Dataset, cp{:});
mxC = 0;  mxG = 0;
for k = 1:numel(e28)
    c1 = e28(k).betaRecCurveMed;  c2 = p28(k).betaRecCurveMed;  g1 = e28(k).betaGenStar;  g2 = p28(k).betaGenStar;
    if ~isequal(isnan(c1), isnan(c2)) || ~isequal(isnan(g1), isnan(g2)), error('acc28:order', 'trial %s: serial and parfor NaN patterns differ', e28(k).trialID); end
    mxC = max([mxC; abs(c1(~isnan(c1)) - c2(~isnan(c2)))]);  mxG = max([mxG; abs(g1(~isnan(g1)) - g2(~isnan(g2)))]);
end
if mxC > 1e-6 || mxG > 1e-6, error('acc28:order', 'serial and parfor differ beyond 1e-6 (curves %.3g, beta_gen* %.3g)', mxC, mxG); end
fprintf('5. serial and parfor agree (curves %.3g, beta_gen* %.3g)\nALL CHECKS PASSED\n', mxC, mxG);
end

function needOk(ok, share, label)
if nnz(ok) < share * numel(ok), error('acc28:ran', '%s: only %d of %d trials ran', label, nnz(ok), numel(ok)); end
end

function eqAll(a, b, fields, label)
if numel(a) ~= numel(b), error('acc28:n', '%s: trial counts differ', label); end
for k = 1:numel(a)
    for f = fields
        if ~isequaln(a(k).(f), b(k).(f)), error('acc28:equal', '%s: trial %s field %s differs', label, a(k).trialID, f); end
    end
end
end
