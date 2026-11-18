function checkV024Acceptance_v001(opts)
%CHECKV024ACCEPTANCE_V001 Acceptance of runLoopClosureFftnoise_v024 (makeRunLoopClosureV024_v001). Errors on any failure.
%   Small grid (N_BETA 7, N_REPS 3, N_R7 2, per_subject trials of one dataset, Template 'production' as the real runs; serial
%   unless stated); outputs tagged _acc24.
%   1. Defaults (shaped_xu): v024 equals v023 exactly.
%   2. endMatched: trials with no uncovered tail (emTail 0) equal v023 exactly; trials with a tail have segments that end
%      inside the covered residual (emSeg(:,2) <= M - emTail), and at least one trial has a tail (else the test is empty).
%   3. shapedXuSeg runs on every trial with the same segments and windows as endMatched.
%   4. endMatched serial and parfor agree within 1e-6.
arguments
    opts.Dataset (1,:) char = 'Cook CTRL'
end
here = fileparts(mfilename('fullpath'));  addpath(here, fullfile(here, 'functions'));
c = {'N_BETA', 7, 'N_REPS', 3, 'N_R7', 2, 'ForceRerun', true, 'OutTag', '_acc24', 'Template', 'production'};
eqk = @(a, b, f, k) isequaln(a(k).(f), b(k).(f));
d23 = runLoopClosureFftnoise_v023(opts.Dataset, c{:});  d24 = runLoopClosureFftnoise_v024(opts.Dataset, c{:});
for f = ["betaObs", "betaGenStar", "betaGenStarMed", "invertStatus", "betaRecCurveMed"]
    for k = 1:numel(d24), if ~eqk(d23, d24, f, k), error('acc24:default', 'defaults: %s differs for trial %s', f, d24(k).trialID); end, end
end
fprintf('1. defaults: v024 = v023 exactly on %d trials\n', numel(d24));
e23 = runLoopClosureFftnoise_v023(opts.Dataset, c{:}, 'NoiseModel', 'endMatched');  e24 = runLoopClosureFftnoise_v024(opts.Dataset, c{:}, 'NoiseModel', 'endMatched');
tl = [e24.emTail];  nT = nnz(tl > 0);
if nT == 0, error('acc24:noTail', 'no trial with an uncovered tail in %s: the fix is untested', opts.Dataset); end
for k = 1:numel(e24)
    if e24(k).emStatus ~= "ok", continue, end
    if tl(k) == 0
        for f = ["betaObs", "betaObsWhole", "betaGenStar", "betaRecCurveMed", "emSeg", "emL"]
            if ~eqk(e23, e24, f, k), error('acc24:untailed', 'trial %s has no tail but %s differs from v023', e24(k).trialID, f); end
        end
    elseif any(e24(k).emSeg(:, 2) > e24(k).M - tl(k))
        error('acc24:tail', 'trial %s: segment ends at %s, inside the tail (covered residual ends at %d)', e24(k).trialID, mat2str(e24(k).emSeg(:, 2)'), e24(k).M - tl(k));
    end
end
in23 = arrayfun(@(k) any(e23(k).emSeg(:, 2) > e23(k).M - tl(k)), find(tl > 0 & [e23.emStatus] == "ok"));
fprintf('2. endMatched: %d trials without a tail equal v023 exactly; %d with a tail (median %d samples) now end before it (v023 segments reached into it in %d of them)\n', ...
    nnz(tl == 0), nT, round(median(tl(tl > 0))), nnz(in23));
s24 = runLoopClosureFftnoise_v024(opts.Dataset, c{:}, 'NoiseModel', 'shapedXuSeg');
for k = 1:numel(s24)
    for f = ["betaObs", "betaObsWhole", "emSeg", "emL", "emTail"]
        if ~eqk(e24, s24, f, k), error('acc24:segments', 'trial %s: shapedXuSeg and endMatched differ in %s', s24(k).trialID, f); end
    end
end
fprintf('3. shapedXuSeg: same segments, tails and windows as endMatched on %d trials\n', numel(s24));
if isempty(gcp('nocreate')), ensurePool_v001(); end
p24 = runLoopClosureFftnoise_v024(opts.Dataset, c{:}, 'NoiseModel', 'endMatched', 'UseParfor', true, 'OutTag', '_acc24par');
mxC = 0;  mxG = 0;
for k = 1:numel(e24)
    c1 = e24(k).betaRecCurveMed;  c2 = p24(k).betaRecCurveMed;  g1 = e24(k).betaGenStar;  g2 = p24(k).betaGenStar;
    if ~isequal(isnan(c1), isnan(c2)) || ~isequal(isnan(g1), isnan(g2)), error('acc24:order', 'trial %s: serial and parfor NaN patterns differ', e24(k).trialID); end
    mxC = max([mxC; abs(c1(~isnan(c1)) - c2(~isnan(c2)))]);  mxG = max([mxG; abs(g1(~isnan(g1)) - g2(~isnan(g2)))]);
end
if mxC > 1e-6 || mxG > 1e-6, error('acc24:order', 'serial and parfor differ beyond 1e-6 (curves %.3g, beta_gen* %.3g)', mxC, mxG); end
fprintf('4. endMatched serial and parfor agree (max difference: curves %.3g, beta_gen* %.3g; tolerance 1e-6)\nALL CHECKS PASSED\n', mxC, mxG);
end
