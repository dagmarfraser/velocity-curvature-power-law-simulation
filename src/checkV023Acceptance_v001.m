function checkV023Acceptance_v001(opts)
%CHECKV023ACCEPTANCE_V001 Acceptance of runLoopClosureFftnoise_v023 (makeRunLoopClosureV023_v001). Errors on any failure.
%   Small grid (N_BETA 7, N_REPS 3, N_R7 2, per_subject trials of one dataset, serial unless stated); outputs tagged _acc23.
%   1. Defaults (shaped_xu): v023 equals v022 exactly (same machine, both serial).
%   2. endMatched: v023 equals v022 exactly (betaObs, betaObsWhole, betaGenStar, betaGenStarWhole, betaRecCurveMed, emL).
%   3. shapedXuSeg: every trial carries the endMatched fields; same segments, lengths and window beta_obs as endMatched
%      (identical betaObs, betaObsWhole, emL, emSeg); forward maps differ from endMatched (the noise differs).
%   4. shapedXuSeg serial and parfor agree within 1e-6 (curves and beta_gen*; NaN patterns exact), as checkV022Acceptance_v002.
arguments
    opts.Dataset (1,:) char = 'Cook CTRL'
end
here = fileparts(mfilename('fullpath'));  addpath(here, fullfile(here, 'functions'));
c = {'N_BETA', 7, 'N_REPS', 3, 'N_R7', 2, 'ForceRerun', true, 'OutTag', '_acc23'};
eqf = @(a, b, f) all(arrayfun(@(k) isequaln(a(k).(f), b(k).(f)), 1:numel(a)));
d22 = runLoopClosureFftnoise_v022(opts.Dataset, c{:});  d23 = runLoopClosureFftnoise_v023(opts.Dataset, c{:});
for f = ["betaObs", "betaGenStar", "betaGenStarMed", "invertStatus", "betaRecCurveMed"]
    if ~eqf(d22, d23, f), error('acc23:default', 'defaults: field %s differs between v022 and v023', f); end
end
fprintf('1. defaults: v023 = v022 exactly on %d trials\n', numel(d23));
e22 = runLoopClosureFftnoise_v022(opts.Dataset, c{:}, 'NoiseModel', 'endMatched');  e23 = runLoopClosureFftnoise_v023(opts.Dataset, c{:}, 'NoiseModel', 'endMatched');
for f = ["betaObs", "betaObsWhole", "betaGenStar", "betaGenStarWhole", "betaRecCurveMed", "emL"]
    if ~eqf(e22, e23, f), error('acc23:endMatched', 'endMatched: field %s differs between v022 and v023', f); end
end
fprintf('2. endMatched: v023 = v022 exactly on %d trials\n', numel(e23));
s23 = runLoopClosureFftnoise_v023(opts.Dataset, c{:}, 'NoiseModel', 'shapedXuSeg');
need = string(fieldnames(e23))';  miss = need(~isfield(s23, need));
if ~isempty(miss), error('acc23:fields', 'shapedXuSeg lacks fields: %s', strjoin(miss, ", ")); end
for f = ["betaObs", "betaObsWhole", "emL", "emSeg"]
    if ~eqf(e23, s23, f), error('acc23:segments', 'shapedXuSeg and endMatched differ in %s (they must share segments and windows)', f); end
end
dC = arrayfun(@(k) max(abs(e23(k).betaRecCurveMed(:) - s23(k).betaRecCurveMed(:)), [], 'omitnan'), 1:numel(e23));
if all(dC == 0 | isnan(dC)), error('acc23:same', 'shapedXuSeg forward maps equal endMatched: the shaped_xu noise is not being used'); end
ran = [s23.emStatus] == "ok";
fprintf('3. shapedXuSeg: %d of %d trials ran; same segments and windows as endMatched; forward maps differ (median max |d| %.3f)\n', nnz(ran), numel(s23), median(dC, 'omitnan'));
fprintf('   betaGenStarMed endMatched %s\n   betaGenStarMed shapedXuSeg %s\n', mat2str(round([e23.betaGenStarMed], 3)), mat2str(round([s23.betaGenStarMed], 3)));
if isempty(gcp('nocreate')), ensurePool_v001(); end
p23 = runLoopClosureFftnoise_v023(opts.Dataset, c{:}, 'NoiseModel', 'shapedXuSeg', 'UseParfor', true, 'OutTag', '_acc23par');
mxC = 0;  mxG = 0;
for k = 1:numel(s23)
    c1 = s23(k).betaRecCurveMed;  c2 = p23(k).betaRecCurveMed;  g1 = s23(k).betaGenStar;  g2 = p23(k).betaGenStar;
    if ~isequal(isnan(c1), isnan(c2)) || ~isequal(isnan(g1), isnan(g2)), error('acc23:order', 'trial %s: serial and parfor NaN patterns differ', s23(k).trialID); end
    mxC = max([mxC; abs(c1(~isnan(c1)) - c2(~isnan(c2)))]);  mxG = max([mxG; abs(g1(~isnan(g1)) - g2(~isnan(g2)))]);
end
if mxC > 1e-6 || mxG > 1e-6, error('acc23:order', 'shapedXuSeg serial and parfor differ beyond 1e-6 (curves %.3g, beta_gen* %.3g)', mxC, mxG); end
fprintf('4. shapedXuSeg serial and parfor agree (max difference: curves %.3g, beta_gen* %.3g; tolerance 1e-6)\nALL CHECKS PASSED\n', mxC, mxG);
end
