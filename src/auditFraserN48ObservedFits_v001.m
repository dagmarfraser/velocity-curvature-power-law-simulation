function A = auditFraserN48ObservedFits_v001(opts)
%AUDITFRASERN48OBSERVEDFITS_V001  Step 2 audit: the observed-fit block on the N=48 import.
%
%   Re-runs the observed betaObs block of runLoopClosureFftnoise_v019
%   (processOneTrial_local L930-949: differentiateKinematicsEBR, then
%   curvatureKinematicEBR, then regressDataEBR for OLS, LMLS, IRLS with
%   limitBreak = 0 and the OLS fit seeding LMLS/IRLS), with three changes
%   that do not alter any value:
%     - the empty catch records the error identifier and message,
%     - regressDataEBR's disp('Returning NaN ...') is captured with evalc and
%       counted as a non-convergence NaN,
%     - each pipeline's outcome is classed as ok / nonConverged / error.
%
%   Checks:
%     1. Trials importDB_fraser_v002 left unchanged (no gap, no direct touch)
%        must reproduce the v015 corpus's betaObs BIT-EXACTLY
%        (the betaObs reproducibility claim, Finding #204).
%     2. Trials v002 changed (split or direct-touch dropped): betaObs
%        change against v015, per pipeline.
%     3. Failure counts per pipeline and reason, over all N=48 trials.
%
%   Writes results/auditFraserN48ObservedFits_v001.mat.
%
%   Fraser, D.S. (2026)  v001

arguments
    opts.EdgeClip (1,1) double = 20      % runner default (v019 L174)
    opts.Save     (1,1) logical = true
end

srcDir = fileparts(mfilename('fullpath'));
addpath(srcDir); addpath(genpath(fullfile(srcDir, 'functions')));

FS = 240; EC = opts.EdgeClip;
derivConfigs   = struct('filterType', {2, 6}, 'filterParams', {[2 10 1], [4 17]});   % v015/v019 defaults
regressTypes   = [3 4 5];
labels = ["BWFD-OLS","SG-OLS","BWFD-LMLS","SG-LMLS","BWFD-IRLS","SG-IRLS"];

trials = importDB_fraser_v002(Verbose=false);
nT = numel(trials);
betaObs = NaN(nT, 6); vgfObs = NaN(nT, 6);
status  = strings(nT, 6);            % ok | nonConverged | error
errMsg  = strings(nT, 6);

tStart = tic;
parfor k = 1:nT
    tr = trials(k);
    b6 = NaN(2,3); v6 = NaN(2,3); st = strings(2,3); em = strings(2,3);
    for dIdx = 1:2
        [dx, dy] = differentiateKinematicsEBR(tr.x, tr.y, ...
            derivConfigs(dIdx).filterType, derivConfigs(dIdx).filterParams, FS);
        vx = dx(EC:end-EC, 2); vy = dy(EC:end-EC, 2);
        ax = dx(EC:end-EC, 3); ay = dy(EC:end-EC, 3);
        sp = hypot(vx, vy);
        kp = curvatureKinematicEBR(vx, vy, ax, ay);
        lm = [1, -1/3];
        for rIdx = 1:3
            try
                [b, v, out] = fitCapture_local(sp, kp, regressTypes(rIdx), lm);
                b6(dIdx, rIdx) = b; v6(dIdx, rIdx) = v;
                if contains(out, "Returning NaN"), st(dIdx, rIdx) = "nonConverged";
                elseif isnan(b),                   st(dIdx, rIdx) = "nanOther";
                else,                              st(dIdx, rIdx) = "ok";
                end
                if rIdx == 1, lm = [v, b]; end
            catch ME
                st(dIdx, rIdx) = "error";
                em(dIdx, rIdx) = string(ME.identifier) + ": " + string(ME.message);
            end
        end
    end
    betaObs(k, :) = b6(:)'; vgfObs(k, :) = v6(:)';
    status(k, :) = st(:)'; errMsg(k, :) = em(:)';
end
fprintf('[auditFraserN48ObservedFits] %d trials, %.0f s\n', nT, toc(tStart));

%% Failure counts
fprintf('\nOutcome per pipeline (all %d N=48 trials):\n', nT);
for p = 1:6
    fprintf('  %-10s ok %4d  nonConverged %3d  nanOther %3d  error %3d\n', labels(p), ...
        sum(status(:,p)=="ok"), sum(status(:,p)=="nonConverged"), ...
        sum(status(:,p)=="nanOther"), sum(status(:,p)=="error"));
end
u = unique(errMsg(errMsg ~= ""));
if ~isempty(u), fprintf('Error messages caught:\n'); fprintf('  %s\n', u); end

%% Compare with v015
S = load(fullfile(srcDir, 'loopClosureResults_Fraser_all_shaped_xu_v015.mat'), 'results');
R = S.results;
[inV15, loc] = ismember(string({trials.trialID})', string({R.trialID})');
b15 = NaN(nT, 6);
b15(inV15, :) = cell2mat(arrayfun(@(r) r.betaObs(:)', R(loc(inV15)), 'UniformOutput', false)');
changed = ([trials.nGaps] > 0 | [trials.nDirectTouchDropped] > 0)';

clean = inV15 & ~changed;
eqMask = (betaObs == b15) | (isnan(betaObs) & isnan(b15));
nBitExact = sum(all(eqMask(clean, :), 2));
fprintf('\nCheck 1: unchanged trials in v015: %d; betaObs bit-exact on all 6 pipelines: %d\n', ...
    sum(clean), nBitExact);
if nBitExact < sum(clean)
    bad = find(clean & ~all(eqMask, 2));
    fprintf('  NOT bit-exact: %d trials, max |diff| %.3g\n', numel(bad), ...
        max(abs(betaObs(bad,:) - b15(bad,:)), [], 'all', 'omitnan'));
end

ch = inV15 & changed;
d = betaObs(ch, :) - b15(ch, :);
fprintf('\nCheck 2: changed trials in v015: %d. betaObs(v002) - betaObs(v015):\n', sum(ch));
for p = 1:6
    dp = d(:, p);
    fprintf('  %-10s median %+.4f  IQR [%+.4f %+.4f]  |d|>0.03 (MDC): %d  NaN flips: %d\n', labels(p), ...
        median(dp, 'omitnan'), prctile(dp, 25), prctile(dp, 75), sum(abs(dp) > 0.03), ...
        sum(isnan(betaObs(ch,p)) ~= isnan(b15(ch,p))));
end
fprintf('New trials (ID401-403, not in v015): %d\n', sum(~inV15));

A = struct('trialID', {string({trials.trialID})'}, 'betaObs', betaObs, 'vgfObs', vgfObs, ...
    'status', status, 'errMsg', errMsg, 'betaObsV015', b15, 'inV015', inV15, ...
    'changed', changed, 'labels', labels, 'edgeClip', EC, 'runDate', string(datetime('now')));
if opts.Save
    out = fullfile(fileparts(srcDir), 'results', 'auditFraserN48ObservedFits_v001.mat');
    save(out, 'A');
    fprintf('Saved %s\n', out);
end
end

function [b, v, out] = fitCapture_local(sp, kp, regressType, lm)
% Same call as the runner (limitBreak = 0); evalc only captures the disp text.
    b = NaN; v = NaN; %#ok<NASGU>
    out = evalc('[b, v] = regressDataEBR(sp, kp, regressType, lm, 0, 0);');
end
