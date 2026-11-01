%% checkSGCurvatureResidue_v001.m
% Floating-point curvature residue in the SG pipelines (found by obsBetaRateCutoff_v001's anchor A, 2026-10-24).
% With quantised positions, where a coordinate is constant over a whole SG window the curvature is rounding residue
% of zero (below RESID), its log is an extreme point, and SG fits then depend on the machine's rounding.
% (1) Machine, Fraser: iMac beta_obs (stage S1, A240p) against BlueBEAR beta_obs (v015), trials with and without residue.
% (2) Every dataset at the production settings (as runLoopClosureFftnoise_v015: SG filterType 6 [4 17], edge clip 20):
%     trials with residue; for them, the SG fits with the residue samples (as published) and without them (the
%     samples exact arithmetic would make log(0)), and the published v015 value.
% (3) Fraser cutoff sweep of stage S1: residue trials per setting and the SG-IRLS change when residue samples are removed.
% Reads results/obsBetaRateCutoff_v001.mat and src/loopClosureResults_<dataset>_all_shaped_xu_v015.mat.
% Writes results/checkSGCurvatureResidue_v001.mat and a dated log in logs/. Run from src/.

%% CONFIG
ROOT  = fileparts(fileparts(mfilename("fullpath")));  SRC = fullfile(ROOT, "src");
IN_S1 = fullfile(ROOT, "results", "obsBetaRateCutoff_v001.mat");
OUT   = fullfile(ROOT, "results", "checkSGCurvatureResidue_v001.mat");
DS = struct("name", {"Fraser" "Cook_CTRL" "Cook_ASD" "Hickman_PLAC" "Hickman_HALO" "Zarandi" "Dhieb"}, ...
    "fs", {240 133 133 133 133 100 100}, ...
    "imp", {@() importDB_fraser_v001(Verbose=false), @() importDB_cook_v002(Group="CTRL", Tasks=7, Verbose=false), ...
            @() importDB_cook_v002(Group="ASD", Tasks=7, Verbose=false), @() importDB_hickman_v003(Study=2, Group="PLAC", Verbose=false), ...
            @() importDB_hickman_v003(Study=2, Group="HALO", Verbose=false), @() importDB_zarandi_v001(Verbose=false), ...
            @() importDB_dhieb_v001(Verbose=false)});            % importers and rates as runLoopClosureFftnoise_v015's switch
EDGE = 20;  REG = [3 4 5];  TOL = 1e-6;  SEM = 0.0108;

%% Setup
logFile = fullfile(ROOT, "logs", "checkSGCurvatureResidue_v001_" + string(datetime("now", "Format", "yyyyMMdd_HHmmss")) + ".log");
LG = fopen(logFile, "w");  if LG < 0, error("sgr:log", "%s", "Cannot open log: " + logFile); end
say_local(LG, "checkSGCurvatureResidue_v001 run %s\n", string(datetime("now")));
addpath(SRC);  addpath(genpath(fullfile(SRC, "functions")));  addpath(genpath(fullfile(SRC, "req")));
if ~isfile(IN_S1), error("sgr:input", "%s", "Missing: " + IN_S1); end
S1 = load(IN_S1, "B", "RS", "RESID", "SET", "T");  RESID = S1.RESID;  ids = S1.T.trialID;   % v001 saves the trial IDs in T

%% (1) Machine, Fraser
V = load(fullfile(SRC, "loopClosureResults_Fraser_all_shaped_xu_v015.mat"), "results").results;
if ~isequal(ids, string({V.trialID})'), error("sgr:ids", "%s", "S1 and v015 trial order differ"); end
s = find([S1.SET.name] == "A240p");  d = abs(squeeze(S1.B(s, :, :)) - vertcat(V.betaObs));  rz = S1.RS(s, :)' > 0;  M = table();
for g = [false true]
    k = rz == g;
    M = [M; table(g, sum(k), max(d(k, :), [], 1, "omitnan"), sum(any(d(k, :) > TOL, 2)), 'VariableNames', ["residue" "nTrials" "maxDiff" "nBeyondTol"])]; %#ok<AGROW>
end
say_local(LG, "(1) Fraser, |beta_obs iMac (S1 A240p) - beta_obs BlueBEAR (v015)|, pipelines BWFD-OLS SG-OLS BWFD-LMLS SG-LMLS BWFD-IRLS SG-IRLS:\n%s", evalc("disp(M)"));

%% (2) Every dataset at the production settings
P = table();  perTrial = struct();
for q = 1:numel(DS)
    V = load(fullfile(SRC, "loopClosureResults_" + DS(q).name + "_all_shaped_xu_v015.mat"), "results").results;
    tr = DS(q).imp();  [ok, loc] = ismember(string({V.trialID}), string({tr.trialID}));
    if ~all(ok), error("sgr:import", "%s", sprintf("%s: %d v015 trials not imported", DS(q).name, sum(~ok))); end
    tr = tr(loc);  n = numel(tr);  nres = zeros(n, 1);  keep = NaN(n, 3);  drop = NaN(n, 3);  fs = DS(q).fs;
    parfor k = 1:n
        [sp, kp] = sgSpeedCurv_local(double(tr(k).x), double(tr(k).y), [4 17], fs, EDGE);
        nres(k) = sum(kp < RESID);
        if nres(k) > 0, keep(k, :) = fits_local(sp, kp, REG); g = kp >= RESID; drop(k, :) = fits_local(sp(g), kp(g), REG); end
    end
    pub = vertcat(V.betaObs);  pub = pub(:, [2 4 6]);  r = nres > 0;
    if any(r) && median(abs(keep(r, :) - pub(r, :)), "all", "omitnan") > 0.01, error("sgr:keep", "%s", DS(q).name + ": refits far from v015 (wrong trials or settings?)"); end
    ch = drop(r, :) - keep(r, :);  mc = abs(keep(r, :) - pub(r, :));
    P = [P; table(string(DS(q).name), n, sum(r), 100 * mean(r), colstat_local(abs(ch), "max"), colstat_local(abs(ch), "median"), ...
        colstat_local(mc, "max"), sum(abs(ch(:, 3)) >= SEM), ...
        'VariableNames', ["dataset" "nTrials" "nResidue" "pctResidue" "maxAbsChangeIfDropped" "medianAbsChangeIfDropped" "maxDiffFromV015" "nSGIRLSchangeAtLeastSEM"])]; %#ok<AGROW>
    perTrial.(DS(q).name) = table(string({V(r).trialID})', nres(r), keep(r, :), drop(r, :), pub(r, :), 'VariableNames', ["trialID" "nResidueSamples" "withResidue" "withoutResidue" "v015"]);
    say_local(LG, "%-13s %4d trials at %3d Hz: residue in %d\n", DS(q).name, n, fs, sum(r));
end
say_local(LG, "\n(2) Production settings; SG-OLS, SG-LMLS, SG-IRLS columns; 'changeIfDropped' = refit without the residue samples:\n%s", evalc("disp(P)"));

%% (3) Fraser cutoff sweep of stage S1 (240 Hz, clip 41)
trF = importDB_fraser_v001(Shapes=3, Visits=[1 2], Verbose=false);  [~, loc] = ismember(ids, string({trF.trialID}));  trF = trF(loc);
W = table();
for s = find(startsWith([S1.SET.name], "S240"))
    ref = S1.SET(s).ref;  k0 = find(S1.RS(s, :) > 0);  dIR = NaN(numel(k0), 1);  aIR = NaN(numel(k0), 1);  trs = trF(k0);
    parfor j = 1:numel(k0)
        [sp, kp] = sgSpeedCurv_local(double(trs(j).x), double(trs(j).y), [4 ref], 240, 41);
        g = kp >= RESID;  a = fits_local(sp, kp, REG);  b = fits_local(sp(g), kp(g), REG);  aIR(j) = a(3);  dIR(j) = b(3) - a(3);
    end
    if max(abs(aIR - S1.B(s, k0, 6)'), [], "omitnan") > TOL, error("sgr:sweep", "%s", S1.SET(s).name + ": refits with residue do not reproduce S1"); end
    W = [W; table(string(S1.SET(s).name), numel(k0), median(abs(dIR), "omitnan"), max(abs(dIR), [], "omitnan"), sum(abs(dIR) >= SEM), ...
        'VariableNames', ["setting" "nResidue" "medianAbsChangeSGIRLS" "maxAbsChangeSGIRLS" "nAtLeastSEM"])]; %#ok<AGROW>
end
say_local(LG, "\n(3) Fraser sweep (S1): SG-IRLS change when residue samples are removed:\n%s", evalc("disp(W)"));
save(OUT, "M", "P", "W", "perTrial", "RESID", "TOL", "SEM");
say_local(LG, "Saved: %s\nLog: %s\n", OUT, logFile);  fclose(LG);

%% =========================================================================
function [sp, kp] = sgSpeedCurv_local(x, y, prm, fs, clip)
% Speed and curvature as runLoopClosureFftnoise_v015 processOneTrial_local (SG, filterType 6).
    [dx, dy] = differentiateKinematicsEBR(x, y, 6, prm, fs);
    vx = dx(clip:end-clip, 2);  vy = dy(clip:end-clip, 2);  ax = dx(clip:end-clip, 3);  ay = dy(clip:end-clip, 3);
    sp = hypot(vx, vy);  kp = curvatureKinematicEBR(vx, vy, ax, ay);
end

function b = fits_local(sp, kp, REG)
% OLS, then LMLS and IRLS started from the OLS values, as processOneTrial_local.
    b = NaN(1, 3);  lm = [1, -1/3];
    for r = 1:3
        try
            [bb, v] = regressDataEBR(sp, kp, REG(r), lm, 0, 0);  b(r) = bb;
            if r == 1, lm = [v, bb]; end
        catch
        end
    end
end

function v = colstat_local(A, f)
% Column statistic over rows, NaN(1, 3) when there are no rows.
    if isempty(A), v = NaN(1, 3); elseif f == "max", v = max(A, [], 1, "omitnan"); else, v = median(A, 1, "omitnan"); end
end

function say_local(fid, fmt, varargin)
    fprintf(fmt, varargin{:});  fprintf(fid, fmt, varargin{:});
end
