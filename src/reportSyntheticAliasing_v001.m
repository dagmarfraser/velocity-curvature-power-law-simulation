%% reportSyntheticAliasing_v001.m
% TODO_PreSubmission B1: does synthetic power-law noise, subsampled without an anti-aliasing filter, show the alpha drift
% seen in the real trials (+0.187 per octave of rate, Finding #260)? Read-only on results/validateXuNoiseSamplingRate_v001.mat
% (Finding #23): alpha_gen 0:5 at sigma 2 mm, 30 replicates per cell, fs 240/120/60 Hz, fit band [1 Hz, fs/2.2].
%   A  noise at 240 Hz, subsampled (no AA filter) to each rate      B  noise generated independently at each rate
%   C  as A, added to an analytical ellipse
% Per replicate, slope = least-squares slope of fitted alpha on log2(fs) (alpha change per octave of rate). Intervals are
% 95% percentile bootstrap over replicates. The slope is also interpolated linearly in alpha_gen to the median observed
% alpha at 240 Hz (obsBetaRateCutoff_v002 T.alpha), the case to compare with the real +0.187 (R1).
% The band differs from production (min(20, fs/2-1) Hz in characteriseBiologicalNoise_v004), so agreement is a mechanism
% check, not a reproduction. Run from src/. Writes results/reportSyntheticAliasing_v001.mat and a log.

%% CONFIG
ROOT = fileparts(fileparts(mfilename("fullpath")));
IN   = fullfile(ROOT, "results", "validateXuNoiseSamplingRate_v001.mat");
OBS  = fullfile(ROOT, "results", "obsBetaRateCutoff_v002.mat");
OUT  = fullfile(ROOT, "results", "reportSyntheticAliasing_v001.mat");
N_BOOT = 2000;  SEED = 20261029;  COND = ["A" "B" "C"];

%% Setup
logFile = fullfile(ROOT, "logs", "reportSyntheticAliasing_v001_" + string(datetime("now", "Format", "yyyyMMdd_HHmmss")) + ".log");
LG = fopen(logFile, "w");  if LG < 0, error("alias:log", "%s", "Cannot open log: " + logFile); end
closeLog = onCleanup(@() fclose(LG));
for f = [IN OBS], if ~isfile(f), fail_local(LG, "alias:input", "Missing: " + f); end, end
V = load(IN);  O = load(OBS, "T", "r1");
nA = numel(V.alpha_gen);  nF = numel(V.fs_exam);  nR = V.nRep;
if ~isequal(V.fs_exam(:)', [240 120 60]), fail_local(LG, "alias:fs", "fs_exam is not [240 120 60]: " + mat2str(V.fs_exam)); end
for c = COND
    if ~isequal(size(V.("a" + c + "_alpha")), [nR nA nF]) || ~isequal(size(V.("a" + c + "_sigma")), [nR nA nF])
        fail_local(LG, "alias:size", "Condition " + c + " array has an unexpected size"); end
    if any(isnan(V.("a" + c + "_alpha")), "all"), fail_local(LG, "alias:nan", "Condition " + c + " contains NaN fits"); end
end
aObs = median(O.T.alpha(:, 1), "omitnan");  r1 = O.r1.mean;
if aObs < V.alpha_gen(1) || aObs > V.alpha_gen(end), fail_local(LG, "alias:range", "Observed alpha lies outside the simulated range"); end
lo = min(find(V.alpha_gen <= aObs, 1, "last"), nA - 1);
w  = (aObs - V.alpha_gen(lo)) / (V.alpha_gen(lo + 1) - V.alpha_gen(lo));
X  = [ones(nF, 1) log2(V.fs_exam(:))];  rng(SEED);
say_local(LG, "reportSyntheticAliasing_v001 %s\n%d replicates, alpha_gen %s, fs %s.\nObserved median alpha(240) %.3f (interpolating alpha_gen %g to %g, weight %.3f); observed R1 slope %.3f.\n", ...
    string(datetime("now")), nR, mat2str(V.alpha_gen), mat2str(V.fs_exam), aObs, V.alpha_gen(lo), V.alpha_gen(lo + 1), w, r1);

%% Slopes per condition
R = struct();
for c = COND
    al = V.("a" + c + "_alpha");  sg = V.("a" + c + "_sigma");
    sl = NaN(nR, nA);                                          % per-replicate slope, replicates x alpha levels
    for r = 1:nR, for a = 1:nA, b = X \ squeeze(al(r, a, :));  sl(r, a) = b(2); end, end
    est = mean(sl, 1)';  ci = NaN(nA, 2);
    for a = 1:nA, ci(a, :) = prctile(mean(reshape(sl(randi(nR, nR * N_BOOT, 1), a), nR, N_BOOT), 1), [2.5 97.5]); end
    bsLo = mean(reshape(sl(randi(nR, nR * N_BOOT, 1), lo), nR, N_BOOT), 1);
    bsHi = mean(reshape(sl(randi(nR, nR * N_BOOT, 1), lo + 1), nR, N_BOOT), 1);
    bsI  = (1 - w) * bsLo + w * bsHi;
    R.(c) = struct("alphaGen", V.alpha_gen(:), "slope", est, "lo95", ci(:, 1), "hi95", ci(:, 2), ...
        "bias240", squeeze(mean(al(:, :, 1), 1))' - V.alpha_gen(:), "sigmaRatio", squeeze(mean(sg, 1)) ./ V.sigma_gen, ...
        "slopeAtObs", (1 - w) * est(lo) + w * est(lo + 1), "slopeAtObsCI", prctile(bsI, [2.5 97.5]));
    say_local(LG, "\nCONDITION %s: alpha change per octave of rate, by alpha_gen\n%s", c, ...
        evalc("disp(table(R.(c).alphaGen, R.(c).slope, R.(c).lo95, R.(c).hi95, R.(c).bias240, 'VariableNames', [""alphaGen"" ""slope"" ""lo95"" ""hi95"" ""biasAt240""]))"));
    say_local(LG, "  sigma ratio to truth at 240/120/60 Hz, alpha_gen %g to %g: %s\n", V.alpha_gen(lo), V.alpha_gen(lo + 1), ...
        mat2str(R.(c).sigmaRatio([lo lo + 1], :), 4));
    say_local(LG, "  slope interpolated to observed alpha %.2f: %.3f [%.3f, %.3f]   (real trials: %.3f)\n", aObs, ...
        R.(c).slopeAtObs, R.(c).slopeAtObsCI, r1);
end
save(OUT, "R", "aObs", "r1", "lo", "w", "N_BOOT", "SEED");
say_local(LG, "\nSaved %s\nLog %s\n", OUT, logFile);

%% Local functions
function fail_local(fid, id, msg)
fprintf(fid, "FAIL %s: %s\n", id, msg);  error(id, "%s", msg);
end

function say_local(fid, fmt, varargin)
msg = sprintf(fmt, varargin{:});  fprintf("%s", msg);  fprintf(fid, "%s", msg);
end
