%% reportSamplingRateEffects_v001.m
% TODO_PreSubmission B1 / Rerun D16: the numbers behind the sampling-rate paragraph of section 2.3, read from the saved
% results/obsBetaRateCutoff_v002.mat only (S1: clean subsampling x(1:D:end), no filter). Nothing is recomputed upstream.
% Production arms A240p, A120p, A60p. Person-level percentile bootstrap (persons resampled with replacement, all their
% trials kept), 95% intervals. Reports:
%   levels   median alpha, sigma (mm) and beta_obs (six pipelines) at 240, 120 and 60 Hz;
%   changes  paired per-trial change against 240 Hz: alpha (difference), sigma (ratio), beta_obs (difference), and the
%            share of trials whose |change in beta_obs| reaches the D16 tolerance (SPEC_D16 R3, 0.0108);
%   clean    the same on trials with no SG curvature residue (RS = 0) in all three production arms.
% Anchors (fail loud): the saved R1 slope, the alpha medians, the Q1 beta_obs medians (BWFD-OLS, SG-OLS, SG-IRLS) must
% equal Finding #260's published values to their published digits, and the 240 Hz sigma median must equal the v015
% median quoted in TODO_PreSubmission (1.71 mm). Sigma at 120 and 60 Hz has no published value to anchor to; it rests on
% obsBetaRateCutoff_v002's own anchor B (240 Hz noise reproduces noiseCharacterisation_fraser.mat).
% Run from src/. Reads results/obsBetaRateCutoff_v002.mat. Writes results/reportSamplingRateEffects_v001.mat and a log.

%% CONFIG
ROOT   = fileparts(fileparts(mfilename("fullpath")));
IN     = fullfile(ROOT, "results", "obsBetaRateCutoff_v002.mat");
OUT    = fullfile(ROOT, "results", "reportSamplingRateEffects_v001.mat");
PROD   = ["A240p" "A120p" "A60p"];  FS = [240 120 60];
N_BOOT = 2000;  SEED = 20261029;  TOL_D16 = 0.0108;
% Published (Finding #260; TODO_PreSubmission): alpha and sigma to 2 d.p., beta_obs and R1 to 3 d.p.
PUB_R1 = 0.187;  PUB_ALPHA = [4.38 4.25 3.96];  PUB_SIGMA240 = 1.71;  TOL2 = 0.00501;  TOL3 = 0.000501;
PUB_P  = [1 2 6];                              % rows below: BWFD-OLS, SG-OLS, SG-IRLS (PIPE order); columns 240, 120, 60 Hz
PUB_B  = [0.187 0.187 0.190; 0.186 0.186 0.171; 0.181 0.182 0.158];

%% Setup
logFile = fullfile(ROOT, "logs", "reportSamplingRateEffects_v001_" + string(datetime("now", "Format", "yyyyMMdd_HHmmss")) + ".log");
LG = fopen(logFile, "w");  if LG < 0, error("rate:log", "%s", "Cannot open log: " + logFile); end
closeLog = onCleanup(@() fclose(LG));
say_local(LG, "reportSamplingRateEffects_v001 run %s\n", string(datetime("now")));
if ~isfile(IN), fail_local(LG, "rate:input", "Missing: " + IN); end
S = load(IN, "B", "RS", "SET", "PIPE", "T", "r1");
need = ["B" "RS" "SET" "PIPE" "T" "r1"];  miss = need(~isfield(S, need));
if ~isempty(miss), fail_local(LG, "rate:fields", "Saved file lacks: " + strjoin(miss, ", ")); end
T = S.T;  nT = height(T);  PIPE = S.PIPE;  pn = replace(PIPE, "-", "_");
[~, sIdx] = ismember(PROD, [S.SET.name]);
if any(sIdx == 0), fail_local(LG, "rate:sets", "A production arm is missing from the saved settings"); end
if ~isequal(size(S.B), [numel(S.SET) nT 6]), fail_local(LG, "rate:size", "B has size " + mat2str(size(S.B))); end
A = T.alpha;  Sg = T.sigmaMM;
if ~isequal(size(A), [nT 3]) || ~isequal(size(Sg), [nT 3]), fail_local(LG, "rate:size", "alpha or sigma is not trials x 3 rates"); end
Bp = permute(S.B(sIdx, :, :), [2 3 1]);        % trials x pipelines x rates
[up, ~, pIdx] = unique(string(T.person));  byP = accumarray(pIdx, (1:nT)', [], @(v) {v});
say_local(LG, "%d trials, %d persons; %d bootstrap resamples, seed %d.\n", nT, numel(up), N_BOOT, SEED);

%% Anchors against the published values
gotB = squeeze(median(Bp(:, PUB_P, :), 1, "omitnan"));  gotA = median(A, 1, "omitnan");  gotS = median(Sg(:, 1), "omitnan");
if any(abs(gotB - PUB_B) > TOL3, "all") || any(abs(gotA - PUB_ALPHA) > TOL2) || abs(gotS - PUB_SIGMA240) > TOL2 || abs(S.r1.mean - PUB_R1) > TOL3
    fail_local(LG, "rate:anchor", sprintf("Saved results do not reproduce the published values (alpha %s, sigma240 %.4f, R1 %.4f)", mat2str(gotA, 4), gotS, S.r1.mean));
end
say_local(LG, "Anchors passed: alpha medians, beta_obs medians (3 pipelines), sigma at 240 Hz and the R1 slope match Finding #260 / TODO_PreSubmission.\n");

%% Levels and paired changes (all trials), then the residue-free subset
rng(SEED);
keep = all(S.RS(sIdx, :) == 0, 1)';  nKeep = sum(keep);
byPc = cellfun(@(v) v(keep(v)), byP, "UniformOutput", false);
med  = @(M) median(M, 1, "omitnan");
share = @(M) sum(abs(M) >= TOL_D16, 1) ./ sum(isfinite(M), 1);
lev  = [A, Sg, reshape(Bp, nT, 18)];
levN = [("alpha_" + FS), ("sigmaMM_" + FS), bNames_local(pn, FS, 1:3)];
dB   = reshape(Bp(:, :, 2:3) - Bp(:, :, 1), nT, 12);
chg  = [A(:, 2:3) - A(:, 1), Sg(:, 2:3) ./ Sg(:, 1), dB];
chgN = ["dAlpha_" + FS(2:3), "sigmaRatio_" + FS(2:3), "d" + bNames_local(pn, FS, 2:3)];
Lev = summ_local(lev, byP, N_BOOT, med, levN);   Chg = summ_local(chg, byP, N_BOOT, med, chgN);
Shr = summ_local(dB, byP, N_BOOT, share, "shareAtTol_" + chgN(5:end));
LevC = summ_local(lev, byPc, N_BOOT, med, levN); ChgC = summ_local(chg, byPc, N_BOOT, med, chgN);
ShrC = summ_local(dB, byPc, N_BOOT, share, "shareAtTol_" + chgN(5:end));
say_local(LG, "\nLEVELS, medians over trials (all %d trials):\n%s", nT, evalc("disp(Lev)"));
say_local(LG, "\nPAIRED CHANGES against 240 Hz, median per-trial change (sigma as a ratio):\n%s", evalc("disp(Chg)"));
say_local(LG, "\nSHARE of trials with |change in beta_obs| >= %.4f:\n%s", TOL_D16, evalc("disp(Shr)"));
say_local(LG, "\nRESIDUE-FREE SUBSET (RS = 0 in all three arms; %d of %d trials):\n%s%s%s", nKeep, nT, evalc("disp(LevC)"), evalc("disp(ChgC)"), evalc("disp(ShrC)"));
save(OUT, "Lev", "Chg", "Shr", "LevC", "ChgC", "ShrC", "nKeep", "nT", "SEED", "N_BOOT", "TOL_D16", "PROD", "FS");
say_local(LG, "\nSaved %s\nLog %s\n", OUT, logFile);

%% Local functions
function names = bNames_local(pn, FS, qs)   % beta_obs column names, pipeline fastest within rate, as reshape(Bp, nT, []) orders them
names = strings(1, numel(pn) * numel(qs));  k = 0;
for q = qs, for p = 1:numel(pn), k = k + 1;  names(k) = "betaObs_" + pn(p) + "_" + FS(q); end, end
end

function t = summ_local(X, byP, nBoot, fn, names)   % estimate and person-bootstrap 95% interval of fn(X), column by column
rows = vertcat(byP{:});  est = fn(X(rows, :));  nP = numel(byP);  bs = NaN(nBoot, numel(est));
for b = 1:nBoot, pick = randi(nP, nP, 1);  bs(b, :) = fn(X(vertcat(byP{pick}), :)); end
ci = prctile(bs, [2.5 97.5], 1);  n = sum(isfinite(X(rows, :)), 1);
t = table(names(:), est(:), ci(1, :)', ci(2, :)', n(:), 'VariableNames', ["quantity" "estimate" "lo95" "hi95" "n"]);
end

function fail_local(fid, id, msg)
fprintf(fid, "FAIL %s: %s\n", id, msg);  error(id, "%s", msg);
end

function say_local(fid, fmt, varargin)
msg = sprintf(fmt, varargin{:});  fprintf("%s", msg);  fprintf(fid, "%s", msg);
end
