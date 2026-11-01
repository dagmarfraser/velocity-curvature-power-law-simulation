%% reportSamplingRateLmm_v001.m
% TODO_PreSubmission B1 (closes it): the sampling-rate terms of the Stage 1 LMM and SEM adequacy by fs x pipeline, on the
% registered grid. Read-only. Run from src/.
%   LMM   results/L9ModelInspect_v002.mat (out.fixed: name, est, SE; N = 17,458,535). Model: deltaBeta ~ betaGenerated *
%         VGF * samplingRate * filterType * regressionType * noiseMagnitude * noiseColor + (1|paramComboID). betaGenerated,
%         VGF, noiseMagnitude and noiseColor were z-scored; samplingRate was NOT (Finding on Dev-2), so fs terms are per Hz
%         and the implied slopes hold at the grid means of the z-scored predictors. Reference: BWFD, OLS. filterType_6 = SG,
%         regressionType_4 = LMLS, regressionType_5 = IRLS. The saved file holds standard errors but no covariance, so an
%         SE is given for single terms only; pipeline slopes are sums of terms and carry no SE.
%   SEM   src/perCoordinateSEM_v2_001.mat (coordTable: one row per pipeline x betaGen x VGF x fs x alpha x sigma, SEM over 5
%         repetitions). Strata as stored: adequate (SEM <= 0.0108), marginal (<= 0.0200), inadequate. The share is over
%         grid coordinates, so it depends on the registered alpha 0-6 and sigma 0-20 mm ranges.
% Not regenerated here: the fs-invariant edge clip (Finding #188), which discards 15.4/7.7/3.8% of samples at 60/120/240 Hz
% and is confounded with every fs term below.
% Writes results/reportSamplingRateLmm_v001.mat and a dated log.

%% CONFIG
ROOT = fileparts(fileparts(mfilename("fullpath")));
LMM  = fullfile(ROOT, "results", "L9ModelInspect_v002.mat");
SEMF = fullfile(ROOT, "src", "perCoordinateSEM_v2_001.mat");
OUT  = fullfile(ROOT, "results", "reportSamplingRateLmm_v001.mat");
FS   = [60 120 240];  PIPE = ["BWFD-OLS" "SG-OLS" "BWFD-LMLS" "SG-LMLS" "BWFD-IRLS" "SG-IRLS"];
SG   = [0 1 0 1 0 1];  REG = [0 0 4 4 5 5];  NOBS = 17458535;

%% Setup
logFile = fullfile(ROOT, "logs", "reportSamplingRateLmm_v001_" + string(datetime("now", "Format", "yyyyMMdd_HHmmss")) + ".log");
LG = fopen(logFile, "w");  if LG < 0, error("lmm:log", "%s", "Cannot open log: " + logFile); end
closeLog = onCleanup(@() fclose(LG));
for f = [LMM SEMF], if ~isfile(f), fail_local(LG, "lmm:input", "Missing: " + f); end, end
say_local(LG, "reportSamplingRateLmm_v001 %s\n", string(datetime("now")));

%% LMM: fs terms
L = load(LMM, "out");  fx = L.out.fixed;
if L.out.nObs ~= NOBS || ~contains(string(L.out.verdict), "clean"), fail_local(LG, "lmm:model", "Inspect file is not the clean N = 17,458,535 fit"); end
nm = string(fx.name);
get = @(t) lookup_local(LG, nm, fx, t);
[b0, s0] = get("samplingRate");  [bF, sF] = get("samplingRate:filterType_6");
[bL, sL] = get("samplingRate:regressionType_4");  [bI, sI] = get("samplingRate:regressionType_5");
[bFL] = get("samplingRate:filterType_6:regressionType_4");  [bFI] = get("samplingRate:filterType_6:regressionType_5");
slope = b0 + SG * bF + (REG == 4) * bL + (REG == 5) * bI + SG .* (REG == 4) * bFL + SG .* (REG == 5) * bFI;
lmmTerms = table(["samplingRate" "samplingRate:filterType_6" "samplingRate:regressionType_4" "samplingRate:regressionType_5"]', ...
    [b0 bF bL bI]', [s0 sF sL sI]', 'VariableNames', ["term" "estimatePerHz" "SE"]);
lmmTerms.perHalvingFrom240 = lmmTerms.estimatePerHz * 120;  lmmTerms.span60to240 = lmmTerms.estimatePerHz * 180;
impl = table(PIPE(:), slope(:), slope(:) * 180, 'VariableNames', ["pipeline" "slopePerHz" "impliedChange60to240"]);
say_local(LG, "\nLMM sampling-rate terms (reference BWFD-OLS; deltaBeta per Hz; others at grid mean):\n%s", evalc("disp(lmmTerms)"));
say_local(LG, "\nImplied fs slope by pipeline (sum of terms, no SE: covariance not saved):\n%s", evalc("disp(impl)"));

%% SEM adequacy by fs x pipeline
C = load(SEMF, "coordTable");  C = C.coordTable;
need = ["pipeline" "fs" "sem" "semStratum"];  miss = need(~ismember(need, string(C.Properties.VariableNames)));
if ~isempty(miss), fail_local(LG, "lmm:sem", "coordTable lacks: " + strjoin(miss, ", ")); end
pp = string(C.pipeline);  st = string(C.semStratum);
if ~isempty(setdiff(unique(pp), PIPE')) || ~isempty(setdiff(unique(C.fs), FS')), fail_local(LG, "lmm:grid", "Unexpected pipelines or rates"); end
rows = table();
for p = PIPE, for f = FS
    k = pp == p & C.fs == f;  n = sum(k);
    rows = [rows; table(p, f, n, mean(st(k) == "adequate"), mean(st(k) == "marginal"), mean(st(k) == "inadequate"), median(C.sem(k), "omitnan"), ...
        'VariableNames', ["pipeline" "fs" "nCoordinates" "shareAdequate" "shareMarginal" "shareInadequate" "medianSEM"])]; %#ok<AGROW>
end, end
say_local(LG, "\nSEM strata by fs x pipeline (share of grid coordinates):\n%s", evalc("disp(rows)"));
byFs = table(FS(:), arrayfun(@(f) mean(st(C.fs == f) == "adequate"), FS(:)), 'VariableNames', ["fs" "shareAdequate"]);
say_local(LG, "\nAll pipelines pooled:\n%s", evalc("disp(byFs)"));
save(OUT, "lmmTerms", "impl", "rows", "byFs");
say_local(LG, "\nSaved %s\nLog %s\n", OUT, logFile);

%% Local functions
function [b, s] = lookup_local(fid, nm, fx, term)
k = find(nm == term);
if numel(k) ~= 1, fail_local(fid, "lmm:term", "Term not found exactly once: " + term); end
b = fx.est(k);  s = fx.SE(k);
end

function fail_local(fid, id, msg)
fprintf(fid, "FAIL %s: %s\n", id, msg);  error(id, "%s", msg);
end

function say_local(fid, fmt, varargin)
msg = sprintf(fmt, varargin{:});  fprintf("%s", msg);  fprintf(fid, "%s", msg);
end
