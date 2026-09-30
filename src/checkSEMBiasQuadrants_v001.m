%% checkSEMBiasQuadrants_v001.m
% PR7 (coherence pass v003_v001): the registered dual assessment (prereg §3.2) classes a
% pipeline configuration by systematic bias and SEM jointly ("High bias, low SEM", "Low bias,
% high SEM", "High bias, high SEM"). Finding #3 gave the joint split (about 49% precise but
% biased, about 31% precise and accurate) on the pre-#186 basis (threshold 0.011, the 79.8%
% mixed denominator). This script restates it on Finding #186's basis:
%   denominator  the evaluable coordinates of perCoordinateSEM_v2_001 (3,493,466);
%   precision    SEM < MDC/2.77, from semAdequacyThreshold_v001 (0.0108303);
%   bias         |mean(beta_rec - beta_gen)| > MDC (0.03), the threshold
%                visualiseSEMBiasQuadrants_v001/v2_001 used, stated here explicitly.
% Reported pooled and per pipeline, and on the registered noise range (alpha <= 3) as well
% as the full extended grid.
% SELF-CHECK (Fail Loud): the adequate count reproduces #186 (2,872,808 of 3,493,466).
% Reads:  src/perCoordinateSEM_v2_001.mat (coordTable)
% Writes: results/checkSEMBiasQuadrants_v001.mat
% USAGE:  from the project root: checkSEMBiasQuadrants_v001
% Fraser, D.S. (2026)  v001

%% CONFIG
ROOT = fileparts(fileparts(mfilename("fullpath")));
addpath(fullfile(ROOT, "src")); addpath(genpath(fullfile(ROOT, "src", "functions")));
MDC = 0.03;
SEM_ADEQ = semAdequacyThreshold_v001();
EXPECT = [2872808 3493466];                            % Finding #186
ALPHA_REG = 3;                                         % prereg §6.1/§8.1 noise range
OUT_MAT = fullfile(ROOT, "results", "checkSEMBiasQuadrants_v001.mat");
f = fullfile(ROOT, "src", "perCoordinateSEM_v2_001.mat");
if ~isfile(f), error("quad:input", "%s", "FAILED PATH: " + f); end
T = load(f, "coordTable").coordTable;

%% Self-check against #186
nAdeq = nnz(T.sem < SEM_ADEQ);
if ~isequal([nAdeq height(T)], EXPECT)
    error("quad:anchor", "%s", sprintf("adequate %d of %d; #186 says %d of %d", nAdeq, height(T), EXPECT));
end
fprintf("SELF-CHECK passed: %d of %d evaluable coordinates adequate (%.2f%%), threshold %.7f.\n", ...
    nAdeq, height(T), 100 * nAdeq / height(T), SEM_ADEQ);
nNaNBias = nnz(~isfinite(T.meanBias));
if nNaNBias > 0
    warning("quad:nanBias", "%s", sprintf("%d coordinates have no finite mean bias; counted in no quadrant", nNaNBias));
end

%% Quadrants
Q = @(M) quad_local(M, MDC, SEM_ADEQ);
R = [Q(T); Q(T(T.alpha <= ALPHA_REG, :))];
R.basis = ["all evaluable coordinates"; "registered noise range, alpha <= 3"];
pipes = unique(T.pipeline, "stable");
P = table();
for p = pipes'
    r = Q(T(T.pipeline == p, :));  r.pipeline = string(p);  P = [P; r]; %#ok<AGROW>
end
fprintf("\nPercent of coordinates (precise = SEM < MDC/2.77; biased = |mean bias| > MDC = %.2f):\n", MDC);
disp(R(:, ["basis" "n" "pctPreciseAccurate" "pctPreciseBiased" "pctImpreciseAccurate" "pctImpreciseBiased" "pctBiasedGivenPrecise"]))
fprintf("Per pipeline, all evaluable coordinates:\n");
disp(P(:, ["pipeline" "n" "pctPreciseAccurate" "pctPreciseBiased" "pctImpreciseAccurate" "pctImpreciseBiased" "pctBiasedGivenPrecise"]))

runDate = string(datetime("now", "Format", "yyyy-MM-dd"));
save(OUT_MAT, "R", "P", "MDC", "SEM_ADEQ", "runDate");
fprintf("Saved: %s\n", OUT_MAT);

%% =========================================================================
function r = quad_local(M, mdc, semT)
    ok = isfinite(M.meanBias) & isfinite(M.sem);
    b = abs(M.meanBias(ok)) > mdc;  s = M.sem(ok) < semT;  n = nnz(ok);
    r = table(n, 100 * mean(s & ~b), 100 * mean(s & b), 100 * mean(~s & ~b), 100 * mean(~s & b), 100 * mean(b(s)), ...
        'VariableNames', ["n" "pctPreciseAccurate" "pctPreciseBiased" "pctImpreciseAccurate" "pctImpreciseBiased" "pctBiasedGivenPrecise"]);
end
