%% checkEmpiricalTempoVsGridVGF_v002.m
% Where does each empirical trial sit on the registered grid's VGF axis, at its OWN
% generator exponent?  Successor to v001 (Finding #240), which converted every trial's
% f0 with one constant, K_conv = 175.2636 evaluated at beta = 1/3. That constant is 5.28%
% below the grid's true tempo constant K(1/3) = 185.02 (testGridKConv_v001), and K itself
% falls steeply with beta (d ln K / d beta about -5.3), so a single beta was the larger
% problem. v002 places each trial exactly:
%       VGF_match_i = f0_i * K(beta_i),   inside the grid iff VGF_LO <= VGF_match <= VGF_HI,
% with beta_i the trial's own exponent (gridKConv_v001 requires beta; no default).
%
% Exponent per trial (declared choices, Dagmar 2026-09-29):
%   PRIMARY      betaGenStarMed  the constellation median of the six pointwise-inversion
%                                estimates (generator scale, at the trial's own tempo).
%                                A trial with no finite value is NOT imputed: it is counted
%                                ("no exponent") and excluded from every share below.
%   SHIFT        the primary exponent moved by -/+ SHIFT (0.03, one MDC), same trials.
%   REF-SUBSET   beta = 1/3 on exactly the trials that have a primary exponent, so the
%                selection effect (trials without an exponent are the slow, non-invertible
%                ones, Finding #227) is separable from the beta effect.
%   OBS          median(betaObs) over the six pipelines. NOT a generator exponent: betaObs
%                is pulled towards the noise exponent (Part 1B), so this only shows how much
%                the placement depends on beta. Negative values are outside K's domain and are
%                excluded and counted, not clipped.
%   REFERENCE    beta = 1/3 for every trial, the corrected constant. Reported so the old
%                v001 convention can be compared like for like; it is a stated convention,
%                not an estimate.
% Zarandi and Dhieb are outside the validity domain (#242): treated by the same rule and
% flagged (their beta_gen* is not a generator estimate; the placement needs only an exponent).
%
% SELF-CHECKS (Fail Loud): K(1/3) reproduces 185.0242; the five in-domain trial counts equal
%   Table 6a's all-trial Ns (2829, 94, 102, 359, 338; total 3722, coherence item A4.3).
% Reads:  src/loopClosureResults_<dataset>_all_shaped_xu_v015.mat (results struct array).
% Writes: results/checkEmpiricalTempoVsGridVGF_v002.mat  (results, perTrial, constants).
% USAGE:  from the project root: checkEmpiricalTempoVsGridVGF_v002
% Fraser, D.S. (2026)  v002

%% CONFIG
ROOT      = fileparts(fileparts(mfilename("fullpath")));
DATASETS  = ["Fraser","Cook_CTRL","Cook_ASD","Hickman_PLAC","Hickman_HALO","Zarandi","Dhieb"];
IN_DOMAIN = [true true true true true false false];              % #242
N_TABLE6A = [2829 94 102 359 338 NaN NaN];                       % all-trial Ns; out-of-domain not asserted
BETA_REF  = 1/3;
SHIFT     = 0.03;                                                % one MDC
VGF_LO    = exp(4.5);   VGF_HI = exp(5.8);                        % registered grid VGF range
K13_PUB   = 185.0242;   K13_TOL = 1e-3;
OUT_MAT   = fullfile(ROOT, "results", "checkEmpiricalTempoVsGridVGF_v002.mat");
addpath(ROOT); addpath(fullfile(ROOT, "src")); addpath(genpath(fullfile(ROOT, "src", "functions")));

%% Self-check: the constant
K13 = gridKConv_v001(BETA_REF);
if abs(K13 - K13_PUB) > K13_TOL
    error("tempoVsGrid2:K13", "%s", sprintf("K(1/3) %.4f does not reproduce %.4f", K13, K13_PUB));
end

%% Per-trial placement
perTrial = table();
for d = 1:numel(DATASETS)
    f = fullfile(ROOT, "src", "loopClosureResults_" + DATASETS(d) + "_all_shaped_xu_v015.mat");
    if ~isfile(f), error("tempoVsGrid2:NoCorpus", "%s", "FAILED PATH: " + f); end
    R = load(f, "results").results;
    n = numel(R);
    if ~isnan(N_TABLE6A(d)) && n ~= N_TABLE6A(d)
        error("tempoVsGrid2:Count", "%s", sprintf("%s: %d trials, Table 6a says %d", DATASETS(d), n, N_TABLE6A(d)));
    end
    f0   = arrayfun(@(r) double(r.f0), R(:));
    bGen = arrayfun(@(r) double(r.betaGenStarMed), R(:));
    bObs = arrayfun(@(r) median(double(r.betaObs), "omitnan"), R(:));     % NaN if all six NaN
    ids  = string({R.trialID}');
    T = table(repmat(DATASETS(d), n, 1), ids, f0, bGen, bObs, ...
        'VariableNames', ["dataset" "trialID" "f0" "betaGenStarMed" "betaObsMed"]);
    perTrial = [perTrial; T]; %#ok<AGROW>
end

% VGF_match under each exponent rule; NaN where the rule gives no usable exponent.
usable = @(b) isfinite(b) & b >= 0;
perTrial.hasF0    = isfinite(perTrial.f0) & perTrial.f0 > 0;
perTrial.vgfPrim  = nan(height(perTrial), 1);
perTrial.vgfSens  = nan(height(perTrial), 1);
perTrial.vgfRef   = nan(height(perTrial), 1);
perTrial.vgfDn    = nan(height(perTrial), 1);
perTrial.vgfUp    = nan(height(perTrial), 1);
ok = perTrial.hasF0 & usable(perTrial.betaGenStarMed);
perTrial.vgfPrim(ok) = perTrial.f0(ok) .* gridKConv_v001(perTrial.betaGenStarMed(ok));
ok = perTrial.hasF0 & usable(perTrial.betaObsMed);
perTrial.vgfSens(ok) = perTrial.f0(ok) .* gridKConv_v001(perTrial.betaObsMed(ok));
ok = perTrial.hasF0;
perTrial.vgfRef(ok)  = perTrial.f0(ok) * K13;
ok = perTrial.hasF0 & usable(perTrial.betaGenStarMed);
perTrial.vgfUp(ok)   = perTrial.f0(ok) .* gridKConv_v001(perTrial.betaGenStarMed(ok) + SHIFT);
okDn = ok & perTrial.betaGenStarMed >= SHIFT;                       % K needs beta >= 0
perTrial.vgfDn(okDn) = perTrial.f0(okDn) .* gridKConv_v001(perTrial.betaGenStarMed(okDn) - SHIFT);
perTrial.vgfRefSub = nan(height(perTrial), 1);
perTrial.vgfRefSub(ok) = perTrial.f0(ok) * K13;

%% Dataset summary
rows = numel(DATASETS);
S = cell(rows, 1);
for d = 1:rows
    s = perTrial.dataset == DATASETS(d);
    P = perTrial(s, :);
    nAll = height(P);
    r = struct("dataset", DATASETS(d), "inDomain", IN_DOMAIN(d), "nTrials", nAll, ...
        "nNoF0", nnz(~P.hasF0), "nNoExpPrim", nnz(P.hasF0 & ~usable(P.betaGenStarMed)), ...
        "nNoExpSens", nnz(P.hasF0 & ~usable(P.betaObsMed)));
    for tag = ["Prim" "Sens" "Ref" "RefSub" "Dn" "Up"]
        v = P.("vgf" + tag);  v = v(isfinite(v));
        r.("n" + tag)         = numel(v);
        r.("medVGF" + tag)    = median(v);
        r.("pctBelow" + tag)  = 100 * mean(v < VGF_LO);
        r.("pctAbove" + tag)  = 100 * mean(v > VGF_HI);
    end
    r.nBelowZeroShift = nnz(P.hasF0 & usable(P.betaGenStarMed) & P.betaGenStarMed < SHIFT);
    r.medF0       = median(P.f0(P.hasF0));
    r.medBetaGen  = median(P.betaGenStarMed, "omitnan");
    r.medBetaObs  = median(P.betaObsMed, "omitnan");
    S{d} = r;
end
results = struct2table(vertcat(S{:}));
if any(results.nNoF0 > 0)
    warning("tempoVsGrid2:NoF0", "%s", sprintf("Trials with no usable f0 (excluded everywhere): %s", ...
        strjoin(DATASETS(results.nNoF0 > 0) + "=" + string(results.nNoF0(results.nNoF0 > 0))', ", ")));
end
in5 = results.inDomain;
fprintf("In-domain trials (Table 6a check): %d = %d\n", sum(results.nTrials(in5)), sum(N_TABLE6A(in5)));

%% Report
f0Lo13 = VGF_LO / K13;  f0Hi13 = VGF_HI / K13;
fprintf("\nK(1/3) = %.4f; grid VGF %.2f-%.2f = tempo %.3f-%.3f Hz at beta = 1/3 (v001 said 0.514-1.885 with 175.2636)\n", ...
    K13, VGF_LO, VGF_HI, f0Lo13, f0Hi13);
fprintf("\nPLACEMENT OF TRIALS ON THE VGF AXIS (percent of trials below / above the grid's VGF range)\n");
fprintf("PRIMARY beta_i = betaGenStarMed; trials without one are excluded (noExp), not imputed.\n");
fprintf("%-13s %3s %5s %5s %6s %7s %8s | %13s | %13s\n", "dataset", "dom", "n", "noExp", "nPrim", "medBeta", "medVGF", "PRIMARY b/a", "REFERENCE b/a");
for d = 1:rows
    q = results(d, :);
    fprintf("%-13s %3d %5d %5d %6d %7.3f %8.1f | %6.1f %6.1f | %6.1f %6.1f\n", q.dataset, q.inDomain, q.nTrials, ...
        q.nNoExpPrim, q.nPrim, q.medBetaGen, q.medVGFPrim, q.pctBelowPrim, q.pctAbovePrim, q.pctBelowRef, q.pctAboveRef);
end
fprintf("\nDECOMPOSITION of %%below (same trials unless stated): REF(all) -> REF(subset) is selection; REF(subset) -> PRIMARY is beta.\n");
fprintf("%-13s %9s %9s %9s | %9s %9s   (shift +/-%.2f)\n", "dataset", "REF(all)", "REF(sub)", "PRIMARY", "beta-0.03", "beta+0.03", SHIFT);
for d = 1:rows
    q = results(d, :);
    fprintf("%-13s %9.1f %9.1f %9.1f | %9.1f %9.1f\n", q.dataset, q.pctBelowRef, q.pctBelowRefSub, q.pctBelowPrim, q.pctBelowDn, q.pctBelowUp);
end
fprintf("\nSame decomposition for %%above:\n");
for d = 1:rows
    q = results(d, :);
    fprintf("%-13s %9.1f %9.1f %9.1f | %9.1f %9.1f\n", q.dataset, q.pctAboveRef, q.pctAboveRefSub, q.pctAbovePrim, q.pctAboveDn, q.pctAboveUp);
end
fprintf("\nOBS (betaObs is NOT a generator exponent; shows beta dependence only): nObs, medBetaObs, medVGF, %%below, %%above\n");
for d = 1:rows
    q = results(d, :);
    fprintf("%-13s %5d %7.3f %8.1f %7.1f %7.1f\n", q.dataset, q.nSens, q.medBetaObs, q.medVGFSens, q.pctBelowSens, q.pctAboveSens);
end

%% Save
constants = struct("K13", K13, "VGF_LO", VGF_LO, "VGF_HI", VGF_HI, "BETA_REF", BETA_REF, ...
    "runDate", string(datetime("now", "Format", "yyyy-MM-dd")));
save(OUT_MAT, "results", "perTrial", "constants");
fprintf("\nSaved %s\n", OUT_MAT);
