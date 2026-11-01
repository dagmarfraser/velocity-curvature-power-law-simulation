%% checkSGResidueSensitivity_v001.m
% Sol's closure 1 (2026-10-27): do the headline outputs move when the Savitzky-Golay curvature-residue trials
% (checkSGCurvatureResidue_v001: Fraser 61, Dhieb 13, Hickman PLAC 4, Cook ASD 1) are left out?
% Variants of the pooled within/between tempo decomposition (analyseTempoDecomposition_v002 L81-89, L119-124, unchanged):
%   V0  saved table, Dhieb in             regression anchor: must reproduce the saved D1 (lfw, lfb) to 1e-4.
%   V1  Dhieb out (validity domain, #242) the manuscript's headline domain.
%   V2  Dhieb out, residue trials out     lfMean, lfb, lsa, al re-derived on the subset, as in checkTempoDecompositionDomain_v001.
% Also reported per dataset and pipeline, with and without the residue trials: median beta_obs, median beta_gen*,
% share invertible. Both pipelines lose the same trials so the trial set is common.
% Reads:  src/loopClosureResults_<dataset>_all_shaped_xu_v015.mat, results/tempoDecomposition_v002.mat,
%         results/checkSGCurvatureResidue_v001.mat
% Writes: results/checkSGResidueSensitivity_v001.mat. USAGE: run from src/.
% Fraser, D.S. (2026)  v001

%% CONFIG
ROOT = fileparts(fileparts(mfilename("fullpath")));
OUT_MAT = fullfile(ROOT, "results", "checkSGResidueSensitivity_v001.mat");
DATASETS = ["Fraser" "Cook_CTRL" "Cook_ASD" "Hickman_PLAC" "Hickman_HALO" "Dhieb"];
PL = ["BWFD-OLS" "SG-OLS" "BWFD-LMLS" "SG-LMLS" "BWFD-IRLS" "SG-IRLS"];   % v015 runner order
PIPES = ["SG-IRLS" "BWFD-OLS"];  PAIRED = ["Hickman_PLAC" "Hickman_HALO"];
TOL = 1e-4;
IN_TD = fullfile(ROOT, "results", "tempoDecomposition_v002.mat");
IN_RS = fullfile(ROOT, "results", "checkSGCurvatureResidue_v001.mat");
for f = [IN_TD IN_RS], if ~isfile(f), error("sgSens:input", "%s", "FAILED PATH: " + f); end, end
SV = load(IN_TD, "D1");  RS = load(IN_RS, "perTrial");
resid = strings(0, 1);  for fn = string(fieldnames(RS.perTrial))', resid = [resid; RS.perTrial.(fn).trialID]; end %#ok<AGROW>

%% Trial table, as analyseTempoDecomposition_v002 plus trialID
rows = {};
for d = DATASETS
    f = fullfile(ROOT, "src", "loopClosureResults_" + d + "_all_shaped_xu_v015.mat");
    if ~isfile(f), error("sgSens:input", "%s", "Missing: " + f); end
    R = load(f, "results").results;
    if ~isfield(R, "trialID"), error("sgSens:id", "%s", "No trialID in " + f); end
    person = d;  if ismember(d, PAIRED), person = "Hickman"; end
    for i = 1:numel(R)
        r = R(i);
        for p = PIPES
            k = find(PL == p);
            rows{end+1, 1} = struct("dataset", d, "trialID", string(r.trialID), "session", d + "_" + string(r.subjectID), ...
                "person", person + "_" + string(r.subjectID), "pipeline", p, "f0", r.f0, "sigma", r.sigmaMM, ...
                "a", r.a_mm, "alpha", r.alphaIRA, "betaObs", r.betaObs(k), "betaGen", r.betaGenStar(k), ...
                "invertible", r.invertStatus(k) == "rise"); %#ok<SAGROW>
        end
    end
end
T0 = struct2table(vertcat(rows{:}));
T0 = T0(isfinite(T0.f0) & T0.f0 > 0 & isfinite(T0.sigma) & isfinite(T0.a) & T0.a > 0 & isfinite(T0.alpha), :);
T0.isResid = ismember(T0.trialID, resid);
nMiss = sum(~ismember(resid, T0.trialID));
fprintf("Residue trials listed %d; present in the table %d (%d absent: dropped by the finite filters or outside the domain)\n", numel(resid), numel(resid) - nMiss, nMiss);

%% Variants
V0 = prep_local(T0, PIPES);                                   % saved-table construction
V1 = prep_local(T0(T0.dataset ~= "Dhieb", :), PIPES);
V2 = prep_local(T0(T0.dataset ~= "Dhieb" & ~T0.isResid, :), PIPES);
VNAME = ["V0 Dhieb in (anchor)" "V1 Dhieb out" "V2 Dhieb out, residue trials out"];  V = {V0, V1, V2};
D1 = cell(1, 3);
for v = 1:3, D1{v} = fitD1_local(V{v}, PIPES); end

%% Regression anchor
A = D1{1}(ismember(D1{1}.term, ["lfw" "lfb"]), :);  B = SV.D1(ismember(SV.D1.term, ["lfw" "lfb"]), :);
key = @(X) X.pipeline + "|" + X.outcome + "|" + X.term;
[ok, loc] = ismember(key(B), key(A));
if ~all(ok) || any(abs(A.b(loc) - B.b) > TOL) || any(abs(A.SE(loc) - B.SE) > TOL)
    error("sgSens:anchor", "V0 does not reproduce the saved D1 (lfw, lfb)");
end
fprintf("REGRESSION ANCHOR passed: V0 reproduces the saved D1.\n");

%% Report
SUM = table();
for v = 1:3
    fprintf("\n=== %s: %d trial x pipeline rows ===\n", VNAME(v), height(V{v}));
    for p = PIPES
        w = D1{v}(D1{v}.pipeline == p & D1{v}.term == "lfw", :);
        bo = w.b(w.outcome == "beta_obs");  bg = w.b(w.outcome == "beta_gen*");
        seo = w.SE(w.outcome == "beta_obs");  seg = w.SE(w.outcome == "beta_gen*");
        fprintf("   %-8s within-session slope: beta_obs %+.4f (SE %.4f), beta_gen* %+.4f (SE %.4f); reduction %.1f%%\n", p, bo, seo, bg, seg, 100 * (1 - bg / bo));
        SUM = [SUM; table(VNAME(v), p, bo, seo, bg, seg, 100 * (1 - bg / bo), 'VariableNames', ["variant" "pipeline" "bObs" "seObs" "bGen" "seGen" "reductionPct"])]; %#ok<AGROW>
    end
end
disp(SUM)
fprintf("\nChange V1 to V2 (SG-IRLS, BWFD-OLS): slope beta_obs, beta_gen*, reduction (pp)\n");
for p = PIPES
    a = SUM(SUM.variant == VNAME(2) & SUM.pipeline == p, :);  b = SUM(SUM.variant == VNAME(3) & SUM.pipeline == p, :);
    fprintf("   %-8s %+.4f  %+.4f  %+.1f\n", p, b.bObs - a.bObs, b.bGen - a.bGen, b.reductionPct - a.reductionPct);
end

% Per dataset: medians and invertible share with all trials and without residue trials
DS = table();
for d = unique(T0.dataset, "stable")'
    for p = PIPES
        X = T0(T0.dataset == d & T0.pipeline == p, :);  Y = X(~X.isResid, :);
        DS = [DS; table(d, p, height(X), sum(X.isResid), median(X.betaObs, "omitnan"), median(Y.betaObs, "omitnan"), ...
            median(X.betaGen, "omitnan"), median(Y.betaGen, "omitnan"), 100 * mean(X.invertible), 100 * mean(Y.invertible), ...
            'VariableNames', ["dataset" "pipeline" "n" "nResid" "medObsAll" "medObsDrop" "medGenAll" "medGenDrop" "pctInvAll" "pctInvDrop"])]; %#ok<AGROW>
    end
end
DS.dMedObs = DS.medObsDrop - DS.medObsAll;  DS.dMedGen = DS.medGenDrop - DS.medGenAll;  DS.dPctInv = DS.pctInvDrop - DS.pctInvAll;
disp(DS)
fprintf("Largest |change|: median beta_obs %.4f, median beta_gen* %.4f, invertible share %.2f pp\n", ...
    max(abs(DS.dMedObs)), max(abs(DS.dMedGen)), max(abs(DS.dPctInv)));

runDate = string(datetime("now", "Format", "yyyy-MM-dd"));
save(OUT_MAT, "D1", "SUM", "DS", "VNAME", "resid", "runDate");
fprintf("Saved: %s\n", OUT_MAT);

%% =========================================================================
function T = prep_local(T, pipes)
% analyseTempoDecomposition_v002 L62-72, unchanged: centring and z-scoring on the subset.
    T.lf = log2(T.f0);
    G = groupsummary(T(T.pipeline == pipes(1), :), "session", "mean", "lf");
    [~, loc] = ismember(T.session, G.session);  T.lfMean = G.mean_lf(loc);
    T.lfw = T.lf - T.lfMean;  T.lfb = T.lfMean - mean(G.mean_lf);
    z = @(x) (x - mean(x, "omitnan")) ./ std(x, "omitnan");
    T.lsa = z(log2(T.sigma ./ T.a));  T.al = z(T.alpha);
    T.datasetC = removecats(categorical(T.dataset));  T.sessionC = removecats(categorical(T.session));
end

function D1 = fitD1_local(T, pipes)
% analyseTempoDecomposition_v002 L81-89, unchanged.
    terms = ["lfw" "lfb"];  D1 = table();
    for p = pipes
        V = T(T.pipeline == p, :);
        mo = fitlme(V(isfinite(V.betaObs), :), "betaObs ~ 1 + lfw + lfb + lsa + al + datasetC + (1|sessionC)");
        mg = fitlme(V(isfinite(V.betaGen), :), "betaGen ~ 1 + lfw + lfb + lsa + al + datasetC + (1|sessionC)");
        D1 = [D1; tag_local(coef_local(mo, terms), p, "beta_obs", sum(isfinite(V.betaObs))); ...
                  tag_local(coef_local(mg, terms), p, "beta_gen*", sum(isfinite(V.betaGen)))]; %#ok<AGROW>
    end
end

function t = coef_local(mdl, terms)
    c = mdl.Coefficients;  k = ismember(string(c.Name), terms);
    t = table(string(c.Name(k)), round(c.Estimate(k), 4), round(c.SE(k), 4), 'VariableNames', ["term" "b" "SE"]);
end

function t = tag_local(c, p, outcome, n)
    t = [table(repmat(string(p), height(c), 1), repmat(string(outcome), height(c), 1), repmat(n, height(c), 1), ...
        'VariableNames', ["pipeline" "outcome" "nTrials"]), c];
end
