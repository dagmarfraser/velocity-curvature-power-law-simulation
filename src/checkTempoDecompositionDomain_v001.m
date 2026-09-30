%% checkTempoDecompositionDomain_v001.m
% analyseTempoDecomposition_v002 (Session 123) predates Finding #242: its DATASETS (L43) still
% holds Dhieb, so the pooled within/between decomposition (section 1) and the selection GLMM
% (section 2a, the manuscript's "logit +2.28 per doubling") include Dhieb's 90 trials.
% This script refits exactly those pooled models on the #242 validity domain, from the saved
% trial table, the same way Finding #243 re-ran #231/#234:
%   V0  the saved table as is (six datasets)   regression anchor: must reproduce the saved D1 and S2
%                                              (b, SE, df) to their stored precision.
%   V1  Dhieb removed                          lfb re-centred and lsa, al re-z-scored on the subset.
% Sections 1 (per dataset), 2b, 3 and 3b are per dataset or per contrast (Cook, Hickman) and do not
% pool Dhieb; their coefficients on group/drug are invariant to the pooled centring and z-scoring,
% so they are not refitted here.
% Reads:  results/tempoDecomposition_v002.mat
% Writes: results/checkTempoDecompositionDomain_v001.mat
% USAGE:  from the project root: checkTempoDecompositionDomain_v001
% Fraser, D.S. (2026)  v001

%% CONFIG
ROOT   = fileparts(fileparts(mfilename("fullpath")));
IN_MAT = fullfile(ROOT, "results", "tempoDecomposition_v002.mat");
OUT_MAT = fullfile(ROOT, "results", "checkTempoDecompositionDomain_v001.mat");
OUT_DOM = "Dhieb";                                            % Finding #242 (Zarandi was never in)
TOL    = struct("b", 1e-4, "SE", 1e-4);
if ~isfile(IN_MAT), error("tempoDom:input", "%s", "FAILED PATH: " + IN_MAT); end
S = load(IN_MAT, "T", "D1", "S2", "PIPES");
PIPES = S.PIPES;

%% Variants
T0 = S.T;
T1 = T0(T0.dataset ~= OUT_DOM, :);
G = groupsummary(T1(T1.pipeline == PIPES(1), :), "session", "mean", "lf");
T1.lfb = T1.lfMean - mean(G.mean_lf);
z = @(x) (x - mean(x, "omitnan")) ./ std(x, "omitnan");
T1.lsa = z(log2(T1.sigma ./ T1.a));  T1.al = z(T1.alpha);
T1.datasetC = removecats(categorical(string(T1.datasetC)));  T1.sessionC = removecats(T1.sessionC);
V = {T0, T1};  VNAME = ["V0 as saved (Dhieb in)", "V1 #242 domain (Dhieb out)"];
D1 = cell(1, 2);  S2 = cell(1, 2);
for v = 1:2
    [D1{v}, S2{v}] = fitPooled_local(V{v}, PIPES);
end

%% Regression anchor
for pair = {{D1{1}, S.D1, "D1"}, {S2{1}, S.S2, "S2"}}
    [A, B, nm] = pair{1}{:};
    B = B(ismember(B.term, ["lfw" "lfb"]), :);  A = A(ismember(A.term, ["lfw" "lfb"]), :);
    key = @(X) X.pipeline + "|" + X.outcome + "|" + X.term;
    [ok, loc] = ismember(key(B), key(A));
    if ~all(ok), error("tempoDom:anchorKeys", "%s", nm + ": rows missing from V0"); end
    A = A(loc, :);
    if any(abs(A.b - B.b) > TOL.b) || any(abs(A.SE - B.SE) > TOL.SE) || any(A.df ~= B.df) || any(A.nTrials ~= B.nTrials)
        error("tempoDom:anchor", "%s", nm + ": V0 does not reproduce the saved table");
    end
end
fprintf("REGRESSION ANCHOR passed: V0 reproduces the saved D1 and S2 (lfw, lfb).\n");

%% Report
for v = 1:2
    fprintf("\n=== %s: %d trial x pipeline rows, %d sessions ===\n", VNAME(v), height(V{v}), numel(unique(V{v}.session)));
    fprintf("1. Within/between decomposition (per doubling of f0):\n");  disp(D1{v}(ismember(D1{v}.term, ["lfw" "lfb"]), :))
    fprintf("2a. Selection GLMM (logit per doubling):\n");  disp(S2{v})
    for p = PIPES
        w = D1{v}(D1{v}.pipeline == p & D1{v}.term == "lfw", :);
        bo = w.b(w.outcome == "beta_obs");  bg = w.b(w.outcome == "beta_gen*");
        fprintf("   %s within-session slope: beta_obs %+.4f, beta_gen* %+.4f; reduction %.1f%%\n", p, bo, bg, 100 * (1 - bg / bo));
    end
end

runDate = string(datetime("now", "Format", "yyyy-MM-dd"));
save(OUT_MAT, "D1", "S2", "VNAME", "OUT_DOM", "runDate");
fprintf("Saved: %s\n", OUT_MAT);

%% =========================================================================
function [D1, S2] = fitPooled_local(T, pipes)
% analyseTempoDecomposition_v002 L81-89 and L119-124, unchanged.
    terms = ["lfw" "lfb" "lf" "group" "drug" "lsa" "al"];
    D1 = table();  S2 = table();
    for p = pipes
        V = T(T.pipeline == p, :);
        mo = fitlme(V(isfinite(V.betaObs), :), "betaObs ~ 1 + lfw + lfb + lsa + al + datasetC + (1|sessionC)");
        mg = fitlme(V(isfinite(V.betaGen), :), "betaGen ~ 1 + lfw + lfb + lsa + al + datasetC + (1|sessionC)");
        D1 = [D1; tag_local(coef_local(mo, terms), p, "beta_obs", sum(isfinite(V.betaObs))); ...
                  tag_local(coef_local(mg, terms), p, "beta_gen*", sum(isfinite(V.betaGen)))]; %#ok<AGROW>
        ms = fitglme(V, "invertible ~ 1 + lfw + lfb + lsa + al + datasetC + (1|sessionC)", "Distribution", "Binomial", "FitMethod", "Laplace");
        S2 = [S2; tag_local(coef_local(ms, ["lfw" "lfb"]), p, "invertible", height(V))]; %#ok<AGROW>
    end
end

function t = coef_local(mdl, terms)
    c = mdl.Coefficients;  k = ismember(string(c.Name), terms);
    t = table(string(c.Name(k)), round(c.Estimate(k), 4), round(c.SE(k), 4), round(c.tStat(k), 2), c.DF(k), ...
        c.pValue(k), round(c.Lower(k), 4), round(c.Upper(k), 4), 'VariableNames', ["term" "b" "SE" "t" "df" "p" "CIlo" "CIhi"]);
end

function t = tag_local(c, p, outcome, n)
    t = [table(repmat(string(p), height(c), 1), repmat(string(outcome), height(c), 1), repmat(n, height(c), 1), ...
        'VariableNames', ["pipeline" "outcome" "nTrials"]), c];
end
