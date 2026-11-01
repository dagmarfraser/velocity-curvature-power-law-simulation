%% checkTemplateTailOmission_v001.m
% Sol's closure 2 (2026-10-27): does the headline pointwise inversion depend on the template routine's zero-residual tail?
% Under the PRODUCTION template (4 harmonics, 4 cycles per window; loopClosureCore_v004 L77) the window rule of
% templateTailStratify_v002's edge_local is replayed for every trial of every dataset (f0 and trace length as the core uses them);
% trials with zero stretch above DZ_TOL = 0.05 cycle are "tail trials". Then, with and without them:
%   (a) per dataset: invertible share, median beta_gen* (constellation), primary cccMed (checkPart3PrimaryMetric_v001 construction);
%   (b) the pooled within-session tempo slope and its reduction (analyseTempoDecomposition_v002 L81-89), Dhieb out as in the headline.
% Anchor: replayed analysed length equals the core's saved M for every trial (error otherwise).
% Reads src/loopClosureResults_<dataset>_all_shaped_xu_v015.mat and the dataset importers. Writes results/checkTemplateTailOmission_v001.mat.
% Run from src/.   Fraser, D.S. (2026)  v001

%% CONFIG
SRC = fileparts(mfilename("fullpath"));  ROOT = fileparts(SRC);
addpath(SRC);  addpath(genpath(fullfile(SRC, "functions")));  addpath(genpath(fullfile(SRC, "req")));
OUT = fullfile(ROOT, "results", "checkTemplateTailOmission_v001.mat");
NH = 4;  NCW = 4;  DZ_TOL = 0.05;  MDC = 0.03;  DEG = 0.005;
PL = ["BWFD-OLS" "SG-OLS" "BWFD-LMLS" "SG-LMLS" "BWFD-IRLS" "SG-IRLS"];  PIPES = ["SG-IRLS" "BWFD-OLS"];
DS = struct("name", {"Fraser" "Cook_CTRL" "Cook_ASD" "Hickman_PLAC" "Hickman_HALO" "Zarandi" "Dhieb"}, ...
    "fs", {240 133 133 133 133 100 100}, ...
    "imp", {@() importDB_fraser_v001(Verbose=false), @() importDB_cook_v002(Group="CTRL", Tasks=7, Verbose=false), ...
            @() importDB_cook_v002(Group="ASD", Tasks=7, Verbose=false), @() importDB_hickman_v003(Study=2, Group="PLAC", Verbose=false), ...
            @() importDB_hickman_v003(Study=2, Group="HALO", Verbose=false), @() importDB_zarandi_v001(Verbose=false), ...
            @() importDB_dhieb_v001(Verbose=false)});
PAIRED = ["Hickman_PLAC" "Hickman_HALO"];
[pP, pQ] = pairs_local(6);

%% Per dataset
A = table();  rows = {};
for q = 1:numel(DS)
    nm = DS(q).name;  fs = DS(q).fs;
    r = load(fullfile(SRC, "loopClosureResults_" + nm + "_all_shaped_xu_v015.mat"), "results").results;  N = numel(r);
    tr = DS(q).imp();  [ok, loc] = ismember(string({r.trialID}), string({tr.trialID}));
    if ~all(ok), error("tailOm:import", "%s", sprintf("%s: %d trials not imported", nm, sum(~ok))); end
    tr = tr(loc);  f0 = [r.f0]';  Nn = cellfun(@numel, {tr.x})';
    [z, Lc] = edge_local(Nn, f0, fs, NCW);
    Mrep = round(Lc .* fs ./ f0);  Msaved = [r.M]';
    if any(isfinite(Mrep) & abs(Mrep - Msaved) > 1), error("tailOm:anchor", "%s", nm + ": replayed analysed length differs from the core's M"); end
    isT = z > DZ_TOL;  isT(~isfinite(z)) = false;
    inv = arrayfun(@(t) isfinite(t.betaGenStarMed), r)';  gMed = [r.betaGenStarMed]';
    dbO = NaN(N, numel(pP));  dbP = dbO;  okc = false(N, 1);
    for t = 1:N
        bo = r(t).betaObs(:)';  sl = r(t).betaRecSlice;
        if numel(bo) ~= 6 || any(~isfinite(bo)) || isempty(sl) || size(sl, 2) ~= 6 || ~inv(t), continue; end
        bp = median(sl, 1, "omitnan");  if any(~isfinite(bp)), continue; end
        dbO(t, :) = bo(pP) - bo(pQ);  dbP(t, :) = bp(pP) - bp(pQ);  okc(t) = true;
    end
    valid = okc & max(abs(dbO), [], 2) >= DEG & ~all(abs(dbO) < MDC, 2);
    cA = ccc_local(dbP, dbO, valid);  cD = ccc_local(dbP, dbO, valid & ~isT);
    A = [A; table(string(nm), N, nnz(isT), 100 * mean(isT), median(z(isT)), cA, cD, cD - cA, median(gMed, "omitnan"), median(gMed(~isT), "omitnan"), ...
        100 * mean(inv), 100 * mean(inv(~isT)), 100 * mean(inv(isT)), nnz(valid), nnz(valid & ~isT), ...
        'VariableNames', ["dataset" "N" "nTail" "pctTail" "medZtail" "cccMedAll" "cccMedDrop" "dCCC" "medGenAll" "medGenDrop" "pctInvAll" "pctInvDrop" "pctInvTailOnly" "NvalidAll" "NvalidDrop"])]; %#ok<AGROW>
    person = nm;  if ismember(nm, PAIRED), person = "Hickman"; end
    for i = 1:N
        for p = PIPES
            k = find(PL == p);
            rows{end+1, 1} = struct("dataset", string(nm), "isTail", isT(i), "session", string(nm) + "_" + string(r(i).subjectID), ...
                "person", person + "_" + string(r(i).subjectID), "pipeline", p, "f0", r(i).f0, "sigma", r(i).sigmaMM, "a", r(i).a_mm, ...
                "alpha", r(i).alphaIRA, "betaObs", r(i).betaObs(k), "betaGen", r(i).betaGenStar(k)); %#ok<SAGROW>
        end
    end
end
A.dMedGen = A.medGenDrop - A.medGenAll;  A.dPctInv = A.pctInvDrop - A.pctInvAll;
disp(A)

%% (b) pooled tempo slope, Dhieb out and Zarandi out (validity domain), with and without tail trials
T0 = struct2table(vertcat(rows{:}));
T0 = T0(T0.dataset ~= "Zarandi" & T0.dataset ~= "Dhieb" & isfinite(T0.f0) & T0.f0 > 0 & isfinite(T0.sigma) & isfinite(T0.a) & T0.a > 0 & isfinite(T0.alpha), :);
VN = ["domain, all trials" "domain, tail trials out"];  V = {prep_local(T0, PIPES), prep_local(T0(~T0.isTail, :), PIPES)};
SUM = table();
for v = 1:2
    D1 = fitD1_local(V{v}, PIPES);
    for p = PIPES
        w = D1(D1.pipeline == p & D1.term == "lfw", :);
        bo = w.b(w.outcome == "beta_obs");  bg = w.b(w.outcome == "beta_gen*");
        SUM = [SUM; table(VN(v), p, height(V{v}) / 2, bo, bg, 100 * (1 - bg / bo), 'VariableNames', ["variant" "pipeline" "nTrials" "bObs" "bGen" "reductionPct"])]; %#ok<AGROW>
    end
end
disp(SUM)
fprintf("Largest |change|: cccMed %.4f, median beta_gen* %.4f, invertible share %.2f pp\n", max(abs(A.dCCC)), max(abs(A.dMedGen)), max(abs(A.dPctInv)));
runDate = string(datetime("now", "Format", "yyyy-MM-dd"));
save(OUT, "A", "SUM", "runDate");  fprintf("Saved: %s\n", OUT);

%% =========================================================================
function [z, L] = edge_local(Nn, f0, fs, nCW)
% templateTailStratify_v002 edge_local, unchanged: z = zero-residual stretch inside the analysed segment (cycles), L = analysed length (cycles).
    n = numel(Nn);  z = nan(n, 1);  L = z;
    for i = 1:n
        N = Nn(i);  if ~isfinite(f0(i)) || f0(i) <= 0, continue, end
        win0 = round(nCW / f0(i) * fs);  hop = round(win0 / 2);  win = max(min(win0, N), round(2 / f0(i) * fs));
        cl = max(round(win / 4), 1);  cs = cl;  ce = min(max(N - cl, cs + 100), N);
        if win > N, lastEnd = 0; else, lastEnd = (win - 1) + hop * floor((N - win) / hop) + 1; end
        z(i) = max(0, ce - max(lastEnd, cs - 1)) * f0(i) / fs;  L(i) = (ce - cs + 1) * f0(i) / fs;
    end
end

function m = ccc_local(dbP, dbO, valid)
    vi = find(valid);  c = NaN(numel(vi), 1);
    for i = 1:numel(vi), c(i) = linCCC_v001(dbP(vi(i), :)', dbO(vi(i), :)').ccc; end
    m = median(c, "omitnan");
end

function [pP, pQ] = pairs_local(n)
    pP = [];  pQ = [];
    for p = 2:n, for q = 1:p-1, pP(end+1) = p; pQ(end+1) = q; end, end %#ok<AGROW>
end

function T = prep_local(T, pipes)
% analyseTempoDecomposition_v002 L62-72: centring and z-scoring on the subset.
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
    D1 = table();
    for p = pipes
        V = T(T.pipeline == p, :);
        mo = fitlme(V(isfinite(V.betaObs), :), "betaObs ~ 1 + lfw + lfb + lsa + al + datasetC + (1|sessionC)");
        mg = fitlme(V(isfinite(V.betaGen), :), "betaGen ~ 1 + lfw + lfb + lsa + al + datasetC + (1|sessionC)");
        D1 = [D1; tag_local(mo.Coefficients, p, "beta_obs"); tag_local(mg.Coefficients, p, "beta_gen*")]; %#ok<AGROW>
    end
end

function t = tag_local(c, p, outcome)
    k = ismember(string(c.Name), ["lfw" "lfb"]);
    t = table(repmat(string(p), nnz(k), 1), repmat(string(outcome), nnz(k), 1), string(c.Name(k)), c.Estimate(k), ...
        'VariableNames', ["pipeline" "outcome" "term" "b"]);
end
