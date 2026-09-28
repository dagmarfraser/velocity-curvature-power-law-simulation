%% checkPracticalIdentifiability_v003.m
% Is the precision/identifiability dissociation real, or an artefact of how
% precision is measured? v002 (Finding #214 addendum) found 11/42 cells
% trial-level adequate on the RECOVERED scale, but 5 of the 11 are Dhieb
% (rise coverage 6-26%): recovered-scale SEM shrinks where beta_rec stops
% responding to beta_gen, so it rewards flatness. v002's conditioning proxy
% (genW/recW) used v012 CI fields affected by Finding #166, so it is not
% citable either.
%
% v003 reads the v015 corpora (N_REPS=200, same config and numerics as
% v013/v014 apart from run identity, Finding #204), which persist the
% per-trial forward map: betaRecCurveMed / betaRecCurveSD (nP x N_BETA).
% Verified 2026-09-19 on Hickman PLAC: interp1(betaGenVec, Med, betaGenStar)
% returns betaObs (median error 0), CurveSD equals std(betaRecSlice) at the
% slice node, and the slice node is the node nearest the cross-pipeline
% median betaGenStar.
%
% Conditioning, measured directly (no CI fields involved):
%   slope   = secant slope of Med over the grid interval containing the
%             pipeline's own betaGenStar. The inversion is linear
%             interpolation, so 1/slope is the EXACT local derivative of the
%             inverse actually applied, not an approximation to it.
%   ampl    = 1/|slope|
%   sdOp    = CurveSD interpolated at betaGenStar (recovered scale)
%   semGen  = sdOp * ampl  (delta method; generator scale; Inf where slope
%             is 0). This is the column that separates "precise" from
%             "flat": a flat cell has small sdOp but large ampl.
% Caveats: N_BETA=25 (step 0.03125), so slopes are interval averages; the
% delta method is first order and understates error where the map curves
% sharply within an interval. offCurve counts trials whose betaGenStar does
% not reproduce betaObs on Med (inversion off the plain curve).
%
% Replicate: the same verdict split and v002-definition medSEMrec are
% recomputed on v013/v014 (exactly one per dataset; Cook_ASD and Dhieb are
% v014) and compared cell by cell, so run-to-run drift is measured, not
% assumed irrelevant.
%
% Provenance: v003 2026-09-19, Session 113. Forks v002 (Session 112).

%% CONFIG
SRC_DIR   = fileparts(mfilename("fullpath"));
TAGS      = ["Zarandi","Cook_CTRL","Cook_ASD","Dhieb","Hickman_PLAC","Hickman_HALO","Fraser"];
COORD     = [3.18 4.77 100; 4.77 8.15 133; 5.06 7.84 133; 2.50 7.50 100; ...
             5.34 7.17 133; 5.42 7.42 133; 4.289 2.009 240];   % alpha, sigma, fs per TAGS (Fig 0 registry)
PRIMARY   = "v015";
REPLICATE = ["v013","v014"];                        % exactly one must exist per dataset
EXPECTED  = [8 2 32];                               % Finding #171 N=200 split; checked, not assumed
Z90       = 2 * 1.6449;                             % 90% interval -> SD
ZERO_TOL  = 1e-12;
OUT_BASE  = fullfile(SRC_DIR, "practicalIdentifiability_v003");
addpath(genpath(fullfile(SRC_DIR, "functions")));
THR       = semAdequacyThreshold_v001();            % single source, MDC/2.77

%% RUN
semFile = fullfile(SRC_DIR, "perCoordinateSEM_v2_001.mat");
if ~isfile(semFile), error("pid3:missing", "%s", "FAILED PATH: " + semFile); end
G = load(semFile, "coordTable"); G = G.coordTable;
rows = {};
for d = 1:numel(TAGS)
    tag = TAGS(d);
    L  = loadCorpus_local(SRC_DIR, tag, PRIMARY);
    Lr = loadCorpus_local(SRC_DIR, tag, resolveReplicate_local(SRC_DIR, tag, REPLICATE));
    labels = string(L.pipelineLabels); nP = numel(labels);
    if ~isequal(labels, string(Lr.pipelineLabels)) || numel(L.results) ~= numel(Lr.results)
        error("pid3:mismatch", "%s", tag + ": primary and replicate differ in pipelineLabels or trial count.");
    end
    g  = L.betaGenVec(:).';
    P  = perTrial_local(L.results,  g, nP, Z90, ZERO_TOL, true,  tag);
    Pr = perTrial_local(Lr.results, g, nP, Z90, ZERO_TOL, false, tag);
    semGrid = gridSEM_local(G, COORD(d,:), labels, tag);
    for p = 1:nP
        [verdict,  cov]  = verdict_local(P.st(:,p));
        [verdictR, covR] = verdict_local(Pr.st(:,p));
        inv  = P.inv(:,p);  invR = Pr.inv(:,p);
        sg   = P.semGen(inv,p);  a = P.ampl(inv,p);
        mis  = (P.st(:,p) == "rise" & P.slope(:,p) <= 0) | (P.st(:,p) == "desc" & P.slope(:,p) >= 0);
        rows(end+1,:) = {tag, labels(p), verdict, cov, nnz(inv), semGrid(p), ...
            median(P.semRec(inv,p), "omitnan"), median(P.sdOp(inv,p), "omitnan"), ...
            median(sg, "omitnan"), mean(sg < THR), median(a, "omitnan"), prctile_local(a(isfinite(a)), 90), ...
            nnz(mis), nnz(P.offCurve(inv,p)), median(P.semGenCI(inv,p), "omitnan"), nnz(P.zeroCI(inv,p)), ...
            verdictR, covR, median(Pr.semRec(invR,p), "omitnan")}; %#ok<SAGROW>
    end
end
T = cell2table(rows, "VariableNames", ["dataset","pipeline","verdict","riseCov","nInverted","semGrid", ...
    "medSEMrec","medSDop","medSEMgen","fracSEMgenAdequate","medAmpl","p90Ampl", ...
    "nSignMismatch","nOffCurve","medSEMgenCI","nZeroCI","verdictRep","riseCovRep","medSEMrecRep"]);

%% SELF-CHECK + REPORT
split = countVerdicts_local(T.verdict);
fprintf("Primary (%s) split %d/%d/%d (expected %d/%d/%d)\n", PRIMARY, split, EXPECTED);
if ~isequal(split, EXPECTED)
    warning("pid3:split", "%s", "Verdict split departs from Finding #171's N=200 split; resolve before interpreting.");
end
fprintf("Replicate split %d/%d/%d\n", countVerdicts_local(T.verdictRep));
flip = T.verdict ~= T.verdictRep;
fprintf("Verdict flips primary vs replicate: %d/42\n", nnz(flip));
if any(flip), disp(T(flip, ["dataset","pipeline","verdict","riseCov","verdictRep","riseCovRep"])); end
fprintf("medSEMrec max |primary - replicate| = %.4g; adequate %d vs %d\n", ...
    max(abs(T.medSEMrec - T.medSEMrecRep)), nnz(T.medSEMrec < THR), nnz(T.medSEMrecRep < THR));
fprintf("Trials off the Med curve: %d; sign-mismatched slopes: %d; zero-width CIs: %d\n\n", ...
    sum(T.nOffCurve), sum(T.nSignMismatch), sum(T.nZeroCI));

fprintf("Adequacy threshold (MDC/2.77) = %.5f\n", THR);
crossTab_local(T, T.semGrid   < THR, "grid SEM (registered basis)");
crossTab_local(T, T.medSEMrec < THR, "trial, recovered scale (v002 definition)");
crossTab_local(T, T.medSEMgen < THR, "trial, generator scale (delta method)");
crossTab_local(T, T.medSEMgenCI < THR, "trial, generator scale (CI width, cross-check)");

disp(T(:, ["dataset","pipeline","verdict","semGrid","medSEMrec","medSEMgen","medAmpl","p90Ampl","medSEMgenCI"]));
save(OUT_BASE + ".mat", "T", "THR", "PRIMARY", "COORD");
writetable(T, OUT_BASE + ".txt", "Delimiter", "\t");
fprintf("Saved %s.mat / .txt\n", OUT_BASE);

%% ---------------------------------------------------------------------------
function L = loadCorpus_local(srcDir, tag, ver)
f = fullfile(srcDir, sprintf("loopClosureResults_%s_all_shaped_xu_%s.mat", tag, ver));
if ~isfile(f), error("pid3:missing", "%s", "FAILED PATH: " + f); end
L = load(f, "results", "pipelineLabels", "betaGenVec", "config");
if L.config.N_REPS ~= 200
    error("pid3:nreps", "%s", sprintf("%s %s: N_REPS = %d, expected 200.", tag, ver, L.config.N_REPS));
end
end

function ver = resolveReplicate_local(srcDir, tag, candidates)
hit = arrayfun(@(v) isfile(fullfile(srcDir, sprintf("loopClosureResults_%s_all_shaped_xu_%s.mat", tag, v))), candidates);
if nnz(hit) ~= 1
    error("pid3:replicate", "%s", sprintf("%s: %d of {%s} present, expected exactly 1.", tag, nnz(hit), strjoin(candidates, ",")));
end
ver = candidates(hit);
end

function P = perTrial_local(R, g, nP, z90, zeroTol, needCurve, tag)
nT = numel(R); hasCurve = isfield(R, "betaRecCurveMed");
if needCurve && ~hasCurve, error("pid3:curve", "%s", tag + ": primary corpus lacks betaRecCurveMed."); end
P.st = strings(nT, nP);
[P.semRec, P.semGenCI, P.slope, P.sdOp] = deal(NaN(nT, nP));
[P.zeroCI, P.offCurve] = deal(false(nT, nP));
for t = 1:nT
    r = R(t); s = string(r.invertStatus);
    if numel(s) ~= nP, error("pid3:status", "%s", sprintf("%s trial %d: invertStatus not per-pipeline.", tag, t)); end
    P.st(t,:) = s(:).';
    bgs = r.betaGenStar(:).';
    if ~any(isfinite(bgs)), continue, end
    P.semRec(t,:) = (prctile(r.betaRecSlice, 95, 1) - prctile(r.betaRecSlice, 5, 1)) / z90;
    w = (r.ciHi(:) - r.ciLo(:)).';
    P.zeroCI(t,:) = isfinite(w) & w < zeroTol;  w(P.zeroCI(t,:)) = NaN;   % degenerate != precise
    P.semGenCI(t,:) = w / z90;
    if ~hasCurve, continue, end
    M = r.betaRecCurveMed; S = r.betaRecCurveSD;
    if ~isequal(size(M), [nP numel(g)]) || ~isequal(size(S), [nP numel(g)])
        error("pid3:curveSize", "%s", sprintf("%s trial %d: curve is %s, expected [%d %d].", tag, t, mat2str(size(M)), nP, numel(g)));
    end
    for p = find(isfinite(bgs))
        k = find(g <= bgs(p), 1, "last");
        if isempty(k), k = 1; end
        k = min(k, numel(g) - 1);
        P.slope(t,p)    = (M(p,k+1) - M(p,k)) / (g(k+1) - g(k));
        P.sdOp(t,p)     = interp1(g, S(p,:), bgs(p));
        P.offCurve(t,p) = abs(interp1(g, M(p,:), bgs(p)) - r.betaObs(p)) > 1e-6;
    end
end
P.inv    = P.st == "rise" | P.st == "desc";
P.ampl   = 1 ./ abs(P.slope);                 % Inf where slope == 0: locally non-identifiable
P.semGen = P.sdOp .* P.ampl;
end

function semGrid = gridSEM_local(G, coord, labels, tag)
% Canonical grid SEM: mean over ALL betaGen and ALL VGF at the snapped
% (alpha, sigma, fs); same lookup as plotPaperRoadmapSchematic_v001.
uA = unique(G.alpha); uS = unique(G.sigma); uF = unique(G.fs);
[~,ai] = min(abs(uA - coord(1))); [~,si] = min(abs(uS - coord(2))); [~,fi] = min(abs(uF - coord(3)));
atCoord = G.alpha == uA(ai) & G.sigma == uS(si) & G.fs == uF(fi);
semGrid = NaN(1, numel(labels));
for p = 1:numel(labels)
    m = atCoord & G.pipeline == labels(p);
    if ~any(m), error("pid3:grid", "%s", tag + "/" + labels(p) + ": no grid SEM rows at snapped coordinate."); end
    semGrid(p) = mean(G.sem(m), "omitnan");
end
end

function [v, cov] = verdict_local(st)
valid = st ~= "no_beta_obs";
cov = sum(st(valid) == "rise") / max(nnz(valid), 1);
if cov >= 0.95, v = "PASS"; elseif cov >= 0.90, v = "CONDITIONAL"; else, v = "FAIL"; end
end

function n = countVerdicts_local(v)
n = [nnz(v == "PASS"), nnz(v == "CONDITIONAL"), nnz(v == "FAIL")];
end

function crossTab_local(T, adequate, label)
fprintf("\nAdequacy basis: %s  (%d/42 adequate)\n", label, nnz(adequate));
fprintf("%-12s %9s %11s\n", "verdict", "adequate", "inadequate");
for v = ["PASS","CONDITIONAL","FAIL"]
    m = T.verdict == v;
    fprintf("%-12s %9d %11d", v, nnz(m & adequate), nnz(m & ~adequate));
    ds = unique(T.dataset(m & adequate));
    if ~isempty(ds), fprintf("   adequate: %s", strjoin(ds, ", ")); end
    fprintf("\n");
end
end

function q = prctile_local(x, p)
if isempty(x), q = NaN; else, q = prctile(x, p); end
end
