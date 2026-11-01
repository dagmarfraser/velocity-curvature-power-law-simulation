function out = analyseTempoSlopeRecovery_v003(source, opts)
%ANALYSETEMPOSLOPERECOVERY_V002  = v001 plus opts.GateOnly (spec SPEC_T2 10.18): returns the synthetic and empirical-matched gate pass rates for EVERY pipeline
%   (v001 reports PIPES(1) only) using the same trial table and filter, saves analyseTempoGateAll_<arm>.mat, and fits nothing. Without GateOnly it is v001.
%   Tempo-slope decomposition on empirical or synthetic (T2/T2b) trials.
%   out = analyseTempoSlopeRecovery_v001("empirical")   ACCEPTANCE TEST: must reproduce the published
%                                                       within-session slopes (checkTempoDecompositionDomain_v001, V1).
%   out = analyseTempoSlopeRecovery_v001("A")           any simulated arm (A,B1,B2,E,Sm2,...); needs results/tempoSim.
%   Spec: docs/SPEC_T2_TempoSlopeRecovery_v001.md sections 3-4.
%
%   Model, per pipeline (SG-IRLS primary, BWFD-OLS legacy), exactly analyseTempoDecomposition_v002 L81-89:
%     outcome ~ 1 + lfw + lfb + lsa + al + datasetC + (1|sessionC)   outcome = betaObs | betaGen* (rise-gated)
%   Variants (synthetic only; true f0 is known there):
%     V1 as published (ESTIMATED f0, estimated noise covariates)
%     V2 TRUE f0 in place of estimated f0
%     V3 no noise covariates (drops lsa, al)
%   Synthetic arms also report an EMPIRICAL-MATCHED baseline (same trial IDs, real data) and the synthetic
%   vs empirical gate pass rate per dataset (a synthetic rate far above the empirical one = arm too easy).
%   Datasets: the #242 validity domain (Dhieb and Zarandi out), same as the published slope.
arguments
    source (1,1) string
    opts.Tol (1,1) double = 1e-4
    opts.GateOnly (1,1) logical = false
    opts.KeepTable (1,1) logical = false   % keep the synthetic and matched trial tables in out.Tv, out.Te (diagnostics, spec 10.22)
    opts.Tag (1,1) string = ""   % appended to saved file names so a restricted-dataset run never overwrites the standard arm file
    opts.Inversion (1,1) string {mustBeMember(opts.Inversion,["gated","ungated"])} = "gated"   % ungated = gate-off diagnostic (arm G)
    opts.Datasets (1,:) string = ["Fraser","Cook_CTRL","Cook_ASD","Hickman_PLAC","Hickman_HALO"]   % override only for smoke tests
end
root = fileparts(fileparts(mfilename('fullpath')));
srcDir = fullfile(root,'src'); simDir = fullfile(root,'results','tempoSim');
DS  = opts.Datasets;
PL  = ["BWFD-OLS","SG-OLS","BWFD-LMLS","SG-LMLS","BWFD-IRLS","SG-IRLS"];
PIPES = ["SG-IRLS","BWFD-OLS"];
out = struct('source',source);

if source == "empirical"
    T = buildTable_local(DS, "", srcDir, simDir, PL, PIPES, "", "gated");
    R = [fitAll_local(T, PIPES, "V1", "empirical"); fitAll_local(T, PIPES, "V3", "empirical")];   % V3 = no noise covariates (adjustment-set sensitivity, spec 10.10)
    out.slopes = R;
    A = load(fullfile(root,'results','checkTempoDecompositionDomain_v001.mat'),'D1');
    D = A.D1{2};
    bad = false;
    for p = PIPES
        for oc = ["beta_obs","beta_gen*"]
            ref = D.b(D.pipeline==p & D.outcome==oc & D.term=="lfw");
            got = R.b(R.pipeline==p & R.outcome==strrep(strrep(oc,"beta_obs","betaObs"),"beta_gen*","betaGen") & R.variant=="V1");
            fprintf('acceptance %-8s %-9s published %+.4f   here %+.4f\n', p, oc, ref, round(got,4));
            if abs(round(got,4) - ref) > opts.Tol, bad = true; end
        end
    end
    if bad, error('analyseTempoSlope:acceptance','Empirical slopes do not reproduce the published V1 table; do not trust synthetic arms.'); end
    fprintf('ACCEPTANCE passed: published within-session slopes reproduced.\n');
    Wd = table(); for p = PIPES, for oc = ["betaObs","betaGen"]   % adjustment-set sensitivity: with (V1) vs without (V3) the noise covariates lsa, al
        a = R(R.variant=="V1" & R.pipeline==p & R.outcome==oc,:); b = R(R.variant=="V3" & R.pipeline==p & R.outcome==oc,:);
        Wd = [Wd; table(p, oc, a.n, a.b, a.SE, b.n, b.b, b.SE, b.b-a.b, 'VariableNames', ["pipeline","outcome","nV1","bV1","seV1","nV3","bV3","seV3","V3minusV1"])]; end, end %#ok<AGROW>
    out.adjustment = Wd; fprintf('Empirical within-session tempo slope: published (V1) vs WITHOUT noise covariates lsa, al (V3); spec 10.10.\n'); disp(Wd)
    save(fullfile(simDir, 'analyseTempoSlopeRecovery_empirical.mat'), 'out');
    return
end

% ---------------- synthetic arm ----------------
Tv = buildTable_local(DS, source, srcDir, simDir, PL, PIPES, "synthetic", opts.Inversion);
Te = buildTable_local(DS, source, srcDir, simDir, PL, PIPES, "matched", opts.Inversion);
Pm = load(fullfile(simDir, sprintf('tempoSim_%s_%s.mat', source, strrep(DS(1)," ","_"))), 'P'); trueSlope = Pm.P.slope;
if opts.GateOnly
    out.gateAll = gateAll_local(Tv, Te, DS, PIPES);  disp(out.gateAll)
    save(fullfile(simDir, sprintf('analyseTempoGateAll_%s.mat', source)), 'out');  return
end
R = [fitAll_local(Tv, PIPES, "V1", "synthetic"); fitAll_local(Tv, PIPES, "V2", "synthetic"); ...
     fitAll_local(Tv, PIPES, "V3", "synthetic"); fitAll_local(Te, PIPES, "V1", "empirical-matched"); ...
     fitAll_local(Te(ismember(Te.dataset + "|" + Te.trialSrc, unique(Tv.dataset + "|" + Tv.trialSrc)),:), PIPES, "V1", "empirical-sametrials")];   % matched trials that survive in the synthetic table (spec 10.25 R19)
out.slopes = R; out.trueSlope = trueSlope;
if opts.KeepTable, out.Tv = Tv; out.Te = Te; end
% gate pass rate: synthetic vs empirical on the same source trials
G = table();
for d = DS
    s = Tv(Tv.dataset==d & Tv.pipeline==PIPES(1),:); e = Te(Te.dataset==d & Te.pipeline==PIPES(1),:);
    G = [G; table(d, height(s), mean(s.invertible), mean(e.invertible), 'VariableNames', ["dataset","nTrials","passSynthetic","passEmpiricalSame"])]; %#ok<AGROW>
end
out.gate = G;
fprintf('\n=== arm %s [%s inversion] (true within-session tempo slope on beta_gen: %+.4f) ===\n', source, upper(opts.Inversion), trueSlope);
if opts.Inversion == "ungated"   % disclosure: how often the ungated value is forced to a grid edge (observed beta outside the whole forward map)
    fprintf('UNGATED diagnostic: share of synthetic trials with beta_obs outside the forward map (value pinned to nearest edge) = %.3f\n', mean(Tv.atEdge(Tv.pipeline==PIPES(1))));
end
disp(R(R.pipeline==PIPES(1),:)); fprintf('Gate pass rate, SG-IRLS (synthetic vs empirical on the same source trials):\n'); disp(G)
out.verdict = verdict_local(source, R, trueSlope, PIPES(1));
fprintf('%s\n', out.verdict);
if startsWith(source,"E")   % strict estimator-coupling bound (replaces the slope test, which is underpowered when true tempo is constant)
    out.coupling = coupling_local(Tv, Te, PIPES(1)); disp(out.coupling); out.verdict = out.verdict + newline + out.coupling.text;
    fprintf('%s\n', out.coupling.text);
end
tag = ''; if opts.Inversion == "ungated", tag = '_ungated'; end
tag = [tag char(opts.Tag)];
save(fullfile(simDir, sprintf('analyseTempoSlopeRecovery_%s%s.mat', source, tag)), 'out');
end

% =========================================================================
function T = buildTable_local(DS, arm, srcDir, simDir, PL, PIPES, kind, inv)
% one row per trial x pipeline, as analyseTempoDecomposition_v002 L52-77
rows = {};
for d = DS
    if kind == "" || kind == "matched"
        f = fullfile(srcDir, "loopClosureResults_" + d + "_all_shaped_xu_v015.mat");
    else
        f = fullfile(simDir, "loopClosureResults_SYN-" + arm + "-" + d + "_all_shaped_xu_syn001.mat");
    end
    if ~isfile(f), error('analyseTempoSlope:input','Missing %s',f); end
    R = load(f,'results').results; gv = [];
    if inv == "ungated", gv = load(f,'betaGenVec').betaGenVec; end
    ids = string({R.trialID})';
    if kind == "synthetic"
        tr = load(fullfile(simDir, "tempoSim_" + arm + "_" + d + ".mat"), 'truth').truth;
        [ok, loc] = ismember(ids, tr.trialID); if ~all(ok), error('analyseTempoSlope:truth','%s: trials without truth',d); end
        f0True = tr.f0True(loc);
        sf = fullfile(simDir, "tempoSim_" + arm + "_" + d + ".mat");
        w = whos('-file', sf);
        if any(strcmp({w.name}, 'covariatesEstimated'))
            bio = load(sf, 'bioResults').bioResults;
        else   % older corpus: stored covariates are copies of the source trial's, so re-estimate from the record
            warning('analyseTempoSlope:recomputeCovariates', '%s: re-estimating noise covariates from the stored synthetic trials.', sf);
            S = load(sf, 'trials', 'FS'); trs = S.trials; for q = 1:numel(trs), trs(q).fs = S.FS; end
            bio = characteriseBiologicalNoise_v004(trs, Method="irasa", Verbose=false);
        end
        [okb, lb] = ismember(ids, string(bio.trialID)); if ~all(okb), error('analyseTempoSlope:bio', '%s: trials without noise covariates', d); end
        alphaE = bio.ira_alphaMean(lb); sigmaE = bio.sigmaMean(lb) * load(sf, 'sigToMM').sigToMM;
    elseif kind == "matched"
        srcIDs = extractAfter(string(load(fullfile(simDir,"tempoSim_"+arm+"_"+d+".mat"),'truth').truth.trialID), "SYN-"+arm+"-");
        keep = ismember(ids, srcIDs); R = R(keep); f0True = [R.f0]';
    else
        f0True = [R.f0]';
    end
    if kind ~= "synthetic", alphaE = [R.alphaIRA]'; sigmaE = [R.sigmaMM]'; end   % empirical: stored covariates
    for i = 1:numel(R)
        r = R(i);
        for p = PIPES
            k = find(PL == p); edge = false; bG = r.betaGenStar(k);
            if inv == "ungated", [bG, edge] = ungatedInvert_local(r.betaObs(k), r.betaRecCurveMed(k,:), gv); end
            tsrc = string(r.trialID); if kind == "synthetic", tsrc = extractAfter(tsrc, "SYN-" + arm + "-"); end   % source trial ID (spec 10.25)
            rows{end+1,1} = struct("dataset",d,"session",d+"_"+string(r.subjectID),"pipeline",p,"trialSrc",tsrc,"f0",r.f0, ...
                "f0True",f0True(i),"sigma",sigmaE(i),"a",r.a_mm,"alpha",alphaE(i),"betaObs",r.betaObs(k), ...
                "betaGen",bG,"atEdge",edge,"invertible",r.invertStatus(k)=="rise"); %#ok<AGROW>
        end
    end
end
T = struct2table(vertcat(rows{:}));
T = T(isfinite(T.f0) & T.f0>0 & isfinite(T.sigma) & isfinite(T.a) & T.a>0 & isfinite(T.alpha), :);
end

function G = gateAll_local(Tv, Te, DS, PIPES)
% pass rate (invertStatus == "rise") per dataset x pipeline, synthetic vs empirical on the same source trials
G = table();
for d = DS, for p = PIPES
    s = Tv(Tv.dataset==d & Tv.pipeline==p,:); e = Te(Te.dataset==d & Te.pipeline==p,:);
    G = [G; table(d, p, height(s), height(e), mean(s.invertible), mean(e.invertible), 'VariableNames', ["dataset","pipeline","nSyn","nEmp","passSyn","passEmp"])]; %#ok<AGROW>
end, end
end

function R = fitAll_local(T, PIPES, variant, label)
% variant-specific tempo and covariates, then the published model
z = @(x) (x - mean(x,"omitnan")) ./ std(x,"omitnan");
if variant == "V2", T.lf = log2(T.f0True); else, T.lf = log2(T.f0); end
first = T(T.pipeline==PIPES(1),:);
G = groupsummary(first,"session","mean","lf"); [~,loc] = ismember(T.session,G.session);
T.lfw = T.lf - G.mean_lf(loc); T.lfb = G.mean_lf(loc) - mean(G.mean_lf);
T.lsa = z(log2(T.sigma./T.a)); T.al = z(T.alpha);
T.datasetC = categorical(T.dataset); T.sessionC = categorical(T.session);
rhs = " ~ 1 + lfw + lfb + lsa + al + datasetC + (1|sessionC)";
oneDS = numel(unique(T.dataset)) == 1;   % a single dataset: datasetC is constant and is dropped (spec 10.20)
if oneDS, rhs = " ~ 1 + lfw + lfb + lsa + al + (1|sessionC)"; end
if variant == "V3", rhs = " ~ 1 + lfw + lfb + datasetC + (1|sessionC)"; if oneDS, rhs = " ~ 1 + lfw + lfb + (1|sessionC)"; end, end
R = table();
for p = PIPES
    V = T(T.pipeline==p,:);
    for oc = ["betaObs","betaGen"]
        W = V(isfinite(V.(oc)),:);
        if max(abs(W.lfw)) < 1e-9   % no within-session tempo variation: lfw column is zero, model unidentifiable
            warning('analyseTempoSlope:noLfw','%s %s: lfw has no within-session variation (true tempo constant by design); fit skipped, reported as NaN.', label, variant);
            R = [R; table(string(label), string(variant), p, oc, height(W), NaN, NaN, NaN, NaN, NaN, NaN, 'VariableNames',["source","variant","pipeline","outcome","n","b","SE","CIlo","CIhi","bB","SEB"])]; continue
        end
        m = fitlme(W, oc + rhs); c = m.Coefficients; k = string(c.Name)=="lfw";
        R = [R; table(string(label), string(variant), p, oc, height(W), c.Estimate(k), c.SE(k), c.Lower(k), c.Upper(k), c.Estimate(string(c.Name)=="lfb"), c.SE(string(c.Name)=="lfb"), ...
            'VariableNames',["source","variant","pipeline","outcome","n","b","SE","CIlo","CIhi","bB","SEB"])]; %#ok<AGROW>
    end
end
end

function s = verdict_local(arm, R, trueSlope, p)
g = @(src,var,oc) R(R.source==src & R.variant==var & R.pipeline==p & R.outcome==oc,:);
v1g = g("synthetic","V1","betaGen"); v1o = g("synthetic","V1","betaObs");
v2g = g("synthetic","V2","betaGen"); v3g = g("synthetic","V3","betaGen"); em = g("empirical-matched","V1","betaObs");
s = string(sprintf('VERDICT arm %s (spec section 4; thresholds fixed before running):\n', arm));
s = s + sprintf('  beta_obs slope: synthetic %+.4f vs empirical-matched %+.4f\n', v1o.b, em.b);
switch armClass_local(arm, trueSlope)
    case "null"
        s = s + sprintf('  beta_gen* slope %+.4f (CI %+.4f..%+.4f): %s (null needs |b|<=0.01)\n', v1g.b, v1g.CIlo, v1g.CIhi, nullv_local(v1g));
        s = s + sprintf('  beta_obs slope within 0.02 of +0.084: %s\n', passfail_local(abs(v1o.b-0.084)<=0.02));
    case "recovery"
        r = v1g.b/trueSlope; excl = v1g.CIlo>0 || v1g.CIhi<0;
        s = s + sprintf('  recovered/true = %.2f (CI %.2f..%.2f), CI excludes 0: %d -> %s\n', r, v1g.CIlo/trueSlope, v1g.CIhi/trueSlope, excl, ...
            ternary_local(r>=0.8 && excl, 'RECOVERED (no attenuation)', ternary_local(r<0.5, 'NOT IDENTIFIED (attenuated >50%)', 'PARTIAL')));
    case "E"
        s = s + sprintf('  slope on estimated tempo %+.4f: %s (estimator coupling needs |b|<=0.01 to be negligible)\n', v1g.b, nullv_local(v1g));
end
s = s + sprintf('  published spec vs true f0 (V1-V2) %+.4f, vs no noise covariates (V1-V3) %+.4f: %s\n', v1g.b-v2g.b, v1g.b-v3g.b, ...
    ternary_local(max(abs([v1g.b-v3g.b, v1g.b-v2g.b]),[],'omitnan')>0.01, 'BIAS of published specification (>0.01)', 'no material bias'));
end

% ESTIMATOR-COUPLING BOUND (arm E: true f0 is constant within session, so the within-session spread of the ESTIMATED f0 is pure estimator error e).
% If e covaries with the outcome y, the published within-session slope is biased by about cov_w(e,y)/var_w(lfw_empirical).
% Session-bootstrap CI (resampling whole sessions). Reported for beta_obs and beta_gen*.
function C = coupling_local(Tv, Te, p)
S = Tv(Tv.pipeline==p,:); E = Te(Te.pipeline==p,:);
vEmp = var(wdev_local(log2(E.f0), E.session));
[~,~,gs] = unique(S.session); nS = max(gs); rows = accumarray(gs, (1:height(S))', [], @(v) {v});
rng(1); nb = 500; res = table();
for oc = ["betaObs","betaGen"]
    ok = isfinite(S.(oc));
    est = covW_local(log2(S.f0(ok)), S.(oc)(ok), S.session(ok)) / vEmp;
    bs = NaN(nb,1);
    for b = 1:nb
        pick = randi(nS, nS, 1); idx = []; grp = [];
        for q = 1:nS, v = rows{pick(q)}; v = v(ok(v)); idx = [idx; v]; grp = [grp; q*ones(numel(v),1)]; end %#ok<AGROW>
        bs(b) = covW_local(log2(S.f0(idx)), S.(oc)(idx), grp) / vEmp;
    end
    res = [res; table(oc, est, prctile(bs,2.5), prctile(bs,97.5), sqrt(var(wdev_local(log2(S.f0), S.session))/vEmp), 'VariableNames', ["outcome","bound","lo","hi","sdRatio"])]; %#ok<AGROW>
end
g = res(res.outcome=="betaGen",:);
txt = string(sprintf('  ESTIMATOR COUPLING (strict): possible bias of the published within-session beta_gen* slope = %+.4f (CI %+.4f..%+.4f); estimator-error SD / empirical within-session SD = %.2f. Negligible if bound and CI lie within +-0.01.', g.bound, g.lo, g.hi, g.sdRatio));
C = struct('table', res, 'text', txt);
end

function d = wdev_local(x, g)
[~,~,k] = unique(g); m = accumarray(k, x, [], @mean); d = x - m(k);
end

function c = covW_local(x, y, g)
xd = wdev_local(x, g); yd = wdev_local(y, g); c = sum(xd .* yd) / (numel(x) - 1);
end

% UNGATED beta_gen* (gate-off diagnostic): nearest-curve inversion of beta_obs through the trial's own stored forward map, no monotonicity gate.
function [bg, atEdge] = ungatedInvert_local(bObs, curve, g)
bg = NaN; atEdge = false; ok = isfinite(curve(:)) & isfinite(g(:)); c = curve(ok); x = g(ok); c = c(:)'; x = x(:)';
if ~isfinite(bObs) || numel(c) < 2, return, end
atEdge = bObs < min(c) || bObs > max(c);
[~,i] = min(abs(c - bObs)); bg = x(i);
for j = [i-1 i+1]   % refine linearly toward a neighbour that brackets bObs
    if j >= 1 && j <= numel(c) && (c(i)-bObs)*(c(j)-bObs) < 0 && c(j) ~= c(i), bg = x(i) + (bObs-c(i))*(x(j)-x(i))/(c(j)-c(i)); break, end
end
end

% null-truth arms (incl. reseeds As2.. and generator controls AXb/AXu) vs recovery arms (true slope > 0) vs estimator-coupling arm E
function c = armClass_local(arm, trueSlope)
if startsWith(arm,"E"), c = "E"; elseif trueSlope > 0, c = "recovery"; else, c = "null"; end
end

function s = passfail_local(tf), s = ternary_local(tf,'PASS','FAIL'); end
function v = ternary_local(c,a,b), if c, v = a; else, v = b; end, end

function s = nullv_local(r)
% a null/near-zero claim needs both a small estimate and a CI tight enough to mean something
% 16 Oct fix (labelling only, thresholds unchanged): a CI lying wholly outside +-0.01 is a FAIL however wide it is
if r.CIlo > 0.01 || r.CIhi < -0.01, s = 'FAIL (CI excludes the +-0.01 null band)';
elseif (r.CIhi - r.CIlo)/2 > 0.02, s = 'INCONCLUSIVE (CI half-width > 0.02)';
elseif abs(r.b) <= 0.01, s = 'PASS';
else, s = 'FAIL'; end
end
