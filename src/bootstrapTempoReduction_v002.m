function out = bootstrapTempoReduction_v002(opts)
%BOOTSTRAPTEMPOREDUCTION_V001  Person-level bootstrap of the within-session tempo slopes and of the inversion's reduction of them.
%   out = bootstrapTempoReduction_v002()                    HPC: B = 2000, 72 workers (parpool 'local')
%   out = bootstrapTempoReduction_v002(B=2, NWorkers=0)     smoke test, serial
%   Spec and pre-stated reading rule: docs/SPEC_T2_TempoSlopeRecovery_v001.md section 10.13 (written before the run).
%   Data and model exactly as checkTempoDecompositionDomain_v001 (V1: Dhieb out; lfb re-centred; lsa, al re-z-scored):
%     V1  outcome ~ 1 + lfw + lfb + lsa + al + datasetC + (1|sessionC)     as published
%     V3  outcome ~ 1 + lfw + lfb + datasetC + (1|sessionC)                no noise covariates (spec 10.10)
%   Resampling: persons within dataset, with replacement; a person's sessions and trials travel together (in this table
%   each person has one session). lfb, lsa and al are recomputed on every resample. Replicate b uses rng(Seed + b), so the
%   result does not depend on the number of workers. A replicate that fails to fit is counted, never replaced.
%   Estimands per pipeline: bObsV1 bGenV1 bObsV3 bGenV3 (within-session slope per doubling of f0), redV1 = 1 - bGenV1/bObsV1,
%   redV3, dGen = bGenV3 - bGenV1, dObs = bObsV3 - bObsV1. CIs are percentile 95%.
%   Writes results/tempoSim/bootstrapTempoReduction_v002.mat and .txt.   Fraser, D.S. (2026)  v001
arguments
    opts.B (1,1) double {mustBePositive, mustBeInteger} = 2000
    opts.Seed (1,1) double = 20261017
    opts.NWorkers (1,1) double {mustBeNonnegative, mustBeInteger} = 72
    opts.MaxFailFrac (1,1) double = 0.01
    opts.InMat (1,1) string = ""   % v002: alternative trial table (e.g. end-anchored v019); "" = results/tempoDecomposition_v002.mat
    opts.Tag (1,1) string = ""     % v002: suffix on output names
end
EST = ["bObsV1","bGenV1","bObsV3","bGenV3","redV1","redV3","dGen","dObs"];
root = fileparts(fileparts(mfilename('fullpath')));  simDir = fullfile(root,'results','tempoSim');
inMat = fullfile(root,'results','tempoDecomposition_v002.mat');  if opts.InMat ~= "", inMat = char(opts.InMat); end
if ~isfile(inMat), error('bootTempo:input','%s', "FAILED PATH: " + inMat); end
if ~isfolder(simDir), mkdir(simDir); end
diary(fullfile(simDir,char("bootstrapTempoReduction_v002" + opts.Tag + ".txt"))); cleanup = onCleanup(@() diary('off')); %#ok<NASGU>
fprintf('bootstrapTempoReduction_v002 | %s | B=%d seed=%d workers=%d | root %s\n', string(datetime('now')), opts.B, opts.Seed, opts.NWorkers, root);
S = load(inMat,'T','PIPES');  PIPES = string(S.PIPES);  T0 = S.T(S.T.dataset ~= "Dhieb", :);   % #242 validity domain
T0.dataset = string(T0.dataset);  T0.session = string(T0.session);  T0.person = string(T0.person);  T0.pipeline = string(T0.pipeline);

%% point estimates and regression anchors
est = zeros(numel(PIPES), numel(EST));  Tp = prep_local(T0, PIPES);
for k = 1:numel(PIPES), est(k,:) = fitEst_local(Tp, PIPES(k)); end
if opts.InMat == "", anchor_local(est, PIPES, EST, root); else, warning('bootTempo:noAnchor', 'REGRESSION ANCHOR SKIPPED: alternative table %s; V1 not checked against the published slopes (expected to differ).', opts.InMat); end

%% bootstrap
ds = unique(T0.dataset);  rowsByPerson = cell(numel(ds),1);
for i = 1:numel(ds)
    ri = find(T0.dataset == ds(i));  [~,~,g] = unique(T0.person(ri));  rowsByPerson{i} = accumarray(g, ri, [], @(v) {v});
end
B = opts.B;  M = opts.NWorkers;  seed0 = opts.Seed;  E = nan(B, numel(PIPES), numel(EST));  failed = false(B,1);  msg = strings(B,1);
if M > 0, ensurePool_local(M); end
t0 = tic;
parfor (b = 1:B, M)
    try
        e = oneRep_local(T0, rowsByPerson, seed0 + b, PIPES);
    catch ME
        e = nan(numel(PIPES), numel(EST));  failed(b) = true;  msg(b) = string(ME.message);
    end
    E(b,:,:) = reshape(e, 1, numel(PIPES), numel(EST));
end
fprintf('bootstrap done: %d replicates in %.1f min, %d failed\n', B, toc(t0)/60, nnz(failed));
if mean(failed) > opts.MaxFailFrac
    error('bootTempo:fail','%d of %d replicates failed (limit %.0f%%); first message: %s', nnz(failed), B, 100*opts.MaxFailFrac, msg(find(failed,1)));
end

%% report
R = table();
for k = 1:numel(PIPES)
    for j = 1:numel(EST)
        x = E(~failed,k,j);  ci = prctile(x,[2.5 97.5]);
        R = [R; table(PIPES(k), EST(j), est(k,j), std(x,'omitnan'), ci(1), ci(2), 'VariableNames',["pipeline","estimand","estimate","bootSE","ciLo","ciHi"])]; %#ok<AGROW>
    end
end
disp(R)
g = @(n) R(R.pipeline==PIPES(1) & R.estimand==n,:);
r1 = g("redV1"); d = g("dGen"); b1 = g("bGenV1");
fprintf('\nREADING (SG-IRLS; rules fixed in spec 10.13 before the run):\n');
fprintf('  R1 reduction V1 %.3f (95%% CI %.3f..%.3f): %s\n', r1.estimate, r1.ciLo, r1.ciHi, tern_local(r1.ciLo > 0.5, 'a MAJORITY is supported (lower limit > 0.50)', 'a majority is NOT supported; quote the interval'));
fprintf('  R2 beta_gen* slope V3 - V1 %+.4f (CI %+.4f..%+.4f): %s\n', d.estimate, d.ciLo, d.ciHi, tern_local(d.ciLo >= -0.02 && d.ciHi <= 0.02, 'NOT materially dependent on the noise-covariate adjustment (CI within +-0.02)', 'SENSITIVE to the adjustment (CI leaves +-0.02); report both'));
fprintf('  R3 beta_gen* V1 slope %+.4f (CI %+.4f..%+.4f): CI %s 0\n', b1.estimate, b1.ciLo, b1.ciHi, tern_local(b1.ciLo > 0 || b1.ciHi < 0, 'excludes', 'includes'));
out = struct('R',R,'est',est,'E',E,'failed',failed,'msg',msg,'EST',EST,'PIPES',PIPES,'opts',opts,'runDate',datetime('now'),'matlab',version);
save(fullfile(simDir,char("bootstrapTempoReduction_v002" + opts.Tag + ".mat")),'out');
fprintf('Saved %s\n', fullfile(simDir,char("bootstrapTempoReduction_v002" + opts.Tag + ".mat")));
end

% =========================================================================
function e = oneRep_local(T0, rowsByPerson, seed, PIPES)
% one person-level resample within each dataset, then the same fits as the point estimate
rng(seed,'twister');  idx = [];  lab = [];  n0 = 0;
for i = 1:numel(rowsByPerson)
    R = rowsByPerson{i};  pick = randi(numel(R), numel(R), 1);
    for q = 1:numel(pick), v = R{pick(q)};  idx = [idx; v];  lab = [lab; repmat(n0 + q, numel(v), 1)]; end %#ok<AGROW>
    n0 = n0 + numel(pick);
end
Tb = T0(idx,:);  Tb.session = "b" + string(lab);  Tb = prep_local(Tb, PIPES);
e = zeros(numel(PIPES), 8);
for k = 1:numel(PIPES), e(k,:) = fitEst_local(Tb, PIPES(k)); end
end

function T = prep_local(T, PIPES)
% checkTempoDecompositionDomain_v001 L28-34: lfb re-centred over sessions, lsa and al re-z-scored, factors rebuilt
z = @(x) (x - mean(x,'omitnan')) ./ std(x,'omitnan');
G = groupsummary(T(T.pipeline == PIPES(1),:), 'session', 'mean', 'lf');
T.lfb = T.lfMean - mean(G.mean_lf);  T.lsa = z(log2(T.sigma ./ T.a));  T.al = z(T.alpha);
T.datasetC = categorical(T.dataset);  T.sessionC = categorical(T.session);
end

function e = fitEst_local(T, p)
% V1 then V3, each for betaObs then betaGen; lfw coefficient; then the derived estimands
V = T(T.pipeline == p,:);  s = nan(1,4);  k = 0;
for rhs = [" ~ 1 + lfw + lfb + lsa + al + datasetC + (1|sessionC)", " ~ 1 + lfw + lfb + datasetC + (1|sessionC)"]
    for oc = ["betaObs","betaGen"]
        k = k + 1;  W = V(isfinite(V.(oc)),:);  m = fitlme(W, oc + rhs);  c = m.Coefficients;
        s(k) = c.Estimate(string(c.Name) == "lfw");
    end
end
e = [s, 1 - s(2)/s(1), 1 - s(4)/s(3), s(4) - s(2), s(3) - s(1)];
end

function anchor_local(est, PIPES, EST, root)
% V1 must reproduce the published table (tolerance 1e-4); V3 must reproduce the empirical adjustment run (5e-4) if it is present
A = load(fullfile(root,'results','checkTempoDecompositionDomain_v001.mat'),'D1');  D = A.D1{2};
oc = struct('bObsV1',"beta_obs",'bGenV1',"beta_gen*");
for k = 1:numel(PIPES)
    for n = ["bObsV1","bGenV1"]
        ref = D.b(D.pipeline == PIPES(k) & D.outcome == oc.(n) & D.term == "lfw");  got = est(k, EST == n);
        fprintf('anchor V1 %-8s %-6s published %+.4f   here %+.4f\n', PIPES(k), n, ref, got);
        if abs(round(got,4) - ref) > 1e-4, error('bootTempo:anchor','V1 %s %s does not reproduce the published slope.', PIPES(k), n); end
    end
end
f = fullfile(root,'results','tempoSim','analyseTempoSlopeRecovery_empirical.mat');
if ~isfile(f), warning('bootTempo:noV3anchor','V3 ANCHOR SKIPPED: %s not found; V3 point estimates are not checked against an earlier run.', f); return, end
W = load(f,'out').out.adjustment;
for k = 1:numel(PIPES)
    for n = ["bObsV3","bGenV3"]
        o = ternary_oc_local(n);  ref = W.bV3(string(W.pipeline) == PIPES(k) & string(W.outcome) == o);  got = est(k, EST == n);
        fprintf('anchor V3 %-8s %-6s earlier run %+.4f   here %+.4f\n', PIPES(k), n, ref, got);
        if abs(got - ref) > 5e-4, error('bootTempo:anchorV3','V3 %s %s differs from the empirical adjustment run by more than 5e-4.', PIPES(k), n); end
    end
end
end

function o = ternary_oc_local(n), o = "betaGen"; if n == "bObsV3", o = "betaObs"; end, end
function v = tern_local(c, a, b), if c, v = a; else, v = b; end, end

function p = ensurePool_local(nWorkers)
% as runTempoSimAll_v001: IdleTimeout disabled, fail loud if no pool
p = gcp('nocreate');
if ~isempty(p)
    try, alive = p.Connected && p.NumWorkers > 0; catch, alive = false; end
    if ~alive, try, delete(p); catch, end, p = []; end
end
if isempty(p), p = parpool('local', nWorkers, 'IdleTimeout', Inf); end
if isempty(p) || ~p.Connected, error('bootTempo:noPool','Could not establish a parallel pool.'); end
end
