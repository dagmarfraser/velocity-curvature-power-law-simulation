function out = analyseR31c_v001(nRep)
%% analyseR31c_v001  R31c scenario S2 (spec 10.62): does adjusting for tempo, or matching on it, keep a REAL condition difference?
% Arms B1Xr (true within-session slope +0.0265) and B2Xr (+0.084) from the T2 chain; Fraser = full-corpus arms (suffix F), Hickman PLAC = no suffix.
% Same assignment designs and estimators as analyseR31b_v001 (W: within-session; B: between-session; unadjusted / adjusted / tempo-matched; OLS point estimates).
% Difference from R31b: the TRUE condition effect is read per assignment from the generator (truth.betaTrue), not assumed zero. Spec 10.62 (amended) explains why:
% the B arms put the slope on the WITHIN-session deviation only, so a between-session assignment (design B) has a true effect near 0 whatever the tempo gap.
% Stated expectation (fixed in spec 10.62 before this script was run): unadjusted beta_gen* tracks the true total effect; adjusted and matched estimates track
% the true effect remaining after tempo is equalised (about 0). Reported as estimate minus truth; no pass/fail threshold.
% Fails loud if an arm file is missing (the B arms were run on the HPC and are not in a fresh clone).
if nargin < 1, nRep = 200; end
ROOT = fileparts(fileparts(mfilename('fullpath')));  addpath(genpath(fullfile(ROOT,'src','functions')));
simDir = fullfile(ROOT,'results','tempoSim');
txt = fullfile(simDir,'analyseR31c_v001.txt');  if isfile(txt), delete(txt); end
diary(txt); cleanup = onCleanup(@() diary('off'));
fprintf('analyseR31c_v001  nRep=%d  %s\n', nRep, char(datetime('now','TimeZone','Europe/London','Format','yyyy-MM-dd HH:mm:ss z')));
CASES = {"B1XrF","Fraser"; "B2XrF","Fraser"; "B1Xr","Hickman_PLAC"; "B2Xr","Hickman_PLAC"};  KS = [1 2 4];  DES = ["W","B"];  OC = ["betaObs","betaGen"];
need = strings(0,1);
for c = 1:size(CASES,1)
    need = [need; fullfile(simDir, "tempoSim_" + CASES{c,1} + "_" + CASES{c,2} + ".mat"); ...
        fullfile(simDir, "loopClosureResults_SYN-" + extractBefore(CASES{c,1}+"-", "-") + "-" + CASES{c,2} + "_all_shaped_xu_syn001.mat")]; %#ok<AGROW>
end
miss = need(~arrayfun(@isfile, need));
if ~isempty(miss), error('analyseR31c:input', 'Missing arm files (copy from RDS results/tempoSim, no overwrite):\n%s', strjoin(miss, newline)); end
z = @(x) (x - mean(x,'omitnan')) ./ std(x,'omitnan');
rows = table();  rng(32);
for c = 1:size(CASES,1)
    arm = CASES{c,1};  ds = CASES{c,2};
    o = analyseTempoSlopeRecovery_v003(arm, Datasets=ds, Tag="R31c_"+ds, KeepTable=true);
    T = o.Tv(o.Tv.pipeline=="SG-IRLS",:);
    tr = load(fullfile(simDir,"tempoSim_"+arm+"_"+ds+".mat"),'truth').truth;
    [ok, loc] = ismember(T.trialSrc, extractAfter(string(tr.trialID), "SYN-"+arm+"-"));
    if ~all(ok), error('analyseR31c:truth', '%s %s: %d trials without generator truth', arm, ds, sum(~ok)); end
    bT = tr.betaTrue(loc);
    lf = log2(T.f0);  [~,~,g] = unique(T.session);  nS = max(g);  sm = accumarray(g, lf, [], @mean);
    lfw = lf - sm(g);  lfb = sm(g) - mean(sm);  lsa = z(log2(T.sigma./T.a));  al = z(T.alpha);
    zW = lfw ./ std(lfw);  zS = (sm - mean(sm)) ./ std(sm);  bin = floor(lf/0.1);
    fprintf('\n%s %s: %d trials, %d sessions\n', arm, ds, height(T), nS);
    for d = DES
        for k = KS
            B = nan(nRep, numel(OC), 3);  dl = nan(nRep,1);  tf = nan(nRep,2);
            for r = 1:nRep
                if d == "W", cond = rand(height(T),1) < 1./(1+exp(-k*zW));
                else, cs = rand(nS,1) < 1./(1+exp(-k*zS));  cond = cs(g);  if sum(cs) < 3 || sum(~cs) < 3, continue, end, end
                if sum(cond) < 20 || sum(~cond) < 20, continue, end
                dl(r) = mean(lf(cond)) - mean(lf(~cond));
                keep = matchBins_local(cond, bin);
                if sum(keep & cond) < 10 || sum(keep & ~cond) < 10, continue, end
                tf(r,1) = mean(bT(cond)) - mean(bT(~cond));
                tf(r,2) = mean(bT(keep & cond)) - mean(bT(keep & ~cond));
                for q = 1:numel(OC)
                    y = T.(OC(q));
                    B(r,q,1) = est_local(y, double(cond));
                    B(r,q,2) = est_local(y, [double(cond) lfw lfb lsa al]);
                    B(r,q,3) = est_local(y(keep), double(cond(keep)));
                end
            end
            for q = 1:numel(OC)
                m = squeeze(mean(B(:,q,:),1,'omitnan'))';  s = squeeze(std(B(:,q,:),0,1,'omitnan'))';
                rows = [rows; table(string(arm), string(ds), d, k, mean(dl,'omitnan'), OC(q), mean(tf(:,1),'omitnan'), mean(tf(:,2),'omitnan'), ...
                    m(1), s(1), m(2), s(2), m(3), s(3), sum(isfinite(dl)), ...
                    'VariableNames',["arm","dataset","design","k","deltaLog2f0","outcome","trueEffect","trueEffectMatched","unadj","sdUnadj","adj","sdAdj","matched","sdMatched","nValid"])]; %#ok<AGROW>
            end
        end
    end
end
rows.errUnadj = rows.unadj - rows.trueEffect;  rows.errAdj = rows.adj - rows.trueEffect;  rows.errMatched = rows.matched - rows.trueEffectMatched;
fprintf('\nR31c S2: estimated minus true condition effect, per outcome (true effect read from the generator; matched compared with the true effect in the matched sample).\n');
disp(rows);  out = rows;  save(fullfile(simDir,'analyseR31c_v001.mat'), 'rows');
fprintf('\n=== DONE ===\n');
end

function b = est_local(y, X)
ok = isfinite(y) & all(isfinite(X),2);
if sum(ok) < size(X,2)+5, b = NaN; return, end
c = [ones(sum(ok),1) X(ok,:)] \ y(ok);  b = c(2);
end

function keep = matchBins_local(cond, bin)
% equal counts per tempo bin in the two conditions (random subsample of the larger side)
keep = false(size(cond));
for u = unique(bin)'
    i1 = find(bin==u & cond);  i0 = find(bin==u & ~cond);  n = min(numel(i1), numel(i0));
    if n == 0, continue, end
    keep(i1(randperm(numel(i1), n))) = true;  keep(i0(randperm(numel(i0), n))) = true;
end
end
