function out = analyseR31b_v001(nRep)
%% analyseR31b_v001  R31b (spec 10.58 D): experimenter demonstration. Null-truth synthetic corpus (true beta independent of tempo, spec 10.8).
% Two pseudo-conditions with NO true difference but different tempo distributions; the condition coefficient on beta_obs / beta_gen* is therefore pure bias.
% Design W: conditions within sessions (trial assigned with P = logistic(k * within-session z of log2 f0)). Design B: conditions between sessions (P = logistic(k * z of session-mean log2 f0)).
% Estimators (OLS point estimates, averaged over nRep random assignments; no SEs, so this is bias, not a test): (a) unadjusted cond; (b) adjusted for lfw, lfb, lsa, al;
% (c) tempo-MATCHED: trials resampled to equal counts per 0.1-doubling bin of log2 f0 in each condition, then unadjusted. k = 0 is the control (no tempo difference).
if nargin < 1, nRep = 200; end
ROOT = fileparts(fileparts(mfilename('fullpath')));  addpath(genpath(fullfile(ROOT,'src','functions')));
simDir = fullfile(ROOT,'results','tempoSim');
txt = fullfile(simDir,'analyseR31b_v001.txt');  if isfile(txt), delete(txt); end
diary(txt); cleanup = onCleanup(@() diary('off'));
fprintf('analyseR31b_v001  nRep=%d  %s\n', nRep, char(datetime('now','TimeZone','Europe/London','Format','yyyy-MM-dd HH:mm:ss z')));
CASES = {"AXrF","Fraser"; "AXr","Hickman_PLAC"};  KS = [0 1 2 4];  DES = ["W","B"];  OC = ["betaObs","betaGen"];
z = @(x) (x - mean(x,'omitnan')) ./ std(x,'omitnan');
rows = table();  rng(31);
for c = 1:size(CASES,1)
    arm = CASES{c,1};  ds = CASES{c,2};
    o = analyseTempoSlopeRecovery_v003(arm, Datasets=ds, Tag="R31b_"+ds, KeepTable=true);
    T = o.Tv(o.Tv.pipeline=="SG-IRLS",:);
    lf = log2(T.f0);  [~,~,g] = unique(T.session);  nS = max(g);  sm = accumarray(g, lf, [], @mean);
    lfw = lf - sm(g);  lfb = sm(g) - mean(sm);  lsa = z(log2(T.sigma./T.a));  al = z(T.alpha);
    zW = lfw ./ std(lfw);  zS = (sm - mean(sm)) ./ std(sm);  bin = floor(lf/0.1);
    fprintf('\n%s %s: %d trials, %d sessions, SD lfw %.3f, SD session-mean lf %.3f\n', arm, ds, height(T), nS, std(lfw), std(sm));
    for d = DES
        for k = KS
            B = nan(nRep, numel(OC), 3);  dl = nan(nRep,1);  fk = nan(nRep,1);
            for r = 1:nRep
                if d == "W", cond = rand(height(T),1) < 1./(1+exp(-k*zW));
                else, cs = rand(nS,1) < 1./(1+exp(-k*zS));  cond = cs(g);  if sum(cs) < 3 || sum(~cs) < 3, continue, end, end
                if sum(cond) < 20 || sum(~cond) < 20, continue, end
                dl(r) = mean(lf(cond)) - mean(lf(~cond));
                keep = matchBins_local(cond, bin);  fk(r) = mean(keep);
                for q = 1:numel(OC)
                    y = T.(OC(q));
                    B(r,q,1) = est_local(y, double(cond));
                    B(r,q,2) = est_local(y, [double(cond) lfw lfb lsa al]);
                    B(r,q,3) = est_local(y(keep), double(cond(keep)));
                end
            end
            for q = 1:numel(OC)
                m = squeeze(mean(B(:,q,:),1,'omitnan'))';  s = squeeze(std(B(:,q,:),0,1,'omitnan'))';
                rows = [rows; table(string(ds), d, k, mean(dl,'omitnan'), OC(q), m(1), s(1), m(2), s(2), m(3), s(3), mean(fk,'omitnan'), sum(isfinite(dl)), ...
                    'VariableNames',["dataset","design","k","deltaLog2f0","outcome","biasUnadj","sdUnadj","biasAdj","sdAdj","biasMatched","sdMatched","fracKeptMatched","nValid"])]; %#ok<AGROW>
            end
        end
    end
end
fprintf('\nR31b: condition coefficient when the TRUE condition effect is 0 (= bias). deltaLog2f0 = mean tempo difference between conditions in doublings; sd = across random assignments.\n');
disp(rows);  out = rows;  save(fullfile(simDir,'analyseR31b_v001.mat'), 'rows');
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
