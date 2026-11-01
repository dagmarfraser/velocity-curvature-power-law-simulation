function checkR30Arms_v001()
%CHECKR30ARMS_V001  Smoke test of simulateTempoCorpus_v004 and shapedCommonSurrogate_v001 (spec 10.55 as amended in 10.56, R30).
%   (1) v004 bit-reproduces v003 Xp and Xw and v002 Xr; (2) Xc/Xa/Xs/Xn take bTrue from Xr (Xs may lose trials without a long-enough donor);
%   (3) per-axis SD preserved (corrReal == corrOwn: exact for Xc/Xa/Xs, approximate for Xn); (4) unit test of the Xa surrogate (cross-spectrum kept, SD kept);
%   (5) alpha medians vs Xr, with per-trial gaps; (6) Xs never wraps. Small Hickman PLAC sample in a temp folder.
%   Writes results/tempoSim/checkR30Arms_v001.txt.   Fraser, D.S. (2026)  v001 (revised for 10.56)
srcDir = fileparts(mfilename('fullpath')); addpath(srcDir); addpath(genpath(fullfile(srcDir,'functions'))); addpath(genpath(fullfile(srcDir,'req')));
txt = fullfile(srcDir,'..','results','tempoSim','checkR30Arms_v001.txt'); if isfile(txt), delete(txt); end
diary(txt); c = onCleanup(@() diary('off')); %#ok<NASGU>
fprintf('checkR30Arms_v001 %s\n', char(datetime('now','TimeZone','Europe/London','Format','yyyy-MM-dd HH:mm:ss z')));
% ---- unit test of the Xa surrogate on a coloured correlated pair ----
rs = RandStream('twister','Seed',1); FSu = 133; n = 800;
a = cumsum(cumsum(cumsum(randn(rs,n,1))))*1e-3; a = a - mean(a); b = 0.6*a + cumsum(cumsum(randn(rs,n,1)))*1e-2; b = b - mean(b);
[na, nb] = shapedCommonSurrogate_v001(a, b, FSu, 4.0, 4.0);
Za = fft(a); Zb = fft(b); Zna = fft(na); Znb = fft(nb);
x0 = conj(Za).*Zb; x1 = conj(Zna).*Znb; kk = abs(x0) > 1e-3*max(abs(x0)); pd = angle(x1(kk).*conj(x0(kk)));
fprintf('UNIT Xa: real %d | SD ratio major %.6f minor %.6f | cross-spectrum phase max |err| %.2e rad over %d bins | magnitude ratio spread (std/mean of |x1|/|x0|) %.3f | zero-lag corr own %.3f surrogate %.3f\n', ...
    isreal(na)&&isreal(nb), std(na)/std(a), std(nb)/std(b), max(abs(pd)), nnz(kk), std(abs(x1(kk))./abs(x0(kk)))/mean(abs(x1(kk))./abs(x0(kk))), corr(a,b), corr(na,nb));
% ---- corpus arms ----
od = fullfile(tempname,'r30'); mkdir(od); d2 = fullfile(od,'v2'); d3 = fullfile(od,'v3'); d4 = fullfile(od,'v4'); mkdir(d2); mkdir(d3); mkdir(d4);
ds = "Hickman PLAC"; N = 40; ld = @(o) load(o.files(1));
r2 = ld(simulateTempoCorpus_v002("AXrF",Datasets=ds,MaxTrials=N,OutDir=d2));
p3 = ld(simulateTempoCorpus_v003("AXpF",Datasets=ds,MaxTrials=N,OutDir=d3));
w3 = ld(simulateTempoCorpus_v003("AXwF",Datasets=ds,MaxTrials=N,OutDir=d3));
r4 = ld(simulateTempoCorpus_v004("AXrF",Datasets=ds,MaxTrials=N,OutDir=d4));
p4 = ld(simulateTempoCorpus_v004("AXpF",Datasets=ds,MaxTrials=N,OutDir=d4));
w4 = ld(simulateTempoCorpus_v004("AXwF",Datasets=ds,MaxTrials=N,OutDir=d4));
fprintf('CHECK1 v004 vs earlier: Xr trials identical %d | Xp %d | Xw %d | Xp alpha identical %d | Xw alpha identical %d\n', ...
    isequal(r2.trials,r4.trials), isequal(p3.trials,p4.trials), isequal(w3.trials,w4.trials), isequaln(p3.bioResults.ira_alphaMean,p4.bioResults.ira_alphaMean), isequaln(w3.bioResults.ira_alphaMean,w4.bioResults.ira_alphaMean));
arms = ["AXcF","AXaF","AXsF","AXnF"]; S = struct();
for q = 1:numel(arms)
    S.(char(arms(q))) = ld(simulateTempoCorpus_v004(arms(q),Datasets=ds,MaxTrials=N,OutDir=d4));
end
fprintf('trials: Xr %d | Xp %d | Xc %d | Xa %d | Xs %d | Xn %d  (Xs skips trials with no donor at least as long, so its set can be smaller)\n', numel(r4.trials), numel(p4.trials), numel(S.AXcF.trials), numel(S.AXaF.trials), numel(S.AXsF.trials), numel(S.AXnF.trials));
for q = 1:numel(arms)
    A = S.(char(arms(q))); same = numel(A.trials)==numel(r4.trials) && isequal(A.truth.betaTrue, r4.truth.betaTrue);
    fprintf('CHECK2 %s: bTrue equal to Xr %d | same trial set as Xr %d', arms(q), same, numel(A.trials)==numel(r4.trials));
    if arms(q) == "AXsF", fprintf(' | wrapped %d (must be 0)', sum(A.truth.wrapped)); end
    fprintf('\n');
    fprintf('CHECK3 %s: corrOwn %+.6f corrReal %+.6f | |diff| %.2e (exact for Xc/Xa/Xs, approximate for Xn)\n', arms(q), A.dDiag.corrOwn, A.dDiag.corrReal, abs(A.dDiag.corrOwn-A.dDiag.corrReal));
end
fprintf('CHECK4 median estimated alpha: Xr %.3f | Xp %.3f | Xc %.3f | Xa %.3f | Xs %.3f | Xn %.3f   (criterion: Xc, Xa within 0.05 of Xr; Xs compared with Xr only for information)\n', median(r4.bioResults.ira_alphaMean,'omitnan'), median(p4.bioResults.ira_alphaMean,'omitnan'), ...
    median(S.AXcF.bioResults.ira_alphaMean,'omitnan'), median(S.AXaF.bioResults.ira_alphaMean,'omitnan'), median(S.AXsF.bioResults.ira_alphaMean,'omitnan'), median(S.AXnF.bioResults.ira_alphaMean,'omitnan'));
ar = r4.bioResults.ira_alphaMean;
for q = 1:numel(arms)
    A = S.(char(arms(q))); [tf, loc] = ismember(string(A.truth.trialID), strrep(string(r4.truth.trialID),'AXrF',char(arms(q))));
    d = A.bioResults.ira_alphaMean(tf) - ar(loc(tf));
    fprintf('CHECK4b %s: per-trial alpha minus Xr alpha: median %+.3f, IQR [%+.3f %+.3f], share within 0.5: %.0f%%\n', arms(q), median(d,'omitnan'), prctile(d,25), prctile(d,75), 100*mean(abs(d)<=0.5));
end
fprintf('CHECK5 median estimated sigma: Xr %.4f | Xc %.4f | Xa %.4f | Xs %.4f | Xn %.4f\n', median(r4.bioResults.sigmaMean,'omitnan'), median(S.AXcF.bioResults.sigmaMean,'omitnan'), median(S.AXaF.bioResults.sigmaMean,'omitnan'), median(S.AXsF.bioResults.sigmaMean,'omitnan'), median(S.AXnF.bioResults.sigmaMean,'omitnan'));
fprintf('=== DONE ===\n');
end
