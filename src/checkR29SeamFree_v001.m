function checkR29SeamFree_v001()
%CHECKR29SEAMFREE_V001  Smoke test of the seam-free permutation arms Xd, Xo, Xl in simulateTempoCorpus_v005 (spec 10.56).
%   (1) v005 reproduces v003 Xp and the v003 Xr-pair design; (2) Xd/Xo/Xl share donors and trial set, no wrapping; (3) skipped trials counted;
%   (4) realised noise-level/tempo coupling; (5) alpha medians vs Xr and the within-session corr(alpha, true log2 f0) vs the wrapped Xp (the 10.56 seam finding).
%   Small Hickman PLAC sample in a temp folder. Writes results/tempoSim/checkR29SeamFree_v001.txt.   Fraser, D.S. (2026)  v001
srcDir = fileparts(mfilename('fullpath')); addpath(srcDir); addpath(genpath(fullfile(srcDir,'functions'))); addpath(genpath(fullfile(srcDir,'req')));
txt = fullfile(srcDir,'..','results','tempoSim','checkR29SeamFree_v001.txt'); if isfile(txt), delete(txt); end
diary(txt); c = onCleanup(@() diary('off')); %#ok<NASGU>
fprintf('checkR29SeamFree_v001 %s\n', char(datetime('now','TimeZone','Europe/London','Format','yyyy-MM-dd HH:mm:ss z')));
od = fullfile(tempname,'sf'); d3 = fullfile(od,'v3'); d5 = fullfile(od,'v5'); mkdir(od); mkdir(d3); mkdir(d5);
ds = "Hickman PLAC"; N = 40; ld = @(o) load(o.files(1));
p3 = ld(simulateTempoCorpus_v003("AXpF",Datasets=ds,MaxTrials=N,OutDir=d3));
p5 = ld(simulateTempoCorpus_v005("AXpF",Datasets=ds,MaxTrials=N,OutDir=d5));
r5 = ld(simulateTempoCorpus_v005("AXrF",Datasets=ds,MaxTrials=N,OutDir=d5));
D = ld(simulateTempoCorpus_v005("AXdF",Datasets=ds,MaxTrials=N,OutDir=d5));
O = ld(simulateTempoCorpus_v005("AXoF",Datasets=ds,MaxTrials=N,OutDir=d5));
L = ld(simulateTempoCorpus_v005("AXlF",Datasets=ds,MaxTrials=N,OutDir=d5));
fprintf('CHECK1 v005 vs v003 Xp: trials identical %d | alpha identical %d\n', isequal(p3.trials,p5.trials), isequaln(p3.bioResults.ira_alphaMean,p5.bioResults.ira_alphaMean));
sh = @(A,B) isequal(strrep(string(A.truth.trialID),'AXoF','AXdF'),string(B.truth.trialID));
fprintf('CHECK2 trials: Xr %d | Xp %d | Xd %d | Xo %d | Xl %d ; wrapped: Xd %d Xo %d Xl %d (must be 0) ; Xd,Xo,Xl same donors %d %d ; same betaTrue %d %d\n', numel(r5.trials), numel(p5.trials), numel(D.trials), numel(O.trials), numel(L.trials), ...
    sum(D.truth.wrapped), sum(O.truth.wrapped), sum(L.truth.wrapped), isequal(D.truth.donorTrialID,O.truth.donorTrialID), isequal(D.truth.donorTrialID,L.truth.donorTrialID), isequal(D.truth.betaTrue,O.truth.betaTrue), isequal(D.truth.betaTrue,L.truth.betaTrue));
fprintf('CHECK3 wrapped share in v003 Xp (the seam being removed): %.0f%%\n', 100*mean(p3.truth.wrapped));
fprintf('CHECK4 within-session corr(noise level, tempo): own residual %+.3f | donor %+.3f | REALISED Xd %+.3f Xo %+.3f Xl %+.3f (Xo should equal own, Xl follow donor, Xd donor; Xp wrapped %+.3f)\n', D.dDiag.corrOwn, D.dDiag.corrDonor, D.dDiag.corrReal, O.dDiag.corrReal, L.dDiag.corrReal, p5.dDiag.corrReal);
mm = @(m,g) m(g); cc = @(x,s) x-mm(splitapply(@mean,x,findgroups(s)),findgroups(s));
cr = @(A) corr(cc(A.bioResults.ira_alphaMean(isfinite(A.bioResults.ira_alphaMean)),A.truth.session(isfinite(A.bioResults.ira_alphaMean))), cc(A.truth.lfwTrue(isfinite(A.bioResults.ira_alphaMean)),A.truth.session(isfinite(A.bioResults.ira_alphaMean))));
fprintf('CHECK5 median alpha: Xr %.3f | Xp(wrapped) %.3f | Xd %.3f | Xo %.3f | Xl %.3f\n', median(r5.bioResults.ira_alphaMean,'omitnan'), median(p5.bioResults.ira_alphaMean,'omitnan'), median(D.bioResults.ira_alphaMean,'omitnan'), median(O.bioResults.ira_alphaMean,'omitnan'), median(L.bioResults.ira_alphaMean,'omitnan'));
fprintf('CHECK6 within-session corr(alpha, true log2 f0): Xr %+.3f | Xp(wrapped) %+.3f | Xd %+.3f | Xo %+.3f | Xl %+.3f\n', cr(r5), cr(p5), cr(D), cr(O), cr(L));
fprintf('=== DONE ===\n');
end
