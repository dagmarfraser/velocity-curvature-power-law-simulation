function checkR29Arms_v001()
%CHECKR29ARMS_V001  Smoke test of simulateTempoCorpus_v003 (spec 10.52, R29): Xp unchanged from v002; Xv and Xw share Xp's draws; realised noise-level/tempo coupling; alpha structure.
%   Small Hickman PLAC sample in a temp folder (nothing enters results/). Writes results/tempoSim/checkR29Arms_v001.txt.   Fraser, D.S. (2026)  v001
srcDir = fileparts(mfilename('fullpath')); addpath(srcDir);
txt = fullfile(srcDir,'..','results','tempoSim','checkR29Arms_v001.txt'); if isfile(txt), delete(txt); end
diary(txt); c = onCleanup(@() diary('off')); %#ok<NASGU>
fprintf('checkR29Arms_v001 %s\n', char(datetime('now','TimeZone','Europe/London','Format','yyyy-MM-dd HH:mm:ss z')));
od = fullfile(tempname,'r29'); mkdir(od); mkdir(fullfile(od,'v3'));
ds = "Hickman PLAC"; N = 40;
o2 = simulateTempoCorpus_v002("AXpF",Datasets=ds,MaxTrials=N,OutDir=od);
o3 = simulateTempoCorpus_v003("AXpF",Datasets=ds,MaxTrials=N,OutDir=fullfile(od,'v3'));
ov = simulateTempoCorpus_v003("AXvF",Datasets=ds,MaxTrials=N,OutDir=fullfile(od,'v3'));
ow = simulateTempoCorpus_v003("AXwF",Datasets=ds,MaxTrials=N,OutDir=fullfile(od,'v3'));
A = load(o2.files(1)); B = load(o3.files(1)); V = load(ov.files(1)); W = load(ow.files(1));
fprintf('bioResults fields: %s\n', strjoin(string(B.bioResults.Properties.VariableNames), ' '));
fprintf('CHECK1 Xp v002 vs v003: trials identical %d | betaTrue identical %d | alpha identical %d\n', isequal(A.trials,B.trials), isequal(A.truth.betaTrue,B.truth.betaTrue), isequaln(A.bioResults.ira_alphaMean,B.bioResults.ira_alphaMean));
fprintf('CHECK2 same betaTrue as Xp: Xv %d Xw %d | same donors: Xv %d Xw %d\n', ...
    isequal(V.truth.betaTrue,B.truth.betaTrue), isequal(W.truth.betaTrue,B.truth.betaTrue), isequal(V.truth.donorTrialID,B.truth.donorTrialID), isequal(W.truth.donorTrialID,B.truth.donorTrialID));
fprintf('CHECK3 realised within-session corr(noise level, tempo): Xp %+.3f | Xv %+.3f | Xw %+.3f   (own residual %+.3f, donor %+.3f)\n', B.dDiag.corrReal, V.dDiag.corrReal, W.dDiag.corrReal, B.dDiag.corrOwn, B.dDiag.corrDonor);
a = @(S) S.bioResults.ira_alphaMean;
fprintf('CHECK4 median estimated alpha: Xp %.3f | Xv %.3f | Xw %.3f ; corr(alpha Xp, alpha Xv) %.3f, corr(alpha Xp, alpha Xw) %.3f\n', median(a(B),'omitnan'), median(a(V),'omitnan'), median(a(W),'omitnan'), ...
    corr(a(B),a(V),'rows','complete'), corr(a(B),a(W),'rows','complete'));
fprintf('trials %d (Xp) %d (Xv) %d (Xw)\n', numel(B.trials), numel(V.trials), numel(W.trials));
fprintf('=== DONE ===\n');
end
