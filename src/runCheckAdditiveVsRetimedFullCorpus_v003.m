%% runCheckAdditiveVsRetimedFullCorpus_v003
% Full-corpus generative-model check (Findings #222/#223, D20):
%   1. RUN      checkAdditiveVsRetimed_v002 over every trial, all seven arms
%               -> results/checkAdditiveVsRetimed_v002_<stamp>.mat
%   2. ANALYSE  analyseCheckAdditiveVsRetimed_v001 on the table THIS run wrote
%               -> results/checkAdditiveVsRetimed_subjectSummary_v002_<stamp>.mat
% To re-analyse without re-running, call step 2 directly on a saved table:
%   analyseCheckAdditiveVsRetimed_v001("results/checkAdditiveVsRetimed_v002_<stamp>.mat")
% Cost: ~4 min on the M5, ~2 h on the iMac (pool = cores-1).
%
% Fraser, D.S. (2026) v003 (v002 with analysis moved to its own function;
% CONFIG and run unchanged)

%% CONFIG
CFG.Datasets       = ["Fraser","Cook CTRL","Cook ASD","Hickman PLAC", ...
                      "Hickman HALO","Zarandi","Dhieb"];
CFG.TrialSelection = "all";
CFG.NSurr          = 10;
CFG.EdgeClip       = 20;
CFG.RngSeed        = 1729;
CFG.SepK           = 2;
CFG.CycleBandMult  = 3;
CFG.UseParfor      = true;
CFG.MinSubj        = 6;

%% RUN
srcDir = fileparts(mfilename("fullpath"));
cd(srcDir);
tStart = datetime("now");
fprintf("%d cores | start %s\n", feature("numcores"), string(tStart));
checkAdditiveVsRetimed_v002(Datasets=CFG.Datasets, TrialSelection=CFG.TrialSelection, ...
    NSurr=CFG.NSurr, EdgeClip=CFG.EdgeClip, RngSeed=CFG.RngSeed, SepK=CFG.SepK, ...
    CycleBandMult=CFG.CycleBandMult, UseParfor=CFG.UseParfor);
fprintf("Run time %.1f min\n", minutes(datetime("now") - tStart));

%% LOCATE the table this run wrote (newest v002 table, created after tStart)
F = dir(fullfile(srcDir, "results", "checkAdditiveVsRetimed_v002_*.mat"));
F = F([F.datenum] >= datenum(tStart) - 1/86400);  % 1 s tolerance
if numel(F) ~= 1
    error("runCheckAdditiveVsRetimed_v003:TableNotFound", "%s", sprintf( ...
        "Expected exactly one trial table written since %s, found %d.", string(tStart), numel(F)));
end
trialFile = string(fullfile(F.folder, F.name));

%% ANALYSE
Ssum = analyseCheckAdditiveVsRetimed_v001(trialFile, MinSubj=CFG.MinSubj);
