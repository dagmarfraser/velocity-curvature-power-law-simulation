function runT2H1FraserHPC_v002()
%RUNT2H1FRASERHPC_V002  v001 plus the two RECOVERY arms B1XrF and B2XrF (true within-session slopes 0.0265 and 0.084; simulateTempoCorpus_v005 cases B1, B2), 2026-11-05.
%   Why: runT2FraserAll_v001 reads AXrF, B2XrF, AXrs2F, B1XrF; the N=45 results include all four (R14/R15 recovery readings), v001 produced only the two null arms, and the
%   recovery readings (B2XrF PARTIAL, B1XrF NOT IDENTIFIED at N=45) cannot be reproduced without B1XrF and B2XrF. AXrF and AXrs2F already exist in the copy and are skipped; only B1XrF and B2XrF compute (about 2 h each).
%   Same settings as v001 (N_REPS 50, Fraser only, all sessions). Log: src/logs/H1Fraser_<stamp>.log. After it finishes, rerun runT2FraserAll_v001.
%RUNT2H1FRASERHPC_V001  T2 stage H1 for Fraser N=48: the null-truth arms AXrF and AXrs2F (spec SPEC_T2 10.8, 10.24), one restartable job.
%   Written 2026-11-04 because stage H1 had no recorded launcher (runT2Diagnosis_v001: "the exact arguments that made these files are not
%   recorded in the code"). The arguments here are read from the N=45 arm files, not guessed: loopClosureResults_SYN-AXrF-Fraser config has
%   N_REPS 50, N_BETA 25, N_R7 10, seed 1729, shaped_xu, 72 workers; tempoSim_AXrF_Fraser has maxTrials = Inf; arm names ending in F keep ALL
%   Fraser sessions (simulateTempoCorpus_v005 line 49). It runs runTempoSimAll_v005 with Arms = AXrF, AXrs2F and Datasets = Fraser only
%   (Hickman PLAC's AXr / AXrs2 arms exist and are not touched).
%   Open this file in the MATLAB GUI on the HPC, from the _N48 copy's src/ (sinteractive, 72 cores) and press Run. Estimated about 2 h per arm
%   (R29 timing, not timed for this arm). Restartable: an arm whose output exists is skipped.
%   Output: results/tempoSim/tempoSim_AXrF_Fraser.mat, loopClosureResults_SYN-AXrF-Fraser_all_shaped_xu_syn001.mat (and AXrs2F), summary
%   results/tempoSim/runTempoSimAll_H1_FraserN48_summary.mat. Log: src/logs/H1Fraser_<stamp>.log.
%   Not included: runT2FraserOnly_v001's un-suffixed Fraser arms (AXr, AXrs2, Sm1Xr, Sm2Xr, Sp1Xr, B1Xr, B2Xr on the legacy retained set): no such
%   Fraser files exist in the N=45 results either; that script is the record of an earlier analysis.   Fraser, D.S. (2026)  v001
srcDir = fileparts(mfilename('fullpath')); addpath(srcDir);
need = ["runTempoSimAll_v005.m","simulateTempoCorpus_v005.m","shapedCommonSurrogate_v001.m","makeRunLoopClosureSynthetic_v001.m","makeRunLoopClosureV015_v001.m","runLoopClosureFftnoise_v013.m"];
for f = need, if ~isfile(fullfile(srcDir,f)), error('runT2H1FraserHPC:missing','Missing %s',fullfile(srcDir,f)); end, end
f = fullfile(srcDir, "loopClosureResults_Fraser_all_shaped_xu_v015.mat");
if ~isfile(f), error('runT2H1FraserHPC:noCorpus','Missing v015 corpus %s',f); end
if ~contains(string(srcDir), "_N48"), error('runT2H1FraserHPC:notN48', '%s', "Not the N=48 copy: " + string(srcDir) + " (open the launcher from <project>_N48/src)."); end
logDir = fullfile(srcDir,'logs'); if ~isfolder(logDir), mkdir(logDir); end
diary(fullfile(logDir, ['H1FraserB_' char(datetime('now','Format','yyyyMMdd_HHmmss')) '.log']));
cleanup = onCleanup(@() diary('off')); %#ok<NASGU>
[~, host] = system('hostname');
fprintf('T2 H1 Fraser v002 (null + recovery arms) HPC | %s | host %s | MATLAB %s | cores %d | src %s\n', string(datetime('now','TimeZone','Europe/London','Format','yyyy-MM-dd HH:mm:ss z')), strtrim(host), version, feature('numcores'), srcDir);
for a = ["AXrF" "AXrs2F" "B1XrF" "B2XrF"]
    for g = ["tempoSim_" + a + "_Fraser.mat", "loopClosureResults_SYN-" + a + "-Fraser_all_shaped_xu_syn001.mat"]
        p = fullfile(srcDir, '..', 'results', 'tempoSim', g);  if isfile(p), fprintf('  existing: %s\n', g); end
    end
end
nw = min(72, feature('numcores'));
runTempoSimAll_v005(Arms=["AXrF" "AXrs2F" "B1XrF" "B2XrF"], Datasets="Fraser", N_REPS=50, NWorkers=nw, SummaryFile="runTempoSimAll_H1_FraserN48_B_summary.mat");
for a = ["AXrF" "AXrs2F" "B1XrF" "B2XrF"]
    p = fullfile(srcDir, '..', 'results', 'tempoSim', "loopClosureResults_SYN-" + a + "-Fraser_all_shaped_xu_syn001.mat");
    if ~isfile(p), error('runT2H1FraserHPC:noOutput', '%s', "Expected output missing: " + p); end
end
fprintf('T2 H1 Fraser v002 finished | %s\n=== COMPLETED OK ===\n', string(datetime('now','TimeZone','Europe/London','Format','yyyy-MM-dd HH:mm:ss z')));
end
