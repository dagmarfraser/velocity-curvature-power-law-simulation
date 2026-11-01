function runOvernightHPC_v001()
%RUNOVERNIGHTHPC_V001  R30 (spec 10.55 as amended in 10.56) and the R29 seam-free set (10.56) on the HPC, one restartable job.
%   Arms AXcF, AXaF, AXsF, AXnF (R30) then AXdF, AXoF, AXlF (R29 seam-free) on Hickman PLAC first, then Fraser (all sessions).
%   Open this file in the MATLAB GUI on the HPC and press Run. Estimated 13-15 h in all (about 2 h per arm x both datasets on 72 workers, from R29; not timed). Restartable.
%   Log: src/logs/Overnight_<stamp>.log. Afterwards the Mac runs analyseR30_v001 and analyseR29SeamFree_v001.   Fraser, D.S. (2026)  v001
srcDir = fileparts(mfilename('fullpath')); addpath(srcDir);
need = ["runTempoSimAll_v005.m","simulateTempoCorpus_v005.m","shapedCommonSurrogate_v001.m","makeRunLoopClosureSynthetic_v001.m","makeRunLoopClosureV015_v001.m","runLoopClosureFftnoise_v013.m"];
for f = need, if ~isfile(fullfile(srcDir,f)), error('runOvernightHPC:missing','Missing %s',fullfile(srcDir,f)); end, end
for ds = ["Fraser","Hickman_PLAC"]
    f = fullfile(srcDir, "loopClosureResults_"+ds+"_all_shaped_xu_v015.mat");
    if ~isfile(f), error('runOvernightHPC:noCorpus','Missing v015 corpus %s',f); end
end
logDir = fullfile(srcDir,'logs'); if ~isfolder(logDir), mkdir(logDir); end
diary(fullfile(logDir, ['Overnight_' char(datetime('now','Format','yyyyMMdd_HHmmss')) '.log']));
cleanup = onCleanup(@() diary('off')); %#ok<NASGU>
[~, host] = system('hostname');
fprintf('Overnight R30 + R29 seam-free HPC | %s | host %s | MATLAB %s | cores %d | src %s\n', string(datetime('now','TimeZone','Europe/London','Format','yyyy-MM-dd HH:mm:ss z')), strtrim(host), version, feature('numcores'), srcDir);
nw = min(72, feature('numcores'));
runTempoSimAll_v005(NWorkers=nw);
fprintf('Overnight HPC finished | %s\n=== COMPLETED OK ===\n', string(datetime('now','TimeZone','Europe/London','Format','yyyy-MM-dd HH:mm:ss z')));
end
