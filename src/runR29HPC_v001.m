function runR29HPC_v001()
%RUNR29HPC_V001  R29 on the HPC (spec SPEC_T2 10.52): simulate arms AXvF and AXwF (Hickman PLAC, Fraser all sessions) and run them through the runner.
%   Open this file in the MATLAB GUI on the HPC and press Run. About 2.2 h per arm on 72 workers (estimate, not timed for these datasets). Restartable.
%   Log: src/logs/R29_<stamp>.log. Afterwards the Mac runs analyseR29_v001.   Fraser, D.S. (2026)  v001
srcDir = fileparts(mfilename('fullpath')); addpath(srcDir);
need = ["runTempoSimAll_v003.m","simulateTempoCorpus_v003.m","makeRunLoopClosureSynthetic_v001.m","makeRunLoopClosureV015_v001.m","runLoopClosureFftnoise_v013.m"];
for f = need, if ~isfile(fullfile(srcDir,f)), error('runR29HPC:missing','Missing %s',fullfile(srcDir,f)); end, end
for ds = ["Fraser","Hickman_PLAC"]
    f = fullfile(srcDir, "loopClosureResults_"+ds+"_all_shaped_xu_v015.mat");
    if ~isfile(f), error('runR29HPC:noCorpus','Missing v015 corpus %s',f); end
end
logDir = fullfile(srcDir,'logs'); if ~isfolder(logDir), mkdir(logDir); end
diary(fullfile(logDir, ['R29_' char(datetime('now','Format','yyyyMMdd_HHmmss')) '.log']));
cleanup = onCleanup(@() diary('off')); %#ok<NASGU>
[~, host] = system('hostname');
fprintf('R29 HPC | %s | host %s | MATLAB %s | cores %d | src %s\n', string(datetime('now','TimeZone','Europe/London','Format','yyyy-MM-dd HH:mm:ss z')), strtrim(host), version, feature('numcores'), srcDir);
nw = min(72, feature('numcores'));
runTempoSimAll_v003(NWorkers=nw);
fprintf('R29 HPC finished | %s\n=== COMPLETED OK ===\n', string(datetime('now','TimeZone','Europe/London','Format','yyyy-MM-dd HH:mm:ss z')));
end
