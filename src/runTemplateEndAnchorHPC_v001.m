function runTemplateEndAnchorHPC_v001()
%RUNTEMPLATEENDANCHORHPC_V001  BlueBEAR launcher for templateEndAnchorTest_v001 (end-anchored last window, causal test of the template tail; SESSION_LOG Addendum 47). EXPLORATORY.
%   Part 1: simulation, 3 datasets x 432 trials x 8 arms (same seed/trials as templateVariantTest_v004). Part 2: real trials, 5 datasets x 300 random trials x 2 arms.
%   Open in the MATLAB GUI on the HPC (sinteractive, 72 cores) and press Run. No arguments. Resumable (chunks results/tea1_*_chunks); pauses safely at the deadline, re-run to continue.
%   Output results/templateEndAnchorTest_v001.mat, log src/logs/TEA_<stamp>.log.   Fraser, D.S. (2026)  v001
WALL_HOURS = 11;  NWORKERS = 72;
srcDir = fileparts(mfilename('fullpath'));  addpath(srcDir);
need = ["templateEndAnchorTest_v001.m" "loopClosureCore_v004.m" "synthWorld_v001.m"];
for f = need, if ~isfile(fullfile(srcDir, f)), error('runTea:missing', 'Missing %s in %s', f, srcDir); end, end
logDir = fullfile(srcDir, 'logs');  if ~isfolder(logDir), mkdir(logDir); end
diary(fullfile(logDir, ['TEA_' char(datetime('now', 'Format', 'yyyyMMdd_HHmmss')) '.log']));
cleanup = onCleanup(@() diary('off')); %#ok<NASGU>
p0 = gcp('nocreate');  if ~isempty(p0), delete(p0); end   % a pool left from an earlier run keeps its workers' cached functions
clear functions;  rehash;
cleanupPool = onCleanup(@() delete(gcp('nocreate'))); %#ok<NASGU>
[~, host] = system('hostname');
fprintf('End-anchored template test HPC | %s | host %s | MATLAB %s | cores %d\n', stamp_local(), strtrim(host), version, feature('numcores'));
R = templateEndAnchorTest_v001(NWorkers=NWORKERS, Deadline=datetime('now') + hours(WALL_HOURS));
if isempty(R), fprintf('PAUSED at the deadline. Open this file and Run again to resume.\n');  return, end
fprintf('End-anchored template test finished | %s\n=== COMPLETED OK ===\n', stamp_local());
end

function s = stamp_local()
    s = string(datetime('now', 'TimeZone', 'Europe/London', 'Format', 'yyyy-MM-dd HH:mm:ss z'));
end
