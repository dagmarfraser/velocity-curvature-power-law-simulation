function runMatchedControlSeedsHPC_v002(master, nReps)
% RUNMATCHEDCONTROLSEEDSHPC_V002  v001 with the pool opened at numcores - 1 workers, the size checkAdditiveVsRetimed_v003 asks for, so it never
% deletes and restarts the pool inside a cell (v001 opened the default 72, v003 then restarted at 71 and a worker was killed for memory during
% the restart, BlueBEAR 2026-11-07 16:28). Also prints the pool size and the node memory.
% (v001 header follows.) BlueBEAR launcher for positiveControlMatchedByDataset_v006 (2026-11-07). Runs master seed MASTER for replicates
% 1:NREPS in chunks of 25, one output file per chunk (resumable: a chunk whose file exists is skipped). Planned sessions: master 1 x 200
% replicates (primary), masters 2 and 3 x 50 (stability check). Run in sinteractive from src/ of the _N48 tree (Fraser = N=48; the other
% datasets do not depend on N). Prints COMPLETED OK at the end.   Fraser, D.S. (2026)  v001
if nargin < 2, nReps = 200; end
here = fileparts(mfilename("fullpath"));  cd(here);  logDir = fullfile(here, "logs");  if ~isfolder(logDir), mkdir(logDir); end
diary(fullfile(logDir, sprintf("MCS2_m%d_%s.log", master, string(datetime("now", "Format", "yyyyMMdd_HHmmss")))));  c = onCleanup(@() diary("off"));
fprintf("runMatchedControlSeedsHPC_v002 master %d reps 1:%d | %s\n", master, nReps, string(datetime("now")));
nW = max(1, feature("numcores") - 1);  p = gcp("nocreate");
if ~isempty(p) && p.NumWorkers ~= nW, delete(p); p = []; end
if isempty(p), parpool("Processes", nW); end
fprintf("Pool %d workers on %d cores\n", gcp("nocreate").NumWorkers, feature("numcores"));  [~, m] = system("free -g | head -2");  fprintf("%s\n", m);
for r0 = 1:25:nReps
    reps = r0:min(r0 + 24, nReps);
    tag = sprintf("_m%d_r%d-%d", master, reps(1), reps(end));
    done = dir(fullfile(here, "..", "results", "positiveControlMatchedByDataset_v006" + tag + "_*.mat"));
    if ~isempty(done), fprintf("chunk %s exists, skipped\n", tag); continue, end
    t0 = tic;  positiveControlMatchedByDataset_v006(Master=master, Reps=reps);
    fprintf("chunk %s done in %.1f min | %s\n", tag, toc(t0)/60, string(datetime("now")));
end
fprintf("=== COMPLETED OK ===\n");
end
