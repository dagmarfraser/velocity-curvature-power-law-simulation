function runV022HPC_v002(datasets, deadline, opts)
%RUNV022HPC_V002  Batch (sbatch) launcher for runLoopClosureFftnoise_v022, NoiseModel 'endMatched', on BlueBEAR (2026-11-15).
%   v002 (2026-11-15): Acceptance mode calls checkV022Acceptance_v002 (serial vs parfor within 1e-6; v001 failed on BlueBEAR
%   on 4e-8 rounding differences). Nothing else changed.
%   Called by runV022HPC_v002.sh through `matlab -batch`; no GUI, no clicking. Also runs by hand in an sinteractive session.
%   runV022HPC_v002("Cook CTRL", "2026-11-16T07:00:00")         one dataset, pause cleanly before this time
%   runV022HPC_v002(["Hickman PLAC" "Hickman HALO" "Cook CTRL" "Cook ASD"], deadline)
%   runV022HPC_v002("FraserN48", deadline, Acceptance=true)       the small acceptance grid only (checkV022Acceptance_v002)
%   Settings (decisions 2026-11-14/15, TODO v009 N-1 to N-7): NoiseModel 'endMatched', Template 'production' (the v015
%   reporting basis), N_BETA 25, N_REPS 200, N_R7 10, RngSeed 1729, TrialSelection 'all', ChunkSize 64 with checkpoints, so a
%   job that pauses at the deadline (or is killed) resumes from its finished chunks when resubmitted unchanged.
%   Parallel pool as BEAR Apps recommends for MATLAB: parcluster('local') with JobStorageLocation = $TMPDIR (so concurrent jobs
%   cannot collide), sized to the cores Slurm gave the job (SLURM_CPUS_PER_TASK, else SLURM_TASKS_PER_NODE, else all cores).
%   The pool is opened once here and reused by everything below (CLAUDE.md: never restart a pool inside a called function).
%   Each dataset's own failure is caught, logged with its full report, and the next dataset still runs; the function errors at
%   the end if any dataset failed, so `matlab -batch` exits non-zero and sacct shows FAILED.   Fraser, D.S. (2026)  v002
arguments
    datasets (1,:) string
    deadline (1,1) string = ""                 % ISO time, e.g. from the job's end time; "" = no deadline
    opts.Acceptance (1,1) logical = false
    opts.ChunkSize (1,1) double {mustBeInteger, mustBePositive} = 64
    opts.Template (1,1) string {mustBeMember(opts.Template, ["production", "endAnchored"])} = "production"
end
here = fileparts(mfilename('fullpath'));  addpath(here, genpath(fullfile(here, 'functions')));
fprintf('runV022HPC_v002 on %s, MATLAB %s, start %s\n', strtrim(evalc('system(''hostname'');')), version, char(datetime('now')));
if opts.Acceptance
    openPool();  checkV022Acceptance_v002(Dataset=char(datasets(1)));  return
end
dl = NaT;  if deadline ~= "", dl = datetime(deadline, 'InputFormat', "yyyy-MM-dd'T'HH:mm:ss") - minutes(20); end   % 20 min margin for the save
fprintf('Deadline for starting a new chunk: %s\n', char(dl));
openPool();
failed = strings(0);
for d = datasets
    fprintf('\n===== %s: start %s =====\n', d, char(datetime('now')));
    try
        runLoopClosureFftnoise_v022(char(d), 'NoiseModel', 'endMatched', 'Template', char(opts.Template), 'TrialSelection', 'all', ...
            'N_BETA', 25, 'N_REPS', 200, 'N_R7', 10, 'RngSeed', 1729, 'UseParfor', true, 'ChunkSize', opts.ChunkSize, 'Deadline', dl);
        fprintf('===== %s: returned %s =====\n', d, char(datetime('now')));
    catch ME
        failed(end + 1) = d; %#ok<AGROW>
        fprintf(2, '===== %s: FAILED %s =====\n%s\n', d, char(datetime('now')), getReport(ME, 'extended', 'hyperlinks', 'off'));
    end
end
if ~isempty(failed), error('runV022HPC:failed', 'datasets failed: %s (see the log above)', strjoin(failed, ", ")); end
fprintf('\nAll requested datasets returned (a PAUSED message above means resubmit to continue). End %s\n', char(datetime('now')));
end

function openPool()
p = gcp('nocreate');
if ~isempty(p), fprintf('Reusing pool of %d workers\n', p.NumWorkers);  return, end
n = str2double(getenv('SLURM_CPUS_PER_TASK'));
if ~isfinite(n), n = str2double(getenv('SLURM_TASKS_PER_NODE')); end      % "72" or "72(x1)": str2double gives NaN for the latter
if ~isfinite(n), n = feature('numcores'); end
pc = parcluster('local');  t = getenv('TMPDIR');
if ~isempty(t), pc.JobStorageLocation = t; end
pc.NumWorkers = n;  parpool(pc, n);
fprintf('Opened pool of %d workers (JobStorageLocation %s)\n', n, pc.JobStorageLocation);
end
