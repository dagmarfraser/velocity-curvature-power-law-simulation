function runV023HPC_v001(datasets, deadline, opts)
%RUNV023HPC_V001  Batch (sbatch) launcher for runLoopClosureFftnoise_v023 on BlueBEAR (2026-11-15). Supersedes runV022HPC_v002.
%   Called by runV023HPC_v001.sh through `matlab -batch`; also runs by hand. New in this version (Dagmar, 2026-11-15 09:10):
%     NoiseModel  "endMatched" (primary) or "shapedXuSeg" (the modernised shapeXu sensitivity arm: same segments, length,
%                 windows, two inversions and seeds; only the noise construction differs; makeRunLoopClosureV023_v001)
%     RngSeed / NReps  for the seed-stability reruns (README_SeedPractice rule 7). Any non-default seed or N_REPS adds
%                 "_seed<S>_r<N>" to the output name (through the runner's OutTag), so a rerun never collides with, skips
%                 or overwrites the primary file, and gets its own checkpoint folder.
%   runV023HPC_v001(["Hickman PLAC" "Hickman HALO" "Cook CTRL" "Cook ASD" "Zarandi" "Dhieb"], "2026-11-17T09:00:00")
%   runV023HPC_v001("FraserN48", deadline, NoiseModel="shapedXuSeg")
%   runV023HPC_v001("Cook CTRL", deadline, RngSeed=271828, NReps=50)     -> ..._endMatched_tplProd_seed271828_r50_v023.mat
%   runV023HPC_v001("Cook CTRL", "", Acceptance=true)                    -> checkV023Acceptance_v001
%   Fixed settings: Template 'production', N_BETA 25, N_R7 10, TrialSelection 'all', ChunkSize 64, pause 20 min before the
%   deadline (the job's end time), pool opened once the BEAR way (parcluster('local'), JobStorageLocation = $TMPDIR, sized to
%   the job's cores). A dataset's failure is logged and the rest still run; the function errors at the end if any failed.
arguments
    datasets (1,:) string
    deadline (1,1) string = ""
    opts.NoiseModel (1,1) string {mustBeMember(opts.NoiseModel, ["endMatched", "shapedXuSeg"])} = "endMatched"
    opts.RngSeed (1,1) double {mustBeInteger, mustBeNonnegative} = 1729
    opts.NReps (1,1) double {mustBeInteger, mustBePositive, mustBeLessThan(opts.NReps, 1000)} = 200
    opts.Acceptance (1,1) logical = false
    opts.ChunkSize (1,1) double {mustBeInteger, mustBePositive} = 64
end
here = fileparts(mfilename('fullpath'));  addpath(here, genpath(fullfile(here, 'functions')));
[~, host] = system('hostname');
fprintf('runV023HPC_v001 on %s, MATLAB %s, start %s\n', strtrim(host), version, char(datetime('now')));
openPool();
if opts.Acceptance, checkV023Acceptance_v001(Dataset=char(datasets(1)));  return, end
outTag = "";  if opts.RngSeed ~= 1729 || opts.NReps ~= 200, outTag = "_seed" + opts.RngSeed + "_r" + opts.NReps; end
dl = NaT;  if deadline ~= "", dl = datetime(deadline, 'InputFormat', "yyyy-MM-dd'T'HH:mm:ss") - minutes(20); end
fprintf('Model %s, seed %d, N_REPS %d, output tag "%s", deadline for starting a chunk %s\n', opts.NoiseModel, opts.RngSeed, opts.NReps, outTag, char(dl));
failed = strings(0);
for d = datasets
    fprintf('\n===== %s: start %s =====\n', d, char(datetime('now')));
    try
        runLoopClosureFftnoise_v023(char(d), 'NoiseModel', char(opts.NoiseModel), 'Template', 'production', 'TrialSelection', 'all', ...
            'N_BETA', 25, 'N_REPS', opts.NReps, 'N_R7', 10, 'RngSeed', opts.RngSeed, 'UseParfor', true, 'ChunkSize', opts.ChunkSize, ...
            'Deadline', dl, 'OutTag', char(outTag));
        fprintf('===== %s: returned %s =====\n', d, char(datetime('now')));
    catch ME
        failed(end + 1) = d; %#ok<AGROW>
        fprintf(2, '===== %s: FAILED %s =====\n%s\n', d, char(datetime('now')), getReport(ME, 'extended', 'hyperlinks', 'off'));
    end
end
if ~isempty(failed), error('runV023HPC:failed', 'datasets failed: %s (see the log above)', strjoin(failed, ", ")); end
fprintf('\nAll requested datasets returned (a PAUSED message above means resubmit to continue). End %s\n', char(datetime('now')));
end

function openPool()
p = gcp('nocreate');
if ~isempty(p), fprintf('Reusing pool of %d workers\n', p.NumWorkers);  return, end
n = str2double(getenv('SLURM_CPUS_PER_TASK'));
if ~isfinite(n), n = str2double(getenv('SLURM_TASKS_PER_NODE')); end
if ~isfinite(n), n = feature('numcores'); end
pc = parcluster('local');  t = getenv('TMPDIR');
if ~isempty(t), pc.JobStorageLocation = t; end
pc.NumWorkers = n;  parpool(pc, n);
fprintf('Opened pool of %d workers (JobStorageLocation %s)\n', n, pc.JobStorageLocation);
end
