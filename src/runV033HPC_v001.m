function runV033HPC_v001(datasets, deadline, opts)
%RUNV028HPC_V003  Batch (sbatch) launcher for runLoopClosureFftnoise_v033 on BlueBEAR (2026-11-15). Supersedes runV024HPC_v001.
%   v001 (2026-11-18): runs runLoopClosureFftnoise_v033, which also accepts --harmonics 0 (template strips no harmonic; constant and trend only). Otherwise identical to runV028HPC_v003. Output names: ..._nH<K>_v033.mat.
%   Decisions (Dagmar, 2026-11-15): end-anchored template, 4 harmonics; end-matched noise with the detectability segment rule
%   as primary; the fixed 1% budget as a sensitivity arm; shapedXuSeg as the noise sensitivity arm; seed-stability reruns.
%   v003 (2026-11-16): adds option Harmonics ("4" default, or "fc" = floor(BW cutoff / f0) harmonics per trial, or an integer) for the harmonic-count sensitivity (TODO N-22); non-default values add _nH<value> to the output names (runner v025+), so nothing is overwritten. Nothing else changed.
%   v002 (2026-11-15 19:25): fixes the failure printout that crashed AFTER each dataset's results were saved (v001 used extractBefore(..., 60)
%   on failure reasons shorter than 60 characters, so every finished dataset was reported FAILED and the job exited 1; results were intact).
%   runV033HPC_v001(datasets, deadline)                                     primary: endMatched, SegmentRule "detect"
%   runV033HPC_v001(datasets, deadline, Rule="budget")                      1% budget arm (SegmentBudget 0.01)
%   runV033HPC_v001(datasets, deadline, NoiseModel="shapedXuSeg")           shapeXu arm (detect)
%   runV033HPC_v001(datasets, deadline, RngSeed=271828, NReps=50)           seed rerun (adds "_seed<S>_r<N>")
%   runV033HPC_v001("Cook CTRL", "", Acceptance=true)                       checkV033EqualsV028_v001
%   Fixed: Template 'endAnchored', NHarmonics 4, N_BETA 25, N_R7 10, TrialSelection 'all', ChunkSize 72 (one trial per worker on
%   a 72-core node), pause 20 min before the deadline, pool opened once the BEAR way (parcluster('local'), JobStorageLocation =
%   $TMPDIR, sized to the job's cores). Output names carry the runner's tags (_segdetect or _bud1, _nH..., model) plus the
%   seed tag, so no run collides with or overwrites another. A dataset's failure is logged; the function errors at the end.
arguments
    datasets (1,:) string
    deadline (1,1) string = ""
    opts.NoiseModel (1,1) string {mustBeMember(opts.NoiseModel, ["endMatched", "shapedXuSeg"])} = "endMatched"
    opts.Rule (1,1) string {mustBeMember(opts.Rule, ["detect", "budget"])} = "detect"
    opts.Budget (1,1) double {mustBePositive} = 0.01
    opts.RngSeed (1,1) double {mustBeInteger, mustBeNonnegative} = 1729
    opts.NReps (1,1) double {mustBeInteger, mustBePositive, mustBeLessThan(opts.NReps, 1000)} = 200
    opts.Harmonics (1,1) string = "4"
    opts.Acceptance (1,1) logical = false
    opts.ChunkSize (1,1) double {mustBeInteger, mustBePositive} = 72
end
here = fileparts(mfilename('fullpath'));  addpath(here, genpath(fullfile(here, 'functions')));
[~, host] = system('hostname');
fprintf('runV033HPC_v001 on %s, MATLAB %s, start %s\n', strtrim(host), version, char(datetime('now')));
openPool();
if opts.Acceptance, checkV033EqualsV028_v001(Dataset=char(datasets(1)));  return, end
outTag = "";  if opts.RngSeed ~= 1729 || opts.NReps ~= 200, outTag = "_seed" + opts.RngSeed + "_r" + opts.NReps; end
nh = 4;  if opts.Harmonics == "fc", nh = 'fc'; elseif opts.Harmonics ~= "4", nh = str2double(opts.Harmonics); end
seg = {'SegmentRule', 'detect'};  if opts.Rule == "budget", seg = {'SegmentBudget', opts.Budget}; end
dl = NaT;  if deadline ~= "", dl = datetime(deadline, 'InputFormat', "yyyy-MM-dd'T'HH:mm:ss") - minutes(20); end
fprintf('Harmonics %s. Model %s, rule %s, seed %d, N_REPS %d, tag "%s", deadline for starting a chunk %s\n', opts.Harmonics, opts.NoiseModel, opts.Rule, opts.RngSeed, opts.NReps, outTag, char(dl));
failed = strings(0);
for d = datasets
    fprintf('\n===== %s: start %s =====\n', d, char(datetime('now')));
    try
        R = runLoopClosureFftnoise_v033(char(d), 'NoiseModel', char(opts.NoiseModel), 'Template', 'endAnchored', 'NHarmonics', nh, seg{:}, ...
            'TrialSelection', 'all', 'N_BETA', 25, 'N_REPS', opts.NReps, 'N_R7', 10, 'RngSeed', opts.RngSeed, 'UseParfor', true, ...
            'ChunkSize', opts.ChunkSize, 'Deadline', dl, 'OutTag', char(outTag));
        if ~isempty(R)
            st = [R.emStatus];  fprintf('===== %s: %d trials, %d ran, failures by reason:\n', d, numel(R), nnz(st == "ok"));
            [u, ~, k] = unique(extractBefore(pad(st(st ~= "ok"), 60), 60));   % v002: padded, so reasons shorter than 60 characters no longer crash the printout  for i = 1:numel(u), fprintf('   %4d  %s\n', nnz(k == i), u(i)); end
        end
        fprintf('===== %s: returned %s =====\n', d, char(datetime('now')));
    catch ME
        failed(end + 1) = d; %#ok<AGROW>
        fprintf(2, '===== %s: FAILED %s =====\n%s\n', d, char(datetime('now')), getReport(ME, 'extended', 'hyperlinks', 'off'));
    end
end
if ~isempty(failed), error('runV028HPC:failed', 'datasets failed: %s (see the log above)', strjoin(failed, ", ")); end
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
