function results = runToolLoopClosure_v002(bidsRoot, options)
% RUNTOOLLOOPCLOSURE_V002  Batch Tier 1 processing across a Motion-BIDS tree.
%
% Discovers every *_motion.tsv file under bidsRoot (sub-*/motion/), reads
% each trial's SamplingFrequency from its companion *_motion.json sidecar,
% and calls processTrialLoopClosure_v001 once per trial. Function/
% dispatcher split (plain per-trial worker, for/parfor chosen by
% UseParfor) rather than parfor baked in: MATLAB Online does not support
% the local parfor profile without a separately-configured Cloud Center
% Cluster, so this tool's own default is UseParfor=false, unlike
% runLoopClosureFftnoise_v009's HPC-oriented default of relying on an
% explicit flag from a BlueBEAR caller. Same worker function either way --
% only the loop construct differs.
%
% v002 (2026-09-20): exposes the session label already present in every
%   BIDS filename this tool's own ingest writes
%   (convertFraserToMotionBIDS_v001.m's "..._ses-<sesLabel>_..." stem) but
%   previously discarded by this file's own parsing regex, which only
%   captured sub- and run-. Needed so a within-subject paired design
%   (e.g. two drug sessions per subject) can be identified downstream by
%   tier3DesignPrecision_v001.m without re-deriving it from the raw file
%   path. Nothing else changed: trial discovery, sidecar reading, the
%   fs/unit contract, and Tier 1 dispatch (function/dispatcher split,
%   UseParfor) are byte-identical to v001 -- diff this file against v001
%   to confirm.
%
% UNIT ASSUMPTION: motion.tsv values are already millimetres, per this
% tool's own ingest convention (e.g. convertFraserToMotionBIDS_v001.m's
% PX_PER_MM conversion). This function does not re-convert units, and
% does not verify the assumption -- a device-specific converter that
% wrote raw pixels into motion.tsv would silently produce wrong results.
%
% FAILURE HANDLING: one trial's failure (too short, IRASA non-finite,
% regression failure on all pipelines) does not abort the batch. That
% trial's row in the results table is marked ok=false with the caught
% error's message; the value fields are left empty/NaN, not silently
% zero-filled or dropped from the table -- fail loud, but keep going.
%
% SYNTAX:
%   results = runToolLoopClosure_v001(bidsRoot)
%   results = runToolLoopClosure_v001(bidsRoot, UseParfor=true, N_REPS=10)
%
% INPUTS:
%   bidsRoot - Path to a Motion-BIDS dataset root (the folder containing
%              dataset_description.json and sub-<label>/ folders), e.g.
%              tool/examples/fraser_sub-103.
%   Name-value options:
%     UseParfor       false (default here; true is fine on BlueBEAR/desktop
%                     MATLAB with a local pool, not MATLAB Online without
%                     a Cloud Center Cluster)
%     ReportPipeline  "SG-IRLS" (passed through to processTrialLoopClosure_v001)
%     N_BETA, BetaRange, N_REPS, EdgeClip, MonoSlopeTol, MonoSmoothWidth,
%     MonoMinSegWidth, MonoTrimToleranceK  -- all passed through to
%     processTrialLoopClosure_v001; see that function's own help for
%     defaults. (Renamed 2026-08-15 from MonoSmoothWindow/MonoMinSegLength
%     to match processTrialLoopClosure_v001's own gate swap to
%     findMonotonicSegments -- this file's passthrough was previously
%     broken, still using the old names; fixed same date.)
%
% OUTPUT:
%   results - Nx1 struct array, one entry per discovered trial:
%     .motionFile   - full path to the source *_motion.tsv
%     .subLabel, .runLabel  - parsed from the BIDS filename entities
%     .sesLabel     - parsed from the BIDS filename's ses- entity (v002);
%                     empty string "" if the filename carries no ses-
%                     entity at all, not an error -- ses- is optional in
%                     Motion-BIDS and plenty of valid datasets omit it
%     .fs           - sampling frequency read from the sidecar JSON
%     .ok           - logical, whether processTrialLoopClosure_v001 succeeded
%     .errorMessage - char, the caught error's message if ok==false, else ''
%     .trial        - the full processTrialLoopClosure_v001 output struct
%                     if ok==true, else []
%
% See also: processTrialLoopClosure_v001, convertFraserToMotionBIDS_v001
%
% Fraser, D.S. (2026)

arguments
    bidsRoot (1,1) string {mustBeFolder}
    options.UseParfor (1,1) logical = false
    options.ReportPipeline (1,1) string = "SG-IRLS"
    options.N_BETA (1,1) double = 25
    options.BetaRange (1,2) double = [0 0.75]
    options.N_REPS (1,1) double = 20
    options.EdgeClip (1,1) double = 20
    options.MonoSlopeTol (1,1) double = 0.05
    options.MonoSmoothWidth (1,1) double = 0.0625
    options.MonoMinSegWidth (1,1) double = 0.0625
    options.MonoTrimToleranceK (1,1) double = 1.0
end

%% Discover trials
motionFiles = dir(fullfile(bidsRoot, 'sub-*', 'motion', '*_motion.tsv'));
if isempty(motionFiles)
    error('runToolLoopClosure_v001:NoTrialsFound', '%s', ...
        sprintf('No *_motion.tsv files found under %s/sub-*/motion/. Check bidsRoot.', bidsRoot));
end
nTrials = numel(motionFiles);
fprintf('runToolLoopClosure_v001: found %d trials under %s\n', nTrials, bidsRoot);

filePaths = cell(nTrials, 1);
subLabels = strings(nTrials, 1);
sesLabels = strings(nTrials, 1);
runLabels = strings(nTrials, 1);
fsVec     = NaN(nTrials, 1);

for i = 1:nTrials
    fp = fullfile(motionFiles(i).folder, motionFiles(i).name);
    filePaths{i} = fp;

    tok = regexp(motionFiles(i).name, 'sub-([^_]+).*run-([^_]+)_motion\.tsv', 'tokens', 'once');
    if ~isempty(tok)
        subLabels(i) = tok{1};
        runLabels(i) = tok{2};
    end
    sesTok = regexp(motionFiles(i).name, 'ses-([^_]+)', 'tokens', 'once');
    if ~isempty(sesTok)
        sesLabels(i) = sesTok{1};
    end

    jsonPath = strrep(fp, '_motion.tsv', '_motion.json');
    if ~isfile(jsonPath)
        error('runToolLoopClosure_v001:NoSidecar', '%s', ...
            sprintf('No sidecar JSON found for %s -- cannot read SamplingFrequency. Fix the dataset, do not guess fs.', fp));
    end
    sidecar = jsondecode(fileread(jsonPath));
    if ~isfield(sidecar, 'SamplingFrequency')
        error('runToolLoopClosure_v001:NoSamplingFrequency', '%s', ...
            sprintf('%s has no SamplingFrequency field.', jsonPath));
    end
    fsVec(i) = sidecar.SamplingFrequency;
end

%% Process, function/dispatcher split
resultsCells = cell(nTrials, 1);

if options.UseParfor
    parfor i = 1:nTrials
        resultsCells{i} = processOneFile_local(filePaths{i}, subLabels(i), sesLabels(i), runLabels(i), fsVec(i), options);
    end
else
    tRun = tic;
    for i = 1:nTrials
        resultsCells{i} = processOneFile_local(filePaths{i}, subLabels(i), sesLabels(i), runLabels(i), fsVec(i), options);
        fprintf('[%3d/%3d] sub-%s run-%s  ok=%d  (%.0fs elapsed)\n', ...
            i, nTrials, subLabels(i), runLabels(i), resultsCells{i}.ok, toc(tRun));
    end
end

results = [resultsCells{:}]';

nOk = sum([results.ok]);
fprintf('\nDone: %d/%d trials processed successfully.\n', nOk, nTrials);
if nOk < nTrials
    fprintf('Failed trials:\n');
    for i = 1:nTrials
        if ~results(i).ok
            fprintf('  sub-%s run-%s: %s\n', results(i).subLabel, results(i).runLabel, results(i).errorMessage);
        end
    end
end

end

% --------------------------------------------------------------------------
function r = processOneFile_local(fp, subLabel, sesLabel, runLabel, fs, options)
    r.motionFile = fp;
    r.subLabel   = subLabel;
    r.sesLabel   = sesLabel;
    r.runLabel   = runLabel;
    r.fs         = fs;
    r.ok         = false;
    r.errorMessage = '';
    r.trial      = [];

    try
        xy = readmatrix(fp, 'FileType', 'text');
        if size(xy, 2) < 2
            error('runToolLoopClosure_v001:BadMotionFile', '%s', ...
                sprintf('%s has fewer than 2 columns -- expected x,y.', fp));
        end
        trialResult = processTrialLoopClosure_v001(xy(:,1), xy(:,2), fs, ...
            'ReportPipeline', options.ReportPipeline, ...
            'N_BETA', options.N_BETA, 'BetaRange', options.BetaRange, ...
            'N_REPS', options.N_REPS, 'EdgeClip', options.EdgeClip, ...
            'MonoSlopeTol', options.MonoSlopeTol, ...
            'MonoSmoothWidth', options.MonoSmoothWidth, ...
            'MonoMinSegWidth', options.MonoMinSegWidth, ...
            'MonoTrimToleranceK', options.MonoTrimToleranceK);
        r.trial = trialResult;
        r.ok = true;
    catch ME
        r.errorMessage = ME.message;
    end
end
