% runVGFSEMCentroidPipeline_v001.m  Orchestrates the two VGF-SEM scripts
% in dependency order, on BlueBEAR via sinteractive.
%
% Runs, in this order:
%   1. computePerCoordinateSEM_VGF_v001()  -> perCoordinateSEM_VGF_v001.mat
%   2. checkVGFSEMCentroid_v001            -> console table (captured below)
%
% All console output is captured via diary() to
% vgfSEMCentroidRun_log_<timestamp>.txt in this script's own directory.
%
% BUG FOUND AND FIXED, 2026-09-13 (first real run, RDS/BlueBEAR): the
% original version called Stage 2 via run(fullfile(srcDir,
% 'checkVGFSEMCentroid_v001.m')) directly from this script's own
% top-level body. Because this file is itself a script, that call shares
% the base workspace -- so checkVGFSEMCentroid_v001.m's own clearvars (at
% its own top, correct for its standalone use) wiped every variable here,
% including the onCleanup object guarding the diary. Clearing an
% onCleanup object fires its cleanup function immediately, so diary('off')
% fired the instant Stage 2 started -- before Stage 2 printed a single
% line. Confirmed from the actual log file: it cut off mid-Stage-2, the
% entire results table only ever reached the screen, never disk.
%
% Fix: Stage 2 is now invoked through the local function runIsolated()
% below. Local functions get their own workspace even when called from a
% script's top-level body, so a clearvars inside the target script can
% only clear that local function's own (empty) workspace, not this
% script's. Confirmed this is the actual mechanism, not assumed: Stage 1
% never had this problem because computePerCoordinateSEM_VGF_v001 is a
% FUNCTION file (arguments block, own workspace already), only Stage 2
% (a script calling a script) was ever at risk.
%
% Fail Loud: every precondition is checked and reported BEFORE any
% compute starts, rather than letting Stage 2 fail after Stage 1's
% (slower) database query has already run. A try/catch around both
% stages guarantees diary('off') runs on any failure path too, rather
% than depending on an onCleanup object's survival (the same class of
% fragility that caused the bug above).
%
% USAGE (BlueBEAR, sinteractive session, from the RDS src/ directory):
%   >> runVGFSEMCentroidPipeline_v001
%
% Fraser, D.S. (2026)  v002 -- diary-truncation fix

clearvars
srcDir = fileparts(mfilename('fullpath'));
if isempty(srcDir), srcDir = pwd; end
cd(srcDir);
addpath(genpath(fullfile(srcDir, 'functions')));

fprintf('=== VGF SEM/CENTROID PIPELINE RUNNER ===\n');
fprintf('Working directory: %s\n\n', srcDir);

%% ---------------------------------------------------------------------
%% Precondition checks -- Fail Loud, all checked before anything runs
%% ---------------------------------------------------------------------
required = {
    fullfile(srcDir, 'computePerCoordinateSEM_VGF_v001.m'), 'VGF per-coordinate script';
    fullfile(srcDir, 'checkVGFSEMCentroid_v001.m'),          'VGF centroid-lookup script';
    fullfile(srcDir, 'functions', 'semAdequacyThreshold_v001.m'), 'SEM threshold helper';
    fullfile(srcDir, 'perCoordinateSEM_v2_001.mat'),         'canonical beta per-coordinate table';
    fullfile(srcDir, '..', 'results', 'powerlaw_debug_v058.db'), 'v058 simulation database';
};

fprintf('--- Precondition check ---\n');
missing = false(size(required,1),1);
for i = 1:size(required,1)
    ok = isfile(required{i,1});
    missing(i) = ~ok;
    fprintf('  [%s] %-40s %s\n', string(ok), required{i,2}, required{i,1});
end

if any(missing)
    fprintf('\n');
    error('runVGFSEMCentroidPipeline:MissingPrereqs', '%s', sprintf([...
        '%d precondition(s) missing (see [false] rows above). ' ...
        'Nothing was run. Sync the missing file(s) onto RDS before retrying.'], ...
        sum(missing)));
end
fprintf('  All preconditions present.\n\n');

%% ---------------------------------------------------------------------
%% Diary: capture everything from here on to a timestamped log file
%% ---------------------------------------------------------------------
logFile = fullfile(srcDir, sprintf('vgfSEMCentroidRun_log_%s.txt', ...
    datestr(now, 'yyyymmdd_HHMMSS'))); %#ok<TNOW1,DATST>
diary(logFile);
diary on

fprintf('Log file: %s\n', logFile);
fprintf('Run started: %s\n\n', datestr(now)); %#ok<TNOW1,DATST>

try
    %% -------------------------------------------------------------
    %% Stage 1: VGF per-coordinate SEM/bias (a FUNCTION -- own
    %% workspace already, no isolation trick needed)
    %% -------------------------------------------------------------
    fprintf('=== STAGE 1: computePerCoordinateSEM_VGF_v001 ===\n');
    stage1Start = tic;
    computePerCoordinateSEM_VGF_v001();
    fprintf('Stage 1 wall time: %.1f s\n\n', toc(stage1Start));

    vgfMat = fullfile(srcDir, 'perCoordinateSEM_VGF_v001.mat');
    if ~isfile(vgfMat)
        error('runVGFSEMCentroidPipeline:Stage1Failed', '%s', sprintf(...
            'Stage 1 completed without error but %s was not created. Stopping before Stage 2.', ...
            vgfMat));
    end

    %% -------------------------------------------------------------
    %% Stage 2: centroid lookup, gated by beta SEM adequacy.
    %% checkVGFSEMCentroid_v001 is a SCRIPT (clearvars at its own top,
    %% by design, for standalone use) -- routed through runIsolated()
    %% below so its clearvars cannot reach this workspace. See header
    %% comment for the bug this fixes.
    %% -------------------------------------------------------------
    fprintf('=== STAGE 2: checkVGFSEMCentroid_v001 ===\n');
    stage2Start = tic;
    runIsolated(fullfile(srcDir, 'checkVGFSEMCentroid_v001.m'));
    fprintf('Stage 2 wall time: %.1f s\n\n', toc(stage2Start));

    fprintf('=== PIPELINE COMPLETE ===\n');
    fprintf('Run finished: %s\n', datestr(now)); %#ok<TNOW1,DATST>
    fprintf('Outputs on disk: %s, %s\n', vgfMat, logFile);

    diary off

catch ME
    fprintf('\n=== PIPELINE FAILED ===\n%s\n', getReport(ME, 'extended'));
    diary off
    rethrow(ME);
end

%% ------------------------------------------------------------------
function runIsolated(scriptPath)
% Executes scriptPath inside THIS function's own workspace. A clearvars
% at the top of scriptPath clears this function's (otherwise-empty)
% workspace only -- it cannot reach the calling script's workspace, so
% the diary guard and any other caller-side variables survive.
    run(scriptPath);
end
