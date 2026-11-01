function runEndAnchorHeadlineHPC_v001()
%RUNENDANCHORHEADLINEHPC_V001  Sol's closure: the five core datasets at N_REPS = 200 with the END-ANCHORED template (shaped_xu), for 12 h sessions.
%   Open this file in the MATLAB GUI on the HPC (sinteractive) and press Run. No arguments, no setup.
%   Runs runLoopClosureFftnoise_v019 (= v017 + one final template window ending at the last sample, so the zero-residual tail cannot occur).
%   Checkpoint every CHUNK trials, resume from them, clean stop (PAUSED, not an error) before any chunk that would end after WALL_HOURS.
%   If a session ends, press Run again in the next one: nothing finished is redone. Expected about 9 h in total, Fraser dominating.
%   Outputs: src/loopClosureResults_<dataset>_all_shaped_xu_v019.mat. Nothing existing is replaced (v015 corpora untouched). Log: src/logs/EndAnchor_<stamp>.log.
%   Then run analyseEndAnchorHeadline_v001 on the Mac.   Fraser, D.S. (2026)  v001

%% CONFIG
MODEL      = "shaped_xu";
DATASETS   = ["Cook CTRL" "Cook ASD" "Hickman PLAC" "Hickman HALO" "Fraser"];   % cheapest first, Fraser last
NWORKERS   = 72;
CHUNK      = 360;
WALL_HOURS = 11;
RUN        = struct('N_BETA', 25, 'N_REPS', 200, 'N_R7', 10, 'RngSeed', 1729);   % as v015

%% Preflight (fail loud before hours are spent)
srcDir = fileparts(mfilename('fullpath'));  addpath(srcDir);
need = ["runLoopClosureFftnoise_v017.m" "makeRunLoopClosureV019_v001.m" "noiseCharacterisation_cook.mat" "noiseCharacterisation_cookASD.mat" ...
        "noiseCharacterisation_hickmanPLAC.mat" "noiseCharacterisation_hickmanHALO.mat" "noiseCharacterisation_fraser.mat"];
for f = need, if ~isfile(fullfile(srcDir, f)), error('runEndAnchor:missing', 'Missing %s', fullfile(srcDir, f)); end, end
if isempty(dir(fullfile(srcDir, 'functions', 'generateLoopClosureNoise_v003.m')))
    error('runEndAnchor:noGenerator', 'Missing %s', fullfile(srcDir, 'functions', 'generateLoopClosureNoise_v003.m'));
end
if ~isfile(fullfile(srcDir, 'runLoopClosureFftnoise_v019.m')), makeRunLoopClosureV019_v001('SrcDir', srcDir); end
if ~contains(fileread(fullfile(srcDir, 'runLoopClosureFftnoise_v019.m')), 'END_ANCHORED_TEMPLATE'), error('runEndAnchor:notAnchored', 'v019 lacks the end-anchored template'); end
DEADLINE = datetime('now') + hours(WALL_HOURS);
logDir = fullfile(srcDir, 'logs');  if ~isfolder(logDir), mkdir(logDir); end
diary(fullfile(logDir, ['EndAnchor_' char(datetime('now', 'Format', 'yyyyMMdd_HHmmss')) '.log']));
cleanup = onCleanup(@() diary('off')); %#ok<NASGU>
[~, host] = system('hostname');
fprintf('End-anchored headline v001 | %s | host %s | MATLAB %s | cores %d | src %s\n', stamp_local(), strtrim(host), version, feature('numcores'), srcDir);
fprintf('Chunk %d trials | stop starting chunks that would end after %s (WALL_HOURS %g)\n', CHUNK, char(DEADLINE), WALL_HOURS);
for d = dir(fullfile(srcDir, 'loopClosureResults_*_v019*'))', fprintf('  existing: %s\n', d.name); end

%% Run
args = namedargs2cell(RUN);
for ds = DATASETS
    pool_local(min(NWORKERS, feature('numcores')));
    fprintf('\n=== %s | %s | %s ===\n', ds, MODEL, stamp_local());
    r = runLoopClosureFftnoise_v019(char(ds), 'TrialSelection', 'all', 'NoiseModel', char(MODEL), 'UseParfor', true, ...
        args{:}, 'ChunkSize', CHUNK, 'Deadline', DEADLINE);
    if isempty(r)
        fprintf('\nEND-ANCHOR PAUSED in %s (not an error) | %s | press Run again in the next session\n', ds, stamp_local());
        return
    end
end

%% Completion check
for ds = DATASETS
    f = fullfile(srcDir, sprintf('loopClosureResults_%s_all_%s_v019.mat', strrep(ds, ' ', '_'), MODEL));
    if ~isfile(f), error('runEndAnchor:noOutput', 'Expected output missing: %s', f); end
    w = whos('-file', f);
    if ~all(ismember({'results', 'config'}, {w.name})), error('runEndAnchor:noResults', '%s lacks results/config', f); end
    S = load(f, 'config');
    if ~isfield(S.config, 'END_ANCHORED_TEMPLATE') || ~S.config.END_ANCHORED_TEMPLATE || S.config.N_REPS ~= 200, error('runEndAnchor:config', '%s config not end-anchored N_REPS=200', f); end
    fprintf('OK  %s  (%.1f MB)\n', f, sum([w.bytes]) / 1e6);
end
fprintf('End-anchored headline finished | %s\n=== COMPLETED OK ===\n', stamp_local());
end

function p = pool_local(nw)
    p = gcp('nocreate');
    if ~isempty(p)
        try, alive = p.Connected && p.NumWorkers > 0; catch, alive = false; end
        if ~alive, try, delete(p); catch, end, p = []; end
    end
    if isempty(p), p = parpool('local', nw, 'IdleTimeout', Inf); elseif ~isinf(p.IdleTimeout), p.IdleTimeout = Inf; end
    if isempty(p) || ~p.Connected, error('runEndAnchor:noPool', 'Could not establish a parallel pool'); end
end

function s = stamp_local()
    s = string(datetime('now', 'TimeZone', 'Europe/London', 'Format', 'yyyy-MM-dd HH:mm:ss z'));
end
