function runFraserN48HPC_v002(templates)
%RUNFRASERN48HPC_V002  v001 with the template as an argument, so the two templates run as two parallel sessions (2026-11-03).
%   Use the one-press wrappers: runFraserN48HPC_production_v001 in one sinteractive session and
%   runFraserN48HPC_endAnchored_v001 in another. Each is about 6.5 h on 72 cores (the v019 Fraser run took 391 min for
%   2829 trials across two sessions), so each should finish inside one 11 h window. With no argument both run in turn, as v001.
%   The two runs write different output files, checkpoint folders and logs, so they cannot collide.
%   v001's header follows.
%RUNFRASERN48HPC_V001  The N=48 Fraser corpus on BlueBEAR at N_REPS = 200, under BOTH templates, for 12 h sessions (2026-11-03).
%   Open this file in the MATLAB GUI on the HPC (sinteractive) and press Run. No arguments, no setup.
%   Runs runLoopClosureFftnoise_v020 on dataset 'FraserN48' (importDB_fraser_v002: 2877 Ellipse trials, 48 participants):
%     1. Template = 'production'   (v015's windows: the reporting basis, Z12)  -> loopClosureResults_FraserN48_all_shaped_xu_tplProd_v020.mat
%     2. Template = 'endAnchored'  (v019's windows: the D30 sensitivity)        -> loopClosureResults_FraserN48_all_shaped_xu_v020.mat
%   Same hardware, settings and seed as the other datasets' v015 and v019 runs (N_BETA 25, N_REPS 200, N_R7 10, seed 1729).
%   Checkpoint every CHUNK trials, resume from them, clean stop (PAUSED, not an error) before any chunk that would end after
%   WALL_HOURS. If a session ends, press Run again in the next one: nothing finished is redone. Expect about two sessions
%   (Fraser dominated the 9 h five-dataset v019 run; this is two Fraser-sized runs).
%   Nothing existing is replaced: the 'Fraser' v015/v019 corpora are untouched. Log: src/logs/FraserN48_<stamp>.log.
%   Needs on RDS src/: runLoopClosureFftnoise_v020.m, noiseCharacterisation_fraserN48.mat, functions/importDB_fraser_v002.m,
%   and data/fraser/automatic/ with the ID401-403 folders.   Fraser, D.S. (2026)  v001

arguments
    templates (1,:) string {mustBeMember(templates, ["production" "endAnchored"])} = ["production" "endAnchored"]
end

%% CONFIG
MODEL      = "shaped_xu";
TEMPLATES  = templates;   % v002: from the argument (reporting basis first by default)
NWORKERS   = 72;
CHUNK      = 360;
WALL_HOURS = 11;
N_TRIALS   = 2877;   % importDB_fraser_v002 on 2026-11-03 (2880 Ellipse trials, 3 rejected with reasons)
RUN        = struct('N_BETA', 25, 'N_REPS', 200, 'N_R7', 10, 'RngSeed', 1729);   % as v015 and v019

%% Preflight (fail loud before hours are spent)
srcDir = fileparts(mfilename('fullpath'));  addpath(srcDir);  addpath(genpath(fullfile(srcDir, 'functions')));
need = ["runLoopClosureFftnoise_v020.m" "noiseCharacterisation_fraserN48.mat" "functions/importDB_fraser_v002.m" ...
        "functions/generateLoopClosureNoise_v003.m"];
for f = need, if ~isfile(fullfile(srcDir, f)), error('runFraserN48HPC_v002:missing', 'Missing %s', fullfile(srcDir, f)); end, end
if ~contains(fileread(fullfile(srcDir, 'runLoopClosureFftnoise_v020.m')), "case 'FraserN48'")
    error('runFraserN48HPC_v002:oldRunner', '%s', 'runLoopClosureFftnoise_v020.m has no FraserN48 case: copy the current file to RDS.');
end
trials = importDB_fraser_v002(Verbose=false);
if numel(trials) ~= N_TRIALS || numel(unique(string({trials.subjectID}))) ~= 48
    error('runFraserN48HPC_v002:corpus', 'importDB_fraser_v002 gives %d trials from %d participants here; expected %d from 48. Check data/fraser/automatic on RDS.', ...
        numel(trials), numel(unique(string({trials.subjectID}))), N_TRIALS);
end
bio = load(fullfile(srcDir, 'noiseCharacterisation_fraserN48.mat'), 'bioResults').bioResults;
if height(bio) ~= N_TRIALS, error('runFraserN48HPC_v002:noiseFile', 'noiseCharacterisation_fraserN48.mat has %d rows; expected %d.', height(bio), N_TRIALS); end
clear trials bio
DEADLINE = datetime('now') + hours(WALL_HOURS);
logDir = fullfile(srcDir, 'logs');  if ~isfolder(logDir), mkdir(logDir); end
diary(fullfile(logDir, ['FraserN48_' char(strjoin(TEMPLATES, '+')) '_' char(datetime('now', 'Format', 'yyyyMMdd_HHmmss')) '.log']));
cleanup = onCleanup(@() diary('off')); %#ok<NASGU>
[~, host] = system('hostname');
fprintf('FraserN48 HPC v002 [%s] | %s | host %s | MATLAB %s | cores %d | src %s\n', strjoin(TEMPLATES, ' + '), stamp_local(), strtrim(host), version, feature('numcores'), srcDir);
fprintf('Chunk %d trials | stop starting chunks that would end after %s (WALL_HOURS %g)\n', CHUNK, char(DEADLINE), WALL_HOURS);
for d = dir(fullfile(srcDir, 'loopClosureResults_FraserN48*'))', fprintf('  existing: %s\n', d.name); end

%% Run
args = namedargs2cell(RUN);
for tpl = TEMPLATES
    pool_local(min(NWORKERS, feature('numcores')));
    fprintf('\n=== FraserN48 | %s | template %s | %s ===\n', MODEL, tpl, stamp_local());
    r = runLoopClosureFftnoise_v020('FraserN48', 'TrialSelection', 'all', 'NoiseModel', char(MODEL), 'UseParfor', true, ...
        args{:}, 'Template', char(tpl), 'ChunkSize', CHUNK, 'Deadline', DEADLINE);
    if isempty(r)
        fprintf('\nFRASERN48 PAUSED in template %s (not an error) | %s | press Run again in the next session\n', tpl, stamp_local());
        return
    end
end

%% Completion check
tags = struct('production', "_tplProd", 'endAnchored', "");
for tpl = TEMPLATES
    f = fullfile(srcDir, sprintf('loopClosureResults_FraserN48_all_%s%s_v020.mat', MODEL, tags.(tpl)));
    if ~isfile(f), error('runFraserN48HPC_v002:noOutput', 'Expected output missing: %s', f); end
    w = whos('-file', f);
    if ~all(ismember({'results', 'config'}, {w.name})), error('runFraserN48HPC_v002:noResults', '%s lacks results/config', f); end
    S = load(f, 'config', 'results');
    if ~strcmp(S.config.TEMPLATE, tpl) || S.config.N_REPS ~= 200 || numel(S.results) ~= N_TRIALS
        error('runFraserN48HPC_v002:config', '%s: template %s, N_REPS %d, %d trials; expected %s, 200, %d', f, ...
            S.config.TEMPLATE, S.config.N_REPS, numel(S.results), tpl, N_TRIALS);
    end
    nErr = sum([S.results.nEbrErr]);
    fprintf('OK  %s  (%.1f MB; regressDataEBR errors counted: %d in %d trials)\n', f, sum([w.bytes]) / 1e6, nErr, sum([S.results.nEbrErr] > 0));
end
fprintf('FraserN48 finished | %s\n=== COMPLETED OK ===\n', stamp_local());
end

function p = pool_local(nw)
    p = gcp('nocreate');
    if ~isempty(p)
        try, alive = p.Connected && p.NumWorkers > 0; catch, alive = false; end
        if ~alive, try, delete(p); catch, end, p = []; end
    end
    if isempty(p), p = parpool('local', nw, 'IdleTimeout', Inf); elseif ~isinf(p.IdleTimeout), p.IdleTimeout = Inf; end
    if isempty(p) || ~p.Connected, error('runFraserN48HPC_v002:noPool', 'Could not establish a parallel pool'); end
end

function s = stamp_local()
    s = string(datetime('now', 'TimeZone', 'Europe/London', 'Format', 'yyyy-MM-dd HH:mm:ss z'));
end
