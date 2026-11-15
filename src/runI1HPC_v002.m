function runI1HPC_v002()
%RUNI1HPC_V002  v001 for the N=48 copy (2026-11-04): calls loopClosureLOPOBaselines_v002 (the Finding #243 anchor is exact for the four unchanged datasets, and Fraser N=48 is reported against N=45); writes results/lopoBaselines_v002.mat.
%RUNI1HPC_V001  D-I1 on the HPC: the out-of-sample leave-one-pipeline-out figure against four baselines.
%   Open this file in the MATLAB GUI on the HPC (sinteractive) and press Run. No arguments, no pool, a few minutes.
%   Reads the five validity-domain v015 corpora in src/, checks itself against Finding #243 (20,514 cells, 84.53% within MDC),
%   then scores the framework against a naive mean, a participant-cross-validated linear fit, a mismatched map and a
%   wrong-colour map on the same cells. Writes results/lopoBaselines_v001.mat. Log: src/logs/I1_<stamp>.log.
%   Safe to run beside the other jobs: it is light and writes only its own files.   Fraser, D.S. (2026)  v001
srcDir = fileparts(mfilename('fullpath'));  addpath(srcDir);
ds = ["Cook_CTRL" "Cook_ASD" "Hickman_PLAC" "Hickman_HALO" "Fraser"];
if ~isfile(fullfile(srcDir, 'loopClosureLOPOBaselines_v002.m')), error('runI1HPC_v002:missing', 'Missing loopClosureLOPOBaselines_v002.m in %s', srcDir); end
for d = ds
    f = fullfile(srcDir, sprintf('loopClosureResults_%s_all_shaped_xu_v015.mat', d));
    if ~isfile(f), error('runI1HPC_v002:noCorpus', 'Missing v015 corpus %s', f); end
end
logDir = fullfile(srcDir, 'logs');  if ~isfolder(logDir), mkdir(logDir); end
diary(fullfile(logDir, ['I1v2_' char(datetime('now', 'Format', 'yyyyMMdd_HHmmss')) '.log']));
cleanup = onCleanup(@() diary('off')); %#ok<NASGU>
[~, host] = system('hostname');
fprintf('I1 LOPO baselines HPC v002 | %s | host %s | MATLAB %s | src %s\n', stamp_local(), strtrim(host), version, srcDir);
out = loopClosureLOPOBaselines_v002();
o = fullfile(srcDir, '..', 'results', 'lopoBaselines_v002.mat');
if ~isfile(o), error('runI1HPC_v002:noOutput', 'Expected output missing: %s', o); end
fprintf('Beats all four baselines: %s\n', mat2str(out.beats));
fprintf('I1 finished | %s\n=== COMPLETED OK ===\n', stamp_local());
end

function s = stamp_local()
    s = string(datetime('now', 'TimeZone', 'Europe/London', 'Format', 'yyyy-MM-dd HH:mm:ss z'));
end
