function runD16HPC_v001()
%RUND16HPC_V001  D-D16 stage S2 on the HPC (docs/SPEC_D16_SamplingRateCutoff_v001.md, section 3): forward maps and
%   beta_gen* for the Fraser trials at other sampling rates and filter cutoffs, with runner v016.
%   Open this file in the MATLAB GUI on the HPC (sinteractive) and press Run. Restartable: the runner skips any
%   output that exists. About 8 h per 240 Hz setting on 72 workers (scaled from v015, not timed); the subsampled
%   arms (clean subsampling, no filter) are cheaper. Log: src/logs/D16_<stamp>.log. Afterwards the Mac runs the S2 analysis (rule R2).
%   Steps: (1) generate runLoopClosureFftnoise_v016 from v015 if absent; (2) acceptance test: v016 at the
%   defaults on Zarandi reproduces the Zarandi v015 corpus (beta_obs exactly, the stochastic fields within
%   runLoopClosureAll_v015's tolerances, Finding #204); (3) the Fraser settings below, cheapest first.
%   Fraser, D.S. (2026)  v001

%% CONFIG
SETTINGS = struct( ...                                   % name, Subsample (clean stride), BWCutoff, SGRef, EdgeClip
    "name", {"A60p" "A120p" "S240fc10" "S240fc05" "S240fc07" "S240fc14" "S240fc20"}, ...
    "dec",  {4 2 1 1 1 1 1}, "fc", {10 10 10 5 7 14 20}, "ref", {17 17 17 34 24 12 8.5}, "clip", {20 20 41 41 41 41 41});
CLIP_ARMS = struct("name", {"A60d" "A120d"}, "dec", {4 2}, "fc", {10 10}, "ref", {17 17}, "clip", {5 10});
RUN_CLIP_ARMS = false;   % rule R3 not triggered: logs/obsBetaRateCutoff_v002_20261024_131133.log, median |p - d| at most 0.0033 < 0.0108
RUN = struct('TrialSelection', 'all', 'NoiseModel', 'shaped_xu', 'N_BETA', 25, 'N_REPS', 200, 'N_R7', 10, 'RngSeed', 1729);  % as runLoopClosureAll_v015
TOL_OBS = 1e-6;  TOL_MEDIAN = 0.005;  TOL_NANFRAC = 0.05;   % acceptance (stochastic tolerances as runLoopClosureAll_v015)

%% Setup
srcDir = fileparts(mfilename('fullpath'));  addpath(srcDir);
need = ["runLoopClosureFftnoise_v015.m" "makeRunLoopClosureV016_v001.m" "noiseCharacterisation_fraser.mat" ...
        "noiseCharacterisation_zarandi.mat" "loopClosureResults_Zarandi_all_shaped_xu_v015.mat" "loopClosureResults_Fraser_all_shaped_xu_v015.mat"];
for f = need, if ~isfile(fullfile(srcDir, f)), error('runD16HPC:missing', 'Missing %s', fullfile(srcDir, f)); end, end
logDir = fullfile(srcDir, 'logs');  if ~isfolder(logDir), mkdir(logDir); end
diary(fullfile(logDir, ['D16_' char(datetime('now', 'Format', 'yyyyMMdd_HHmmss')) '.log']));
cleanup = onCleanup(@() diary('off')); %#ok<NASGU>
[~, host] = system('hostname');
fprintf('D16 S2 HPC | %s | host %s | MATLAB %s | cores %d | src %s\n', stamp_local(), strtrim(host), version, feature('numcores'), srcDir);
if ~isfile(fullfile(srcDir, 'runLoopClosureFftnoise_v016.m')), makeRunLoopClosureV016_v001('SrcDir', srcDir); end
pool_local(min(72, feature('numcores')));

%% Acceptance: v016 at the defaults is v015 (Zarandi)
t0 = tic;
args = namedargs2cell(RUN);
r16 = runLoopClosureFftnoise_v016('Zarandi', 'UseParfor', true, args{:});
r15 = load(fullfile(srcDir, 'loopClosureResults_Zarandi_all_shaped_xu_v015.mat'), 'results').results;
fprintf('Acceptance (Zarandi, %.1f min): %s\n', toc(t0) / 60, accept_local(r15, r16, TOL_OBS, TOL_MEDIAN, TOL_NANFRAC));

%% Fraser settings
S = SETTINGS;  if RUN_CLIP_ARMS, S = [S CLIP_ARMS]; end
for k = 1:numel(S)
    pool_local(min(72, feature('numcores')));  t0 = tic;
    fprintf('\n=== %s: subsample stride %d, BW %g Hz, SG ref %g, clip %d | %s ===\n', S(k).name, S(k).dec, S(k).fc, S(k).ref, S(k).clip, stamp_local());
    r = runLoopClosureFftnoise_v016('Fraser', 'UseParfor', true, args{:}, ...
        'Subsample', S(k).dec, 'BWCutoff', S(k).fc, 'SGRef', S(k).ref, 'EdgeClip', S(k).clip);
    g = arrayfun(@(q) q.betaGenStar(6), r);  o = arrayfun(@(q) q.betaObs(6), r);
    fprintf('%s done in %.1f h: %d trials; SG-IRLS median beta_obs %.4f, beta_gen* %.4f (finite %d)\n', ...
        S(k).name, toc(t0) / 3600, numel(r), median(o, 'omitnan'), median(g, 'omitnan'), sum(isfinite(g)));
end
fprintf('D16 S2 finished | %s\n=== COMPLETED OK ===\n', stamp_local());
end

%% =========================================================================
function msg = accept_local(r15, r16, tolObs, tolMed, tolNaN)
% As runLoopClosureAll_v015 verify_local: beta_obs is deterministic (must match; NaN pattern exactly); beta_gen*,
% ciLo, ciHi come from the synthetic forward map and are checked against tolerances.
    if numel(r15) ~= numel(r16), error('runD16HPC:count', 'v016 %d trials, v015 %d', numel(r16), numel(r15)); end
    if ~isfield(r16, 'betaRecCurveMed'), error('runD16HPC:noCurve', 'v016 output lacks betaRecCurveMed'); end
    if ~isequal(string({r15.trialID}), string({r16.trialID})), error('runD16HPC:ids', 'Trial order differs from v015'); end
    A = vertcat(r15.betaObs);  B = vertcat(r16.betaObs);
    if ~isequal(isfinite(A), isfinite(B)), error('runD16HPC:obsNaN', 'beta_obs NaN pattern differs from v015'); end
    dObs = max(abs(A(isfinite(A)) - B(isfinite(A))));
    if dObs > tolObs, error('runD16HPC:obsDiff', 'beta_obs differs from v015 by %.3g (tolerance %.0g)', dObs, tolObs); end
    d = [];  flip = 0;  n = 0;
    for f = ["betaGenStar" "ciLo" "ciHi"]
        a = vertcat(r15.(f));  b = vertcat(r16.(f));  n = n + numel(a);  flip = flip + nnz(isfinite(a) ~= isfinite(b));
        m = isfinite(a) & isfinite(b);  d = [d; abs(a(m) - b(m))]; %#ok<AGROW>
    end
    if median(d) > tolMed, error('runD16HPC:stochDrift', 'median |diff| %.4f above %.4f', median(d), tolMed); end
    if flip / n > tolNaN, error('runD16HPC:nanDrift', '%.1f%% of cells flipped NaN state', 100 * flip / n); end
    msg = sprintf('PASSED: beta_obs max diff %.3g; stochastic median %.4f, p90 %.4f, max %.4f; NaN flips %d of %d', ...
        dObs, median(d), quantile(d, 0.9), max(d), flip, n);
end

function p = pool_local(nw)
% A live pool with IdleTimeout disabled (as runLoopClosureAll_v015); fail loud rather than run serially.
    p = gcp('nocreate');
    if ~isempty(p)
        try, alive = p.Connected && p.NumWorkers > 0; catch, alive = false; end
        if ~alive, try, delete(p); catch, end, p = []; end
    end
    if isempty(p), p = parpool('local', nw, 'IdleTimeout', Inf); elseif ~isinf(p.IdleTimeout), p.IdleTimeout = Inf; end
    if isempty(p) || ~p.Connected, error('runD16HPC:noPool', 'Could not establish a parallel pool'); end
end

function s = stamp_local()
    s = string(datetime('now', 'TimeZone', 'Europe/London', 'Format', 'yyyy-MM-dd HH:mm:ss z'));
end
