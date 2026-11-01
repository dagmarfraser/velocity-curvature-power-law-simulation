%% checkV016Local_v001.m
% D-D16 stage S2, local checks of runner v016 before it goes to the HPC (docs/SPEC_D16_SamplingRateCutoff_v001.md).
% (1) Generates src/runLoopClosureFftnoise_v016.m from v015 (makeRunLoopClosureV016_v001) if absent.
% (2) Defaults: v016 and v015 on the same Fraser trials (TrialID override, serial, same seed, small N_REPS) must
%     agree exactly in every saved field (beta_obs, beta_gen*, CIs, status, the forward map). Serial runs share
%     the RNG stream, so here even the stochastic fields must be identical; a difference is a generator bug.
% (3) Options: v016 with Subsample 2 (A120p), Subsample 4 + EdgeClip 5 (A60d) and BWCutoff 5, SGRef 34, EdgeClip 41
%     (S240fc05) must reproduce stage S1's beta_obs for those settings (results/obsBetaRateCutoff_v002.mat), so
%     S1 and S2 apply the same subsampling, filters and clip.
% Runs in a scratch folder under tempdir (outputs and logs go there, not to src/), with links to src/functions,
% src/req and the noise file. Trials: the first two of the v015 Fraser set without SG curvature residue (S1 RS).
% Writes a dated log in logs/. Run from src/ after obsBetaRateCutoff_v002 has finished.

%% CONFIG
ROOT  = fileparts(fileparts(mfilename("fullpath")));  SRC = fullfile(ROOT, "src");
IN_S1 = fullfile(ROOT, "results", "obsBetaRateCutoff_v002.mat");
SMALL = {'TrialSelection', 'all', 'NoiseModel', 'shaped_xu', 'N_BETA', 5, 'N_REPS', 2, 'N_R7', 2, 'RngSeed', 1729, 'UseParfor', false};
OPTS  = struct("name", {"A120p" "A60d" "S240fc05"}, "args", {{'Subsample', 2}, {'Subsample', 4, 'EdgeClip', 5}, {'BWCutoff', 5, 'SGRef', 34, 'EdgeClip', 41}});
N_TRIALS = 2;  TOL_S1 = 1e-12;

%% Setup
logFile = fullfile(ROOT, "logs", "checkV016Local_v001_" + string(datetime("now", "Format", "yyyyMMdd_HHmmss")) + ".log");
LG = fopen(logFile, "w");  if LG < 0, error("v16:log", "%s", "Cannot open log: " + logFile); end
say_local(LG, "checkV016Local_v001 run %s\n", string(datetime("now")));
if ~isfile(IN_S1), error("v16:input", "%s", "Missing (run obsBetaRateCutoff_v002 first): " + IN_S1); end
if ~isfile(fullfile(SRC, "runLoopClosureFftnoise_v016.m")), makeRunLoopClosureV016_v001("SrcDir", SRC); end
S1 = load(IN_S1, "B", "RS", "SET", "ids");
pick = find(S1.RS(1, :) == 0, N_TRIALS);  tids = S1.ids(pick);
say_local(LG, "Trials: %s\n", strjoin(tids, ", "));

TMP = fullfile(tempdir, "d16smoke_" + string(datetime("now", "Format", "yyyyMMdd_HHmmss")));  mkdir(TMP);  mkdir(fullfile(TMP, "logs"));
for f = ["runLoopClosureFftnoise_v015.m" "runLoopClosureFftnoise_v016.m"], copyfile(fullfile(SRC, f), TMP); end
for f = ["functions" "req" "noiseCharacterisation_fraser.mat"]
    [st, out] = system(sprintf('ln -s "%s" "%s"', fullfile(SRC, f), fullfile(TMP, f)));
    if st ~= 0, error("v16:link", "%s", "Cannot link " + f + ": " + out); end
end
here = pwd;  back = onCleanup(@() cd(here));  cd(TMP);
if ~contains(string(which("runLoopClosureFftnoise_v016")), extractAfter(TMP, tempdir)), error("v16:path", "%s", "v016 does not resolve to the scratch copy"); end

%% (2) Defaults: v016 is v015
for t = tids'
    a = runLoopClosureFftnoise_v015('Fraser', SMALL{:}, 'TrialID', char(t));
    b = runLoopClosureFftnoise_v016('Fraser', SMALL{:}, 'TrialID', char(t));
    fa = setdiff(fieldnames(a), "runDate");  fb = setdiff(fieldnames(b), "runDate");
    if ~isequal(fa, fb), error("v16:fields", "%s", "Saved fields differ for " + t); end
    bad = fa(~cellfun(@(f) isequaln(a.(f), b.(f)), fa));
    if ~isempty(bad), error("v16:defaults", "%s", "v016 at the defaults differs from v015 for " + t + " in: " + strjoin(string(bad), ", ")); end
    say_local(LG, "Defaults, %s: v016 = v015 in all %d fields (beta_gen* SG-IRLS %.4f).\n", t, numel(fa), b.betaGenStar(6));
end

%% (3) Options reproduce stage S1's beta_obs
for o = OPTS
    s = find([S1.SET.name] == o.name);
    for j = 1:numel(tids)
        b = runLoopClosureFftnoise_v016('Fraser', SMALL{:}, 'TrialID', char(tids(j)), o.args{:});
        d = max(abs(b.betaObs(:)' - squeeze(S1.B(s, pick(j), :))'));
        if ~(d <= TOL_S1), error("v16:s1", "%s", sprintf("%s, %s: v016 beta_obs differs from S1 by %.3g", o.name, tids(j), d)); end
        say_local(LG, "%-9s %s: v016 beta_obs = S1 (max diff %.2g); beta_gen* SG-IRLS %.4f\n", o.name, tids(j), d, b.betaGenStar(6));
    end
end
say_local(LG, "All checks passed. Scratch folder: %s\nLog: %s\n", TMP, logFile);  fclose(LG);

%% =========================================================================
function say_local(fid, fmt, varargin)
    fprintf(fmt, varargin{:});  fprintf(fid, fmt, varargin{:});
end
