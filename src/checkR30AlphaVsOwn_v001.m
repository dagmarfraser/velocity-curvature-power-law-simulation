%% checkR30AlphaVsOwn_v001.m
% The noise-colour check of spec 10.57 C, from saved files. Each noise replacement's median IRASA
% alpha is set against the median of the trial's own (seam-free) residual run, AXr.
% Why: analyseR30_v001 ran before the AXr run files were on this machine, so its dAlphaVsAXr is NaN,
% and spec 10.57 C read the AXr medians by hand (Fraser 4.3275, PLAC 5.0937). This makes that step
% reproducible. The AXr median is computed as in analyseR30_v001 (alphaIRA over all trials, omitnan);
% arm medians are read as saved, from analyseR30_v001.mat (R30b) and analyseR29SeamFree_v001.mat (corr).
% Anchor (fail loud): the AXr medians reproduce 10.57 C to its 4 decimals.
% Run from src/. Writes results/tempoSim/checkR30AlphaVsOwn_v001.mat and a dated log in logs/.

%% CONFIG
ROOT    = fileparts(fileparts(mfilename("fullpath")));
SIM     = fullfile(ROOT, "results", "tempoSim");
DS      = ["Fraser" "Hickman_PLAC"];
AXR_RUN = ["loopClosureResults_SYN-AXrF-Fraser_all_shaped_xu_syn001.mat" "loopClosureResults_SYN-AXr-Hickman_PLAC_all_shaped_xu_syn001.mat"];
ANCHOR  = [4.3275 5.0937];  TOL_A = 5.01e-5;   % spec 10.57 C, quoted to 4 decimals
CRIT    = 0.05;                                 % spec 10.55 R30b criterion (a judgement call, 10.57 C)
OUT     = fullfile(SIM, "checkR30AlphaVsOwn_v001.mat");

%% Setup
logFile = fullfile(ROOT, "logs", "checkR30AlphaVsOwn_v001_" + string(datetime("now", "Format", "yyyyMMdd_HHmmss")) + ".log");
LG = fopen(logFile, "w");  if LG < 0, error("r30a:log", "%s", "Cannot open log: " + logFile); end
say_local(LG, "checkR30AlphaVsOwn_v001 run %s\n", string(datetime("now")));
R30 = load(fullfile(SIM, "analyseR30_v001.mat"), "R30b").R30b;
SF  = load(fullfile(SIM, "analyseR29SeamFree_v001.mat"), "corr").corr;

%% Own-residual medians and the anchor
medAXr = NaN(1, numel(DS));
for d = 1:numel(DS)
    f = fullfile(SIM, AXR_RUN(d));  if ~isfile(f), error("r30a:input", "%s", "Missing: " + f); end
    r = load(f, "results").results;
    medAXr(d) = median(arrayfun(@(x) x.alphaIRA, r), "omitnan");
    say_local(LG, "%s: own-residual (AXr) median alpha %.5f (%d trials)\n", DS(d), medAXr(d), numel(r));
end
if any(abs(medAXr - ANCHOR) > TOL_A)
    error("r30a:anchor", "%s", sprintf("AXr medians %.5f, %.5f do not reproduce spec 10.57 C (%.4f, %.4f)", medAXr, ANCHOR));
end
say_local(LG, "Anchor passed: own-residual medians reproduce spec 10.57 C.\n");

%% Each arm against its dataset's own residual
A = [R30(:, ["dataset" "arm" "medianAlpha"]); SF(:, ["dataset" "arm" "medianAlpha"])];
[~, di] = ismember(A.dataset, DS);
if any(di == 0), error("r30a:dataset", "%s", "Unknown dataset in the saved tables"); end
A.medAXr = medAXr(di)';
A.dAlpha = A.medianAlpha - A.medAXr;
A.colourPreserved = abs(A.dAlpha) <= CRIT;
A.colourPreserved(A.arm == "AXnF") = false;   % white control: expected not to preserve colour; reported, not tested
say_local(LG, "\nMedian alpha per arm against the own residual (criterion |d| <= %.2f; white control not tested):\n%s", CRIT, evalc("disp(A)"));
save(OUT, "A", "medAXr", "DS", "CRIT");
say_local(LG, "Saved: %s\nLog: %s\n", OUT, logFile);  fclose(LG);

%% =========================================================================
function say_local(fid, fmt, varargin)
    fprintf(fmt, varargin{:});  fprintf(fid, fmt, varargin{:});
end
