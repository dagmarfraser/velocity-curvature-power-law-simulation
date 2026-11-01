%% checkEdgeClipShare_v001.m
% Share of each synthetic trajectory discarded by the harness's fixed edge clip, by sampling rate (Finding #188, now script-backed).
% The clip keeps rows edgeClip:end-edgeClip (Toolchain_func_v032 L299), i.e. N - 2*edgeClip + 1 of N samples. N comes from
% generateSyntheticData_v011 (ellipse, shape 6; canvas, orbitCount and resample as Toolchain_caller_v058) at beta_gen = 1/3, at
% the log-midpoint VGF exp(5.15) and at each of the 14 grid VGF values exp(linspace(4.5, 5.8, 14)). Reads only the cached
% baseline shapes (regen = 0; errors if a cache is missing or changes). Run from src/.
% Writes results/checkEdgeClipShare_v001.mat; log logs/checkEdgeClipShare_v001_<date>.log.   Fraser, D.S. (2026)  v001
srcDir = fileparts(mfilename("fullpath"));  ROOT = fileparts(srcDir);  fDir = fullfile(srcDir, "functions");  addpath(fDir);
L = fullfile(ROOT, "logs");  if ~isfolder(L), mkdir(L); end
diary(fullfile(L, "checkEdgeClipShare_v001_" + string(datetime("now", "Format", "yyyyMMdd_HHmmss")) + ".log"));  cleanupDiary = onCleanup(@() diary("off"));
cs = fileread(fullfile(srcDir, "Toolchain_caller_v058.m"));
get = @(nm) str2double(regexp(cs, "cfg\." + nm + "\s*=\s*(\d+)\s*;", "tokens", "once"));
edgeClip = get("edgeClip");  orbitCount = get("orbitCount");  resampleN = get("resample");
canvas = str2double(regexp(cs, "cfg\.canvas\s*=\s*\[(\d+);\s*(\d+)\]", "tokens", "once"))';
if any(isnan([edgeClip orbitCount resampleN canvas'])), error("edge:cfg", "%s", "could not read edgeClip/orbitCount/resample/canvas from Toolchain_caller_v058.m"); end
SHAPE = 6;  BETA = 1/3;  FS = [60 120 240];  VGF = exp(linspace(4.5, 5.8, 14));  VMID = exp(5.15);
cache = arrayfun(@(f) fullfile(fDir, sprintf("baselineShp%d_%dHz.mat", SHAPE, f)), FS);
if ~all(isfile(cache)), error("edge:cache", "missing baseline cache: %s", strjoin(cache(~isfile(cache)), ", ")); end
t0 = arrayfun(@(c) dir(c).datenum, cache);
vg = [VMID VGF];  N = zeros(numel(FS), numel(vg));
for i = 1:numel(FS), for j = 1:numel(vg)
    xt = generateSyntheticData_v011(SHAPE, canvas, FS(i), BETA, vg(j), orbitCount, resampleN, 0, 0);  N(i, j) = numel(xt);
end, end
if ~isequal(t0, arrayfun(@(c) dir(c).datenum, cache)), error("edge:cacheChanged", "%s", "a baseline cache file changed during the run"); end
kept = N - 2*edgeClip + 1;  share = 1 - kept ./ N;
T = table(FS', N(:, 1), kept(:, 1), 100*share(:, 1), 100*min(share(:, 2:end), [], 2), 100*max(share(:, 2:end), [], 2), ...
    'VariableNames', ["fs" "N_mid" "kept_mid" "pctDiscarded_mid" "pctDiscarded_gridMin" "pctDiscarded_gridMax"]);
fprintf("edgeClip = %d samples; beta_gen = 1/3; mid VGF = exp(5.15) = %.1f\n", edgeClip, VMID);  disp(T);
ratio = share(1, 2:end) ./ share(3, 2:end);  fprintf("60/240 Hz discarded-share ratio across grid VGF: %.3f to %.3f\n", min(ratio), max(ratio));
save(fullfile(ROOT, "results", "checkEdgeClipShare_v001.mat"), "T", "N", "kept", "share", "vg", "FS", "edgeClip", "ratio");
fprintf("=== checkEdgeClipShare_v001 DONE ===\n");
