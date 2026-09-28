%% checkSEMTempoMatched_v001.m
% Does the registered per-dataset SEM adequacy (36/42 cells, SEM < MDC/2.77) depend on
% averaging over grid tempos the dataset never drew at?
% Canonical lookup (plotPaperRoadmapSchematic_v001.m L84-97; knockDownFlags_v001.m):
%   mean SEM over ALL beta_gen x VGF rows at the dataset's snapped (alpha, sigma, fs).
% Tempo-matched lookup: the same mean restricted to (beta_gen, VGF) cells whose REALISED
%   tempo (checkGridTempoAtFold_v004, f0 = 10 orbits / duration) lies inside the dataset's
%   own trial tempo range (v015 per-trial f0, P_LO-P_HI percentiles).
% Also reports each dataset's f0 distribution (Finding #227 / SPEC 4.4 table, now scripted).
% Reads: src/perCoordinateSEM_v2_001.mat, results/gridTempoAtFold_v004.mat,
%        src/loopClosureResults_<dataset>_all_shaped_xu_v015.mat
% Writes: results/semTempoMatched_v001.mat

%% CONFIG
ROOT      = fileparts(fileparts(mfilename("fullpath")));
SEM_MAT   = fullfile(ROOT, "src", "perCoordinateSEM_v2_001.mat");
GRID_MAT  = fullfile(ROOT, "results", "gridTempoAtFold_v004.mat");
OUT_MAT   = fullfile(ROOT, "results", "semTempoMatched_v001.mat");
PIPELINES = ["BWFD-OLS" "SG-OLS" "BWFD-LMLS" "SG-LMLS" "BWFD-IRLS" "SG-IRLS"];
% Registry exactly as plotPaperRoadmapSchematic_v001.m L43-48 (v004 centroids)
REG = table(["Zarandi"; "Cook_CTRL"; "Cook_ASD"; "Dhieb"; "Hickman_PLAC"; "Hickman_HALO"; "Fraser"], ...
            [3.18; 4.77; 5.06; 2.50; 5.34; 5.42; 4.289], [4.77; 8.15; 7.84; 7.50; 7.17; 7.42; 2.009], ...
            [100; 133; 133; 100; 133; 133; 240], 'VariableNames', ["name" "alpha" "sigma" "fs"]);
SEM_ADEQUATE = 0.03 / 2.77;        % MDC / 2.77 = 0.0108
P_LO = 10;  P_HI = 90;             % trial tempo range used for matching
PUBLISHED_SGIRLS = [0.0079 0.0051 0.0046 0.0189 0.0040 0.0040 0.0032];   % self-check, REG order (schematic L105-107)

%% Load
for f = [SEM_MAT GRID_MAT], if ~isfile(f), error("semTempo:input", "%s", "Missing: " + f); end, end
S = load(SEM_MAT, "coordTable");  T = S.coordTable;
for v = ["betaGen" "VGF" "fs" "alpha" "sigma" "sem"], T.(v) = double(T.(v)); end
G = load(GRID_MAT);  N = G.out.perNode;                       % realised f0 per (beta_gen, VGF, fs)
N = renamevars(N(:, ["generated_beta" "vgf_value" "sampling_rate" "f0"]), ...
               ["generated_beta" "vgf_value" "sampling_rate"], ["betaGen" "VGF" "fs"]);
aN = unique(T.alpha);  sN = unique(T.sigma);  fN = unique(T.fs);

%% Per dataset: trial tempo distribution, snapped node, canonical vs tempo-matched SEM
tempo = table();  out = table();
for d = 1:height(REG)
    f = fullfile(ROOT, "src", "loopClosureResults_" + REG.name(d) + "_all_shaped_xu_v015.mat");
    if ~isfile(f), error("semTempo:input", "%s", "Missing: " + f); end
    L = load(f, "results");  f0 = [L.results.f0];  f0 = f0(isfinite(f0));
    lo = prctile(f0, P_LO);  hi = prctile(f0, P_HI);
    tempo = [tempo; table(REG.name(d), numel(f0), median(f0), lo, hi, ...
        'VariableNames', ["dataset" "nTrials" "f0median" "f0_p" + P_LO "f0_p" + P_HI])]; %#ok<AGROW>
    [~, ai] = min(abs(aN - REG.alpha(d)));  [~, si] = min(abs(sN - REG.sigma(d)));  [~, fi] = min(abs(fN - REG.fs(d)));
    sa = aN(ai);  ss = sN(si);  sf = fN(fi);
    tf = N(N.fs == sf, :);
    for p = PIPELINES
        c = T(T.alpha == sa & T.sigma == ss & T.fs == sf & T.pipeline == p, :);
        if isempty(c)
            error("semTempo:noRows", "%s", sprintf("No SEM rows for %s / %s at (%.2f, %.2f, %d)", REG.name(d), p, sa, ss, sf));
        end
        c = outerjoin(c, tf, "Keys", ["betaGen" "VGF" "fs"], "Type", "left", "MergeKeys", true);
        inRange = c.f0 >= lo & c.f0 <= hi;
        semAll = mean(c.sem, "omitnan");  semTM = mean(c.sem(inRange), "omitnan");
        out = [out; table(REG.name(d), p, sa, ss, sf, height(c), sum(isnan(c.f0)), sum(inRange), ...
            semAll, semTM, semAll < SEM_ADEQUATE, semTM < SEM_ADEQUATE, ...
            'VariableNames', ["dataset" "pipeline" "alphaSnap" "sigmaSnap" "fsSnap" "nRows" "nNoTempo" ...
            "nTempoMatched" "semCanonical" "semTempoMatched" "adequateCanonical" "adequateTempoMatched"])]; %#ok<AGROW>
    end
end

%% Self-check against the paper's published SG-IRLS values
chk = out(out.pipeline == "SG-IRLS", :);
for d = 1:height(REG)
    ref = PUBLISHED_SGIRLS(d);
    if isfinite(ref) && abs(chk.semCanonical(d) - ref) > 5e-4
        error("semTempo:repro", "%s", sprintf("%s SG-IRLS canonical SEM %.5f does not reproduce published %.4f", ...
            REG.name(d), chk.semCanonical(d), ref));
    end
end
fprintf("Canonical lookup reproduces the published SG-IRLS SEMs (7 datasets, tolerance 5e-4).\n");

%% Report
tempo{:, 3:end} = round(tempo{:, 3:end}, 3);
fprintf("\nTrial tempo by dataset (v015 per-trial f0, Hz):\n");  disp(tempo)
fprintf("SEM adequacy (SEM < %.4f), 42 cells:\n", SEM_ADEQUATE);
fprintf("  canonical (all beta_gen x VGF):  %d / 42\n", sum(out.adequateCanonical));
fprintf("  tempo-matched (p%d-p%d of trial f0): %d / 42\n", P_LO, P_HI, sum(out.adequateTempoMatched));
dom = out.dataset ~= "Zarandi";
fprintf("  within the validity domain (36 cells): canonical %d, tempo-matched %d\n", ...
    sum(out.adequateCanonical & dom), sum(out.adequateTempoMatched & dom));
chg = out(out.adequateCanonical ~= out.adequateTempoMatched, :);
fprintf("\nCells whose verdict changes: %d\n", height(chg));
if height(chg) > 0, disp(chg(:, ["dataset" "pipeline" "semCanonical" "semTempoMatched" "nTempoMatched"])), end
disp(out(:, ["dataset" "pipeline" "nRows" "nTempoMatched" "semCanonical" "semTempoMatched"]))
fprintf("Ratio tempo-matched / canonical SEM: median %.2f, range %.2f-%.2f\n", ...
    median(out.semTempoMatched ./ out.semCanonical), min(out.semTempoMatched ./ out.semCanonical), ...
    max(out.semTempoMatched ./ out.semCanonical));

if ~isfolder(fileparts(OUT_MAT)), error("semTempo:outDir", "%s", "Missing folder: " + fileparts(OUT_MAT)); end
save(OUT_MAT, "out", "tempo", "REG", "SEM_ADEQUATE", "P_LO", "P_HI", "-v7.3");
fprintf("Saved: %s\n", OUT_MAT);
