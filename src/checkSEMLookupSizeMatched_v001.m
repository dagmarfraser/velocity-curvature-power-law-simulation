%% checkSEMLookupSizeMatched_v001.m
% Coherence item Q5, second half. checkEllipseSizeVsGrid_v001 showed the empirical ellipses
% differ in size from the grid's (semi-axes 49.0 x 22.7 harness-mm): Fraser's are about 0.59 of
% the grid's, Cook's and Hickman's about 1.16 times. The registered per-dataset SEM lookup
% snaps each dataset's sigma in real mm to the grid's sigma in harness-mm, so at the same
% sigma it misstates the dataset's noise relative to its path. Does the registered
% SEM-adequacy classification (36/42 cells) survive matching relative noise instead?
%   V0 CANONICAL      the registered lookup (checkSEMTempoMatched_v001 L46-56; knockDownFlags_v001):
%                     snap (alpha, sigma, fs) to the nearest grid node, mean SEM over all
%                     beta_gen x VGF rows, per pipeline. Regression anchor: reproduces the
%                     published SG-IRLS SEMs (tolerance 5e-4) and 36 of 42 adequate.
%   V1 SIZE-MATCHED   the same with sigma scaled by grid a / empirical a (median semi-axis per
%                     dataset, results/checkEllipseSizeVsGrid_v001.mat), so sigma/a matches.
% Reads:  src/perCoordinateSEM_v2_001.mat, results/checkEllipseSizeVsGrid_v001.mat
% Writes: results/checkSEMLookupSizeMatched_v001.mat
% USAGE:  from the project root: checkSEMLookupSizeMatched_v001
% Fraser, D.S. (2026)  v001

%% CONFIG
ROOT = fileparts(fileparts(mfilename("fullpath")));
addpath(fullfile(ROOT, "src")); addpath(genpath(fullfile(ROOT, "src", "functions")));
SEM_MAT = fullfile(ROOT, "src", "perCoordinateSEM_v2_001.mat");
ELL_MAT = fullfile(ROOT, "results", "checkEllipseSizeVsGrid_v001.mat");
OUT_MAT = fullfile(ROOT, "results", "checkSEMLookupSizeMatched_v001.mat");
PIPELINES = ["BWFD-OLS" "SG-OLS" "BWFD-LMLS" "SG-LMLS" "BWFD-IRLS" "SG-IRLS"];
% Registry exactly as checkSEMTempoMatched_v001 L21-23 (v004 centroids)
REG = table(["Zarandi"; "Cook_CTRL"; "Cook_ASD"; "Dhieb"; "Hickman_PLAC"; "Hickman_HALO"; "Fraser"], ...
            [3.18; 4.77; 5.06; 2.50; 5.34; 5.42; 4.289], [4.77; 8.15; 7.84; 7.50; 7.17; 7.42; 2.009], ...
            [100; 133; 133; 100; 133; 133; 240], 'VariableNames', ["name" "alpha" "sigma" "fs"]);
PUBLISHED_SGIRLS = [0.0079 0.0051 0.0046 0.0189 0.0040 0.0040 0.0032];
N_ADEQ_PUB = 36;
SEM_ADEQ = semAdequacyThreshold_v001();
for f = [SEM_MAT ELL_MAT], if ~isfile(f), error("semSize:input", "%s", "FAILED PATH: " + f); end, end

%% Load
T = load(SEM_MAT, "coordTable").coordTable;
for v = ["alpha" "sigma" "fs" "sem"], T.(v) = double(T.(v)); end
E = load(ELL_MAT, "T", "grid");  ET = E.T;
[ok, loc] = ismember(REG.name, ET.dataset);
if ~all(ok), error("semSize:ell", "%s", "checkEllipseSizeVsGrid_v001 output lacks a dataset"); end
REG.scale = E.grid.aFitMM ./ ET.aMM(loc);                          % grid a / empirical a
REG.sigmaMatched = REG.sigma .* REG.scale;
aN = unique(T.alpha);  sN = unique(T.sigma);  fN = unique(T.fs);

%% Lookups
out = table();
for d = 1:height(REG)
    [~, ai] = min(abs(aN - REG.alpha(d)));  [~, fi] = min(abs(fN - REG.fs(d)));
    [~, s0] = min(abs(sN - REG.sigma(d)));  [~, s1] = min(abs(sN - REG.sigmaMatched(d)));
    for p = PIPELINES
        c0 = T.sem(T.alpha == aN(ai) & T.sigma == sN(s0) & T.fs == fN(fi) & T.pipeline == p);
        c1 = T.sem(T.alpha == aN(ai) & T.sigma == sN(s1) & T.fs == fN(fi) & T.pipeline == p);
        if isempty(c0) || isempty(c1), error("semSize:noRows", "%s", REG.name(d) + " / " + p + ": no SEM rows"); end
        out = [out; table(REG.name(d), p, REG.scale(d), sN(s0), sN(s1), mean(c0, "omitnan"), mean(c1, "omitnan"), ...
            'VariableNames', ["dataset" "pipeline" "scale" "sigmaSnap" "sigmaSnapMatched" "semCanonical" "semMatched"])]; %#ok<AGROW>
    end
end
out.adequateCanonical = out.semCanonical < SEM_ADEQ;
out.adequateMatched   = out.semMatched < SEM_ADEQ;

%% Regression anchor
chk = out(out.pipeline == "SG-IRLS", :);
if any(abs(chk.semCanonical' - PUBLISHED_SGIRLS) > 5e-4) || nnz(out.adequateCanonical) ~= N_ADEQ_PUB
    disp(chk(:, ["dataset" "semCanonical"]));
    error("semSize:anchor", "%s", sprintf("V0 does not reproduce the published SEMs or %d/42 (got %d)", N_ADEQ_PUB, nnz(out.adequateCanonical)));
end
fprintf("REGRESSION ANCHOR passed: canonical lookup reproduces the published SG-IRLS SEMs and %d/42 adequate.\n\n", N_ADEQ_PUB);

%% Report
fprintf("Per dataset: sigma (mm), the size factor, and the snapped grid sigma before and after matching:\n");
disp(REG)
dom = ~ismember(out.dataset, ["Zarandi" "Dhieb"]);
fprintf("Adequate cells (SEM < MDC/2.77): canonical %d/42, size-matched %d/42; in domain %d/30 and %d/30\n", ...
    nnz(out.adequateCanonical), nnz(out.adequateMatched), nnz(out.adequateCanonical & dom), nnz(out.adequateMatched & dom));
chg = out(out.adequateCanonical ~= out.adequateMatched, :);
fprintf("Cells whose verdict changes: %d\n", height(chg));
if height(chg) > 0, disp(chg), end
S = groupsummary(out, "dataset", ["min" "max"], ["semCanonical" "semMatched"]);
disp(S)

runDate = string(datetime("now", "Format", "yyyy-MM-dd"));
save(OUT_MAT, "out", "REG", "SEM_ADEQ", "runDate");
fprintf("Saved: %s\n", OUT_MAT);
