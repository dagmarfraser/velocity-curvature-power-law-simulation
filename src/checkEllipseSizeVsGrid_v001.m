%% checkEllipseSizeVsGrid_v001.m
% Coherence item Q5: the registered SEM lookups (the 36/42 classification, Fig 0, §3c's
% empirical operating points) match each dataset's noise sigma, in real mm, to the grid's
% sigma, in harness-mm, on one ellipse. How much noise matters depends on sigma relative to
% the path's size, so this reports the empirical ellipse sizes beside the grid's and the
% relative noise sigma/a on each side. Pointwise inversion is not affected: each trial's
% forward map is built at that trial's own geometry.
%   grid       gridShapeGeometry_v002: least-squares ellipse semi-axes (px) / pixelScale.
%   empirical  per trial, a_mm and b_mm (semi-axes carried in the v015 corpus; checked against
%              Zarandi's 8 x 2 cm template, Zarandi et al. 2023) and
%              sigmaMM; per dataset, medians over trials.
%   ratio      grid a / empirical a: the factor by which a lookup at the same sigma understates
%              (ratio > 1) or overstates (< 1) the dataset's relative noise.
% Reads:  src/loopClosureResults_<dataset>_all_shaped_xu_v015.mat; the grid shape via
%         gridShapeGeometry_v002 (src/functions/AccuracyMeasure/Sinusoidal_Curves.mat,
%         src/Toolchain_caller_v058.m)
% Writes: results/checkEllipseSizeVsGrid_v001.mat
% USAGE:  from the project root: checkEllipseSizeVsGrid_v001
% Fraser, D.S. (2026)  v001

%% CONFIG
ROOT = fileparts(fileparts(mfilename("fullpath")));
addpath(fullfile(ROOT, "src")); addpath(genpath(fullfile(ROOT, "src", "functions")));
DATASETS = ["Fraser" "Cook_CTRL" "Cook_ASD" "Hickman_PLAC" "Hickman_HALO" "Dhieb" "Zarandi"];
IN_DOM   = [true true true true true false false];               % #242
A_FIT_PUB = [235.06 108.85];                                    % checkGridShapeUnits_v001 (px)
OUT_MAT  = fullfile(ROOT, "results", "checkEllipseSizeVsGrid_v001.mat");

%% Grid
G = gridShapeGeometry_v002(ROOT);
if any(abs([G.aFitPx G.bFitPx] - A_FIT_PUB) > 0.01)
    error("ellSize:grid", "%s", sprintf("grid fit %.2f x %.2f px, expected %.2f x %.2f", G.aFitPx, G.bFitPx, A_FIT_PUB));
end
fprintf("Grid ellipse: semi-axes %.2f x %.2f px = %.1f x %.1f harness-mm (pixelScale %.4f px/mm); a/b %.2f\n", ...
    G.aFitPx, G.bFitPx, G.aFitMM, G.bFitMM, G.pixelScale, G.abFit);

%% Empirical
ZAR_SEMI_A = 80;                                                % Zarandi et al. (2023): 8 x 2 cm semi-axes
T = table();
for k = 1:numel(DATASETS)
    f = fullfile(ROOT, "src", "loopClosureResults_" + DATASETS(k) + "_all_shaped_xu_v015.mat");
    if ~isfile(f), error("ellSize:input", "%s", "FAILED PATH: " + f); end
    R = load(f, "results").results;
    a = arrayfun(@(r) double(r.a_mm), R(:));  b = arrayfun(@(r) double(r.b_mm), R(:));
    s = arrayfun(@(r) double(r.sigmaMM), R(:));
    ok = isfinite(a) & a > 0 & isfinite(b) & b > 0 & isfinite(s);
    if nnz(ok) < numel(R)
        warning("ellSize:missing", "%s", sprintf("%s: %d of %d trials lack a, b or sigma; excluded", DATASETS(k), numel(R) - nnz(ok), numel(R)));
    end
    T = [T; table(DATASETS(k), IN_DOM(k), nnz(ok), median(a(ok)), median(b(ok)), median(a(ok) ./ b(ok)), ...
        median(s(ok)), median(s(ok) ./ a(ok)), median(s(ok)) / G.aFitMM, G.aFitMM / median(a(ok)), ...
        'VariableNames', ["dataset" "inDomain" "n" "aMM" "bMM" "ab" "sigmaMM" "sigmaOverA" "gridSigmaOverA" "gridAoverEmpA"])]; %#ok<AGROW>
end
fprintf("\nPer dataset (medians over trials). gridSigmaOverA: the dataset's sigma on the grid's ellipse.\n");
fprintf("gridAoverEmpA > 1: a lookup at the same sigma understates the dataset's relative noise by that factor.\n");
disp(T)
% Basis check: a_mm must be a semi-axis (Zarandi's template is 80 x 20 mm in semi-axes)
za = T.aMM(T.dataset == "Zarandi");
if abs(za / ZAR_SEMI_A - 1) > 0.25
    error("ellSize:basis", "%s", sprintf("Zarandi median a_mm = %.1f; not a semi-axis of the 80 mm template (full axis would be 160)", za));
end
fprintf("Basis check passed: Zarandi median a_mm %.1f against its 80 mm template semi-axis.\n", za);
d = T(T.inDomain, :);
fprintf("In-domain: semi-axes %.1f-%.1f x %.1f-%.1f mm (grid %.1f x %.1f); sigma/a %.3f-%.3f against %.3f-%.3f on the grid; factor %.2f-%.2f\n", ...
    min(d.aMM), max(d.aMM), min(d.bMM), max(d.bMM), G.aFitMM, G.bFitMM, min(d.sigmaOverA), max(d.sigmaOverA), ...
    min(d.gridSigmaOverA), max(d.gridSigmaOverA), min(d.gridAoverEmpA), max(d.gridAoverEmpA));

runDate = string(datetime("now", "Format", "yyyy-MM-dd"));
grid = struct("aFitPx", G.aFitPx, "bFitPx", G.bFitPx, "aFitMM", G.aFitMM, "bFitMM", G.bFitMM, "pixelScale", G.pixelScale);
save(OUT_MAT, "T", "grid", "runDate");
fprintf("Saved: %s\n", OUT_MAT);
