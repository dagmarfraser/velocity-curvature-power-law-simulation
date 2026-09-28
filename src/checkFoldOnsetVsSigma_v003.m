%% checkFoldOnsetVsSigma_v003.m
% Which feature of the beta_gen -> beta_rec forward map breaks monotonicity
% first as sigma rises from zero, and does either feature have an onset near
% FINDINGS_REFERENCE Key finding #2's "sigma ~0.5 mm (white) to ~1.8 mm (red)"?
%   low-end feature : collapse-and-recovery near the origin (descent at beta_gen <= LOW_MAX)
%   high fold       : bend-back past the peak             (descent at beta_gen >  FOLD_MIN)
% A descent counts only if it exceeds K_SE standard errors of the difference
% between neighbouring replicate means (SE = sem/sqrt(nReps)), plus DY_FLOOR.
% Reads perCoordinateSEM_v2_001.mat (coordTable, written by
% computePerCoordinateSEM_v2_001.m L191). No new simulation. Lives in src/; paths
% resolve from the script location, so it runs from any working folder.
% v002: path resolution only (v001 assumed cwd = project root).
% v003: onset table extended to the empirical colours (alpha 4.2-5.4) out to sigma = 10 mm;
%       alpha matched to the nearest grid node (printed); floor table beta_rec(beta_gen = 0);
%       figure keeps ALPHAS_PLOT; results saved to results/foldOnsetVsSigma_v003.mat.

%% CONFIG
ROOT      = fileparts(fileparts(mfilename("fullpath")));   % src/.. = project root
SRC_MAT   = fullfile(ROOT, "src", "perCoordinateSEM_v2_001.mat");
OUT_PNG   = fullfile(ROOT, "figures", "foldOnsetVsSigma_v003.png");
FS        = 120;            % Hz
ALPHAS    = [0 1 2 3 4.2 4.8 5.0 5.4];  % #2 range (0-3) plus the empirical colours
ALPHAS_PLOT = [0 1 2 3];
SIGMA_MAX = 10;             % mm; empirical sigma reaches ~8
FLOOR_SIGMA = [2 8];        % mm; floor table (Fraser ~2, Cook/Hickman ~7-8)
OUT_MAT   = fullfile(ROOT, "results", "foldOnsetVsSigma_v003.mat");
LOW_MAX   = 0.20;           % beta_gen bound for the low-end feature
FOLD_MIN  = 0.40;           % beta_gen bound for the high fold
K_SE      = 2;
DY_FLOOR  = 1e-6;           % guards sigma = 0, where replicate SE is exactly 0
PLOT_PIPES = ["BWFD-OLS" "SG-IRLS"];

%% Load and validate
if ~isfile(SRC_MAT)
    error("foldOnset:noMat", "%s", "Not found: " + SRC_MAT);
end
S = load(SRC_MAT, "coordTable");
T = S.coordTable;
need = ["pipeline" "betaGen" "VGF" "fs" "alpha" "sigma" "meanBetaRec" "sem" "nReps"];
miss = setdiff(need, string(T.Properties.VariableNames));
if ~isempty(miss)
    error("foldOnset:columns", "%s", "coordTable missing: " + strjoin(miss, ", "));
end
for v = ["betaGen" "VGF" "fs" "alpha" "sigma" "meanBetaRec" "sem" "nReps"]
    T.(v) = double(T.(v));
end
T = T(T.fs == FS, :);
if isempty(T), error("foldOnset:fs", "%s", "No rows at fs = " + FS); end

vgfNodes = unique(T.VGF);
VGF0     = vgfNodes(ceil(numel(vgfNodes)/2));          % median VGF node
T        = T(T.VGF == VGF0, :);
sigNodes = unique(T.sigma);  sigNodes = sigNodes(sigNodes <= SIGMA_MAX);
pipes    = string(categories(removecats(T.pipeline)));
aNodes   = unique(T.alpha);
nodeOf   = @(a) aNodes(find(abs(aNodes - a) == min(abs(aNodes - a)), 1));
fprintf("alpha nodes used: %s\n", strjoin(compose("%.2f->%.2f", ALPHAS(:), arrayfun(nodeOf, ALPHAS(:))), ", "));
fprintf("fs = %d Hz, VGF node = %.1f mm/s, sigma nodes <= %g mm: %s\n", ...
    FS, VGF0, SIGMA_MAX, strjoin(string(sigNodes'), " "));

%% Walk every curve
rows = {};  nMissing = 0;
for p = pipes'
    for a = ALPHAS
        for s = sigNodes'
            c = T(T.pipeline == p & T.alpha == nodeOf(a) & abs(T.sigma - s) < 1e-9, :);
            if height(c) < 3, nMissing = nMissing + 1; continue; end
            c   = sortrows(c, "betaGen");
            b   = c.betaGen;  y = c.meanBetaRec;
            se  = c.sem ./ sqrt(max(c.nReps, 1));
            dy  = diff(y);
            sig = dy < -(K_SE * hypot(se(1:end-1), se(2:end)) + DY_FLOOR);
            bAt = b(2:end);                                 % descent located at its upper node
            [pk, iPk] = max(y);
            rows(end+1, :) = {p, a, s, height(c), any(sig & bAt <= LOW_MAX), ...
                any(sig & bAt > FOLD_MIN), sum(sig), b(iPk), pk, pk - y(end), ...
                y(1), mean(abs(y(b <= 0.3) - b(b <= 0.3)))}; %#ok<SAGROW>
        end
    end
end
if nMissing > 0
    warning("foldOnset:missing", "%d (pipeline, alpha, sigma) curves had < 3 points and were skipped", nMissing);
end
R = cell2table(rows, "VariableNames", ["pipeline" "alpha" "sigma" "nPts" "lowDescent" ...
    "foldDescent" "nDescents" "peakBetaGen" "peakBetaRec" "endDrop" "betaRecAt0" "idErrBelow0p3"]);

%% Onset sigma per pipeline x alpha (first sigma > 0 where each feature appears)
onsetOf = @(x, sg) min([sg(x & sg > 0); NaN]);
fprintf("\nFirst sigma (mm) at which each feature appears (NaN = not within %g mm)\n", SIGMA_MAX);
fprintf("%-10s %5s %10s %10s   %s\n", "pipeline", "alpha", "low-end", "high fold", "sigma=0 descents");
for p = pipes'
    for a = ALPHAS
        r = R(R.pipeline == p & R.alpha == a, :);
        z = r(r.sigma == 0, :);
        fprintf("%-10s %5.1f %10.2f %10.2f   %d\n", p, a, ...
            onsetOf(r.lowDescent, r.sigma), onsetOf(r.foldDescent, r.sigma), ...
            sum(z.nDescents));
    end
end

%% Floor table: beta_rec at beta_gen = 0 (monotone lift of the low end)
fprintf("\nFloor beta_rec(beta_gen = 0), SG-IRLS and BWFD-OLS:\n");
for a = ALPHAS
    for s = FLOOR_SIGMA
        r = R(R.alpha == a & abs(R.sigma - s) < 1e-9 & ismember(R.pipeline, ["SG-IRLS" "BWFD-OLS"]), :);
        fprintf("  alpha %.1f, sigma %g mm: %s\n", a, s, strjoin(compose("%s %.3f", r.pipeline, r.betaRecAt0), ", "));
    end
end

%% Figure: curves by sigma for two pipelines, one row per alpha
figure("Color", "w", "Position", [80 80 900 1100]);
tl = tiledlayout(numel(ALPHAS_PLOT), numel(PLOT_PIPES), "TileSpacing", "compact");
cmap = parula(numel(sigNodes));
for ia = 1:numel(ALPHAS_PLOT)
    for ip = 1:numel(PLOT_PIPES)
        nexttile; hold on;
        plot([0 0.7], [0 0.7], "k--", "HandleVisibility", "off");
        for is = 1:numel(sigNodes)
            c = T(T.pipeline == PLOT_PIPES(ip) & T.alpha == nodeOf(ALPHAS_PLOT(ia)) & ...
                  abs(T.sigma - sigNodes(is)) < 1e-9, :);
            if height(c) < 3, continue; end
            c = sortrows(c, "betaGen");
            plot(c.betaGen, c.meanBetaRec, "-", "Color", cmap(is,:), "LineWidth", 1.2, ...
                "DisplayName", sprintf("\\sigma = %g", sigNodes(is)));
        end
        xline([LOW_MAX FOLD_MIN], ":", "HandleVisibility", "off");
        title(sprintf("%s, \\alpha = %g", PLOT_PIPES(ip), ALPHAS_PLOT(ia)));
        xlim([0 0.7]); ylim([0 0.7]); box on;
        if ia == 1 && ip == numel(PLOT_PIPES), legend("Location", "eastoutside"); end
    end
end
xlabel(tl, "\beta_{gen}"); ylabel(tl, "\beta_{rec} (mean of 5 replicates)");
title(tl, sprintf("Forward map by noise magnitude, fs = %d Hz, VGF node %.0f mm/s", FS, VGF0));
if ~isfolder(fileparts(OUT_PNG)), error("foldOnset:outDir", "%s", "Missing folder: " + fileparts(OUT_PNG)); end
exportgraphics(gcf, OUT_PNG, "Resolution", 150);
fprintf("\nFigure: %s\n", OUT_PNG);
save(OUT_MAT, "R", "ALPHAS", "SIGMA_MAX", "K_SE", "LOW_MAX", "FOLD_MIN", "FS", "VGF0", "-v7.3");
fprintf("Saved: %s\n", OUT_MAT);
