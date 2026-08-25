function plotBetaRecovery_v005(betaGen)
% plotBetaRecovery_v005  Beta recovery vs sigma for all 6 pipelines.
%
% v005: extended to v058 grid (alpha 0-6, sigma 0-20 mm).
% One subplot per alpha level, arranged in colour-band rows:
%   Row 1: alpha 0.0-0.8  (white/near-white)
%   Row 2: alpha 1.0-1.8  (pink)
%   Row 3: alpha 2.0-2.8  (red/brown)
%   Row 4: alpha 3.0-3.8  (dark red)
%   Row 5: alpha 4.0-4.8  (very dark)
%   Row 6: alpha 5.0-5.8  (near-black)
%   Row 7: alpha 6.0      (black) + legend
%
% fs fixed at 120 Hz (mid), VGF fixed at exp(5.1) ~164 (mid of 14 values).
% betaGen snaps to the nearest simulation grid value.
%
% USAGE:
%   plotBetaRecovery_v005           % default beta_gen = 1/3
%   plotBetaRecovery_v005(0.43)     % beta_gen nearest to 0.43
%
% Must be run (or called) from src/.
% Fraser, D.S. (2026)  v005

arguments
    betaGen (1,1) double = 1/3
end

srcDir = fileparts(mfilename("fullpath"));
addpath(genpath(fullfile(srcDir, "functions")));

targetFs  = 120;
targetVGF = exp(5.1);

%% ========== LOAD ==========
matFile = fullfile(srcDir, "perCoordinateSEM_v2_001.mat");
if ~isfile(matFile)
    error("plotBetaRecovery:matNotFound", ...
        "perCoordinateSEM_v2_001.mat not found at %s. Pull from BlueBear RDS first.", ...
        matFile);
end
load(matFile, "coordTable");
T = coordTable;

%% ========== SNAP TO GRID ==========
% Nearest-match for values close to a grid point (handles 1/3 stored as
% 0.333330 due to FP accumulation). Snaps upward only when betaGen falls
% clearly between two grid points.
allBetas = unique(T.betaGen);
allVGFs  = unique(T.VGF);

[minDist, nearestIdx] = min(abs(allBetas - betaGen));
if minDist < 0.01
    useBeta = allBetas(nearestIdx);
else
    above   = allBetas(allBetas >= betaGen);
    useBeta = above(1);
end

[~, vIdx] = min(abs(allVGFs - targetVGF));
useVGF    = allVGFs(vIdx);

mask = abs(T.betaGen - useBeta) < 1e-9 & ...
       T.fs == targetFs                 & ...
       abs(T.VGF - useVGF) < 0.1;
T = T(mask, :);

fprintf("beta_gen = %.6f  (requested %.6f)\n", useBeta, betaGen);
fprintf("VGF      = %.4f  fs = %d Hz  -> %d rows\n\n", useVGF, targetFs, height(T));

if height(T) == 0
    error("plotBetaRecovery:noData", ...
        "No rows after filtering — betaGen=%.4f, fs=%d, VGF=%.2f.", ...
        betaGen, targetFs, targetVGF);
end

%% ========== LAYOUT ==========
nCols     = 5;
rowAlphas = {
    0.0:0.2:0.8, ...
    1.0:0.2:1.8, ...
    2.0:0.2:2.8, ...
    3.0:0.2:3.8, ...
    4.0:0.2:4.8, ...
    5.0:0.2:5.8, ...
    6.0 ...
};
rowLabels = {
    "White  (\alpha = 0 - 0.8)",
    "Pink   (\alpha = 1.0 - 1.8)",
    "Red    (\alpha = 2.0 - 2.8)",
    "Dark   (\alpha = 3.0 - 3.8)",
    "Deeper (\alpha = 4.0 - 4.8)",
    "Near-black (\alpha = 5.0 - 5.8)",
    "Black  (\alpha = 6.0)"
};
nRows = numel(rowAlphas);

%% ========== PIPELINE STYLE ==========
pipeNames  = ["BWFD-OLS", "BWFD-LMLS", "BWFD-IRLS", "SG-OLS", "SG-LMLS", "SG-IRLS"];
pipeColors = [
    0.85 0.33 0.10
    0.93 0.69 0.13
    0.49 0.18 0.56
    0.47 0.67 0.19
    0.30 0.75 0.93
    0.00 0.45 0.74
];
pipeMarkers = ["o", "s", "^", "o", "s", "^"];
pipeStyle   = ["-", "-", "-", "--", "--", "--"];

drawRef13 = abs(useBeta - 1/3) > 0.01;
betaLabel = sprintf("\\beta_{gen} = %.4f", useBeta);

%% ========== FIGURE ==========
fig = figure("Position", [40 40 1550 1400], "Color", "w");
tiledlayout(nRows, nCols, "TileSpacing", "compact", "Padding", "compact");
allA = unique(T.alpha);

for ri = 1:nRows
    alphasInRow = rowAlphas{ri};
    for ci = 1:nCols
        ax = nexttile;

        if ci <= numel(alphasInRow)
            [~, ai] = min(abs(allA - alphasInRow(ci)));
            useAlpha = allA(ai);
            sel      = T(abs(T.alpha - useAlpha) < 1e-9, :);
            hold(ax, "on");

            for pi = 1:numel(pipeNames)
                pData = sortrows(sel(sel.pipeline == pipeNames(pi), :), "sigma");
                if isempty(pData), continue; end
                plot(ax, pData.sigma, pData.meanBetaRec, ...
                    pipeStyle(pi) + pipeMarkers(pi), ...
                    "Color",           pipeColors(pi, :), ...
                    "LineWidth",       1.5, ...
                    "MarkerSize",      4, ...
                    "MarkerFaceColor", pipeColors(pi, :));
            end

            yline(ax, useBeta, "k:",  "LineWidth", 1.1, "Alpha", 0.7);
            if drawRef13
                yline(ax, 1/3, "k--", "LineWidth", 0.8, "Alpha", 0.45);
            end

            xlim(ax, [-0.4 21]);
            ylim(ax, [0 0.72]);
            title(ax, sprintf("\\alpha = %.1f", alphasInRow(ci)), "FontSize", 9);
            xlabel(ax, "\sigma (mm)", "FontSize", 8);
            if ci == 1, ylabel(ax, "\beta_{rec}", "FontSize", 8); end
            set(ax, "FontSize", 8, ...
                "XTick",      [0 5 10 15 20], ...
                "YTick",      uniquetol([0, 1/3, useBeta, 2/3], 0.02), ...
                "TickLength", [0.02 0.02], ...
                "Box",        "on");
            disableDefaultInteractivity(ax);
            grid(ax, "on");

            if ci == 1
                text(ax, -0.32, 0.5, rowLabels{ri}, ...
                    "Units",               "normalized", ...
                    "Rotation",            90, ...
                    "HorizontalAlignment", "center", ...
                    "VerticalAlignment",   "bottom", ...
                    "FontSize",            8, ...
                    "FontWeight",          "bold", ...
                    "Color",               [0.3 0.3 0.3]);
            end

        else
            axis(ax, "off");
            if ci - numel(alphasInRow) == 1
                hold(ax, "on");
                for pi = 1:numel(pipeNames)
                    plot(ax, NaN, NaN, pipeStyle(pi) + pipeMarkers(pi), ...
                        "Color",           pipeColors(pi, :), ...
                        "LineWidth",       1.5, ...
                        "MarkerFaceColor", pipeColors(pi, :), ...
                        "DisplayName",     pipeNames(pi));
                end
                plot(ax, NaN, NaN, "k:", "LineWidth", 1.1, "DisplayName", betaLabel);
                if drawRef13
                    plot(ax, NaN, NaN, "k--", "LineWidth", 0.8, "DisplayName", "\beta = 1/3 ref");
                end
                legend(ax, "show", "Location", "west", "FontSize", 9, "Box", "off");
            end
        end
    end
end

sgtitle(fig, ...
    sprintf("\\beta_{rec} vs \\sigma  |  %s,  fs = %d Hz,  VGF = exp(5.1) \\approx %.0f", ...
        betaLabel, targetFs, useVGF), ...
    "FontSize", 12, "FontWeight", "bold");

%% ========== SAVE ==========
outDir  = fullfile(srcDir, "..", "figures");
if ~exist(outDir, "dir"), mkdir(outDir); end
betaTag = strrep(sprintf("%.4f", useBeta), ".", "p");
outFile = fullfile(outDir, ...
    sprintf("betaRecovery_bg%s_fs%d_midVGF_v005.png", betaTag, targetFs));
exportgraphics(fig, outFile, "Resolution", 200);
fprintf("Saved: %s\n", outFile);

end
