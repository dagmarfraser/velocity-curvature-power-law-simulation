% checkPencilFraction_v001  Do Apple Pencil x,y in the Pilot and Fraser
% StimuliApp enhanced CSVs carry true sub-pixel values, or are they snapped
% to a lattice? Companion to the ShapeMap iPhone finger check in data/field.
%
% Reads raw x,y straight from the -enhanced.csv files (columns x, y), BEFORE
% any importer scaling, negation or de-duplication.
%
% UNITS (decided from the data, fail loud if ambiguous)
%   The device is the iPad Pro 13" M4: 2752 x 2064 px, UIKit scale 2, so the
%   landscape point canvas is 1376 x 1032 pt. Step 1 scans x,y over EVERY
%   enhanced file of the dataset. If any value exceeds the point canvas the
%   units are pixels. If the corpus fills the point canvas (max x > 0.9*1376,
%   max y > 0.9*1032) and never exceeds it, the units are points (pixels would
%   need every touch in 6M+ samples to avoid the right and bottom halves of
%   the screen). Anything else errors.
%
% OUTPUT
%   results/checkPencilFraction_v001.csv  one row per dataset x unit x axis
%   results/checkPencilFraction_v001.png  50-bin histograms of mod(.,1)
%
% Scripts run from src/; paths resolve from this file's location.

clearvars
srcDir  = fileparts(mfilename("fullpath"));
rootDir = fileparts(srcDir);
resDir  = fullfile(rootDir, "results");

datasets   = ["pilot", "fraser"];
PT_W       = 1376;          % iPad Pro 13" M4 landscape width, points
PT_H       = 1032;          % height, points
UI_SCALE   = 2;             % px per pt
N_PER_DS   = 6;             % participant folders sampled for the histograms
TOL        = 1e-3;
N_BINS     = 50;

rows = {};
fig  = figure("Color", "w", "Position", [100 100 1200 700]);
tl   = tiledlayout(fig, numel(datasets), 4, "TileSpacing", "compact");

for iD = 1:numel(datasets)
    ds    = datasets(iD);
    files = dir(fullfile(rootDir, "data", ds, "automatic", "ID*", "*tracingSection-enhanced.csv"));
    if isempty(files)
        error("checkPencilFraction:noFiles", "%s", "No enhanced CSVs found for " + ds);
    end

    %% Step 1: units from the full-corpus range
    mx = -Inf; my = -Inf; mnx = Inf; mny = Inf; nAll = 0;
    for k = 1:numel(files)
        [x, y] = readXY(fullfile(files(k).folder, files(k).name));
        mx = max(mx, max(x)); my = max(my, max(y));
        mnx = min(mnx, min(x)); mny = min(mny, min(y));
        nAll = nAll + numel(x);
    end
    fprintf("%s: %d files, %d samples, x [%.1f, %.1f], y [%.1f, %.1f]\n", ...
        ds, numel(files), nAll, mnx, mx, mny, my);
    if mx > PT_W || my > PT_H
        units = "pixels";
    elseif mx > 0.9*PT_W && my > 0.9*PT_H
        units = "points";
    else
        error("checkPencilFraction:ambiguousUnits", "%s", ...
            sprintf("%s: range x<=%.1f, y<=%.1f neither exceeds nor fills the point canvas", ds, mx, my));
    end
    fprintf("  units: %s\n", units);
    if units == "points"
        toPx = UI_SCALE;
    else
        toPx = 1;
    end

    %% Step 2: fractional parts on a spread of participants
    subj   = unique(string({files.folder}));
    pick   = subj(unique(round(linspace(1, numel(subj), N_PER_DS))));
    xs = []; ys = []; used = strings(0);
    for s = pick
        f = files(string({files.folder}) == s);
        [x, y] = readXY(fullfile(f(1).folder, f(1).name));
        xs = [xs; x]; ys = [ys; y]; %#ok<AGROW>
        used(end+1) = f(1).name; %#ok<AGROW>
    end
    fprintf("  files used (first session of %d participants):\n", numel(pick));
    fprintf("    %s\n", used);

    xPx = xs*toPx; yPx = ys*toPx;
    combos = {"points", xPx/UI_SCALE, yPx/UI_SCALE; "pixels", xPx, yPx};
    iTile = 0;
    for c = 1:2
        for a = 1:2
            v    = combos{c, a+1};
            fr   = mod(v, 1);
            fInt = mean(min(fr, 1-fr) < TOL);
            fHalf = mean(abs(fr - 0.5) < TOL);
            cnt  = histcounts(fr, linspace(0, 1, N_BINS+1));
            share = cnt / sum(cnt);
            ksD  = ksUniform(fr);
            axisName = ["x", "y"];
            rows(end+1, :) = {ds, units, combos{c,1}, axisName(a), numel(v), ...
                fInt, fHalf, max(share), min(share), nnz(cnt), ksD}; %#ok<AGROW>
            iTile = iTile + 1;
            ax = nexttile(tl);
            bar(ax, linspace(0.01, 0.99, N_BINS), share, 1);
            title(ax, sprintf("%s %s, %s", ds, axisName(a), combos{c,1}), "Interpreter", "none");
            xlabel(ax, "mod(value, 1)"); ylabel(ax, "share");
        end
    end
end

T = cell2table(rows, "VariableNames", ["dataset", "rawUnits", "evaluatedIn", "axis", ...
    "N", "fracWithinTolOfInteger", "fracWithinTolOfHalf", "maxBinShare", ...
    "minBinShare", "nonEmptyBins", "ksD_vsUniform"]);
disp(T)
writetable(T, fullfile(resDir, "checkPencilFraction_v001.csv"));
exportgraphics(fig, fullfile(resDir, "checkPencilFraction_v001.png"), "Resolution", 150);

function [x, y] = readXY(fn)
% Raw x,y columns, no scaling, no de-duplication.
opts = detectImportOptions(fn);
opts.SelectedVariableNames = {'x', 'y'};
t = readtable(fn, opts);
x = double(t.x); y = double(t.y);
if any(~isfinite(x)) || any(~isfinite(y))
    error("checkPencilFraction:nonFinite", "%s", "Non-finite x/y in " + fn);
end
end

function D = ksUniform(u)
% One-sample KS statistic against U(0,1), ties included.
u = sort(u(:)); n = numel(u);
D = max(max((1:n)'/n - u), max(u - (0:n-1)'/n));
end
