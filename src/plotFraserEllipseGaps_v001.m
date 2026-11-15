function gapTbl = plotFraserEllipseGaps_v001(opts)
%PLOTFRASERELLIPSEGAPS_V001  Show every Fraser N=48 Ellipse trial with a pen-lift gap.
%
%   Diagnostic for the N=48 import audit (Step 1, 2026-11-03). Reads the raw
%   enhanced CSVs directly (not importDB_fraser_v001), so it shows what the
%   importer's uniform 240 Hz clock hides.
%
%   A trial is flagged when, after sorting by touch_index and clipping
%   EdgeClip samples from each end (as importDB_fraser_v001 does), any step
%   in time_seconds exceeds GapThresh. Expected: 42 trials (audit2.py).
%
%   Each panel: x,y in mm (5.209 UIKit points per mm, y up), Pencil samples
%   grey, the samples either side of each gap joined in red, direct-touch
%   (touch_type 0) samples, if any, as orange dots. Title: participant,
%   session stem, trial, largest gap.
%
%   USAGE
%       gapTbl = plotFraserEllipseGaps_v001();
%       gapTbl = plotFraserEllipseGaps_v001(GapThresh=0.02);
%
%   Writes figures/plotFraserEllipseGaps_v001.png. Returns one row per flagged trial.
%
%   Fraser, D.S. (2026)  v001

arguments
    opts.DataDir   (1,:) char   = defaultDataDir_local()
    opts.Exclude   (1,:) string = ["ID101","ID102","ID214"]
    opts.EdgeClip  (1,1) double = 50
    opts.GapThresh (1,1) double = 0.05      % seconds
    opts.PtPerMm   (1,1) double = 5.209     % UIKit points per mm, M4 13" iPad Pro
    opts.Save      (1,1) logical = true
end

if ~isfolder(opts.DataDir)
    error('plotFraserEllipseGaps:noDir', '%s', ['Data folder not found: ' opts.DataDir]);
end

d = dir(fullfile(opts.DataDir, 'ID*'));
d = d([d.isdir] & ~ismember(string({d.name}), opts.Exclude));
ids = sort(string({d.name}));

rows = {};          % flagged trial records
traj = {};          % their samples, for plotting
nEll = 0;

for s = ids
    subjDir = fullfile(opts.DataDir, s);
    metas = dir(fullfile(subjDir, 'betaTestv002 * tracingSection.csv'));
    for m = 1:numel(metas)
        stem = erase(string(metas(m).name), " tracingSection.csv");
        meta = readtable(fullfile(subjDir, metas(m).name));
        enhFile = fullfile(subjDir, stem + " tracingSection-enhanced.csv");
        if ~isfile(enhFile)
            error('plotFraserEllipseGaps:noEnhanced', '%s', ['Missing: ' char(enhFile)]);
        end
        enh = readtable(enhFile);
        ellTrials = meta.trial(meta.tracingScene_object1_imageNumber == 3);
        for tr = ellTrials'
            T = enh(enh.response_index == tr, :);
            if isempty(T), continue, end           % the known 1-trial-short case
            nEll = nEll + 1;
            T = sortrows(T, 'touch_index');
            N = height(T);
            if N <= 2*opts.EdgeClip, continue, end
            T = T(opts.EdgeClip+1 : N-opts.EdgeClip, :);
            dt = diff(T.time_seconds);
            iGap = find(dt > opts.GapThresh);
            if isempty(iGap), continue, end
            % Longest gap-free section: runs between gaps, in samples and seconds
            edges = [0; iGap; height(T)];
            runLen = diff(edges);
            [longN, j] = max(runLen);
            longS = T.time_seconds(edges(j+1)) - T.time_seconds(edges(j)+1);
            rows(end+1, :) = {s, stem, tr, height(T), max(dt), numel(iGap), ...
                sum(T.touch_type == 0), longN, longS, longN/height(T)}; %#ok<AGROW>
            traj{end+1} = struct('x', T.x / opts.PtPerMm, 'y', -T.y / opts.PtPerMm, ...
                'type', T.touch_type, 'iGap', iGap); %#ok<AGROW>
        end
    end
end

gapTbl = cell2table(rows, 'VariableNames', ...
    {'participant','stem','trial','nSamples','maxGap_s','nGaps','nDirectTouch', ...
     'longestRun_samples','longestRun_s','longestRun_frac'});
gapTbl = sortrows(gapTbl, 'maxGap_s', 'descend');
fprintf('[plotFraserEllipseGaps] %d of %d Ellipse trials (%d participants) have a gap > %.0f ms after a %d-sample clip.\n', ...
    height(gapTbl), nEll, numel(ids), 1000*opts.GapThresh, opts.EdgeClip);

if isempty(rows), return, end

% Plot in descending order of largest gap
[~, ord] = sort(cell2mat(rows(:,5)), 'descend');
nT = numel(ord);
nC = ceil(sqrt(nT * 4/3));
nR = ceil(nT / nC);
fig = figure('Color', 'w', 'Position', [50 50 260*nC 220*nR]);
tl = tiledlayout(fig, nR, nC, 'TileSpacing', 'compact', 'Padding', 'compact');
for k = 1:nT
    r = rows(ord(k), :); q = traj{ord(k)};
    ax = nexttile(tl); hold(ax, 'on');
    pen = q.type == 2;
    plot(ax, q.x(pen), q.y(pen), '.', 'Color', [0.55 0.55 0.55], 'MarkerSize', 3);
    for g = q.iGap'
        plot(ax, q.x(g:g+1), q.y(g:g+1), 'r-', 'LineWidth', 1.5);
        plot(ax, q.x(g:g+1), q.y(g:g+1), 'ro', 'MarkerSize', 4, 'MarkerFaceColor', 'r');
    end
    if any(~pen)
        plot(ax, q.x(~pen), q.y(~pen), '.', 'Color', [1 0.5 0], 'MarkerSize', 8);
    end
    axis(ax, 'equal'); box(ax, 'on'); ax.FontSize = 7;
    title(ax, sprintf('%s %s T%d\nmax gap %.3f s (%d)', r{1}, extractAfter(r{2}, "betaTestv002 "), ...
        r{3}, r{5}, r{6}), 'Interpreter', 'none', 'FontSize', 7);
end
title(tl, sprintf('Fraser N=48 Ellipse trials with a time gap > %.0f ms after %d-sample clip (%d of %d); mm at %.3f pt/mm; red = gap, orange = direct touch', ...
    1000*opts.GapThresh, opts.EdgeClip, nT, nEll, opts.PtPerMm), 'Interpreter', 'none', 'FontSize', 10);

if opts.Save
    figDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'figures');
    out = fullfile(figDir, 'plotFraserEllipseGaps_v001.png');
    exportgraphics(fig, out, 'Resolution', 150);
    fprintf('Saved %s\n', out);
end
end

function dd = defaultDataDir_local()
    dd = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'data', 'fraser', 'automatic');
end
