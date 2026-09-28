function plotPaperRoadmapSchematic_v004()
% PLOTPAPERROADMAPSCHEMATIC_V004  Fig 0 for the v003 skeleton, redrawn from scratch.
%
% v004 (2026-09-28, Session 125). New drawing, same verified Panel B data path.
%   Panel A (new, v003 skeleton spec L292-297): one trial -> the two pulls on
%     the measured exponent, with tempo deciding which dominates -> beta_obs ->
%     pointwise inversion at the trial's own tempo and noise -> beta_gen*, with
%     precision and identifiability as a side branch pointing into Panel B.
%     Stages 2 and 3 are SCHEMATIC (labelled so on the figure): curve shapes are
%     illustrative, anchored only to published constants (filter-pull onset
%     f0 ~ 0.18 f_c, Finding #233; 1/3 as the filter pull's target, #225/#234;
%     gain < 1 limiting per-trial precision, #239). No data are plotted in A.
%   Panel B (same data as v003, redrawn): v003 squeezed SEM and coverage into a
%     jittered 2 x 3 grid; v004 plots both as continuous axes. x = SEM at each
%     dataset's (alpha, sigma, fs) centroid (canonical lookup, unchanged from
%     v003 L126-142); y = per-trial rise coverage (compare42CellVerdict method,
%     unchanged from v003 L166-195), plotted at N_REPS = 200 (v015) only, at
%     Dagmar's request (Session 125: showing both was too busy). N_REPS = 20
%     (v012, first reported) is still loaded, solely to self-check #160's
%     4/6/32. Neither replicate count is pre-registered: pointwise inversion
%     postdates registration.
%   Self-checks are now ERRORS, not warnings (v003 only warned): published
%     SG-IRLS SEMs; 4/6/32 at N_REPS = 20 (#160); 8/2/32 at N_REPS = 200 (#223);
%     36/42 SEM-adequate; 30/36 within the validity domain (Zarandi excluded).
%
% OUTPUT (run from src/):
%   ../figures/paperRoadmap_v004.png (600 dpi), ../figures/paperRoadmap_v004.pdf
%   (vector), ../figures/paperRoadmap_v004_preview.png (150 dpi),
%   ../results/paperRoadmap_v004_panelB.mat (the plotted grids).
%
% Fraser, D.S. (2026)

srcDir  = fileparts(mfilename('fullpath'));
rootDir = fileparts(srcDir);

%% CONFIG ------------------------------------------------------------------
PIPELINES = ["BWFD-OLS","SG-OLS","BWFD-LMLS","SG-LMLS","BWFD-IRLS","SG-IRLS"];
OI = struct('blue',[0 114 178]/255, 'sky',[86 180 233]/255, 'green',[0 158 115]/255, ...
    'orange',[230 159 0]/255, 'verm',[213 94 0]/255, 'purple',[204 121 167]/255, ...
    'grey',[0.62 0.62 0.62], 'ink',[0.13 0.13 0.13]);
REGISTRY = struct( ...
    'name',     {'Fraser','Cook CTRL','Cook ASD','Hickman PLAC','Hickman HALO','Dhieb','Zarandi'}, ...
    'matTag',   {'Fraser','Cook_CTRL','Cook_ASD','Hickman_PLAC','Hickman_HALO','Dhieb','Zarandi'}, ...
    'alpha',    {4.289,    4.77,       5.06,      5.34,          5.42,          2.50,   3.18}, ...
    'sigma',    {2.009,    8.15,       7.84,      7.17,          7.42,          7.50,   4.77}, ...
    'fs',       {240,      133,        133,       133,           133,           100,    100}, ...
    'col',      {OI.green, OI.sky,     OI.blue,   OI.orange,     OI.verm,       OI.purple, OI.grey}, ...
    'inDomain', {true,     true,       true,      true,          true,          true,   false});
PUBLISHED_SG_IRLS = containers.Map( ...
    {'Zarandi','Cook CTRL','Cook ASD','Dhieb','Hickman PLAC','Hickman HALO','Fraser'}, ...
    {0.0079,    0.0051,     0.0046,    0.0189, 0.0040,        0.0040,        0.0032});
MDC = 0.03;  SEM_ADEQUATE = MDC / 2.77;  PASS_T = 0.95;  COND_T = 0.90;
EXPECT20 = [4 6 32];  EXPECT200 = [8 2 32];  EXPECT_ADEQ = 36;  EXPECT_ADEQ_DOM = 30;
FIG_W = 18;  FIG_H = 15.5;   % cm; BRM double-column width

%% Panel B data -------------------------------------------------------------
semGrid = semLookup_local(srcDir, REGISTRY, PIPELINES);
pI = find(PIPELINES == "SG-IRLS", 1);
for d = 1:numel(REGISTRY)
    ref = PUBLISHED_SG_IRLS(REGISTRY(d).name);
    if abs(semGrid(d,pI) - ref) > 5e-4
        error('roadmap:SEMRepro', '%s', sprintf('SG-IRLS SEM for %s is %.5f, not the published %.4f.', ...
            REGISTRY(d).name, semGrid(d,pI), ref));
    end
end
[cov20,  ver20]  = coverage_local(srcDir, REGISTRY, PIPELINES, "v012", PASS_T, COND_T);
[cov200, ver200] = coverage_local(srcDir, REGISTRY, PIPELINES, "v015", PASS_T, COND_T);
split = @(v) [sum(v(:)=="PASS") sum(v(:)=="CONDITIONAL") sum(v(:)=="FAIL")];
if ~isequal(split(ver20), EXPECT20)
    error('roadmap:Split20', '%s', sprintf('N_REPS=20 split %s, expected %s (#160).', mat2str(split(ver20)), mat2str(EXPECT20)));
end
if ~isequal(split(ver200), EXPECT200)
    error('roadmap:Split200', '%s', sprintf('N_REPS=200 split %s, expected %s.', mat2str(split(ver200)), mat2str(EXPECT200)));
end
adeq = semGrid < SEM_ADEQUATE;  inDom = [REGISTRY.inDomain];
if sum(adeq(:)) ~= EXPECT_ADEQ || sum(adeq(inDom,:), 'all') ~= EXPECT_ADEQ_DOM
    error('roadmap:Adequacy', '%s', sprintf('SEM-adequate %d/42 (expected %d), in domain %d/36 (expected %d).', ...
        sum(adeq(:)), EXPECT_ADEQ, sum(adeq(inDom,:),'all'), EXPECT_ADEQ_DOM));
end
fprintf('Self-checks passed: SG-IRLS SEMs; splits %s / %s; adequate %d/42, %d/36 in domain.\n', ...
    mat2str(split(ver20)), mat2str(split(ver200)), sum(adeq(:)), sum(adeq(inDom,:),'all'));
nPass20 = sum(ver20(:)=="PASS");  nPass200 = sum(ver200(:)=="PASS");
nPass200Dom = sum(ver200(inDom,:)=="PASS", 'all');

resDir = fullfile(rootDir, 'results');
if ~isfolder(resDir), error('roadmap:NoResults', '%s', sprintf('FAILED PATH: %s', resDir)); end
names = string({REGISTRY.name});
save(fullfile(resDir, 'paperRoadmap_v004_panelB.mat'), 'semGrid', 'cov20', 'cov200', 'ver20', 'ver200', ...
    'names', 'PIPELINES', 'SEM_ADEQUATE', 'PASS_T', 'COND_T');

%% Figure ---------------------------------------------------------------------
fig = figure('Units','centimeters', 'Position',[2 2 FIG_W FIG_H], 'Color','w', 'Visible','off', ...
    'PaperUnits','centimeters', 'PaperSize',[FIG_W FIG_H], 'PaperPosition',[0 0 FIG_W FIG_H]);
bg = axes(fig, 'Units','centimeters', 'Position',[0 0 FIG_W FIG_H], 'XLim',[0 FIG_W], 'YLim',[0 FIG_H], ...
    'Visible','off');
hold(bg, 'on');
drawPanelA_local(fig, bg, OI, SEM_ADEQUATE);
drawPanelB_local(fig, bg, OI, REGISTRY, PIPELINES, semGrid, cov20, cov200, SEM_ADEQUATE, PASS_T, COND_T, ...
    sum(adeq(:)), sum(adeq(inDom,:),'all'), nPass20, nPass200, nPass200Dom);

figDir = fullfile(rootDir, 'figures');
if ~isfolder(figDir), error('roadmap:NoFigures', '%s', sprintf('FAILED PATH: %s', figDir)); end
exportgraphics(fig, fullfile(figDir, 'paperRoadmap_v004.png'), 'Resolution', 600);
exportgraphics(fig, fullfile(figDir, 'paperRoadmap_v004.pdf'), 'ContentType', 'vector');
exportgraphics(fig, fullfile(figDir, 'paperRoadmap_v004_preview.png'), 'Resolution', 150);
close(fig);
fprintf('Saved paperRoadmap_v004.png/.pdf/_preview.png in %s\n', figDir);
end

%% ======================= data (logic unchanged from v003) ====================
function semGrid = semLookup_local(srcDir, REG, PIPES)
f = fullfile(srcDir, 'perCoordinateSEM_v2_001.mat');
if ~isfile(f), error('roadmap:NoSEM', '%s', sprintf('FAILED PATH: %s', f)); end
T = load(f, 'coordTable').coordTable;
aa = sort(unique(T.alpha)); ss = sort(unique(T.sigma)); ff = sort(unique(T.fs));
semGrid = NaN(numel(REG), numel(PIPES));
for d = 1:numel(REG)
    [~,ai] = min(abs(aa - REG(d).alpha)); [~,si] = min(abs(ss - REG(d).sigma)); [~,fi] = min(abs(ff - REG(d).fs));
    for p = 1:numel(PIPES)
        at = T.alpha==aa(ai) & T.sigma==ss(si) & T.fs==ff(fi) & T.pipeline==PIPES(p);
        if ~any(at)
            error('roadmap:NoSEMRow', '%s', sprintf('No SEM rows for %s / %s.', REG(d).name, PIPES(p)));
        end
        semGrid(d,p) = mean(T.sem(at), 'omitnan');   % mean over ALL betaGen and ALL VGF (canonical)
    end
end
end

function [covGrid, verGrid] = coverage_local(srcDir, REG, PIPES, ver, PASS_T, COND_T)
covGrid = NaN(numel(REG), numel(PIPES)); verGrid = strings(numel(REG), numel(PIPES));
for d = 1:numel(REG)
    f = fullfile(srcDir, sprintf('loopClosureResults_%s_all_shaped_xu_%s.mat', REG(d).matTag, ver));
    if ~isfile(f), error('roadmap:NoCorpus', '%s', sprintf('FAILED PATH: %s', f)); end
    S = load(f, 'results', 'pipelineLabels');
    if ~isfield(S, 'pipelineLabels') || ~isequal(string(S.pipelineLabels(:)).', PIPES)
        error('roadmap:PipeOrder', '%s', sprintf('%s: pipelineLabels missing or not in PIPELINES order.', f));
    end
    for p = 1:numel(PIPES)
        st = arrayfun(@(r) r.invertStatus(p), S.results);
        ok = st ~= "no_beta_obs";
        if ~any(ok), verGrid(d,p) = "no-data"; continue, end
        c = sum(st(ok) == "rise") / sum(ok);
        covGrid(d,p) = c;
        if c >= PASS_T, verGrid(d,p) = "PASS"; elseif c >= COND_T, verGrid(d,p) = "CONDITIONAL"; else, verGrid(d,p) = "FAIL"; end
    end
end
end

%% ============================ Panel A =========================================
function drawPanelA_local(fig, bg, OI, SEMT)
FS = 7; ink = OI.ink; sub = [0.40 0.40 0.40];
text(bg, 0.30, 15.05, 'A', 'FontSize', 11, 'FontWeight', 'bold', 'Color', ink);
text(bg, 0.85, 15.05, 'From one trajectory to a generator estimate', 'FontSize', 8.5, 'FontWeight', 'bold', 'Color', ink);

hdr = {0.55, 'One trial',              'tempo f_0, noise (\alpha, \sigma)';
       4.95, 'Two pulls',              'tempo decides which dominates';
       9.85, 'Pointwise inversion',    'at the trial''s own f_0, \alpha, \sigma';
       14.2, 'Generator estimate',     'what the data can support'};
for k = 1:4
    text(bg, hdr{k,1}, 14.30, sprintf('%d  %s', k, hdr{k,2}), 'FontSize', FS+0.5, 'FontWeight', 'bold', 'Color', ink);
    text(bg, hdr{k,1}+0.30, 13.90, hdr{k,3}, 'FontSize', FS-0.5, 'Color', sub);
end

% --- 1: one trial (ellipse, speed-coloured, noisy samples) ---
a1 = axes(fig, 'Units','centimeters', 'Position',[0.55 10.35 3.3 3.2]); hold(a1,'on');
t = linspace(0, 2*pi, 400); A = 2; B = 1;
x = A*cos(t); y = B*sin(t);
kap = A*B ./ (A^2*sin(t).^2 + B^2*cos(t).^2).^1.5;
v = kap.^(-1/3); v = (v - min(v)) / (max(v) - min(v));
surface(a1, [x; x], [y; y], zeros(2, numel(t)), [v; v], 'EdgeColor','interp', 'FaceColor','none', 'LineWidth', 2.2);
cm = interp1([0 1], [0.80 0.86 0.93; 0.05 0.25 0.55], linspace(0,1,64)); colormap(a1, cm);
rng(7, 'twister'); ts = linspace(0, 2*pi, 46);
plot(a1, A*cos(ts) + 0.07*randn(size(ts)), B*sin(ts) + 0.07*randn(size(ts)), '.', 'Color', [0.45 0.45 0.45], 'MarkerSize', 4);
axis(a1, 'equal'); axis(a1, 'off'); xlim(a1, [-2.3 2.3]); ylim(a1, [-1.6 1.6]);
text(a1, 0, -1.55, 'slow in the bends, fast on the flats', 'FontSize', FS-1, 'Color', sub, 'HorizontalAlignment','center');

% --- 2: two pulls (schematic) ---
a2 = axes(fig, 'Units','centimeters', 'Position',[5.55 10.55 3.25 2.95]); hold(a2,'on');
bG = 0.22; bN = 0.05; xr = logspace(-2, log10(1.5), 300);
wN = 1 ./ (1 + (xr/0.05).^1.5); wF = 1 ./ (1 + (0.45./xr).^3);
bObs = bG + wN*(bN - bG) + wF*(1/3 - bG);
patch(a2, [0.07 0.26 0.26 0.07], [0 0 0.4 0.4], [0.93 0.95 0.93], 'EdgeColor','none');
plot(a2, xr, bG*ones(size(xr)), '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.8);
plot(a2, xr, (1/3)*ones(size(xr)), ':', 'Color', OI.blue, 'LineWidth', 1.1);
plot(a2, xr, bN*ones(size(xr)), ':', 'Color', OI.verm, 'LineWidth', 1.1);
xline(a2, 0.18, '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.6);
plot(a2, xr, bObs, '-', 'Color', ink, 'LineWidth', 1.8);
x0 = 0.013; y0 = interp1(xr, bObs, x0);
plot(a2, [x0 x0], [bG-0.01 y0+0.022], '-', 'Color', OI.verm, 'LineWidth', 1.2);
plot(a2, x0, y0+0.018, 'v', 'MarkerFaceColor', OI.verm, 'MarkerEdgeColor', 'none', 'MarkerSize', 4);
x1 = 1.05; y1 = interp1(xr, bObs, x1);
plot(a2, [x1 x1], [bG+0.01 y1-0.018], '-', 'Color', OI.blue, 'LineWidth', 1.2);
plot(a2, x1, y1-0.016, '^', 'MarkerFaceColor', OI.blue, 'MarkerEdgeColor', 'none', 'MarkerSize', 4);
text(a2, 0.0165, 0.135, 'noise', 'FontSize', FS-1, 'Color', OI.verm);
text(a2, 0.60, 0.385, 'filter', 'FontSize', FS-1, 'Color', OI.blue, 'HorizontalAlignment','center');
text(a2, 0.125, 0.015, 'window', 'FontSize', FS-1.5, 'Color', [0.35 0.5 0.35], 'HorizontalAlignment','center');
text(a2, 1.65, 1/3, '1/3', 'FontSize', FS-1, 'Color', OI.blue);
text(a2, 1.65, bN, '\beta_{noise}', 'FontSize', FS-1, 'Color', OI.verm);
text(a2, 0.0115, 0.238, '\beta_{gen}', 'FontSize', FS-1, 'Color', [0.35 0.35 0.35]);
text(a2, 0.085, 0.262, '\beta_{obs}', 'FontSize', FS, 'Color', ink, 'FontWeight','bold');
text(a2, 0.18, 0.412, '0.18', 'FontSize', FS-1.5, 'Color', [0.55 0.55 0.55], 'HorizontalAlignment','center');
set(a2, 'XScale','log', 'XLim',[0.01 1.5], 'YLim',[0 0.4], 'XTick',[0.01 0.1 1], 'XTickLabel',{'0.01','0.1','1'}, ...
    'YTick',[0 0.2 0.4], 'FontSize', FS-1, 'TickDir','out', 'TickLength',[0.02 0.02], 'Box','off', ...
    'XColor', sub, 'YColor', sub, 'Layer','top');
xlabel(a2, 'tempo  f_0 / f_c', 'FontSize', FS-0.5); ylabel(a2, '\beta', 'FontSize', FS);
text(a2, 1.45, 0.018, 'schematic', 'FontSize', FS-1.5, 'FontAngle','italic', 'Color', sub, 'HorizontalAlignment','right');

% --- 3: pointwise inversion (schematic) ---
a3 = axes(fig, 'Units','centimeters', 'Position',[10.45 10.55 2.95 2.95]); hold(a3,'on');
Tt = 0.20; g = 0.50; bo = 0.25; db = 0.02; gs = Tt + (bo - Tt)/g; dgs = db/g;
xx = [0 0.6];
plot(a3, xx, xx, '--', 'Color', [0.78 0.78 0.78], 'LineWidth', 0.8);
patch(a3, [gs-dgs gs+dgs gs+dgs gs-dgs], [0 0 0.6 0.6], [0.88 0.93 0.97], 'EdgeColor','none', 'FaceAlpha', 0.9);
patch(a3, [0 0.6 0.6 0], [bo-db bo-db bo+db bo+db], [0.92 0.92 0.92], 'EdgeColor','none', 'FaceAlpha', 0.9);
plot(a3, xx, Tt + g*(xx - Tt), '-', 'Color', ink, 'LineWidth', 1.8);
plot(a3, [0 gs], [bo bo], ':', 'Color', ink, 'LineWidth', 1.0);
plot(a3, [gs gs], [bo 0.035], '-', 'Color', OI.blue, 'LineWidth', 1.2);
plot(a3, gs, 0.03, 'v', 'MarkerFaceColor', OI.blue, 'MarkerEdgeColor','none', 'MarkerSize', 4);
plot(a3, gs, bo, 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', ink, 'MarkerSize', 3.5);
text(a3, 0.02, 0.555, 'forward map, slope < 1', 'FontSize', FS-1.5, 'Color', sub);
text(a3, 0.012, bo+0.04, '\beta_{obs}', 'FontSize', FS-1, 'Color', ink);
text(a3, gs+0.035, 0.14, '\beta_{gen}^*', 'FontSize', FS-0.5, 'Color', OI.blue, 'FontWeight','bold');
set(a3, 'XLim',[0 0.6], 'YLim',[0 0.6], 'XTick',[0 0.3 0.6], 'YTick',[0 0.3 0.6], 'FontSize', FS-1, ...
    'TickDir','out', 'TickLength',[0.02 0.02], 'Box','off', 'XColor', sub, 'YColor', sub, 'Layer','top');
xlabel(a3, '\beta_{gen}', 'FontSize', FS); ylabel(a3, '\beta_{rec}', 'FontSize', FS);
text(a3, 0.59, 0.02, 'schematic', 'FontSize', FS-1.5, 'FontAngle','italic', 'Color', sub, 'HorizontalAlignment','right');

% --- 4: generator estimate box ---
rectangle(bg, 'Position',[14.25 10.75 3.35 2.75], 'Curvature', 0.12, 'EdgeColor', ink, 'LineWidth', 0.9, ...
    'FaceColor', [0.97 0.98 1.0]);
text(bg, 15.925, 12.95, '\beta_{gen}^*', 'FontSize', 13, 'Color', OI.blue, 'HorizontalAlignment','center', 'FontWeight','bold');
text(bg, 14.45, 12.10, {'one trial: too imprecise', 'to correct on its own'}, 'FontSize', FS-0.5, 'Color', ink, 'VerticalAlignment','middle');
text(bg, 14.45, 11.30, {'a design: pooled and', 'contrasted (Parts 2-4)'}, 'FontSize', FS-0.5, 'Color', ink, 'VerticalAlignment','middle');

% --- flow arrows ---
arrow_local(bg, 3.95, 12.05, 4.85, 12.05, ink, 1.1);
arrow_local(bg, 9.00, 12.05, 9.85, 12.05, ink, 1.1);
text(bg, 9.42, 12.38, '\beta_{obs}', 'FontSize', FS, 'Color', ink, 'HorizontalAlignment','center');
arrow_local(bg, 13.55, 12.05, 14.15, 12.05, ink, 1.1);

% --- side branch: precision and identifiability (into Panel B) ---
bx = {5.55, 'Precision', sprintf('SEM < MDC/2.77 = %.4f', SEMT), 'Panel B, across';
      9.75, 'Identifiability', 'locally invertible forward map', 'Panel B, up'};
for k = 1:2
    rectangle(bg, 'Position',[bx{k,1} 8.30 3.75 1.05], 'Curvature', 0.15, 'EdgeColor', [0.45 0.45 0.45], ...
        'LineStyle', '--', 'LineWidth', 0.7, 'FaceColor', 'w');
    text(bg, bx{k,1}+0.18, 9.08, bx{k,2}, 'FontSize', FS, 'FontWeight','bold', 'Color', ink);
    text(bg, bx{k,1}+0.18, 8.72, bx{k,3}, 'FontSize', FS-1, 'Color', ink);
    text(bg, bx{k,1}+3.60, 9.08, bx{k,4}, 'FontSize', FS-1.5, 'Color', sub, 'HorizontalAlignment','right', 'FontAngle','italic');
end
plot(bg, [9.42 9.42], [11.85 9.75], '--', 'Color', [0.55 0.55 0.55], 'LineWidth', 0.7);
plot(bg, [7.42 11.62], [9.75 9.75], '--', 'Color', [0.55 0.55 0.55], 'LineWidth', 0.7);
plot(bg, [7.42 7.42], [9.75 9.40], '--', 'Color', [0.55 0.55 0.55], 'LineWidth', 0.7);
plot(bg, [11.62 11.62], [9.75 9.40], '--', 'Color', [0.55 0.55 0.55], 'LineWidth', 0.7);
text(bg, 13.75, 8.82, {'answered separately:', 'a precise estimate can', 'still be uninvertible'}, ...
    'FontSize', FS-1.5, 'Color', sub, 'VerticalAlignment','middle');
end

%% ============================ Panel B =========================================
function drawPanelB_local(fig, bg, OI, REG, PIPES, semGrid, ~, cov200, SEMT, PASS_T, COND_T, ...
    nAdeq, nAdeqDom, ~, nPass200, nPass200Dom)
FS = 7; ink = OI.ink; sub = [0.40 0.40 0.40];
text(bg, 0.30, 7.55, 'B', 'FontSize', 11, 'FontWeight', 'bold', 'Color', ink);
text(bg, 0.85, 7.55, 'Precision and identifiability are separate questions: all 42 dataset \times pipeline cells', ...
    'FontSize', 8.5, 'FontWeight', 'bold', 'Color', ink);

ax = axes(fig, 'Units','centimeters', 'Position',[1.45 1.05 9.9 5.75]); hold(ax,'on');
XL = [0.0022 0.024]; YL = [0 103];
patch(ax, [XL(1) SEMT SEMT XL(1)], [100*PASS_T 100*PASS_T YL(2) YL(2)], [0.90 0.96 0.91], 'EdgeColor','none');
text(ax, XL(1)*1.04, 101.2, 'precise and identifiable', 'FontSize', FS-1.5, 'Color', [0.25 0.45 0.30]);
xline(ax, SEMT, '--', 'Color', [0.35 0.35 0.35], 'LineWidth', 0.8);
yline(ax, 100*PASS_T, ':', 'Color', [0.55 0.55 0.55], 'LineWidth', 0.7);
yline(ax, 100*COND_T, ':', 'Color', [0.55 0.55 0.55], 'LineWidth', 0.7);
text(ax, XL(2), 100*PASS_T+1.6, 'PASS ', 'FontSize', FS-1.5, 'Color', sub, 'HorizontalAlignment','right');
text(ax, XL(2), 100*COND_T+1.6, 'CONDITIONAL ', 'FontSize', FS-1.5, 'Color', sub, 'HorizontalAlignment','right');
text(ax, XL(2), 100*COND_T-2.6, 'FAIL ', 'FontSize', FS-1.5, 'Color', sub, 'HorizontalAlignment','right');
text(ax, 0.0059, 86.5, {'Zarandi: PASS, but', 'outside validity domain'}, 'FontSize', FS-1.5, ...
    'Color', [0.45 0.45 0.45], 'HorizontalAlignment','left');
plot(ax, [0.0068 0.0066], [97 89.8], '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.5);
text(ax, SEMT*1.035, 4, 'MDC/2.77', 'FontSize', FS-1.5, 'Color', sub, 'Rotation', 90);

[mk, filled] = markers_local(PIPES);
order = numel(REG):-1:1;   % draw Zarandi first so in-domain points sit on top
for d = order
    c = REG(d).col;
    for p = 1:numel(PIPES)
        xs = semGrid(d,p); y200 = 100*cov200(d,p);
        if isnan(y200)
            error('roadmap:NoCoverage', '%s', sprintf('%s / %s has no N_REPS=200 coverage.', REG(d).name, PIPES(p)));
        end
        fc = c; if ~filled(p), fc = 'w'; end
        plot(ax, xs, y200, mk(p), 'MarkerSize', 5.5, 'MarkerEdgeColor', c*0.85, 'MarkerFaceColor', fc, 'LineWidth', 0.9);
    end
end
set(ax, 'XScale','log', 'XLim', XL, 'YLim', YL, 'XTick', [0.003 0.005 SEMT 0.02], ...
    'XTickLabel', {'0.003','0.005',sprintf('%.4f',SEMT),'0.020'}, 'YTick', 0:20:100, ...
    'FontSize', FS-0.5, 'TickDir','out', 'TickLength',[0.012 0.012], 'Box','off', 'XColor', sub, 'YColor', sub, 'Layer','top');
xlabel(ax, 'SEM of \beta_{obs} at the dataset''s noise coordinate (log scale)', 'FontSize', FS);
ylabel(ax, 'trials on the rising branch (%, N_{REPS} = 200)', 'FontSize', FS);

% --- key ---
kx = 12.05; ky = 6.75; dy = 0.40;
text(bg, kx, ky+0.05, 'Dataset', 'FontSize', FS, 'FontWeight','bold', 'Color', ink);
for d = 1:numel(REG)
    yy = ky - d*dy;
    rectangle(bg, 'Position', [kx yy-0.11 0.24 0.22], 'FaceColor', REG(d).col, 'EdgeColor', 'none');
    lbl = REG(d).name; if ~REG(d).inDomain, lbl = [lbl ' (outside validity domain; see note)']; end
    text(bg, kx+0.38, yy, lbl, 'FontSize', FS-1, 'Color', ink, 'Interpreter','none');
end
ky2 = ky - (numel(REG)+1)*dy - 0.05;
text(bg, kx, ky2, 'Pipeline', 'FontSize', FS, 'FontWeight','bold', 'Color', ink);
regs = ["OLS","LMLS","IRLS"]; rmk = ['o','s','^'];
for r = 1:3
    xk = kx + 0.15 + (r-1)*1.45;
    plot(bg, xk, ky2-dy, rmk(r), 'MarkerSize', 5, 'MarkerEdgeColor', ink, 'MarkerFaceColor', ink);
    plot(bg, xk+0.35, ky2-dy, rmk(r), 'MarkerSize', 5, 'MarkerEdgeColor', ink, 'MarkerFaceColor', 'w');
    text(bg, xk+0.62, ky2-dy, regs(r), 'FontSize', FS-1, 'Color', ink);
end
text(bg, kx, ky2-2*dy, 'filled SG, open BWFD', 'FontSize', FS-1.5, 'Color', sub);
note = {'\bfWhy Zarandi is set aside.\rm The simulation adds noise', ...
        'to a path timed by the power law. Checked trial by', ...
        'trial (Part 3), Zarandi''s published data behave instead', ...
        'as re-timed paths, and the export looks resampled.', ...
        'Its cells are shown, but its PASSes carry no claim.'};
text(bg, kx, ky2-2*dy-0.30, note, 'FontSize', FS-1.5, 'Color', sub, 'VerticalAlignment','top');

% --- counts (computed, not typed) ---
text(ax, XL(1)*1.04, 2.5, {sprintf('SEM-adequate: %d of 42 cells (%d of 36 in domain)', nAdeq, nAdeqDom), ...
    sprintf('PASS: %d of 42 (%d in domain, all Fraser''s)', nPass200, nPass200Dom)}, ...
    'FontSize', FS-1, 'Color', sub, 'VerticalAlignment','bottom');
end

%% ============================ helpers ========================================
function [mk, filled] = markers_local(PIPES)
mk = repmat('o', 1, numel(PIPES)); filled = false(1, numel(PIPES));
for p = 1:numel(PIPES)
    if endsWith(PIPES(p), "LMLS"), mk(p) = 's'; elseif endsWith(PIPES(p), "IRLS"), mk(p) = '^'; end
    filled(p) = startsWith(PIPES(p), "SG");
end
end

function arrow_local(ax, x1, y1, x2, y2, col, lw)
L = hypot(x2-x1, y2-y1); ux = (x2-x1)/L; uy = (y2-y1)/L; h = 0.20; w = 0.09;
bx = x2 - h*ux; by = y2 - h*uy;
plot(ax, [x1 bx], [y1 by], '-', 'Color', col, 'LineWidth', lw);
patch(ax, [x2 bx - w*uy bx + w*uy], [y2 by + w*ux by - w*ux], col, 'EdgeColor', 'none');
end
