function out = plotConstellationMetrics_v003(opts)
% plotConstellationMetrics_v003  Fig 7a-7d on the v015 basis (replaces plotConstellationMetrics_v002, which read constellationMetrics_v002.mat).
%
% Numbers (cccMed with bootstrap CI, loopCCC, Nvalid, r, MAE) are read from results/constellationMetrics_v005.mat (table T), the file behind Part 3's
% primary table. Panels 7c and 7d need the per-trial contrasts, which v005 does not store: they are rebuilt from loopClosureResults_<ds>_all_shaped_xu_v015.mat
% with v005's validity mask (finite betaObs and slice, finite betaGenStarMed, max|dObs| >= DEG, not all |dObs| < MDC), and the trial count must equal T.Nvalid
% or the function stops. No re-simulation.
%
%   Fig 7a  cccMed (primary, with 95% bootstrap CI) and loopCCC (secondary) per dataset, shaped_xu.
%   Fig 7b  cccMed for xu, shaped_xu and bootstrap (same noise models as the Supplement A three-model table).
%   Fig 7c  Pooled predicted vs observed inter-pipeline delta-beta (in-domain datasets; Dhieb and Zarandi excluded).
%   Fig 7d  Per-dataset tiles, coloured by contrast type, OLS line and slope; titles from T.
%
% USAGE   plotConstellationMetrics_v003                 % save to figures/
%         plotConstellationMetrics_v003(SaveFigs=false) % display only
% Writes figures/fig7{a,b,c,d}_v015.png and results/plotConstellationMetrics_v003.mat (the plotted numbers, for the ledger).
% Fraser, D.S. (2026)  v003
arguments
    opts.SaveFigs (1,1) logical = true
    opts.FigDir   (1,1) string  = ""
    opts.AxLim    (1,1) double  = 0.25
end
srcDir = fileparts(mfilename("fullpath"));  root = fileparts(srcDir);
if opts.FigDir == "", opts.FigDir = fullfile(root, "figures"); end
if opts.SaveFigs && ~isfolder(opts.FigDir), mkdir(opts.FigDir); end
mf = fullfile(root, "results", "constellationMetrics_v005.mat");
if ~isfile(mf), error("pcm3:mat", "%s", "FAILED PATH: " + mf); end
S = load(mf, "results");  T = S.results.T;  MDC = S.results.config.MDC;  DEG = S.results.config.DEG;
need = ["dataset" "noiseModel" "Nvalid" "cccMed" "cccMedLo" "cccMedHi" "loopCCC" "pearsonMed" "maeMed"];
if ~all(ismember(need, string(T.Properties.VariableNames))), error("pcm3:cols", "constellationMetrics_v005 T lacks a required column"); end
T.dataset = string(T.dataset);  T.noiseModel = string(T.noiseModel);
ds7   = ["Fraser" "Cook_CTRL" "Cook_ASD" "Hickman_HALO" "Hickman_PLAC" "Dhieb"];   % bars; Zarandi (1 valid trial) annotated, not drawn
lab   = ["Fraser" "Cook CTRL" "Cook ASD" "Hickman HALO" "Hickman PLAC" "Dhieb †"];
nDs   = numel(ds7);  xPos = 1:nDs;
row   = @(m, d) T(T.noiseModel == m & T.dataset == d, :);
pick  = @(m, f) arrayfun(@(d) row(m, d).(f), ds7);
sx = struct("cccMed", pick("shaped_xu", "cccMed"), "lo", pick("shaped_xu", "cccMedLo"), "hi", pick("shaped_xu", "cccMedHi"), ...
            "loop", pick("shaped_xu", "loopCCC"), "N", pick("shaped_xu", "Nvalid"));
if any(~isfinite(sx.cccMed)), error("pcm3:nan", "a shaped_xu cccMed is not finite"); end
zar = row("shaped_xu", "Zarandi");
colP = [0.18 0.45 0.73];  colS = [0.70 0.80 0.90];  colD = [0.85 0.40 0.20];

%% 7a
f1 = figure("Name", "Fig7a", "Position", [50 50 820 440]);  ax = axes(f1);  hold(ax, "on");  bw = 0.36;
for i = 1:nDs
    cP = colP;  cS = colS;  if i == nDs, cP = colD;  cS = min(colD + 0.25, 1); end
    b1 = bar(ax, xPos(i) - bw/2, sx.cccMed(i), bw, "FaceColor", cP, "EdgeColor", "none");
    b2 = bar(ax, xPos(i) + bw/2, sx.loop(i),   bw, "FaceColor", cS, "EdgeColor", "k", "LineWidth", 0.8);
    errorbar(ax, xPos(i) - bw/2, sx.cccMed(i), sx.cccMed(i) - sx.lo(i), sx.hi(i) - sx.cccMed(i), "k", "LineWidth", 1, "CapSize", 4);
    if i == 1, h1 = b1;  h2 = b2; end
end
for yl = [0.50 0.70 0.90]
    plot(ax, [0.3 nDs+0.7], [yl yl], "k:", "LineWidth", 0.8);  text(ax, nDs + 0.75, yl, sprintf("%.2f", yl), "FontSize", 8, "Color", [.4 .4 .4]);
end
text(ax, 0.4, 1.085, sprintf("Zarandi (outside the validity domain; %d valid trial)\ncccMed not defined; loopCCC = %.3f", zar.Nvalid, zar.loopCCC), ...
    "FontSize", 7.5, "Color", [.4 .4 .4], "VerticalAlignment", "top", "BackgroundColor", [1 1 1 .85], "EdgeColor", [.7 .7 .7]);
set(ax, "XTick", xPos, "XTickLabel", lab, "XTickLabelRotation", 0, "YLim", [0 1.1], "YTick", 0:0.2:1, "FontSize", 9, "Box", "off");
for i = 1:nDs, text(ax, xPos(i), -0.095, sprintf("N = %d", sx.N(i)), "FontSize", 8, "HorizontalAlignment", "center", "Color", [.35 .35 .35], "Clipping", "off"); end
ylabel(ax, "CCC"); title(ax, "Fig 7a: Constellation validity gradient (shaped\_xu, v015, N\_REPS = 200)", "FontSize", 11, "FontWeight", "bold");
legend(ax, [h1 h2], {"cccMed (primary: per-trial contrast CCC; bars = 95% bootstrap CI)", "loopCCC (secondary: pooled across pipelines)"}, "Location", "northeast", "FontSize", 8.5);
ax.Position = [0.08 0.20 0.88 0.72];  annotation(f1, "textbox", [0.08 0.01 0.88 0.05], "String", "† Dhieb is outside the validity domain.  N = trials with a valid contrast pattern.", "EdgeColor", "none", "FontSize", 9, "Color", [.3 .3 .3]);
save1(f1, opts, "fig7a_cccMed_v015.png");

%% 7b
f2 = figure("Name", "Fig7b", "Position", [50 50 820 440]);  ax2 = axes(f2);  hold(ax2, "on");  bwm = 0.26;
mods = ["xu" "shaped_xu" "bootstrap"];  cols = [0.80 0.35 0.35; 0.18 0.45 0.73; 0.25 0.65 0.35];  hh = gobjects(1, 3);  cm = zeros(3, nDs);
for m = 1:3
    cm(m, :) = pick(mods(m), "cccMed");  lo = pick(mods(m), "cccMedLo");  hi = pick(mods(m), "cccMedHi");  xo = xPos + (m - 2) * bwm;
    hh(m) = bar(ax2, xo, cm(m, :), bwm, "FaceColor", cols(m, :), "EdgeColor", "none");
    errorbar(ax2, xo, cm(m, :), cm(m, :) - lo, hi - cm(m, :), "k", "LineStyle", "none", "CapSize", 3);
end
plot(ax2, [0.3 nDs+0.7], [0 0], "k-", "LineWidth", 0.5);
set(ax2, "XTick", xPos, "XTickLabel", lab, "XTickLabelRotation", 15, "YLim", [-0.2 1.1], "FontSize", 10, "Box", "off");
ylabel(ax2, "cccMed (per-trial contrast CCC)"); title(ax2, "Fig 7b: Three noise models (v015 / D-H1 corpora, N\_REPS = 200)", "FontSize", 11, "FontWeight", "bold");
legend(ax2, hh, {"xu (phase without amplitude structure)", "shaped\_xu (primary)", "bootstrap (robustness check)"}, "Location", "northeast", "FontSize", 8.5);
save1(f2, opts, "fig7b_threeModels_v015.png");

%% 7c / 7d: per-trial contrasts rebuilt from the v015 files
[pP, pQ] = pairs_local(6);
pair = struct();  for i = 1:nDs, pair.(ds7(i)) = rebuild_local(srcDir, ds7(i), pP, pQ, MDC, DEG, sx.N(i)); end
inDom = ds7(1:5);  dsCol = lines(nDs);  A = opts.AxLim;
f3 = figure("Name", "Fig7c", "Position", [50 50 520 520]);  ax3 = axes(f3);  hold(ax3, "on");  nOut = 0;  nAll = 0;
for i = 1:5
    p = pair.(inDom(i));  x = p.dO(:);  y = p.dP(:);  nAll = nAll + numel(x);  nOut = nOut + nnz(abs(x) > A | abs(y) > A);
    scatter(ax3, x, y, 2, dsCol(i, :), "filled", "MarkerFaceAlpha", 0.18, "DisplayName", char(lab(i)));
end
plot(ax3, [-A A], [-A A], "k-", "LineWidth", 1.2);  plot(ax3, [-A A], [0 0], "k:");  plot(ax3, [0 0], [-A A], "k:");
set(ax3, "XLim", [-A A], "YLim", [-A A], "FontSize", 10, "Box", "off", "DataAspectRatio", [1 1 1]);
xlabel(ax3, "Observed \Delta\beta (inter-pipeline)"); ylabel(ax3, "Predicted \Delta\beta (shaped\_xu)");
title(ax3, "Fig 7c: Predicted vs observed contrasts (five in-domain datasets)", "FontSize", 10.5, "FontWeight", "bold");
legend(ax3, "show", "Location", "northwest", "FontSize", 8, "Box", "off");
text(ax3, A, -A + 0.012, sprintf("%d of %d points lie outside the axes", nOut, nAll), "FontSize", 8, "HorizontalAlignment", "right", "Color", [.4 .4 .4]);
save1(f3, opts, "fig7c_predVsObs_v015.png");

ctype = zeros(1, 15);  for k = 1:15, ctype(k) = ctype_local(pP(k), pQ(k)); end   % 1 derivation, 2 regression, 3 mixed
ctCol = [0.18 0.45 0.73; 0.85 0.33 0.10; 0.45 0.67 0.19];
ctNm  = {"Derivation (BWFD vs SG, same regression)", "Regression (same derivation)", "Mixed (both differ)"};
f4 = figure("Name", "Fig7d", "Position", [50 50 1200 780]);  tl = tiledlayout(f4, 2, 3, "TileSpacing", "compact", "Padding", "compact");
title(tl, "Fig 7d: Predicted vs observed \Delta\beta by dataset (shaped\_xu, v015)", "FontSize", 11, "FontWeight", "bold");
slopes = NaN(1, nDs);
for i = 1:nDs
    axT = nexttile(tl);  hold(axT, "on");  p = pair.(ds7(i));  r = row("shaped_xu", ds7(i));
    for ct = 1:3
        m = ctype == ct;  h = scatter(axT, reshape(p.dO(:, m), [], 1), reshape(p.dP(:, m), [], 1), 4, ctCol(ct, :), "filled", "MarkerFaceAlpha", 0.30);
        if i == 1, lh(ct) = h; end %#ok<AGROW>
    end
    q = polyfit(p.dO(:), p.dP(:), 1);  slopes(i) = q(1);
    plot(axT, [-A A], polyval(q, [-A A]), "r--", "LineWidth", 1.4);  text(axT, -A*0.95, A*0.88, sprintf("slope=%.2f", q(1)), "FontSize", 8, "Color", [.75 0 0]);
    plot(axT, [-A A], [-A A], "k-", "LineWidth", 1.2);  plot(axT, [-A A], [0 0], "k:");  plot(axT, [0 0], [-A A], "k:");
    set(axT, "XLim", [-A A], "YLim", [-A A], "FontSize", 9, "Box", "off", "DataAspectRatio", [1 1 1]);
    title(axT, sprintf("%s  (N=%d)\ncccMed=%.3f  r(med)=%.3f  MAE(med)=%.4f", lab(i), r.Nvalid, r.cccMed, r.pearsonMed, r.maeMed), "FontSize", 8, "Color", dsCol(i, :));
    if i == 1, legend(axT, lh, ctNm, "Location", "southeast", "FontSize", 7, "Box", "off"); end
    if mod(i - 1, 3) == 0, ylabel(axT, "Predicted \Delta\beta"); end
    if i > 3, xlabel(axT, "Observed \Delta\beta"); end
end
save1(f4, opts, "fig7d_tiled_v015.png");

out = struct("dataset", ds7, "cccMed", sx.cccMed, "cccMedLo", sx.lo, "cccMedHi", sx.hi, "loopCCC", sx.loop, "Nvalid", sx.N, "slope7d", slopes, ...
    "cccMedXu", cm(1, :), "cccMedBootstrap", cm(3, :), "nOutside7c", nOut, "nPoints7c", nAll, "source", "constellationMetrics_v005.mat; loopClosureResults_<ds>_all_shaped_xu_v015.mat");
save(fullfile(root, "results", "plotConstellationMetrics_v003.mat"), "-struct", "out");
fprintf("plotConstellationMetrics_v003 done. 7c: %d of %d points outside +-%.2f\n", nOut, nAll, A);
end

function p = rebuild_local(srcDir, ds, pP, pQ, MDC, DEG, nExpect)
f = fullfile(srcDir, "loopClosureResults_" + ds + "_all_shaped_xu_v015.mat");  if ~isfile(f), error("pcm3:file", "%s", "FAILED PATH: " + f); end
S = load(f, "results");  r = S.results;  N = numel(r);  dO = NaN(N, numel(pP));  dP = dO;  ok = false(N, 1);
for t = 1:N
    bo = r(t).betaObs(:)';  sl = r(t).betaRecSlice;  g = r(t).betaGenStarMed;
    if numel(bo) ~= 6 || any(~isfinite(bo)) || isempty(sl) || size(sl, 2) ~= 6 || ~isfinite(g), continue, end
    bp = median(sl, 1, "omitnan");  if any(~isfinite(bp)), continue, end
    dO(t, :) = bo(pP) - bo(pQ);  dP(t, :) = bp(pP) - bp(pQ);  ok(t) = true;
end
degen = ok & max(abs(dO), [], 2) < DEG;  scenA = ok & all(abs(dO) < MDC, 2);  vi = find(ok & ~degen & ~scenA);
if numel(vi) ~= nExpect, error("pcm3:count", "%s: rebuilt %d valid trials, constellationMetrics_v005 has %d", ds, numel(vi), nExpect); end
p = struct("dO", dO(vi, :), "dP", dP(vi, :));
end

function ct = ctype_local(p, q)
% pipelines 1..6 = BWFD-OLS, SG-OLS, BWFD-LMLS, SG-LMLS, BWFD-IRLS, SG-IRLS: regression = ceil(i/2), derivation = parity of i
sameReg = ceil(p/2) == ceil(q/2);  sameDer = mod(p, 2) == mod(q, 2);
if sameReg, ct = 1; elseif sameDer, ct = 2; else, ct = 3; end
end

function [pP, pQ] = pairs_local(n)
pP = [];  pQ = [];  for p = 2:n, for q = 1:p-1, pP(end+1) = p;  pQ(end+1) = q; end, end %#ok<AGROW>
end

function save1(fig, opts, name)
if opts.SaveFigs, exportgraphics(fig, fullfile(opts.FigDir, name), "Resolution", 300);  fprintf("  Saved %s\n", name);  close(fig); else, drawnow; end
end
