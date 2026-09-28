%% analyseTempoFoldSweepE_v002.m
% Block E analysis (runTempoFoldSweepE_v001 output). Floor statistics exactly as
% attractorLocation_v001 L295-302 (location = mean beta_rec over the 22 beta_gen values,
% spread = max - min, slope = polyfit slope), for three sources:
%   GRID    perCoordinateSEM_v2_001.mat, fs 120, VGF node 8 (the registered §3a floors)
%   REPLAY  v058 generator through the engine (gate: must reproduce GRID)
%   FIXED   exact ellipse at fixed f0 (does the floor survive without the tempo path?)
% Also: realised grid-path tempo at VGF node 8 (V0 v004), sigma = 0 companions, and the
% share of failed fits per cell. regressDataEBR is called with limitBreak = 0, so any fitnlm
% warning (incl. iteration limit) returns NaN: no unconverged estimate is used, but means over
% survivors are selection-conditioned where failures are common (FAIL_MAX).
% Reads: results/tempoFoldSweepE_v001.mat (RDS first via mount if not local),
%        src/perCoordinateSEM_v2_001.mat, results/gridTempoAtFold_v004.mat
% Writes: results/tempoFoldAnalysisE_v002.mat, figures/floorsByTempo_v002.png
% v002: replay-gate expectation from the t distribution with Welch-Satterthwaite df per node
%       (grid SE rests on 5 replicates, so |z| > 3 is expected ~4%, not the normal 0.3% v001
%       printed); floor values and tables otherwise identical to v001.

%% CONFIG
ROOT     = fileparts(fileparts(mfilename("fullpath")));
MOUNT    = "/Volumes/rdsprojects/f/fraserds-mpo-evaluation/2026_prereg/velocity-curvature-power-law-simulation-main/velocity-curvature-power-law-simulation-main";
SEM_MAT  = fullfile(ROOT, "src", "perCoordinateSEM_v2_001.mat");
GRID_MAT = fullfile(ROOT, "results", "gridTempoAtFold_v004.mat");
OUT_MAT  = fullfile(ROOT, "results", "tempoFoldAnalysisE_v002.mat");
OUT_PNG  = fullfile(ROOT, "figures", "floorsByTempo_v002.png");
LABELS   = ["BWFD-OLS" "BWFD-LMLS" "BWFD-IRLS" "SG-OLS" "SG-LMLS" "SG-IRLS"];   % engine order
MS3A     = [0.277 0.188 0.225 0.266 0.178 0.210];   % §3a registered floors, alpha 2, LABELS order
Z_GATE   = 3;                                       % replay gate: |z| per beta_gen node
FAIL_MAX = 10;                                      % % failed fits above which a floor is survivor-conditioned

%% Load
cands = [fullfile(ROOT, "results", "tempoFoldSweepE_v001.mat"); fullfile(MOUNT, "results", "tempoFoldSweepE_v001.mat")];
iE = find(isfile(cands), 1);  if isempty(iE), error("floorsE:input", "%s", "Missing: " + strjoin(cands, " | ")); end
fprintf("Input: %s\n", cands(iE));
E = load(cands(iE));  R = E.R;
for f = [SEM_MAT GRID_MAT], if ~isfile(f), error("floorsE:input", "%s", "Missing: " + f); end, end
S = load(SEM_MAT, "coordTable");  T = S.coordTable;
for v = ["betaGen" "VGF" "fs" "alpha" "sigma" "meanBetaRec" "sem" "nReps"], T.(v) = double(T.(v)); end
T = T(T.fs == E.FS & T.VGF == E.VGF0, :);
G = load(GRID_MAT);  N = G.out.perNode;  N = N(N.sampling_rate == E.FS & N.vgf_value == E.VGF0, :);
fprintf("Grid path at VGF %.1f, fs %d: realised tempo %.2f Hz (beta_gen %.3f) to %.2f Hz (beta_gen %.3f)\n", ...
    E.VGF0, E.FS, min(N.f0), N.generated_beta(N.f0 == min(N.f0)), max(N.f0), N.generated_beta(N.f0 == max(N.f0)));

%% Floor statistics per source x condition x pipeline
F = table();
for a = E.ALPHAS
    for s = [0 E.SIGMAS_MM]
        g = T(T.alpha == a & T.sigma == s, :);
        for p = 1:numel(LABELS)
            gp = sortrows(g(g.pipeline == LABELS(p), :), "betaGen");
            fp = 100 * (1 - sum(gp.nReps) / (5 * height(gp)));   % v058: 5 reps per coordinate
            if height(gp) >= 3, F = [F; row_local("GRID", NaN, a, s, LABELS(p), gp.betaGen, gp.meanBetaRec, fp)]; end %#ok<AGROW>
        end
        for src = ["REPLAY" "FIXED"]
            f0s = unique(R.f0(R.block == src));  if src == "REPLAY", f0s = NaN; end
            for f0 = f0s(:)'
                m = R.block == src & R.sigmaMM == s & ~R.lengthExcluded & (isnan(f0) | R.f0 == f0);
                if s > 0, m = m & R.alpha == a; elseif a ~= E.ALPHAS(1), continue; end   % sigma 0 once
                r = sortrows(R(m, :), "betaGen");
                for p = 1:numel(LABELS)
                    ok = isfinite(r.betaMean(:, p));
                    fp = 100 * sum(r.nFail(:, p)) / sum(r.nReps);
                    if sum(ok) >= 3, F = [F; row_local(src, f0, a, s, LABELS(p), r.betaGen(ok), r.betaMean(ok, p), fp)]; end %#ok<AGROW>
                end
            end
        end
    end
end

%% 1. Which sigma reproduces §3a? (GRID, alpha 2)
fprintf("\n1. Registered floors (GRID, alpha 2, fs %d, VGF %.1f) against §3a:\n", E.FS, E.VGF0);
g2 = F(F.source == "GRID" & F.alpha == 2 & F.sigma > 0, ["sigma" "pipeline" "loc" "spread" "slope"]);
g2.ms3a = MS3A(arrayfun(@(x) find(LABELS == x), g2.pipeline))';
g2.diff = g2.loc - g2.ms3a;  disp(g2)

%% 2. Replay gate: REPLAY vs GRID per beta_gen node, within replicate SE
Gt = table();
for a = E.ALPHAS
    for s = E.SIGMAS_MM
        r = R(R.block == "REPLAY" & R.alpha == a & R.sigmaMM == s, :);
        for p = 1:numel(LABELS)
            gp = T(T.alpha == a & T.sigma == s & T.pipeline == LABELS(p), ["betaGen" "meanBetaRec" "sem" "nReps"]);
            rp = table(r.betaGen, r.betaMean(:, p), r.betaSD(:, p), r.nReps, 'VariableNames', ["betaGen" "eng" "engSD" "engN"]);
            m = innerjoin(rp, gp, "Keys", "betaGen");  m = m(isfinite(m.eng) & isfinite(m.meanBetaRec), :);
            v1 = m.engSD.^2 ./ m.engN;  v2 = m.sem.^2 ./ m.nReps;  se = sqrt(v1 + v2);
            z  = (m.eng - m.meanBetaRec) ./ se;
            df = (v1 + v2).^2 ./ (v1.^2 ./ (m.engN - 1) + v2.^2 ./ (m.nReps - 1));   % Welch-Satterthwaite
            pExp = 2 * tcdf(-Z_GATE, df);
            Gt = [Gt; table(a, s, LABELS(p), height(m), max(abs(m.eng - m.meanBetaRec)), max(abs(z)), ...
                mean(abs(z) > Z_GATE), mean(pExp, "omitnan"), median(df, "omitnan"), 'VariableNames', ...
                ["alpha" "sigma" "pipeline" "n" "maxAbsDiff" "maxAbsZ" "fracZover" "fracZexpected" "medianDf"])]; %#ok<AGROW>
        end
    end
end
fprintf("\n2. Replay gate (engine + v058 generator + sigma x %.1f px vs perCoordinateSEM):\n", E.PIXEL_SCALE);
obs = sum(Gt.fracZover .* Gt.n) / sum(Gt.n);  expd = sum(Gt.fracZexpected .* Gt.n) / sum(Gt.n);
fprintf("   nodes compared %d; max |diff| %.4f; |z| > %d observed %.4f vs expected %.4f under t (median Welch df %.1f)\n", ...
    sum(Gt.n), max(Gt.maxAbsDiff), Z_GATE, obs, expd, median(Gt.medianDf));
k = round(sum(Gt.fracZover .* Gt.n));  pB = 1 - binocdf(k - 1, sum(Gt.n), expd);
fprintf("   one-sided binomial test of excess exceedances: %d of %d nodes, p = %.3f (gate fails if p < .01)\n", k, sum(Gt.n), pB);
if pB < 0.01
    warning("floorsE:replayGate", "%s", "Replay exceeds the grid beyond replicate error: engine and grid disagree.");
end
disp(Gt(Gt.maxAbsZ > Z_GATE, :))

%% 3. Floors at fixed tempo vs the grid path
fprintf("\n3. Floor location (mean beta_rec over beta_gen), by source and tempo:\n");
L = F(F.sigma > 0, :);  L.src = L.source;  k = ~isnan(L.f0);  L.src(k) = L.source(k) + "_" + string(L.f0(k)) + "Hz";
Wl = unstack(L(:, ["alpha" "sigma" "pipeline" "src" "loc"]), "loc", "src", "VariableNamingRule", "preserve");
disp(Wl)
fprintf("Failed fits (%%) by source and tempo:\n");
Wf = unstack(L(:, ["alpha" "sigma" "pipeline" "src" "failPct"]), "failPct", "src", "VariableNamingRule", "preserve");
disp(Wf)
Lm = L;  Lm.loc(Lm.failPct > FAIL_MAX) = NaN;
fprintf("Floor location with cells over %d%% failed fits masked (survivor-conditioned):\n", FAIL_MAX);
Wm = unstack(Lm(:, ["alpha" "sigma" "pipeline" "src" "loc"]), "loc", "src", "VariableNamingRule", "preserve");
disp(Wm)
fprintf("Slope (polyfit) by source and tempo (0 = flat attractor, 1 = identity):\n");
Ws = unstack(L(:, ["alpha" "sigma" "pipeline" "src" "slope"]), "slope", "src", "VariableNamingRule", "preserve");
disp(Ws)
Z0 = F(F.sigma == 0 & F.source ~= "GRID", ["source" "f0" "pipeline" "loc" "spread" "slope"]);
fprintf("sigma = 0 companions:\n");  disp(Z0)

%% Figure: alpha 2, both sigma, location and slope vs f0 with the grid path marked
if ~isfolder(fileparts(OUT_PNG)), error("floorsE:outDir", "%s", "Missing folder: " + fileparts(OUT_PNG)); end
fg = figure("Color", "w", "Position", [80 80 1100 700]);  tl = tiledlayout(2, numel(E.SIGMAS_MM), "TileSpacing", "compact");
cols = lines(numel(LABELS));
for stat = ["loc" "slope"]
    for s = E.SIGMAS_MM
        nexttile; hold on;
        for p = 1:numel(LABELS)
            x = F(F.source == "FIXED" & F.alpha == 2 & F.sigma == s & F.pipeline == LABELS(p), :);
            plot(x.f0, x.(stat), "-o", "Color", cols(p, :), "LineWidth", 1.3, "DisplayName", LABELS(p));
            g = F(F.source == "GRID" & F.alpha == 2 & F.sigma == s & F.pipeline == LABELS(p), :);
            if ~isempty(g), yline(g.(stat), ":", "Color", cols(p, :), "HandleVisibility", "off"); end
        end
        set(gca, "XScale", "log");  box on;  xlabel("f_0 (Hz), fixed tempo");
        if stat == "loc", ylabel("floor location (mean \beta_{rec})"); yline(1/3, "k--", "HandleVisibility", "off");
        else, ylabel("slope d\beta_{rec}/d\beta_{gen}"); end
        title(sprintf("\\alpha = 2, \\sigma = %g mm (dotted: grid path)", s));
        if stat == "loc" && s == E.SIGMAS_MM(1), legend("Location", "best"); end
    end
end
title(tl, "Compression floors: fixed tempo (lines) vs the grid's fixed-VGF path (dotted)");
set(findall(fg, "Type", "axes"), "Toolbar", []);  exportgraphics(fg, OUT_PNG, "Resolution", 200);

if ~isfolder(fileparts(OUT_MAT)), error("floorsE:outDir", "%s", "Missing folder: " + fileparts(OUT_MAT)); end
save(OUT_MAT, "F", "g2", "Gt", "Wl", "Wf", "Wm", "Ws", "Z0", "FAIL_MAX", "-v7.3");
fprintf("\nSaved: %s\nFigure: %s\n", OUT_MAT, OUT_PNG);

%% =========================================================================
function r = row_local(src, f0, a, s, pipe, b, y, failPct)
    p = polyfit(b, y, 1);
    r = table(string(src), f0, a, s, string(pipe), numel(b), mean(y), max(y) - min(y), p(1), failPct, ...
        'VariableNames', ["source" "f0" "alpha" "sigma" "pipeline" "nBeta" "loc" "spread" "slope" "failPct"]);
end

