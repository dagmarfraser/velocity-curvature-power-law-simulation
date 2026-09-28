%% analyseTempoFoldSweep_v001.m
% SPEC_TempoFoldSweep_v001 analysis. Reproduces every table and figure reported from
% V0 and Blocks A-D (2026-09-26) from the saved results; no simulation, no DB access.
% Inputs (results/ on RDS first, then the iMac mount; paths printed):
%   gridTempoAtFold_v004.mat (or v003)  checkGridTempoAtFold_v004  (BlueBEAR)
%   tempoFoldSweep_v001.mat             runTempoFoldSweep_v001     (BlueBEAR)
%   tempoFoldSweepC_v001.mat            runTempoFoldSweepC_v001    (BlueBEAR)
% Outputs: results/tempoFoldAnalysis_v001.mat;
%          figures/zeroNoiseMaps_v001.png, figures/mapTopByTempo_v001.png

%% CONFIG
ROOT   = fileparts(fileparts(mfilename("fullpath")));
MOUNT  = "/Volumes/rdsprojects/f/fraserds-mpo-evaluation/2026_prereg/velocity-curvature-power-law-simulation-main/velocity-curvature-power-law-simulation-main";
OUT_MAT = fullfile(ROOT, "results", "tempoFoldAnalysis_v001.mat");
FIG_Z   = fullfile(ROOT, "figures", "zeroNoiseMaps_v001.png");
FIG_C   = fullfile(ROOT, "figures", "mapTopByTempo_v001.png");
LABELS  = ["BWFD-OLS" "BWFD-LMLS" "BWFD-IRLS" "SG-OLS" "SG-LMLS" "SG-IRLS"];   % engine order
FS_REF  = 120;  VGF_IDX = [1 4 7 10 14];  F0_PLOT = [0.5 1 2 2.5 4];  AB_REF = 2.16;
BETA_TAB = [0 0.1 0.2 1/3 0.5 2/3];  F0_ECC = 2.5;

%% Load
g  = loadFirst_local(["gridTempoAtFold_v004.mat" "gridTempoAtFold_v003.mat"], ROOT, MOUNT);
ab = loadFirst_local("tempoFoldSweep_v001.mat", ROOT, MOUNT);
cd_ = loadFirst_local("tempoFoldSweepC_v001.mat", ROOT, MOUNT);
T = g.out.perCoordinate;  P = g.out.perCurve;  R = ab.R;  S = cd_.S;  P3 = cd_.P3;
F = groupsummary(T, ["generated_beta" "vgf_value" "sampling_rate"], "mean", "f0");
F = renamevars(F, "mean_f0", "f0");
node = @(v, b0) v(abs(v - b0) == min(abs(v - b0)));   % nearest stored node (5 s.f. in DB)

%% 1. V0: realised tempo on the grid
bg = unique(F.generated_beta);  b13 = node(bg, 1/3);  b23 = node(bg, 2/3);
fprintf("\n1. V0 realised tempo (f0 = 10 orbits / duration)\n");
V0tab = table();
for fs = unique(F.sampling_rate)'
    a = F(F.sampling_rate == fs & F.generated_beta == b13, :);  b = F(F.sampling_rate == fs & F.generated_beta == b23, :);
    V0tab = [V0tab; table(fs, min(a.f0), max(a.f0), min(b.f0), max(b.f0), height(b), ...
        'VariableNames', ["fs" "f0_13_min" "f0_13_max" "f0_23_min" "f0_23_max" "nVGF_23"])]; %#ok<AGROW>
end
disp(V0tab)
bp = node(bg, 0.5333);  x = F(F.sampling_rate == FS_REF & F.generated_beta == bp, :);
fprintf("Counterfactual: at fixed beta_gen %.4f (fs %d) tempo spans %.2f-%.2f Hz across VGF (ratio %.2f)\n", ...
    bp, FS_REF, min(x.f0), max(x.f0), max(x.f0)/min(x.f0));
y = sortrows(F(F.sampling_rate == FS_REF & F.vgf_value == min(F.vgf_value), :), "generated_beta");
r = y.f0(2:end) ./ y.f0(1:end-1);  k = y.generated_beta(2:end) > 0.39 & y.generated_beta(2:end) < 0.68;
fprintf("Tempo ratio between adjacent beta_gen nodes (0.40-0.67, slowest VGF): %s\n", strjoin(compose("%.2f", r(k)'), " "));
peakTab = sortrows(P(P.pipeline == "SG-IRLS" & P.sampling_rate == FS_REF, ["vgf_value" "pkBeta" "f0Pk" "drop"]), "vgf_value");
fprintf("SG-IRLS fs %d, fold peak by VGF node:\n", FS_REF);  disp(peakTab)

%% 2. Block B2 vs grid: where the differences exceed 0.01
B = R(R.block == "B2", :);  B2diff = table();
for p = 1:numel(LABELS)
    gp = T(T.sampling_rate == FS_REF & T.pipeline == LABELS(p), ["generated_beta" "vgf_value" "betaRec"]);
    rp = table(B.betaGen, B.vgf, B.f0, B.betaMean(:, p), 'VariableNames', ["generated_beta" "vgf_value" "f0" "eng"]);
    m = innerjoin(rp, gp, "Keys", ["generated_beta" "vgf_value"]);  m.diff = m.eng - m.betaRec;
    m.pipeline = repmat(LABELS(p), height(m), 1);
    B2diff = [B2diff; m(abs(m.diff) > 0.01, :)]; %#ok<AGROW>
end
fprintf("\n2. B2 (analytic ellipse at realised grid tempo) vs grid, |diff| > 0.01: %d cells\n", height(B2diff));
disp(sortrows(B2diff, ["pipeline" "vgf_value" "generated_beta"]))

%% 3. Block A: fixed beta_gen, rising tempo (a/b 2.16, fs 120); eccentricity at 2/3
A = R(R.block == "A" & abs(R.ab - AB_REF) < 1e-9 & ~R.lengthExcluded, :);
f0s = unique(A.f0);  At = array2table(f0s, 'VariableNames', "f0");
for b0 = BETA_TAB
    p = find(LABELS == ifelse_local(b0 == 0, "SG-OLS", "SG-IRLS"));    % IRLS does not converge at beta 0
    v = NaN(numel(f0s), 1);
    for i = 1:numel(f0s)
        w = A.betaMean(abs(A.betaGen - b0) < 1e-9 & A.f0 == f0s(i), p);
        if isscalar(w), v(i) = w; elseif ~isempty(w), error("tempoAnalysis:dup", "duplicate Block A node"); end
    end
    At.(sprintf("b%.3f_%s", b0, erase(LABELS(p), "-"))) = v;
end
fprintf("\n3. Block A, sigma = 0, a/b %.2f, fs %d: beta_rec by tempo at fixed beta_gen\n", AB_REF, FS_REF);  disp(At)
E = R(R.block == "A" & abs(R.betaGen - 2/3) < 1e-9 & R.f0 == F0_ECC & ~R.lengthExcluded, :);
Ecc = sortrows(table(round(E.ab, 2), E.betaMean(:, LABELS == "SG-IRLS"), 'VariableNames', ["ab" "SG_IRLS"]), "ab");
fprintf("beta_gen 2/3 at %.1f Hz by eccentricity (SG-IRLS):\n", F0_ECC);  disp(Ecc)

%% 4. Blocks C/D: noise-created folds and the tempo window
NF = P3(P3.hasDesc & ~P3.hasDesc0, ["dataset" "f0" "pipeline" "floor" "top" "nDesc" "riseCov"]);
fprintf("\n4. Block C noise-created folds (%d of %d curves):\n", height(NF), height(P3));  disp(NF)
C = S(S.block == "C" & S.sigma > 0, :);
Top = groupsummary(C, ["dataset" "f0"], ["min" "max"], "top");
fprintf("Map top beta_rec(0.75), range across six pipelines:\n");  disp(Top(:, ["dataset" "f0" "min_top" "max_top"]))
Dd = S(S.block == "D", ["alpha" "sigma" "pipeline" "top" "nDesc" "lowDesc"]);
fprintf("Block D (1 Hz): curves with a descending run: %d of %d (low-end: %d)\n", sum(Dd.nDesc > 0), height(Dd), sum(Dd.lowDesc));
disp(Dd(Dd.nDesc > 0, :))

%% 5. Figures
checkDir_local(FIG_Z);  checkDir_local(OUT_MAT);
blues = [0.52 0.72 0.92; 0.22 0.54 0.87; 0.09 0.37 0.65; 0.05 0.27 0.49; 0.02 0.17 0.33];
fz = figure("Color", "w", "Position", [80 80 1000 430]);  tl = tiledlayout(1, 2, "TileSpacing", "compact");
nexttile; hold on; vg = unique(T.vgf_value);
for i = 1:numel(VGF_IDX)
    c = sortrows(T(T.pipeline == "SG-IRLS" & T.sampling_rate == FS_REF & T.vgf_value == vg(VGF_IDX(i)), :), "generated_beta");
    plot(c.generated_beta, c.betaRec, "-", "Color", blues(i, :), "LineWidth", 1.5, "DisplayName", sprintf("VGF %.0f", vg(VGF_IDX(i))));
end
plot([0 0.75], [0 0.75], "k:", "DisplayName", "identity");  axis([0 0.75 0 0.75]); axis square; box on;
title("A. v058 grid, \sigma = 0: VGF fixed"); xlabel("\beta_{gen}"); ylabel("\beta_{rec}"); legend("Location", "northwest");
nexttile; hold on;
for i = 1:numel(F0_PLOT)
    c = sortrows(A(A.f0 == F0_PLOT(i), :), "betaGen");
    plot(c.betaGen, c.betaMean(:, LABELS == "SG-IRLS"), "-", "Color", blues(i, :), "LineWidth", 1.5, "DisplayName", sprintf("%.1f Hz", F0_PLOT(i)));
end
plot([0 0.75], [0 0.75], "k:", "DisplayName", "identity");  axis([0 0.75 0 0.75]); axis square; box on;
title(sprintf("B. Ellipse a/b %.2f, \\sigma = 0: tempo fixed", AB_REF)); xlabel("\beta_{gen}"); legend("Location", "northwest");
title(tl, sprintf("SG-IRLS forward maps at zero noise, fs %d Hz", FS_REF));
set(findall(fz, "Type", "axes"), "Toolbar", []);  exportgraphics(fz, FIG_Z, "Resolution", 200);

fc = figure("Color", "w", "Position", [80 80 1000 380]);  dsn = unique(C.dataset);
tl = tiledlayout(1, numel(dsn), "TileSpacing", "compact");
for d = 1:numel(dsn)
    nexttile; hold on;
    t = Top(Top.dataset == dsn(d), :);  z0 = S(S.block == "C" & S.sigma == 0 & S.dataset == dsn(d) & S.pipeline == "SG-IRLS", :);
    fill([t.f0; flipud(t.f0)], [t.min_top; flipud(t.max_top)], [0.22 0.54 0.87], "FaceAlpha", 0.25, "EdgeColor", "none", "DisplayName", "noise, six pipelines");
    plot(z0.f0, z0.top, "k-o", "MarkerSize", 4, "DisplayName", "\sigma = 0, SG-IRLS");
    yline(0.75, ":", "HandleVisibility", "off");  set(gca, "XScale", "log");  ylim([0 0.8]);  box on;
    title(strrep(dsn(d), "_", " ")); xlabel("f_0 (Hz)"); if d == 1, ylabel("\beta_{rec} at \beta_{gen} = 0.75"); legend("Location", "south"); end
end
title(tl, "Map top by tempo: noise compresses slow movements, bandwidth compresses fast ones");
set(findall(fc, "Type", "axes"), "Toolbar", []);  exportgraphics(fc, FIG_C, "Resolution", 200);

save(OUT_MAT, "V0tab", "peakTab", "B2diff", "At", "Ecc", "NF", "Top", "Dd", "-v7.3");
fprintf("\nSaved: %s\nFigures: %s, %s\n", OUT_MAT, FIG_Z, FIG_C);

%% =========================================================================
function s = loadFirst_local(names, root, mount)
    cands = [fullfile(root, "results", names(:)); fullfile(mount, "results", names(:))];
    i = find(isfile(cands), 1);
    if isempty(i), error("tempoAnalysis:input", "%s", "None found: " + strjoin(cands, " | ")); end
    fprintf("Input: %s\n", cands(i));  s = load(cands(i));
end

function v = ifelse_local(c, a, b)
    if c, v = a; else, v = b; end
end

function checkDir_local(f)
    if ~isfolder(fileparts(f)), error("tempoAnalysis:outDir", "%s", "Missing folder: " + fileparts(f)); end
end
