%% checkGridTempoRange_v001.m
% Which tempo range does the v058 grid span at beta_gen = 1/3, and why do two
% numbers appear in the manuscript?  Resolves coherence-pass item A4.4.
%   REALISED  0.49-1.79 Hz (Part 1B; D21): f0 = ORBITS/duration read from the DB
%             (results/gridTempoAtFold_v004.mat, perNode).
%   K_CONV    0.51-1.89 Hz (Part 1 Results (c)): VGF / K_conv, K_conv = 175.26 from
%             the cached baseline ellipse.
% This script (1) reproduces both, per sampling rate, (2) computes the analytic
% path-integral tempo, f0 = VGF / sum(kappa^beta ds), on the generator's own shape,
% and (3) reports which of the two routes the analytic value supports and by how much
% they differ.  It does not decide the wording; it supplies the numbers.
% Self-checks (Fail Loud): the manuscript's 0.49-1.79 (1/3) and 2.76-10.71 (2/3)
% must reproduce at two decimals; K_conv must reproduce 175.2636.
% Reads:  results/gridTempoAtFold_v004.mat (local, else the RDS mount),
%         src/functions/baselineShp6_120Hz.mat, gridShapeGeometry_v002.
% Writes: results/checkGridTempoRange_v001.mat
% USAGE:  from the project root: checkGridTempoRange_v001
% Fraser, D.S. (2026)  v001

%% CONFIG
ROOT   = fileparts(fileparts(mfilename("fullpath")));
V0_CANDS = [fullfile(ROOT, "results", "gridTempoAtFold_v004.mat")
            "/Volumes/rdsprojects/f/fraserds-mpo-evaluation/2026_prereg/velocity-curvature-power-law-simulation-main/velocity-curvature-power-law-simulation-main/results/gridTempoAtFold_v004.mat"];
OUT_MAT   = fullfile(ROOT, "results", "checkGridTempoRange_v001.mat");
BETAS     = [1/3, 2/3];
BETA_TOL  = 1e-3;                 % DB stores beta_gen to 5 s.f. (checkGridTempoAtFold_v004)
PUB       = struct("third", [0.49 1.79], "twoThirds", [2.76 10.71]);   % v004 L1258-1259
KCONV_PUB = 175.2636;
KCONV_TOL = 1e-3;
VGF_LO = exp(4.5);  VGF_HI = exp(5.8);                                % registered grid range

%% 1. Realised tempo from the DB (perNode)
iV = find(isfile(V0_CANDS), 1);
if isempty(iV), error("gridTempoRange:V0", "%s", "V0 output not found: " + strjoin(V0_CANDS, " | ")); end
V = load(V0_CANDS(iV), "out");
F = V.out.perNode;
need = ["generated_beta" "vgf_value" "sampling_rate" "f0"];
if ~all(ismember(need, string(F.Properties.VariableNames)))
    error("gridTempoRange:cols", "%s", "perNode lacks " + strjoin(setdiff(need, string(F.Properties.VariableNames)), ", "));
end
if any(~isfinite(F.f0)), error("gridTempoRange:nan", "%s", sprintf("%d non-finite realised f0 in perNode", nnz(~isfinite(F.f0)))); end

%% 2. K_conv and the analytic path integral
addpath(ROOT); addpath(fullfile(ROOT, "src")); addpath(genpath(fullfile(ROOT, "src", "functions")));
G = gridShapeGeometry_v002(ROOT);
shp = load(fullfile(ROOT, "src", "functions", "baselineShp6_120Hz.mat"), "pathXYresample", "K");
perimPx = sum(sqrt(sum(diff(shp.pathXYresample, 1, 1).^2, 2)));
kap = shp.K(:); kap = kap(kap > 0 & isfinite(kap));
Kconv = perimPx / mean(kap .^ (-1/3));
if abs(Kconv - KCONV_PUB) > KCONV_TOL
    error("gridTempoRange:Kconv", "%s", sprintf("K_conv %.6f does not reproduce %.4f", Kconv, KCONV_PUB));
end
ds  = vecnorm(diff([G.pathPx; G.pathPx(1, :)]), 2, 2);
dsV = (ds + circshift(ds, 1)) / 2;
analytic = @(vgf, b) vgf ./ sum(G.kappaPx .^ b .* dsV);

%% 3. Range at each beta_gen, three routes
rows = table();
for b = BETAS
    sel = abs(F.generated_beta - b) < BETA_TOL;
    if ~any(sel), error("gridTempoRange:node", "%s", sprintf("no perNode rows at beta_gen = %.4f", b)); end
    R = F(sel, :);
    fsAll = unique(R.sampling_rate)';
    for fs = [fsAll, NaN]
        if isnan(fs), r = R; lab = "all fs"; else, r = R(R.sampling_rate == fs, :); lab = string(fs) + " Hz"; end
        lo = min(r.f0); hi = max(r.f0);
        rows = [rows; table(b, lab, height(r), lo, hi, 'VariableNames', ...
            ["beta_gen" "fs" "nNodes" "f0Min_realised" "f0Max_realised"])]; %#ok<AGROW>
    end
end
rows.f0Min_analytic = arrayfun(@(b) analytic(VGF_LO, b), rows.beta_gen);
rows.f0Max_analytic = arrayfun(@(b) analytic(VGF_HI, b), rows.beta_gen);
rows.f0Min_Kconv    = repmat(VGF_LO / Kconv, height(rows), 1);      % K_conv is defined at beta = 1/3 only
rows.f0Max_Kconv    = repmat(VGF_HI / Kconv, height(rows), 1);
rows.f0Min_Kconv(abs(rows.beta_gen - 1/3) > BETA_TOL) = NaN;
rows.f0Max_Kconv(abs(rows.beta_gen - 1/3) > BETA_TOL) = NaN;
disp(rows);

%% 4. Self-checks against the manuscript's published ranges
chk = @(b, pub) rows(rows.fs == "all fs" & abs(rows.beta_gen - b) < BETA_TOL, ["f0Min_realised" "f0Max_realised"]);
c3 = chk(1/3, PUB.third);   c23 = chk(2/3, PUB.twoThirds);
if any(abs(round([c3.f0Min_realised c3.f0Max_realised], 2) - PUB.third) > 0)
    error("gridTempoRange:pub3", "%s", sprintf("realised 1/3 range %.4f-%.4f does not reproduce %.2f-%.2f", ...
        c3.f0Min_realised, c3.f0Max_realised, PUB.third));
end
if any(abs(round([c23.f0Min_realised c23.f0Max_realised], 2) - PUB.twoThirds) > 0)
    error("gridTempoRange:pub23", "%s", sprintf("realised 2/3 range %.4f-%.4f does not reproduce %.2f-%.2f", ...
        c23.f0Min_realised, c23.f0Max_realised, PUB.twoThirds));
end

%% 5. What differs, and by how much (beta_gen = 1/3, all fs)
sumry = struct();
sumry.realised = [c3.f0Min_realised c3.f0Max_realised];
sumry.analytic = [analytic(VGF_LO, 1/3) analytic(VGF_HI, 1/3)];
sumry.Kconv    = [VGF_LO VGF_HI] / Kconv;
sumry.relDiff_Kconv_vs_realised = sumry.Kconv ./ sumry.realised - 1;
sumry.relDiff_analytic_vs_realised = sumry.analytic ./ sumry.realised - 1;
sumry.relDiff_analytic_vs_Kconv = sumry.analytic ./ sumry.Kconv - 1;
fprintf("\nbeta_gen = 1/3, VGF %.3f-%.3f (K_conv %.4f):\n", VGF_LO, VGF_HI, Kconv);
fprintf("  realised (DB)          %.4f - %.4f Hz\n", sumry.realised);
fprintf("  analytic path integral %.4f - %.4f Hz  (%+.2f%%, %+.2f%% vs realised)\n", sumry.analytic, 100 * sumry.relDiff_analytic_vs_realised);
fprintf("  VGF / K_conv           %.4f - %.4f Hz  (%+.2f%%, %+.2f%% vs realised)\n", sumry.Kconv, 100 * sumry.relDiff_Kconv_vs_realised);
fprintf("  K_conv vs analytic     %+.2f%%, %+.2f%%\n", 100 * (sumry.Kconv ./ sumry.analytic - 1));

%% 6. Save
runDate = string(datetime("now", "Format", "yyyy-MM-dd"));
save(OUT_MAT, "rows", "sumry", "Kconv", "runDate");
fprintf("Saved %s\n", OUT_MAT);
