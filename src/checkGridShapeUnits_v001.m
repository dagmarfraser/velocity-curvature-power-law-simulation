%% checkGridShapeUnits_v001.m
% What are the v058 grid's geometry and units, and where did "235 x 109" come from?
%   1. Extents of the nu = 2 shape the generator uses, in px and in harness mm (4.8 px/mm).
%   2. V0 analytic tempo check (checkGridTempoAtFold_v004 L154-159) redone three ways:
%      true shape (closed integral of kappa^beta ds on the path), ellipse at the extents,
%      ellipse at 235 x 109 (the value carried from PaperDraft v002's General Discussion).
%   3. Two ellipse-equivalents, to see whether 235 is one of them.
%   4. VGF units: kappa is in 1/px, so VGF_mm = VGF_px * pixelScale^(beta - 1).
% Reads results/gridTempoAtFold_v004.mat (local, else the RDS mount).
% Writes results/gridShapeUnits_v001.mat.
% Fraser, D.S. (2026)  v001

%% CONFIG
ROOT   = fileparts(fileparts(mfilename("fullpath")));
V0_CANDS = [fullfile(ROOT, "results", "gridTempoAtFold_v004.mat")
            "/Volumes/rdsprojects/f/fraserds-mpo-evaluation/2026_prereg/velocity-curvature-power-law-simulation-main/velocity-curvature-power-law-simulation-main/results/gridTempoAtFold_v004.mat"];
OUT_MAT = fullfile(ROOT, "results", "gridShapeUnits_v001.mat");
OLD_AB  = [235 109];                     % checkGridTempoAtFold_v004 L36; runTempoFoldSweep_v001 L28;
                                         % runTempoFoldSweepE_v001 L28

%% 1. Geometry
addpath(genpath(fullfile(ROOT, "src", "functions")));
G = gridShapeGeometry_v001(ROOT);
fprintf("Shape: %s (caches: %s)\n", G.source, strjoin(G.cacheCheck, ", "));
fprintf("Semi-extents %.2f x %.2f px, a/b %.4f; at %.2f px/mm: %.2f x %.2f mm; perimeter %.1f px\n", ...
    G.aPx, G.bPx, G.ab, G.pixelScale, G.aMM, G.bMM, G.perimPx);
fprintf("Carried value %g x %g px: a/b %.4f; %.2f x %.2f mm\n", OLD_AB, OLD_AB(1)/OLD_AB(2), OLD_AB / G.pixelScale);

%% 2. V0 analytic check, three geometries
iV = find(isfile(V0_CANDS), 1);
if isempty(iV), error("gridUnits:V0", "%s", "V0 output not found: " + strjoin(V0_CANDS, " | ")); end
V = load(V0_CANDS(iV), "out");  F = V.out.perNode;
F = F(F.generated_beta > 0 & isfinite(F.f0), :);
ds = vecnorm(diff([G.pathPx; G.pathPx(1, :)]), 2, 2);          % segment i -> i+1
dsV = (ds + circshift(ds, 1)) / 2;                               % length attributed to vertex i
shapeTempo   = @(vgf, b) vgf / sum(G.kappaPx .^ b .* dsV);
ellTempo     = @(vgf, b, a, bb) vgf / integral(@(p) (a*bb ./ (a^2*sin(p).^2 + bb^2*cos(p).^2).^1.5).^b ...
                   .* sqrt(a^2*sin(p).^2 + bb^2*cos(p).^2), 0, 2*pi, "RelTol", 1e-10);
geoms = ["shape" "extents" "carried235"];
E = table();
for g = geoms
    switch g
        case "shape",      f = arrayfun(shapeTempo, F.vgf_value, F.generated_beta);
        case "extents",    f = arrayfun(@(v, b) ellTempo(v, b, G.aPx, G.bPx), F.vgf_value, F.generated_beta);
        case "carried235", f = arrayfun(@(v, b) ellTempo(v, b, OLD_AB(1), OLD_AB(2)), F.vgf_value, F.generated_beta);
    end
    r = abs(F.f0 ./ f - 1);
    E = [E; table(g, median(r), max(r), 'VariableNames', ["geometry" "medRelErr" "maxRelErr"])]; %#ok<AGROW>
end
fprintf("\nV0 analytic tempo check over %d coordinates (beta_gen > 0), |realised/analytic - 1|:\n", height(F));
disp(E);

%% 3. Ellipse-equivalents
cfit = fitEllipseAxes_local(G.pathPx);
perimE = @(a, b) integral(@(p) sqrt(a^2*sin(p).^2 + b^2*cos(p).^2), 0, 2*pi);
aPer = fzero(@(a) perimE(a, G.bPx) - G.perimPx, G.aPx);
fprintf("Least-squares ellipse fit: %.2f x %.2f px;  equal-perimeter a at b = %.2f: %.2f px\n", ...
    cfit, G.bPx, aPer);

%% 4. VGF units
fprintf("\nVGF is not mm/s: v = VGF*kappa^-beta with kappa in 1/px, so VGF_mm = VGF_px * %.1f^(beta - 1).\n", G.pixelScale);
for b = [0 1/3 2/3]
    fprintf("  beta %.3f: grid VGF 90.0-330.3 -> %.1f-%.1f (harness-mm units)\n", b, [90.017 330.3] * G.pixelScale^(b - 1));
end

%% Save
if ~isfolder(fileparts(OUT_MAT)), error("gridUnits:outDir", "%s", "Missing folder: " + fileparts(OUT_MAT)); end
out = struct("geometry", rmfield(G, ["pathPx" "kappaPx"]), "v0Check", E, "ellipseFit", cfit, ...
             "equalPerimA", aPer, "carried", OLD_AB, "v0File", V0_CANDS(iV), "runDate", string(datetime("now")));
save(OUT_MAT, "out");
fprintf("\nSaved %s\n", OUT_MAT);

function ax = fitEllipseAxes_local(P)
% Axis-aligned least-squares ellipse about the centroid: (x/a)^2 + (y/b)^2 = 1.
    c = mean(P, 1);  X = P - c;
    w = [X(:, 1).^2, X(:, 2).^2] \ ones(size(X, 1), 1);
    if any(w <= 0), error("gridUnits:fit", "%s", "Ellipse fit not positive definite"); end
    ax = 1 ./ sqrt(w');
end
