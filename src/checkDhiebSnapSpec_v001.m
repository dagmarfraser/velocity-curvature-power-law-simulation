%% checkDhiebSnapSpec_v001  Deep Review follow-up: Dhieb's noise centroid (mean IRASA alpha, mean sigma in mm) against the grid node it snaps to for the SEM check.
% Centroid from results/checkVGFSEMCentroid_v004.mat (table `centroids`); nodes = the alpha and sigma values of the SEM grid (perCoordinateSEM_v2_001.mat), nearest in each, as checkVGFSEMCentroid_v004 does.
% Run from src/. Writes results/checkDhiebSnapSpec_v001.mat.   Fraser, D.S. (2026)  v001
ROOT = fileparts(fileparts(mfilename("fullpath")));
f1 = fullfile(ROOT, "results", "checkVGFSEMCentroid_v004.mat");  f2 = fullfile(ROOT, "src", "perCoordinateSEM_v2_001.mat");
for f = [f1 f2], if ~isfile(f), error("dsn:file", "%s", "FAILED PATH: " + f); end, end
C = load(f1, "centroids").centroids;  c = C(string(C.label) == "Dhieb", :);  if height(c) ~= 1, error("dsn:row", "expected one Dhieb row"); end
G = load(f2, "coordTable").coordTable;  A = sort(unique(G.alpha));  S = sort(unique(G.sigma));
[~, ai] = min(abs(A - c.alpha));  [~, si] = min(abs(S - c.sigmaMM));
out = struct("centroidAlpha", c.alpha, "centroidSigmaMM", c.sigmaMM, "nodeAlpha", A(ai), "nodeSigmaMM", S(si), "alphaStep", mean(diff(A)), "sigmaStep", mean(diff(S)));
disp(out);  save(fullfile(ROOT, "results", "checkDhiebSnapSpec_v001.mat"), "-struct", "out");
