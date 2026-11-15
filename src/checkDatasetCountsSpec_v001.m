%% checkDatasetCountsSpec_v001  Deep Review GD18: dataset counts used in section 3b, from results/checkVGFSEMCentroid_v004.mat (table `centroids`).
% nDatasets (all with a centroid), nInDomain (inDom true), nAlphaGt5 (centroid alpha > 5), and the names of each. Fails loud on an unexpected table shape.
% Run from src/. Writes results/checkDatasetCountsSpec_v001.mat.   Fraser, D.S. (2026)  v001
ROOT = fileparts(fileparts(mfilename("fullpath")));
f = fullfile(ROOT, "results", "checkVGFSEMCentroid_v004.mat");  if ~isfile(f), error("dc:file", "%s", "FAILED PATH: " + f); end
C = load(f, "centroids").centroids;
if ~all(ismember(["label" "inDom" "alpha"], string(C.Properties.VariableNames))), error("dc:cols", "centroids lacks label/inDom/alpha"); end
out = struct("nDatasets", height(C), "nInDomain", sum(C.inDom), "nAlphaGt5", sum(C.alpha > 5), ...
    "inDomainNames", string(C.label(C.inDom))', "alphaGt5Names", string(C.label(C.alpha > 5))');
disp(out);
save(fullfile(ROOT, "results", "checkDatasetCountsSpec_v001.mat"), "-struct", "out");
