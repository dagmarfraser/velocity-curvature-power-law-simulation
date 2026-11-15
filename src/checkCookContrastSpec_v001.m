%% checkCookContrastSpec_v001  Deep Review GD17 item 7c: direction and pipeline order of the Cook ASD-CTRL contrasts.
% Reads the saved results behind the two Cook passages and asserts the direction conventions the manuscript now states:
%   tempo paragraph  results/tempoDecomposition_v002.mat (A3: contrast "Cook ASD-CTRL"; pipelines SG-IRLS, BWFD-OLS) and results/checkCookMatchBins_v001.mat
%                    (group = ASD indicator, so the coefficient is ASD minus CTRL; PIPES order SG-IRLS then BWFD-OLS)
%   General Discussion src/compareCookCtrlAsd_v002.mat (diffPoint = pointEstCtrl - pointEstAsd, i.e. CTRL minus ASD; seDiff = SD of the bootstrap difference)
% Run from src/. Writes results/checkCookContrastSpec_v001.mat.   Fraser, D.S. (2026)  v001
ROOT = fileparts(fileparts(mfilename("fullpath")));
fl = @(p) fullfile(ROOT, p);  for p = ["results/tempoDecomposition_v002.mat" "results/checkCookMatchBins_v001.mat" "src/compareCookCtrlAsd_v002.mat"]
    if ~isfile(fl(p)), error("cookSpec:file", "%s", "FAILED PATH: " + fl(p)); end
end
A = load(fl("results/tempoDecomposition_v002.mat"), "A3").A3;  A = A(A.contrast == "Cook ASD-CTRL", :);
M = load(fl("results/checkCookMatchBins_v001.mat"), "M").M;
C = load(fl("src/compareCookCtrlAsd_v002.mat"), "out").out.median;
g = @(p, o) A(A.pipeline == p & A.outcome == o, :);
R = struct();
R.rawObs   = [g("SG-IRLS", "betaObs").b_only   g("BWFD-OLS", "betaObs").b_only];        % +0.039, +0.046
R.tempoObs = [g("SG-IRLS", "betaObs").b_tempo  g("BWFD-OLS", "betaObs").b_tempo];       % +0.004, +0.020
mb = @(bw, p, o) M.b_matched(M.binWidth == bw & M.pipeline == p & M.outcome == o);
R.halfOctObs = [mb(0.5, "SG-IRLS", "betaObs") mb(0.5, "BWFD-OLS", "betaObs")];  R.halfOctGen = [mb(0.5, "SG-IRLS", "betaGen") mb(0.5, "BWFD-OLS", "betaGen")];
R.tenthObs   = [mb(0.1, "SG-IRLS", "betaObs") mb(0.1, "BWFD-OLS", "betaObs")];  R.tenthGen   = [mb(0.1, "SG-IRLS", "betaGen") mb(0.1, "BWFD-OLS", "betaGen")];
assert(R.rawObs(1) > 0 && R.rawObs(2) > 0, "cookSpec:dir", "Raw ASD-CTRL beta_obs contrast is not positive for both pipelines");
assert(abs(C.diffPoint - (C.pointEstCtrl - C.pointEstAsd)) < 1e-12, "cookSpec:gd", "compareCookCtrlAsd diff is not CTRL minus ASD");
R.gdCtrlMinusAsd = C.diffPoint;  R.gdAsdMinusCtrl = -C.diffPoint;  R.gdSeDiff = C.seDiff;  R.gdSeCtrl = C.seCtrl;  R.gdSeAsd = C.seAsd;
R.gdSeHypot = hypot(C.seCtrl, C.seAsd);
disp(R);  save(fl("results/checkCookContrastSpec_v001.mat"), "-struct", "R");
