%% checkLmmSpec_v001  Deep Review GD17 item 6: coding of the Stage 1 LMM (z-scoring, dummy coding, the 192 count, significance counts).
% Parses ModelAdequacy_Stage1_KitchenSink_v2_001.m (a missing pattern is an error) and reads the fitted model's slim extract (results/L9ModelSlim_v001.mat:
% coefficient names, estimates, SEs, VariableInfo). Records: which predictors were z-scored (samplingRate was NOT), the nominal reference levels, the
% coefficient count by order (1 + 8 + 183 = 192), the number of random-effect groups, and Wald-z significance counts (the file holds no covariance, so
% p is the Wald normal approximation, which is what a t with ~17 million df is).
% Run from src/. Writes results/checkLmmSpec_v001.mat.   Fraser, D.S. (2026)  v001
ROOT = fileparts(fileparts(mfilename("fullpath")));
ks = fileread(fullfile(ROOT, "src", "ModelAdequacy_Stage1_KitchenSink_v2_001.m"));
need = @(txt, pat, what) assert(~isempty(regexp(txt, pat, "once")), "lmmSpec:parse", "Pattern not found: %s", what);
need(ks, "contVars = \{'betaGenerated', 'VGF', 'noiseMagnitude', 'noiseColor'\};", "z-scored predictors");
need(ks, "tableTrue\.\(v\) = \(tableTrue\.\(v\) - mu\) / sg;", "z-score = (x - mean)/SD");
need(ks, "deltaBeta ~ betaGenerated \* VGF \* samplingRate \* filterType \* regressionType \* noiseMagnitude \* noiseColor \+ \(1\|paramComboID\)", "formula");
need(ks, "tableTrue\.filterType = nominal\(tableTrue\.filterType\);", "filterType nominal");
need(ks, "tableTrue\.regressionType = nominal\(tableTrue\.regressionType\);", "regressionType nominal");
need(ks, "'FitMethod', 'REML'", "REML");
assert(isempty(regexp(ks, "DummyVarCoding", "once")), "lmmSpec:dummy", "DummyVarCoding is set explicitly; the default-coding statement needs revisiting");
f = fullfile(ROOT, "results", "L9ModelSlim_v001.mat"); if ~isfile(f), error("lmmSpec:slim", "%s", "FAILED PATH: " + f); end
S = load(f); cn = string(S.cn(:)); fx = S.fx(:); se = S.se(:); V = S.V;
order = count(cn, ":") + 1;  order(cn == "(Intercept)") = 0;
R = struct();
R.nCoef = numel(cn); R.nIntercept = nnz(order == 0); R.nMain = nnz(order == 1); R.nInter = nnz(order >= 2);
assert(R.nCoef == 192 && R.nIntercept == 1 && R.nMain == 8 && R.nInter == 183, "lmmSpec:count", "Coefficient count is not 1 + 8 + 183 = 192");
R.mainNames = cn(order == 1);
rng_ = @(n) V.Range{string(V.Properties.RowNames) == n};
R.fsRange = rng_("samplingRate");  assert(isequal(R.fsRange(:)', [60 240]), "lmmSpec:fs", "samplingRate is not on its raw 60-240 Hz scale");
for n = ["betaGenerated" "VGF" "noiseMagnitude" "noiseColor"]
    r = rng_(n); assert(r(1) < 0 && r(2) > 0, "lmmSpec:z", "%s is not centred", n);
end
R.filterLevels = double(string(rng_("filterType")));  R.regLevels = double(string(rng_("regressionType")));
assert(isequal(R.filterLevels(:)', [2 6]) && isequal(R.regLevels(:)', [3 4 5]), "lmmSpec:levels", "Unexpected categorical levels");
assert(all(ismember(["filterType_6" "regressionType_4" "regressionType_5"], cn)), "lmmSpec:dummies", "Expected reference-coded dummies absent");
R.nGroups = numel(rng_("paramComboID"));
z = fx ./ se;  p = 2*normcdf(-abs(z));
R.mainSigAll = all(p(order == 1) < 0.05); R.mainSig = nnz(p(order == 1) < 0.05);
R.interSig = nnz(p(order >= 2) < 0.05);
disp(R); save(fullfile(ROOT, "results", "checkLmmSpec_v001.mat"), "-struct", "R");
