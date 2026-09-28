%% compareL9CSVToModel_v001.m
% Finding #237: how far does src/L9_coefficients_v004.csv (extractL9Results_v004) depart
% from the fitted L9 Stage 1 model (results/L9ModelInspect_v002.mat, 192 fixed effects,
% verdict clean)? And do the terms that v002/v003 Part 1 Results (a) and Finding #43 cite
% (sourced from LMM_coefficients_top_L9_20260525_093556.txt, read from lme) match the
% model? Order-free keys: tokens sorted, so "a:b" and "b:a" are the same term.
% Run from src/:  compareL9CSVToModel_v001
% Writes results/L9CSVvsModel_v001.mat.
% Fraser, D.S. (2026)  v001

%% CONFIG
ROOT    = fileparts(fileparts(mfilename("fullpath")));
IN_MAT  = fullfile(ROOT, "results", "L9ModelInspect_v002.mat");
CSV     = fullfile(ROOT, "src", "L9_coefficients_v004.csv");
OUT_MAT = fullfile(ROOT, "results", "L9CSVvsModel_v001.mat");
TOL     = 1e-9;                                   % CSV text round-off
CITED = table( ...                                % term, value as printed in the draft
    ["betaGenerated"; "VGF"; "filterType_6"; "betaGenerated:noiseMagnitude"; ...
     "betaGenerated:noiseColor"; "betaGenerated:noiseColor:noiseMagnitude"; ...
     "noiseColor:regressionType_4"; "filterType_6:regressionType_4"; "filterType_6:regressionType_5"], ...
    [0.077; 0.007; 0.009; 0.025; -0.024; -0.024; -0.023; NaN; NaN], ...
    ["β=+0.077"; "β=+0.007"; "β=+0.009"; "β=+0.025"; "β=-0.024"; "β=-0.024"; "β=-0.023"; "p=.33"; "p=.25"], ...
    'VariableNames', ["key" "drafted" "draftedText"]);

%% Load
for f = [IN_MAT CSV]
    if ~isfile(f), error("cmpL9:input", "%s", "Missing: " + f); end
end
if ~isfolder(fileparts(OUT_MAT)), error("cmpL9:outDir", "%s", "Missing folder: " + fileparts(OUT_MAT)); end
I = load(IN_MAT, "out");  F = I.out.fixed;
if height(F) ~= 192, error("cmpL9:model", "%s", sprintf("Model table has %d rows, expected 192", height(F))); end
C = readtable(CSV, "TextType", "string");
C.key = canon_local(C.Name);

%% 1. Row-level agreement
[ok, loc] = ismember(C.key, F.key);
if ~all(ok), error("cmpL9:keys", "%s", "CSV keys absent from the model: " + strjoin(C.Name(~ok), ", ")); end
C.trueEst = F.est(loc);  C.trueSE = F.SE(loc);  C.d = C.Estimate - C.trueEst;
C.sameName = C.Name == F.name(loc);  C.sameEst = abs(C.d) <= TOL;
C.signFlip = sign(C.Estimate) ~= sign(C.trueEst) & ~C.sameEst;
anyVerbatim = sum(ismember(round(C.Estimate, 12), round(F.est, 12)));
[uk, ~, j] = unique(C.key);  nDupKeys = sum(accumarray(j, 1) > 1);
fprintf("CSV %d rows vs model %d coefficients; %d unique keys (%d duplicated); %d model terms absent from CSV\n", ...
    height(C), height(F), numel(uk), nDupKeys, numel(setdiff(F.key, C.key)));
fprintf("Same name and estimate: %d | same name, different estimate: %d | token order differs: %d\n", ...
    sum(C.sameName & C.sameEst), sum(C.sameName & ~C.sameEst), sum(~C.sameName));
fprintf("CSV estimates equal to the model's for the same term: %d | occurring anywhere in the model: %d | sign flips: %d\n", ...
    sum(C.sameEst), anyVerbatim, sum(C.signFlip));
isMain = ~contains(C.Name, ":");
fprintf("Main effects (incl. intercept) matching: %d/%d | interactions matching: %d/%d\n", ...
    sum(C.sameEst & isMain), sum(isMain), sum(C.sameEst & ~isMain), sum(~isMain));

%% 2. Terms cited in the draft and #43
[ok2, l2] = ismember(CITED.key, F.key);
if ~all(ok2), error("cmpL9:cited", "%s", "Cited terms absent from the model: " + strjoin(CITED.key(~ok2), ", ")); end
CITED.modelEst = F.est(l2);  CITED.modelP = 2 * normcdf(-abs(F.est(l2) ./ F.SE(l2)));
CITED.csvEst = NaN(height(CITED), 1);
for i = 1:height(CITED)
    c = C.Estimate(C.key == CITED.key(i));
    if isscalar(c), CITED.csvEst(i) = c; elseif numel(c) > 1, CITED.csvEst(i) = NaN; end   % absent or duplicated
end
CITED.draftOK = abs(round(CITED.modelEst, 3) - CITED.drafted) < 5e-4 | isnan(CITED.drafted);
fprintf("\nTerms cited in Part 1 Results (a) / #43 (draft value vs fitted model; CSV for comparison):\n");
disp(CITED(:, ["key" "draftedText" "modelEst" "modelP" "csvEst" "draftOK"]));

%% Save
out = struct("rows", C, "cited", CITED, "nDupKeys", nDupKeys, "anyVerbatim", anyVerbatim, ...
    "absentFromCSV", setdiff(F.key, C.key), "source", [string(IN_MAT) string(CSV)], "runDate", string(datetime("now")));
save(OUT_MAT, "out");
fprintf("Saved %s\n", OUT_MAT);

function k = canon_local(names)
    names = string(names);  k = strings(size(names));
    for i = 1:numel(names)
        if names(i) == "(Intercept)", k(i) = names(i); continue, end
        k(i) = strjoin(sort(split(names(i), ":")), ":");
    end
end
