%% inspectL9Model_v001.m
% Integrity check of the fitted L9 Stage 1 LMM (ModelAdequacy_Stage1_KitchenSink_v2_001).
% L9_coefficients_v004.csv (= lme.Coefficients, extractL9Results_v004 L89-109) has 161
% rows, not the registered 192: one exact duplicate row, four terms twice with different
% estimates, 36 registered terms absent (checkLMMTempoConfound_v001 gate, 2026-09-27).
% This script asks the fitted object, not the CSV:
%   1. Formula as fitted; fixed-effect terms in the formula and the coefficient count they
%      imply (each term's count = product over categorical factors of (levels - 1)).
%   2. Coefficients actually estimated; name uniqueness; order-free key coverage of the
%      registered 192 (2^6 subsets of betaGenerated, VGF, samplingRate, filterType_6,
%      noiseMagnitude, noiseColor, times regressionType none/4/5).
%   Verdict: formula implies 192 but fewer estimated -> columns dropped at fit (warnings
%   were off, Stage 1 L955); formula itself short -> formula; names clean here but not in
%   the CSV -> extraction.
% Writes results/L9ModelInspect_v001.mat: fixed effects (name, key, est, SE), formula,
% scaling. That mat replaces the CSV as the P-1 gate target.
% Run on BlueBEAR from src/:  cd src; inspectL9Model_v001  (paths resolve from this file;
% usage line corrected in place 2026-09-28)
% Fraser, D.S. (2026)  v001

%% CONFIG
ROOT   = fileparts(fileparts(mfilename("fullpath")));
SRC    = fullfile(ROOT, "src");
CHK    = fullfile(SRC, "extractL9_checkpoint_v004.mat");      % slim cache holding lme
CSV    = fullfile(SRC, "L9_coefficients_v004.csv");
OUT    = fullfile(ROOT, "results", "L9ModelInspect_v001.mat");
TOKENS = ["betaGenerated" "VGF" "samplingRate" "filterType_6" "noiseMagnitude" "noiseColor"];
REGS   = ["" "regressionType_4" "regressionType_5"];

%% Load the fitted model
if ~isfolder(fileparts(OUT)), error("inspL9:outDir", "%s", "Missing folder: " + fileparts(OUT)); end
scaling = [];
if isfile(CHK)
    C = load(CHK, "lme");  lme = C.lme;  srcFile = CHK;
else
    m = dir(fullfile(SRC, "stage1_results_L9_*.mat"));
    if isempty(m), error("inspL9:noModel", "%s", "Neither " + CHK + " nor stage1_results_L9_*.mat found"); end
    [~, k] = max([m.datenum]);  srcFile = string(fullfile(m(k).folder, m(k).name));
    fprintf("Checkpoint absent; loading 'results' from %s (large) ...\n", srcFile);
    S = load(srcFile, "results");  lme = S.results.model;
    if isfield(S.results, "predictor_scaling"), scaling = S.results.predictor_scaling; end
end
fprintf("Model from: %s\nClass: %s   N obs: %d\n", srcFile, class(lme), lme.NumObservations);
fTxt = strtrim(string(evalc("disp(lme.Formula)")));
fprintf("Formula: %s\n", fTxt);

%% 1. Terms in the fitted formula and the coefficients they imply
fe = lme.Formula.FELinearFormula;
termNames = string(fe.TermNames(:));
V  = lme.VariableInfo;
vn = string(fe.VariableNames(:));                     % columns of fe.Terms
Tm = fe.Terms;                                        % terms x variables (exponents)
if size(Tm, 2) ~= numel(vn), error("inspL9:terms", "%s", "Terms matrix and VariableNames disagree"); end
[ok, iv] = ismember(vn, string(V.Properties.RowNames));
if ~all(ok), error("inspL9:vars", "%s", "Formula variables missing from VariableInfo: " + strjoin(vn(~ok), ", ")); end
isCat = V.IsCategorical(iv);  nLev = ones(numel(vn), 1);
for i = find(isCat)'
    nLev(i) = numel(V.Range{iv(i)});
end
implied = zeros(size(Tm, 1), 1);
for t = 1:size(Tm, 1)
    implied(t) = prod(max(nLev(Tm(t, :) > 0 & isCat') - 1, 1), "all");   % "all": empty -> 1
end
fprintf("\n1. Formula: %d fixed-effect terms (registered 7-way factorial: 128), implying %d coefficients (registered: 192)\n", ...
    numel(termNames), sum(implied));
fprintf("   Categorical: %s\n", strjoin(compose("%s (%d levels)", vn(isCat), nLev(isCat)), ", "));

%% 2. Coefficients actually estimated
cn  = string(lme.CoefficientNames(:));
fx  = fixedEffects(lme);
se  = sqrt(diag(lme.CoefficientCovariance));
key = canon_local(cn);
[uk, ~, j] = unique(key);  dup = uk(accumarray(j, 1) > 1);
fprintf("2. Estimated: %d coefficients; %d unique names; %d unique order-free keys\n", ...
    numel(cn), numel(unique(cn)), numel(uk));
for d = dup'
    fprintf("   DUPLICATE KEY %s:\n", d);  disp(table(cn(key == d), fx(key == d), se(key == d), 'VariableNames', ["name" "est" "SE"]));
end
reg = registered_local(TOKENS, REGS);
miss = setdiff(reg, key);  extra = setdiff(key, reg);
fprintf("   Registered keys present: %d/%d; absent: %d; unexpected: %d\n", ...
    numel(intersect(reg, key)), numel(reg), numel(miss), numel(extra));
if ~isempty(miss), fprintf("   ABSENT:\n");  fprintf("     %s\n", miss); end
if ~isempty(extra), fprintf("   UNEXPECTED:\n");  fprintf("     %s\n", extra); end

%% 3. Against the CSV
if isfile(CSV)
    L = readtable(CSV, "TextType", "string");
    same = height(L) == numel(cn) && all(L.Name == cn) && max(abs(L.Estimate - fx)) < 1e-12;
    fprintf("3. CSV (%d rows) identical to the fitted object's names and estimates: %d\n", height(L), same);
else
    fprintf("3. CSV not found at %s; comparison skipped\n", CSV);
end

%% Verdict (printed, not decided silently)
if numel(unique(cn)) < numel(cn)
    v = "fitted object itself has repeated coefficient names";
elseif sum(implied) > numel(cn)
    v = "formula implies more coefficients than were estimated: columns dropped at fit";
elseif sum(implied) < 192
    v = "formula as fitted is not the registered 7-way factorial";
else
    v = "fitted object is clean; the CSV anomaly arose downstream";
end
fprintf("\nVERDICT: %s\n", v);

%% Save
out = struct("fixed", table(cn, key, fx, se, 'VariableNames', ["name" "key" "est" "SE"]), ...
    "formula", fTxt, "termNames", termNames, "impliedCoef", sum(implied), ...
    "absent", miss, "unexpected", extra, "duplicates", dup, "verdict", v, "scaling", scaling, ...
    "source", string(srcFile), "nObs", lme.NumObservations, "runDate", string(datetime("now")));
save(OUT, "out");
fprintf("Saved %s\n", OUT);

%% Local functions
function k = canon_local(names)
    k = strings(size(names));
    for i = 1:numel(names)
        if names(i) == "(Intercept)", k(i) = names(i); continue, end
        k(i) = strjoin(sort(split(names(i), ":")), ":");
    end
end

function k = registered_local(tokens, regs)
    k = strings(0, 1);
    for m = 0:2^numel(tokens) - 1
        t = tokens(bitget(m, 1:numel(tokens)) == 1);
        for r = regs
            s = [t r];  s = s(s ~= "");
            if isempty(s), k(end + 1, 1) = "(Intercept)"; else, k(end + 1, 1) = strjoin(sort(s), ":"); end %#ok<AGROW>
        end
    end
end
