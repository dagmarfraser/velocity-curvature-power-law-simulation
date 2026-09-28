%% inspectL9Model_v002.m
% Integrity check of the fitted L9 Stage 1 LMM, run locally from a slim extract.
% v002: v001 loaded lme on BlueBEAR (10 h single core) and stalled after the load, at or
% after its L43 (disp of lme.Formula). The small pieces were saved in the same session:
%   cn  coefficient names          fx  fixedEffects(lme)       se  sqrt(diag(CoefficientCovariance))
%   tn  FELinearFormula.TermNames  Tm  FELinearFormula.Terms   vn  FELinearFormula.VariableNames
%   V   lme.VariableInfo restricted to vn
% -> results/L9ModelSlim_v001.mat (provenance: interactive extraction from
%    src/extractL9_checkpoint_v004.mat, 2026-09-28; commands in SESSION_LOG).
% Checks as v001: terms in the formula and the coefficients they imply (registered
% 128 terms, 192 coefficients); coefficients estimated, name uniqueness, order-free
% coverage of the registered 192; identity with L9_coefficients_v004.csv; verdict.
% Writes results/L9ModelInspect_v002.mat, the P-1 gate target in place of the CSV.
% Run from src/:  inspectL9Model_v002  (paths resolve from this file).
% Fraser, D.S. (2026)  v002

%% CONFIG
ROOT   = fileparts(fileparts(mfilename("fullpath")));
SLIM   = fullfile(ROOT, "results", "L9ModelSlim_v001.mat");
CSV    = fullfile(ROOT, "src", "L9_coefficients_v004.csv");
OUT    = fullfile(ROOT, "results", "L9ModelInspect_v002.mat");
TOKENS = ["betaGenerated" "VGF" "samplingRate" "filterType_6" "noiseMagnitude" "noiseColor"];
REGS   = ["" "regressionType_4" "regressionType_5"];
NOBS   = 17458535;                                    % Finding #189; printed by v001 on BlueBEAR

%% Load the slim extract
if ~isfile(SLIM), error("inspL9:slim", "%s", "Missing: " + SLIM); end
if ~isfolder(fileparts(OUT)), error("inspL9:outDir", "%s", "Missing folder: " + fileparts(OUT)); end
S = load(SLIM);
need = ["cn" "fx" "se" "tn" "Tm" "vn" "V"];
if ~all(isfield(S, need)), error("inspL9:fields", "%s", "Slim extract lacks: " + strjoin(need(~isfield(S, need)), ", ")); end
cn = S.cn(:);  fx = S.fx(:);  se = S.se(:);  termNames = S.tn(:);  Tm = S.Tm;  vn = S.vn(:);  V = S.V;
if numel(fx) ~= numel(cn) || numel(se) ~= numel(cn), error("inspL9:sizes", "%s", "cn, fx, se lengths differ"); end
fprintf("Slim extract: %s\n", SLIM);

%% 1. Terms in the fitted formula and the coefficients they imply
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
    "termNames", termNames, "impliedCoef", sum(implied), ...
    "absent", miss, "unexpected", extra, "duplicates", dup, "verdict", v, ...
    "source", string(SLIM), "nObs", NOBS, "runDate", string(datetime("now")));
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
