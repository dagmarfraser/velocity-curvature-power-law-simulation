%% checkGainTargetLMMDomain_v001.m
% The target-vs-tempo LMM of Finding #234 on the five-dataset validity domain (Finding #242).
% extractGainTarget_v002.m (L41, L157) excluded Zarandi only (D20), so its "in-domain" LMM
% (1,711 trials, 243 subjects; log f0 b = 0.061, SE 0.0025, t(1702) = 24.9) still held Dhieb.
% Since Session 125, Dhieb is outside the domain (ground B). This refits the identical model,
% from the same saved per-trial table, with Dhieb also excluded. No curve is re-extracted.
% Reads results/gainTarget_v002.mat (P: every curve; m3: the published model).
% Self-check (error): the Zarandi-only refit reproduces m3's n and fixed effects exactly.
% Writes results/gainTargetLMMDomain_v001.mat. Run from src/.

%% CONFIG
ROOT    = fileparts(fileparts(mfilename("fullpath")));
IN_MAT  = fullfile(ROOT, "results", "gainTarget_v002.mat");
OUT_MAT = fullfile(ROOT, "results", "gainTargetLMMDomain_v001.mat");
BASES   = struct("name", {"domain6_F234", "domain5_F242"}, "out", {"Zarandi", ["Zarandi" "Dhieb"]});
TERMS   = ["lf" "lsa" "al"];
FORMULA = "target ~ 1 + lf + lsa + al + datasetC + (1|subjC)";     % extractGainTarget_v002 L161

%% Load
if ~isfile(IN_MAT), error("gainTargetLMM:input", "%s", "Missing: " + IN_MAT); end
S = load(IN_MAT, "P", "m3");
V = S.P(S.P.source == "v015", :);
z = @(x) (x - mean(x, "omitnan")) ./ std(x, "omitnan");

%% Fit on each basis (z-scoring within the fitted subset, as v002 L158-159)
M = cell(numel(BASES), 1);  t = table();
for i = 1:numel(BASES)
    D = V(V.reliable & V.pipeline == "SG-IRLS" & ~ismember(V.dataset, BASES(i).out) & isfinite(V.a) & V.a > 0, :);
    D.lf = z(log(D.f0));  D.lsa = z(log(D.sigma ./ D.a));  D.al = z(D.alpha);
    D.datasetC = categorical(D.dataset);  D.subjC = categorical(D.subject);
    M{i} = fitlme(D, FORMULA);
    cm = M{i}.Coefficients;  k = ismember(string(cm.Name), TERMS);
    t = [t; table(repmat(string(BASES(i).name), sum(k), 1), repmat(height(D), sum(k), 1), ...
        repmat(numel(unique(D.subjC)), sum(k), 1), string(cm.Name(k)), cm.Estimate(k), cm.SE(k), ...
        cm.tStat(k), cm.DF(k), cm.pValue(k), cm.Lower(k), cm.Upper(k), 'VariableNames', ...
        ["basis" "nTrials" "nSubjects" "term" "b" "SE" "t" "df" "p" "CIlo" "CIhi"])]; %#ok<AGROW>
    fprintf("%s: %d trials, %d subjects; datasets: %s\n", BASES(i).name, height(D), ...
        numel(unique(D.subjC)), strjoin(unique(D.dataset)', ", "));
end

%% Self-check: the Zarandi-only refit is the published model
e0 = S.m3.Coefficients.Estimate;  e1 = M{1}.Coefficients.Estimate;
if M{1}.NumObservations ~= S.m3.NumObservations || numel(e0) ~= numel(e1) || max(abs(e0 - e1)) > 1e-10
    error("gainTargetLMM:repro", "%s", "Zarandi-only refit does not reproduce the saved m3 (Finding #234)");
end
fprintf("Basis domain6_F234 reproduces the published model m3 exactly (n = %d).\n\n", S.m3.NumObservations);

disp(t)
if ~isfolder(fileparts(OUT_MAT)), error("gainTargetLMM:outDir", "%s", "Missing folder: " + fileparts(OUT_MAT)); end
save(OUT_MAT, "t", "M", "BASES", "FORMULA", "-v7.3");
fprintf("Saved: %s\n", OUT_MAT);
