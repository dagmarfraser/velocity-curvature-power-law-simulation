%% crossReferenceSimexMonotonicity_v002.m
% Findings #219/#220 re-joined against the rebuilt monotonicity classification (Finding #244).
% v001 joined each SIMEX corpus to src/results/queryMonotonicSlopeAtEmpirical_v003.mat, whose
% per-trial invertibility placed each trial on the grid's VGF axis with K_conv = 175.2636
% (5.28% low, blind to beta). v002 joins every corpus to three classifications on identical
% SIMEX rows, so any change is attributable:
%   V0 RETIRED   v003's perTrial                       regression anchor: must reproduce #219/#220's
%                                                      Cohen's d per reporting group (3 d.p.).
%   V1 CONSTANT  queryMonotonicSlopeAtEmpirical_v005, variant 2 (K(1/3))
%   V2 PRIMARY   queryMonotonicSlopeAtEmpirical_v005, variant 3 (K(beta_i), own exponent per trial)
% V1/V2 rows for trials that could not be placed (no exponent; v005's `placed` = false) are
% excluded and counted, never scored as non-invertible. Invertible = invRise, the definition
% #220 validated against #219. Join mechanics unchanged from v001 (trialIdx = bioResults row
% order -> trialID; join on trialID + pipeline; SIMEX successes only).
% Reads:  results/simexBetaCorpus_*.mat (listed in CORPORA), noiseCharacterisation_*.mat,
%         src/results/queryMonotonicSlopeAtEmpirical_v003.mat, results/queryMonotonicSlopeAtEmpirical_v005.mat
% Writes: results/crossReferenceSimexMonotonicity_v002.mat
% USAGE:  from the project root: crossReferenceSimexMonotonicity_v002
% Fraser, D.S. (2026)  v002. Exploratory, not preregistered.

%% CONFIG
ROOT = fileparts(fileparts(mfilename("fullpath")));
addpath(fullfile(ROOT, "src")); addpath(genpath(fullfile(ROOT, "src", "functions")));
% corpus file, reporting group label, simex datasets in it, published d under v001 (#219/#220)
CORPORA = {
    "simexBetaCorpus_hickman_20260921_002210.mat",              "Hickman (pre-fix, cited)", "hickman", 0.270
    "simexBetaCorpus_hickmanPLAC-hickmanHALO_20260921_123014.mat", "Hickman (post-fix)",    "hickman", 0.265
    "simexBetaCorpus_fraser_20260921_112738.mat",               "Fraser",                   "fraser",  0.110
    "simexBetaCorpus_cookCTRL-cookASD_20260921_103217.mat",     "Cook CTRL+ASD",            "cook",    0.039
    "simexBetaCorpus_zarandi-dhieb_20260921_110230.mat",        "Dhieb",                    "dhieb",  -0.021
    "simexBetaCorpus_zarandi-dhieb_20260921_110230.mat",        "Zarandi",                  "zarandi", -0.235 };
REG = struct('simexDataset', {'fraser','zarandi','cook','cook','dhieb','hickman','hickman'}, ...
    'simexGroup', {'','','ctrl','asd','','plac','halo'}, ...
    'monoName', {"Fraser","Zarandi","Cook_CTRL","Cook_ASD","Dhieb","Hickman_PLAC","Hickman_HALO"}, ...
    'noiseFile', {'noiseCharacterisation_fraser.mat','noiseCharacterisation_zarandi.mat','noiseCharacterisation_cook.mat', ...
    'noiseCharacterisation_cookASD.mat','noiseCharacterisation_dhieb.mat','noiseCharacterisation_hickmanPLAC.mat','noiseCharacterisation_hickmanHALO.mat'});
MONO_V003 = fullfile(ROOT, "src", "results", "queryMonotonicSlopeAtEmpirical_v003.mat");
MONO_V005 = fullfile(ROOT, "results", "queryMonotonicSlopeAtEmpirical_v005.mat");
OUT_MAT   = fullfile(ROOT, "results", "crossReferenceSimexMonotonicity_v002.mat");
VNAME = ["V0 RETIRED", "V1 CONSTANT", "V2 PRIMARY"];
for f = [MONO_V003 MONO_V005], if ~isfile(f), error("simexX2:mono", "%s", "FAILED PATH: " + f); end, end

%% Monotonicity sources
P0 = load(MONO_V003, "perTrial").perTrial;  P0.placed = true(height(P0), 1);
Q  = load(MONO_V005, "perTrial").perTrial;
if numel(Q) ~= 3, error("simexX2:v005", "%s", "v005 perTrial must hold three variants"); end
MONO = {P0, Q{2}, Q{3}};

%% trialIdx -> trialID per monoName (bioResults row order, as v001)
idMap = containers.Map();
for r = 1:numel(REG)
    nf = fullfile(ROOT, "src", REG(r).noiseFile);
    if ~isfile(nf), error("simexX2:noise", "%s", "FAILED PATH: " + nf); end
    idMap(char(REG(r).monoName)) = string(load(nf, "bioResults").bioResults.trialID);
end

%% Join and score
res = table();  allJoined = cell(size(CORPORA, 1), 3);
for c = 1:size(CORPORA, 1)
    cf = fullfile(ROOT, "results", CORPORA{c, 1});
    if ~isfile(cf), error("simexX2:corpus", "%s", "FAILED PATH: " + cf); end
    S = load(cf, "summaryTable").summaryTable;
    if ~ismember("dataset", S.Properties.VariableNames), S.dataset = repmat({'hickman'}, height(S), 1); end
    S = S(S.success & strcmpi(string(S.dataset), CORPORA{c, 3}), :);
    regs = REG(strcmpi({REG.simexDataset}, CORPORA{c, 3}));
    for v = 1:3
        J = {};  nUnplaced = 0;
        for r = 1:numel(regs)
            M = MONO{v}(strcmp(string(MONO{v}.dataset), regs(r).monoName), :);
            ids = idMap(char(regs(r).monoName));
            if isempty(M) || max(M.trialIdx) > numel(ids)
                error("simexX2:map", "%s", sprintf("%s %s: no rows or trialIdx beyond bioResults", VNAME(v), regs(r).monoName));
            end
            M.trialID = ids(M.trialIdx);
            Ssub = S(strcmpi(string(S.group), regs(r).simexGroup) | regs(r).simexGroup == "", :);
            j = innerjoin(M(:, ["trialID" "pipeline" "invRise" "placed"]), Ssub, "Keys", ["trialID" "pipeline"], ...
                "RightVariables", ["betaSimexCI_lo" "betaSimexCI_hi"]);
            nUnplaced = nUnplaced + nnz(~j.placed);
            J{end+1} = j(j.placed, :); %#ok<SAGROW>
        end
        J = vertcat(J{:});  J.ciWidth = J.betaSimexCI_hi - J.betaSimexCI_lo;
        allJoined{c, v} = J;
        s = scoreJoin_local(J);
        res = [res; [table(string(CORPORA{c, 2}), VNAME(v), height(S), nUnplaced, CORPORA{c, 4}, ...
            'VariableNames', ["group" "variant" "nSimex" "nUnplaced" "dPublished"]), s]]; %#ok<AGROW>
    end
end

%% Regression anchor
A = res(res.variant == "V0 RETIRED", :);
if any(abs(round(A.d, 3) - A.dPublished) > 0.0005)
    disp(A(:, ["group" "d" "dPublished"]));
    error("simexX2:anchor", "%s", "V0 does not reproduce the published Cohen's d (#219/#220)");
end
fprintf("REGRESSION ANCHOR passed: V0 reproduces #219/#220's d for all %d groups.\n\n", height(A));
disp(res(:, ["group" "variant" "nUnplaced" "nInv" "nNon" "meanInv" "meanNon" "d" "pWelch" "pWilcoxon" "pctInvOverNonMed" "pctNonOverInvMed"]))

runDate = string(datetime("now", "Format", "yyyy-MM-dd"));
save(OUT_MAT, "res", "allJoined", "CORPORA", "VNAME", "runDate", "-v7.3");
fprintf("Saved: %s\n", OUT_MAT);

%% =========================================================================
function s = scoreJoin_local(J)
% v001 L181-214 for the invRise definition.
    inv = J.invRise;  a = J.ciWidth(inv);  b = J.ciWidth(~inv);
    a = a(isfinite(a));  b = b(isfinite(b));
    if isempty(a) || isempty(b), error("simexX2:emptyClass", "%s", "an invertibility class is empty"); end
    sp = sqrt(((numel(a) - 1) * var(a) + (numel(b) - 1) * var(b)) / (numel(a) + numel(b) - 2));
    [~, pT] = ttest2(b, a, "Vartype", "unequal");
    s = table(numel(a), numel(b), mean(a), mean(b), (mean(b) - mean(a)) / sp, pT, ranksum(b, a), ...
        100 * mean(a > median(b)), 100 * mean(b > median(a)), 'VariableNames', ...
        ["nInv" "nNon" "meanInv" "meanNon" "d" "pWelch" "pWilcoxon" "pctInvOverNonMed" "pctNonOverInvMed"]);
end
