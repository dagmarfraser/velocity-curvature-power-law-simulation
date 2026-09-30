%% checkCookTempoInteraction_v001.m
% Coherence item Q4 (v004 Part 4, the Cook group contrast). Two checks on Finding #235's
% attribution, from analyseTempoDecomposition_v002's saved trial table:
%   (1) Missing statistics. The tempo contrasts are reported as "+0.50 log2 units" (Cook,
%       autistic minus non-autistic) and "+0.03 log2 units, p = .50" (Hickman, haloperidol
%       minus placebo) without test, unit of analysis or df. Printed here from the saved
%       models (section 3's "log2 f0 (tempo)" rows), both pipelines.
%   (2) Separability. The Shapley split (97% tempo, 3% noise) treats tempo and noise as
%       additive. Finding #7 found the groups' noise colour responds to tempo differently
%       (ANCOVA f0 x group). Does the attribution survive a tempo x group term?
%         V0  analyseTempoDecomposition_v002's own models   regression anchor: must reproduce
%             the saved A3 rows for Cook (b_only, b_tempo, b_noise, b_both, shares; 4 d.p.).
%         V1  the same with lf centred on the Cook pooled mean and lf x group added to the
%             tempo block (Shapley split over blocks {lf, lf x group} and {lsa, al}); the group
%             coefficient is then the group difference at the pooled mean tempo.
%       Also reported: the lf x group coefficient itself (does beta's tempo slope differ by group?)
%       and the same interaction on alpha (Finding #7's premise, trial level).
% Unit of analysis: trial, subject random intercept, as analyseTempoDecomposition_v002.
% Reads:  results/tempoDecomposition_v002.mat
% Writes: results/checkCookTempoInteraction_v001.mat
% USAGE:  from the project root: checkCookTempoInteraction_v001
% Fraser, D.S. (2026)  v001

%% CONFIG
ROOT   = fileparts(fileparts(mfilename("fullpath")));
IN_MAT = fullfile(ROOT, "results", "tempoDecomposition_v002.mat");
OUT_MAT = fullfile(ROOT, "results", "checkCookTempoInteraction_v001.mat");
TOL = 1e-4;                                            % A3 was saved rounded to 4 d.p.
if ~isfile(IN_MAT), error("cookInt:input", "%s", "FAILED PATH: " + IN_MAT); end
S = load(IN_MAT, "T", "C3", "A3", "PIPES");

%% (1) Tempo contrasts with their statistics
tc = S.C3(S.C3.outcome == "log2 f0 (tempo)", :);
fprintf("(1) Tempo contrasts (log2 f0 ~ group/drug + (1|subject); trial as unit):\n");
disp(tc)

%% (2) Separability
R = table();  I = table();
for p = S.PIPES
    V = S.T(S.T.pipeline == p & ismember(S.T.dataset, ["Cook_CTRL" "Cook_ASD"]), :);
    V.group = double(V.dataset == "Cook_ASD");
    V.lfc = V.lf - mean(V.lf);
    for out = ["betaObs" "betaGen"]
        W = V(isfinite(V.(out)), :);
        g = @(rhs) coef_local(fitlme(W, out + " ~ 1 + group" + rhs + " + (1|personC)"), "group");
        % V0: as analyseTempoDecomposition_v002 L185-189
        c0 = g("");  ct = g(" + lf");  cn = g(" + lsa + al");  cb = g(" + lf + lsa + al");
        [shT0, shN0] = shapley_local(c0.b, ct.b, cn.b, cb.b);
        % V1: tempo block with the interaction, lf centred
        dt = g(" + lfc + lfc:group");  db = g(" + lfc + lfc:group + lsa + al");
        [shT1, shN1] = shapley_local(c0.b, dt.b, cn.b, db.b);
        R = [R; table(p, out, height(W), c0.b, ct.b, cn.b, cb.b, shT0, shN0, dt.b, db.b, db.CIlo, db.CIhi, shT1, shN1, ...
            'VariableNames', ["pipeline" "outcome" "nTrials" "b_only" "b_tempo" "b_noise" "b_both" "tempoShare" "noiseShare" ...
            "b_tempoInt" "b_bothInt" "CIlo_bothInt" "CIhi_bothInt" "tempoShareInt" "noiseShareInt"])]; %#ok<AGROW>
        mi = fitlme(W, out + " ~ 1 + group + lfc + lfc:group + lsa + al + (1|personC)");
        I = [I; [table(p, out, 'VariableNames', ["pipeline" "outcome"]), coef_local(mi, ["group:lfc" "lfc:group"])]]; %#ok<AGROW>
    end
    W = V(isfinite(V.alpha), :);                             % Finding #7's premise, trial level
    ma = fitlme(W, "alpha ~ 1 + group + lfc + lfc:group + (1|personC)");
    I = [I; [table(p, "alpha (noise colour)", 'VariableNames', ["pipeline" "outcome"]), coef_local(ma, ["group:lfc" "lfc:group"])]]; %#ok<AGROW>
end

%% Regression anchor: V0 must reproduce the saved A3 (Cook rows)
A = S.A3(startsWith(S.A3.contrast, "Cook"), :);
for r = 1:height(A)
    j = find(R.pipeline == A.pipeline(r) & R.outcome == A.outcome(r));
    got  = [R.b_only(j) R.b_tempo(j) R.b_noise(j) R.b_both(j) R.tempoShare(j) R.noiseShare(j)];
    want = [A.b_only(r) A.b_tempo(r) A.b_noise(r) A.b_both(r) A.tempoShare(r) A.noiseShare(r)];
    if numel(j) ~= 1 || any(abs(got - want) > TOL)
        error("cookInt:anchor", "%s", sprintf("%s %s: V0 %s does not reproduce A3 %s", A.pipeline(r), A.outcome(r), mat2str(got, 4), mat2str(want, 4)));
    end
end
fprintf("REGRESSION ANCHOR passed: V0 reproduces the saved Cook attribution (A3, %d rows).\n\n", height(A));

fprintf("(2) Group coefficient (autistic minus non-autistic): alone, + tempo, + noise, + both; V1 adds lf x group\n");
fprintf("    (lf centred, so the V1 group coefficient is the difference at the pooled mean tempo).\n");
R{:, 4:end} = round(R{:, 4:end}, 4);
disp(R)
fprintf("Interaction terms (group x centred log2 f0):\n");
disp(I)

runDate = string(datetime("now", "Format", "yyyy-MM-dd"));
save(OUT_MAT, "tc", "R", "I", "runDate");
fprintf("Saved: %s\n", OUT_MAT);

%% =========================================================================
function [shT, shN] = shapley_local(c0, ct, cn, cb)
% analyseTempoDecomposition_v002 L187-189: two-term Shapley split of the attenuation.
    tot = c0 - cb;
    shT = ((c0 - ct) + (cn - cb)) / 2;  shN = ((c0 - cn) + (ct - cb)) / 2;
    if abs(shT + shN - tot) > 1e-12, error("cookInt:shapley", "%s", "shares do not sum to the attenuation"); end
end

function t = coef_local(mdl, terms)
% analyseTempoDecomposition_v002 L270-274 (rounded to 4 d.p. as saved there).
    c = mdl.Coefficients;  k = ismember(string(c.Name), terms);
    if ~any(k), error("cookInt:term", "%s", "term not in model: " + strjoin(string(terms), ", ")); end
    t = table(string(c.Name(k)), round(c.Estimate(k), 4), round(c.SE(k), 4), round(c.tStat(k), 2), c.DF(k), ...
        c.pValue(k), round(c.Lower(k), 4), round(c.Upper(k), 4), 'VariableNames', ["term" "b" "SE" "t" "df" "p" "CIlo" "CIhi"]);
end
