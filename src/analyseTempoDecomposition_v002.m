%% analyseTempoDecomposition_v002.m
% Tempo enters the data two ways (Session 123 discussion, claude.md v003 row item (f)):
%   MEASUREMENT: a trial's tempo sets where its forward map pulls beta_obs (Finding #234).
%                Pointwise inversion builds each map at the trial's own f0, so beta_gen*
%                should be free of it.
%   BEHAVIOUR:   faster drawing may genuinely generate a different exponent, and tempo can
%                be part of a condition (bradykinesia; haloperidol) -> a mediator, not a
%                nuisance, whose effect must not be adjusted away silently.
% This script separates them with existing results (v015 per-trial corpora; no simulation):
%   1. Within/between decomposition. log2 f0 is split into a within-session deviation
%      (trial minus the subject-session mean: population, device and paradigm held
%      constant) and a between-session term. beta_obs and beta_gen* are modelled
%      separately. Predictions: beta_obs slopes with within-session tempo (the sliding
%      target); beta_gen* is flat if inversion removes the measurement effect. A residual
%      beta_gen* slope is behavioural or model error.
%   2. Selection. beta_gen* exists only for invertible trials (NaN otherwise), and
%      invertibility rises with tempo (#227), so dataset beta_gen* rests on a tempo-selected
%      subset. Quantified per dataset.
%   3. Clinical contrasts with and without tempo:
%      - Cook CTRL vs ASD: independent groups; tempo treated as a confounder of measurement.
%      - Hickman PLAC vs HALO: paired crossover, subjects in both arms; tempo treated as a
%        possible MEDIATOR. The total effect is drug only; the direct effect adds tempo; the
%        difference is the tempo-carried part (linear, no-interaction assumptions).
% Units: tempo terms are per doubling of f0 (log2); noise covariates are z-scored.
% Unit of analysis: trial, with a subject(-session) random intercept.
% v002 (adds section 3b; sections 1-3 unchanged from v001):
%   - attribution: each contrast refitted as group/drug only, + tempo, + noise, + both. The
%     attenuation is split into tempo and noise shares by averaging each covariate's
%     contribution over both entry orders (a two-term Shapley split).
%   - Cook, tempo-matched: coarsened exact matching on log2 f0 (bins MATCH_BIN wide); only
%     bins with both groups are kept, and trials are weighted so each group's tempo
%     distribution equals the pooled one. Refitted on beta_gen* (invertible trials) and
%     beta_obs (all trials). This separates selection (section 2b: inverted ASD trials are
%     faster than inverted CTRL trials) from everything else.
% Reads:  src/loopClosureResults_<dataset>_all_shaped_xu_v015.mat
% Writes: results/tempoDecomposition_v002.mat, figures/tempoDecomposition_v002.png

%% CONFIG
ROOT     = fileparts(fileparts(mfilename("fullpath")));
OUT_MAT  = fullfile(ROOT, "results", "tempoDecomposition_v002.mat");
OUT_PNG  = fullfile(ROOT, "figures", "tempoDecomposition_v002.png");
MATCH_BIN = 0.5;                                                               % log2 f0 bin width for tempo matching
DATASETS = ["Fraser" "Cook_CTRL" "Cook_ASD" "Hickman_PLAC" "Hickman_HALO" "Dhieb"];   % validity domain (Zarandi out, D20)
PL       = ["BWFD-OLS" "SG-OLS" "BWFD-LMLS" "SG-LMLS" "BWFD-IRLS" "SG-IRLS"];         % v015 runner order (v012 L118)
PIPES    = ["SG-IRLS" "BWFD-OLS"];                                                     % primary, legacy
PAIRED   = ["Hickman_PLAC" "Hickman_HALO"];                                            % crossover arms, shared subject IDs
TERMS    = ["lfw" "lfb" "lf" "group" "drug" "lsa" "al"];
for f = [OUT_MAT OUT_PNG], if ~isfolder(fileparts(f)), error("tempoDecomp:outDir", "%s", "Missing folder: " + fileparts(f)); end, end

%% Trial table (one row per trial x pipeline)
rows = {};
for d = DATASETS
    f = fullfile(ROOT, "src", "loopClosureResults_" + d + "_all_shaped_xu_v015.mat");
    if ~isfile(f), error("tempoDecomp:input", "%s", "Missing: " + f); end
    R = load(f, "results").results;
    person = d;  if ismember(d, PAIRED), person = "Hickman"; end                 % same person across arms
    for i = 1:numel(R)
        r = R(i);
        for p = PIPES
            k = find(PL == p);
            rows{end+1, 1} = struct("dataset", d, "session", d + "_" + string(r.subjectID), ...
                "person", person + "_" + string(r.subjectID), "pipeline", p, "f0", r.f0, "sigma", r.sigmaMM, ...
                "a", r.a_mm, "alpha", r.alphaIRA, "betaObs", r.betaObs(k), "betaGen", r.betaGenStar(k), ...
                "invertible", r.invertStatus(k) == "rise"); %#ok<SAGROW>
        end
    end
end
T = struct2table(vertcat(rows{:}));
T = T(isfinite(T.f0) & T.f0 > 0 & isfinite(T.sigma) & isfinite(T.a) & T.a > 0 & isfinite(T.alpha), :);
T.lf  = log2(T.f0);
G = groupsummary(T(T.pipeline == PIPES(1), :), "session", "mean", "lf");         % session means (pipeline-free)
[~, loc] = ismember(T.session, G.session);  T.lfMean = G.mean_lf(loc);
T.lfw = T.lf - T.lfMean;                                                        % within-session deviation
T.lfb = T.lfMean - mean(G.mean_lf);                                             % between-session (grand-centred)
z = @(x) (x - mean(x, "omitnan")) ./ std(x, "omitnan");
T.lsa = z(log2(T.sigma ./ T.a));  T.al = z(T.alpha);
T.datasetC = categorical(T.dataset);  T.sessionC = categorical(T.session);  T.personC = categorical(T.person);
fprintf("Trials x pipelines: %d (%d datasets, %d sessions)\n", height(T), numel(DATASETS), numel(unique(T.session)));

%% 1. Within/between decomposition, pooled in-domain, per pipeline
D1 = table();  M1 = struct();
for p = PIPES
    V = T(T.pipeline == p, :);
    mo = fitlme(V(isfinite(V.betaObs), :), "betaObs ~ 1 + lfw + lfb + lsa + al + datasetC + (1|sessionC)");
    mg = fitlme(V(isfinite(V.betaGen), :), "betaGen ~ 1 + lfw + lfb + lsa + al + datasetC + (1|sessionC)");
    D1 = [D1; tagged_local(coef_local(mo, TERMS), p, "beta_obs", sum(isfinite(V.betaObs))); ...
               tagged_local(coef_local(mg, TERMS), p, "beta_gen*", sum(isfinite(V.betaGen)))]; %#ok<AGROW>
    M1.(erase(p, "-")) = struct("obs", mo, "gen", mg);
end
fprintf("\n1. Tempo decomposition, in-domain pooled (dataset fixed, session random intercept).\n");
fprintf("   lfw = within-session change per doubling of tempo; lfb = between-session.\n");
fprintf("   Prediction: beta_obs slopes with lfw (sliding target); beta_gen* flat if inversion removes it.\n");
disp(D1(ismember(D1.term, ["lfw" "lfb"]), :))

% per dataset, within-session slope only
P1 = table();
for p = PIPES
    for d = DATASETS
        V = T(T.pipeline == p & T.dataset == d, :);
        for out = ["betaObs" "betaGen"]
            W = V(isfinite(V.(out)), :);
            nMulti = sum(groupcounts(W.session) >= 2);                         % sessions that can show within-session tempo change
            if height(W) < 20 || nMulti < 5 || std(W.lfw) == 0
                fprintf("  [skipped] %s %s %s: no estimable within-session slope (%d trials; %d sessions with >= 2 trials)\n", ...
                    p, d, out, height(W), nMulti);
                P1 = [P1; [table(p, d, out, height(W), numel(unique(W.session)), 'VariableNames', ["pipeline" "dataset" "outcome" "nTrials" "nSessions"]), ...
                    table(NaN, NaN, NaN, NaN, NaN, NaN, NaN, 'VariableNames', ["b" "SE" "t" "df" "p" "CIlo" "CIhi"])]]; %#ok<AGROW>
                continue
            end
            m = fitlme(W, out + " ~ 1 + lfw + lfb + lsa + al + (1|sessionC)");
            c = coef_local(m, "lfw");
            P1 = [P1; [table(p, d, out, height(W), numel(unique(W.session)), 'VariableNames', ["pipeline" "dataset" "outcome" "nTrials" "nSessions"]), c(:, 2:end)]]; %#ok<AGROW>
        end
    end
end
fprintf("Within-session tempo slope per dataset (per doubling of f0):\n");  disp(P1)

%% 2. Selection: which trials carry a beta_gen*?
S2 = table();
for p = PIPES
    V = T(T.pipeline == p, :);
    ms = fitglme(V, "invertible ~ 1 + lfw + lfb + lsa + al + datasetC + (1|sessionC)", "Distribution", "Binomial", "FitMethod", "Laplace");
    S2 = [S2; tagged_local(coef_local(ms, ["lfw" "lfb"]), p, "invertible", height(V))]; %#ok<AGROW>
end
fprintf("\n2a. Selection: invertibility against within- and between-session tempo (binomial GLMM, logit scale):\n");  disp(S2)
Sel = table();
for p = PIPES
    for d = DATASETS
        V = T(T.pipeline == p & T.dataset == d, :);  I = V(V.invertible, :);
        Sel = [Sel; table(p, d, height(V), mean(V.invertible), median(V.f0), median(I.f0), ...
            median(V.betaObs, "omitnan"), median(I.betaObs, "omitnan"), median(I.betaGen, "omitnan"), ...
            'VariableNames', ["pipeline" "dataset" "nTrials" "shareInvertible" "f0_all" "f0_invertible" ...
            "betaObs_all" "betaObs_invertible" "betaGen"])]; %#ok<AGROW>
    end
end
Sel{:, 4:end} = round(Sel{:, 4:end}, 4);
fprintf("2b. The subset that carries beta_gen*, per dataset (medians):\n");  disp(Sel)

%% 3. Clinical contrasts, with and without tempo
C3 = table();
for p = PIPES
    % Cook: independent groups
    V = T(T.pipeline == p & ismember(T.dataset, ["Cook_CTRL" "Cook_ASD"]), :);  V.group = double(V.dataset == "Cook_ASD");
    mt = fitlme(V, "lf ~ 1 + group + (1|personC)");
    C3 = [C3; label_local(coef_local(mt, "group"), p, "Cook ASD-CTRL", "log2 f0 (tempo)", "group only")]; %#ok<AGROW>
    for out = ["betaObs" "betaGen"]
        W = V(isfinite(V.(out)), :);
        m0 = fitlme(W, out + " ~ 1 + group + (1|personC)");
        m1 = fitlme(W, out + " ~ 1 + group + lf + lsa + al + (1|personC)");
        C3 = [C3; label_local(coef_local(m0, "group"), p, "Cook ASD-CTRL", out, "group only"); ...
                  label_local(coef_local(m1, "group"), p, "Cook ASD-CTRL", out, "+ tempo, noise")]; %#ok<AGROW>
    end
    % Hickman: paired crossover, subjects present in both arms
    V = T(T.pipeline == p & ismember(T.dataset, PAIRED), :);
    both = intersect(unique(V.person(V.dataset == PAIRED(1))), unique(V.person(V.dataset == PAIRED(2))));
    V = V(ismember(V.person, both), :);  V.drug = double(V.dataset == "Hickman_HALO");
    mt = fitlme(V, "lf ~ 1 + drug + (1|personC)");
    C3 = [C3; label_local(coef_local(mt, "drug"), p, sprintf("Hickman HALO-PLAC (paired, n = %d)", numel(both)), "log2 f0 (tempo)", "drug only")]; %#ok<AGROW>
    for out = ["betaObs" "betaGen"]
        W = V(isfinite(V.(out)), :);
        m0 = fitlme(W, out + " ~ 1 + drug + (1|personC)");                    % total effect
        m1 = fitlme(W, out + " ~ 1 + drug + lf + lsa + al + (1|personC)");     % direct effect (tempo as mediator)
        C3 = [C3; label_local(coef_local(m0, "drug"), p, sprintf("Hickman HALO-PLAC (paired, n = %d)", numel(both)), out, "total (drug only)"); ...
                  label_local(coef_local(m1, "drug"), p, sprintf("Hickman HALO-PLAC (paired, n = %d)", numel(both)), out, "direct (+ tempo, noise)")]; %#ok<AGROW>
    end
end
fprintf("\n3. Clinical contrasts. Cook: tempo as a confounder of measurement. Hickman: tempo as a possible\n");
fprintf("   mediator; total minus direct = the tempo-carried part (linear, no-interaction assumptions).\n");
disp(C3)

%% 3b. Attribution (tempo vs noise) and the tempo-matched Cook contrast
A3 = table();  Mt = table();
for p = PIPES
    specs = { ...
        "Cook ASD-CTRL", T(T.pipeline == p & ismember(T.dataset, ["Cook_CTRL" "Cook_ASD"]), :), "group", "Cook_ASD"; ...
        "Hickman HALO-PLAC (paired)", T(T.pipeline == p & ismember(T.dataset, PAIRED), :), "drug", "Hickman_HALO"};
    for s = 1:size(specs, 1)
        V = specs{s, 2};  term = specs{s, 3};  V.(term) = double(V.dataset == specs{s, 4});
        if term == "drug"
            both = intersect(unique(V.person(V.dataset == PAIRED(1))), unique(V.person(V.dataset == PAIRED(2))));
            V = V(ismember(V.person, both), :);
        end
        for out = ["betaObs" "betaGen"]
            W = V(isfinite(V.(out)), :);
            f = @(rhs) coef_local(fitlme(W, out + " ~ 1 + " + term + rhs + " + (1|personC)"), term);
            c0 = f("");  ct = f(" + lf");  cn = f(" + lsa + al");  cb = f(" + lf + lsa + al");
            tot = c0.b - cb.b;                                                    % total attenuation
            shT = ((c0.b - ct.b) + (cn.b - cb.b)) / 2;                            % Shapley share: tempo
            shN = ((c0.b - cn.b) + (ct.b - cb.b)) / 2;                            % Shapley share: noise
            A3 = [A3; table(p, string(specs{s, 1}), out, height(W), c0.b, ct.b, cn.b, cb.b, c0.CIlo, c0.CIhi, cb.CIlo, cb.CIhi, ...
                shT, shN, tot, 'VariableNames', ["pipeline" "contrast" "outcome" "nTrials" "b_only" "b_tempo" ...
                "b_noise" "b_both" "CIlo_only" "CIhi_only" "CIlo_both" "CIhi_both" "tempoShare" "noiseShare" "totalAtten"])]; %#ok<AGROW>
        end
    end
    % Cook, tempo-matched (coarsened exact matching on log2 f0)
    V = T(T.pipeline == p & ismember(T.dataset, ["Cook_CTRL" "Cook_ASD"]), :);  V.group = double(V.dataset == "Cook_ASD");
    for out = ["betaGen" "betaObs"]
        W = V(isfinite(V.(out)), :);
        W.bin = floor(W.lf / MATCH_BIN);
        cnt = groupsummary(W, ["bin" "group"]);
        ok = intersect(cnt.bin(cnt.group == 0), cnt.bin(cnt.group == 1));
        W = W(ismember(W.bin, ok), :);
        if numel(unique(W.personC(W.group == 0))) < 3 || numel(unique(W.personC(W.group == 1))) < 3
            fprintf("  [skipped] %s tempo-matched %s: fewer than 3 subjects per group on common tempo support\n", p, out);
            continue
        end
        w = zeros(height(W), 1);                                           % CEM weights: pooled bin share / group bin share
        for g = 0:1
            Wg = W.group == g;  nG = sum(Wg);
            for b = ok(:)'
                k = Wg & W.bin == b;  w(k) = (sum(W.bin == b) / height(W)) / (sum(k) / nG);
            end
        end
        W.w = w;
        mU = fitlme(W, out + " ~ 1 + group + (1|personC)");
        mW = fitlme(W, out + " ~ 1 + group + (1|personC)", "Weights", W.w);
        cu = coef_local(mU, "group");  cw = coef_local(mW, "group");
        Mt = [Mt; table(p, out, height(W), numel(ok), numel(unique(W.personC(W.group == 0))), numel(unique(W.personC(W.group == 1))), ...
            median(W.f0(W.group == 0)), median(W.f0(W.group == 1)), cu.b, cu.CIlo, cu.CIhi, cw.b, cw.CIlo, cw.CIhi, cw.p, ...
            'VariableNames', ["pipeline" "outcome" "nTrials" "nBins" "nSubjCTRL" "nSubjASD" "f0medCTRL" "f0medASD" ...
            "b_commonSupport" "CIlo_cs" "CIhi_cs" "b_matched" "CIlo_m" "CIhi_m" "p_matched"])]; %#ok<AGROW>
    end
end
A3{:, 5:end} = round(A3{:, 5:end}, 4);  Mt{:, 7:end} = round(Mt{:, 7:end}, 4);
fprintf("\n3b. Attribution: group/drug coefficient alone, + tempo, + noise, + both; Shapley split of the attenuation.\n");
disp(A3)
fprintf("Cook ASD-CTRL on common tempo support, unweighted and CEM-weighted (bins %.1f log2 wide):\n", MATCH_BIN);
disp(Mt)

%% Figure (SG-IRLS)
p = PIPES(1);  V = T(T.pipeline == p, :);
fg = figure("Color", "w", "Position", [60 60 1300 420]);  tl = tiledlayout(1, 3, "TileSpacing", "compact");
nexttile; hold on;
edges = -1.5:0.5:1.5;  ctr = edges(1:end-1) + 0.25;
for out = ["betaObs" "betaGen"]
    W = V(isfinite(V.(out)), :);
    S = groupsummary(W, "session", "mean", out);  [~, loc] = ismember(W.session, S.session);
    dev = W.(out) - S.("mean_" + out)(loc);                                % within-session deviation
    b = discretize(W.lfw, edges);  mu = NaN(size(ctr));  se = mu;
    for k = 1:numel(ctr), x = dev(b == k); if numel(x) >= 10, mu(k) = mean(x); se(k) = std(x) / sqrt(numel(x)); end, end
    errorbar(ctr, mu, 1.96 * se, "-o", "LineWidth", 1.4, "DisplayName", ifelse_local(out == "betaObs", "\beta_{obs}", "\beta_{gen}*"));
end
yline(0, ":", "HandleVisibility", "off");  box on;  legend("Location", "northwest");
xlabel("within-session tempo (log_2 f_0 deviation)");  ylabel("within-session deviation of \beta");
title("Within-subject tempo, pooled in-domain (" + p + ")");
nexttile; hold on;
x = P1(P1.pipeline == p, :);  ds = unique(x.dataset, "stable");
for out = ["betaObs" "betaGen"]
    y = x(x.outcome == out, :);  [~, pos] = ismember(y.dataset, ds);
    off = ifelse_local(out == "betaObs", -0.12, 0.12);
    errorbar(pos + off, y.b, y.b - y.CIlo, y.CIhi - y.b, "o", "LineWidth", 1.3, "DisplayName", ifelse_local(out == "betaObs", "\beta_{obs}", "\beta_{gen}*"));
end
yline(0, ":", "HandleVisibility", "off");  xticks(1:numel(ds));  xticklabels(strrep(ds, "_", " "));  box on;
set(gca, "TickLabelInterpreter", "none");  legend("Location", "best");
ylabel("slope per doubling of tempo (95% CI)");  title("Within-session tempo slope, per dataset");
nexttile; hold on;
y = C3(C3.pipeline == p & C3.outcome ~= "log2 f0 (tempo)", :);
lab = regexprep(y.contrast, " \(paired, n = \d+\)", "") + ": " + ...
      strrep(strrep(y.outcome, "betaObs", "obs"), "betaGen", "gen*") + ", " + y.model;
errorbar(y.b, 1:height(y), y.b - y.CIlo, y.CIhi - y.b, "horizontal", "o", "LineWidth", 1.3);
xline(0, ":");  yticks(1:height(y));  yticklabels(lab);  set(gca, "TickLabelInterpreter", "none", "YDir", "reverse");  box on;
xlabel("group / drug effect on \beta (95% CI)");  title("Clinical contrasts, with and without tempo");
title(tl, "Tempo as measurement vs behaviour (v015 per-trial corpora)");
set(findall(fg, "Type", "axes"), "Toolbar", []);  exportgraphics(fg, OUT_PNG, "Resolution", 200);

save(OUT_MAT, "T", "D1", "P1", "S2", "Sel", "C3", "A3", "Mt", "M1", "DATASETS", "PIPES", "MATCH_BIN", "-v7.3");
fprintf("\nSaved: %s\nFigure: %s\n", OUT_MAT, OUT_PNG);

%% =========================================================================
function t = coef_local(mdl, terms)
    c = mdl.Coefficients;  k = ismember(string(c.Name), terms);
    t = table(string(c.Name(k)), round(c.Estimate(k), 4), round(c.SE(k), 4), round(c.tStat(k), 2), c.DF(k), ...
        c.pValue(k), round(c.Lower(k), 4), round(c.Upper(k), 4), 'VariableNames', ["term" "b" "SE" "t" "df" "p" "CIlo" "CIhi"]);
end

function t = tagged_local(c, p, outcome, n)
    t = [table(repmat(string(p), height(c), 1), repmat(string(outcome), height(c), 1), repmat(n, height(c), 1), ...
        'VariableNames', ["pipeline" "outcome" "nTrials"]), c];
end

function t = label_local(c, p, contrast, outcome, model)
    t = [table(string(p), string(contrast), string(outcome), string(model), 'VariableNames', ["pipeline" "contrast" "outcome" "model"]), c(:, 2:end)];
end

function v = ifelse_local(c, a, b)
    if c, v = a; else, v = b; end
end
