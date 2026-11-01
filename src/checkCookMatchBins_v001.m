%% checkCookMatchBins_v001.m
% A7 (TODO_PreSubmission_v001; Rerun D15). v004 Part 4 reports the Cook ASD-CTRL contrast "at matched tempo" with
% half-octave (0.5 log2) bins; the design advice (PATCH_DesignRecs_TempoMatching_v002) says 0.1-doubling bins.
% Re-runs analyseTempoDecomposition_v002's coarsened exact matching (L195-223) at several bin widths from the saved
% trial table, so only the bin width changes. A bin width whose common tempo support leaves fewer than 3 subjects per
% group is recorded as SKIPPED with its reason, never dropped silently.
%   Regression anchor: the 0.5 rows must reproduce the saved Mt (4 d.p.) or the script stops.
% Reads:  results/tempoDecomposition_v002.mat (T, Mt)
% Writes: results/checkCookMatchBins_v001.mat ; logs/checkCookMatchBins_v001_<yyyymmdd>.txt
% USAGE:  from the project root: checkCookMatchBins_v001
% Fraser, D.S. (2026)  v001

%% CONFIG
ROOT = fileparts(fileparts(mfilename("fullpath")));
IN_MAT = fullfile(ROOT, "results", "tempoDecomposition_v002.mat");
OUT_MAT = fullfile(ROOT, "results", "checkCookMatchBins_v001.mat");
LOGDIR = fullfile(ROOT, "logs");
BINS = [0.5 0.25 0.1];                                   % log2 f0 bin widths; 0.5 is the published one
PIPES = ["SG-IRLS" "BWFD-OLS"];  TOL = 1e-4;
if ~isfile(IN_MAT), error("cookBins:input", "%s", "FAILED PATH: " + IN_MAT); end
if ~isfolder(LOGDIR), mkdir(LOGDIR); end
diary(fullfile(LOGDIR, "checkCookMatchBins_v001_" + string(datetime("now", "Format", "yyyyMMdd")) + ".txt"));  cl = onCleanup(@() diary("off"));
fprintf("checkCookMatchBins_v001  %s\n", string(datetime("now", "Format", "yyyy-MM-dd HH:mm:ss")));
S = load(IN_MAT, "T", "Mt");

%% Matching at each bin width
M = table();
for bw = BINS
    for p = PIPES
        V = S.T(S.T.pipeline == p & ismember(S.T.dataset, ["Cook_CTRL" "Cook_ASD"]), :);  V.group = double(V.dataset == "Cook_ASD");
        for out = ["betaGen" "betaObs"]
            W = V(isfinite(V.(out)), :);  n0 = height(W);
            W.bin = floor(W.lf / bw);
            cnt = groupsummary(W, ["bin" "group"]);
            ok = intersect(cnt.bin(cnt.group == 0), cnt.bin(cnt.group == 1));
            W = W(ismember(W.bin, ok), :);
            nC = numel(unique(W.personC(W.group == 0)));  nA = numel(unique(W.personC(W.group == 1)));
            row = table(bw, p, out, n0, height(W), numel(ok), nC, nA, 'VariableNames', ["binWidth" "pipeline" "outcome" "nAll" "nKept" "nBins" "nSubjCTRL" "nSubjASD"]);
            if nC < 3 || nA < 3
                fprintf("  [SKIPPED] bin %.2f %s %s: %d CTRL and %d ASD subjects on common tempo support (need 3 each); %d of %d trials kept\n", bw, p, out, nC, nA, height(W), n0);
                row.status = "SKIPPED";  M = [M; [row, nanrow_local()]];  continue %#ok<AGROW>
            end
            w = zeros(height(W), 1);
            for g = 0:1
                Wg = W.group == g;  nG = sum(Wg);
                for b = ok(:)'
                    k = Wg & W.bin == b;  w(k) = (sum(W.bin == b) / height(W)) / (sum(k) / nG);
                end
            end
            W.w = w;
            cu = coef_local(fitlme(W, out + " ~ 1 + group + (1|personC)"), "group");
            cw = coef_local(fitlme(W, out + " ~ 1 + group + (1|personC)", "Weights", W.w), "group");
            row.status = "ok";
            M = [M; [row, table(cu.b, cu.CIlo, cu.CIhi, cw.b, cw.CIlo, cw.CIhi, cw.p, 'VariableNames', ["b_commonSupport" "CIlo_cs" "CIhi_cs" "b_matched" "CIlo_m" "CIhi_m" "p_matched"])]]; %#ok<AGROW>
        end
    end
end

%% Regression anchor: the 0.5 rows must reproduce the saved Mt
A = M(M.binWidth == 0.5, :);
for r = 1:height(S.Mt)
    j = find(A.pipeline == S.Mt.pipeline(r) & A.outcome == S.Mt.outcome(r));
    got = [A.nKept(j) A.nBins(j) A.nSubjCTRL(j) A.nSubjASD(j) round(A.b_commonSupport(j), 4) round(A.b_matched(j), 4) round(A.p_matched(j), 4)];
    want = [S.Mt.nTrials(r) S.Mt.nBins(r) S.Mt.nSubjCTRL(r) S.Mt.nSubjASD(r) S.Mt.b_commonSupport(r) S.Mt.b_matched(r) S.Mt.p_matched(r)];
    if numel(j) ~= 1 || any(abs(got - want) > TOL)
        error("cookBins:anchor", "%s", sprintf("%s %s: 0.5-bin rerun %s does not reproduce saved Mt %s", S.Mt.pipeline(r), S.Mt.outcome(r), mat2str(got, 4), mat2str(want, 4)));
    end
end
fprintf("REGRESSION ANCHOR passed: the 0.5 rows reproduce the saved Cook matched-tempo table (%d rows).\n\n", height(S.Mt));

M.fracKept = M.nKept ./ M.nAll;
num = varfun(@(x) round(x, 4), M, "InputVariables", ["b_commonSupport" "CIlo_cs" "CIhi_cs" "b_matched" "CIlo_m" "CIhi_m" "p_matched" "fracKept"]);
num.Properties.VariableNames = erase(num.Properties.VariableNames, "Fun_");
Mshow = [M(:, ["binWidth" "pipeline" "outcome" "status" "nAll" "nKept" "nBins" "nSubjCTRL" "nSubjASD"]), num];
fprintf("Cook ASD minus CTRL at matched tempo, by bin width (log2 f0):\n");  disp(Mshow)
runDate = string(datetime("now", "Format", "yyyy-MM-dd"));
save(OUT_MAT, "M", "Mshow", "BINS", "runDate");
fprintf("Saved: %s\n", OUT_MAT);

%% =========================================================================
function t = nanrow_local()
    t = table(NaN, NaN, NaN, NaN, NaN, NaN, NaN, 'VariableNames', ["b_commonSupport" "CIlo_cs" "CIhi_cs" "b_matched" "CIlo_m" "CIhi_m" "p_matched"]);
end

function t = coef_local(mdl, terms)
% analyseTempoDecomposition_v002 L270-274.
    c = mdl.Coefficients;  k = ismember(string(c.Name), terms);
    if ~any(k), error("cookBins:term", "%s", "term not in model: " + strjoin(string(terms), ", ")); end
    t = table(string(c.Name(k)), round(c.Estimate(k), 4), round(c.SE(k), 4), round(c.tStat(k), 2), c.DF(k), c.pValue(k), round(c.Lower(k), 4), round(c.Upper(k), 4), ...
        'VariableNames', ["term" "b" "SE" "t" "df" "p" "CIlo" "CIhi"]);
end
