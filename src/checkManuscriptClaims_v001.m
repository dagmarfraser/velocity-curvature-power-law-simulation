%% checkManuscriptClaims_v001.m
% Machine check that numbers printed in the manuscript (and Finding #252) equal what the saved
% result files give today. Two parts.
%  Part A (static audit): every versioned name in the manuscript (x_v###) is classified as a script,
%         function, document, figure, output file, declared RDS-only file or UNRESOLVED; scripts are checked against the .gitignore
%         allowlist by name (informational: a directory rule may also cover a script).
%  Part B (claims ledger): each ledger row rebuilds a literal string from a saved .mat (a handle
%         that formats the value exactly as the prose prints it) and requires that literal to occur
%         in the named document after the named anchor. MISMATCH: the prose and the result file disagree,
%         or the number was edited. ANCHOR LOST: the phrase that locates the claim was reworded or
%         removed (update the anchor or drop the row). Either status ends the run with an error, after
%         every claim has been reported. The manuscript is the highest PaperDraft_IdentifiabilitySkeleton_v*.md.
% Why not Quarto/R Markdown: results come from long simulations saved as .mat; the weak link is
% number-to-source traceability, which this ledger gives without a second language.
% How to add a claim: append one line to claims_local() with an id, the document key, an anchor
% string that precedes the claim, and a handle M -> string. M.sd2 and M.unc hold the loaded mats.
% Whitespace in the prose is normalised (line breaks in markdown), the minus sign is U+2212 where the
% manuscript uses it (mn_local), and p-values print without a leading zero (pnl_local).
% Seeded (2026-10-11) with Finding #252's claims and the older "VERIFIED against" tags (T1-T3, L1-L9, F1-F5, V1-V10,
% E1-E7, S1-S5, A1-A2, G1-G5, K1-K8, W1, N1-N6, D1-D2, C1-C5, R1, X1-X5, P1-P5, Q1-Q2, D3-D4). The 15.9 GB vgfLMM RDS is read through
% Finding #217's small cached result (src/checkD10VGFCoeffCI_v002_results.mat); anything a script prints but does not save is not checkable.
% Fail loud: a missing document, mat, anchor or a stale mat (older than its script) is an error or a
% printed warning; any unmatched claim ends the run with an error after the report is saved.
% Run from src/. Reads docs/ and results/; writes results/checkManuscriptClaims_v001.mat and
% results/checkManuscriptClaims_v001.csv

%% CONFIG
ROOT  = fileparts(fileparts(mfilename("fullpath")));
msList = dir(fullfile(ROOT, "docs", "PaperDraft_IdentifiabilitySkeleton_v*.md"));
if isempty(msList), error("claims:doc", "%s", "No docs/PaperDraft_IdentifiabilitySkeleton_v*.md found"); end
msName = string(sort({msList.name}));  msName = msName(end);          % highest version number wins
DOCS  = struct("ms", msName, "fr", "FINDINGS_REFERENCE.md");
fprintf("Manuscript checked: %s (latest version in docs/)\n", DOCS.ms);
MATS  = struct("sd2", "checkSelectionDefect_v002", "unc", "checkDatasetMedianUncertainty_v001", "pb", "checkPart2VerdictBasis_v002", "pm", "checkPeakAtTopMapped_v001", "em", "checkEmptyMaps_v001", ...
    "sq", "checkSEMBiasQuadrants_v001", "ab", "checkAbstractBetaGenRange_v001", "tv", "checkEmpiricalTempoVsGridVGF_v002", ...
    "ct", "checkCookTempoInteraction_v001", ...
    "el", "checkEllipseSizeVsGrid_v001", "sl", "checkSEMLookupSizeMatched_v001", "td", "checkTempoDecompositionDomain_v001", ...
    "vc", "checkVGFSEMCentroid_v005", "gt", "checkGridTempoRange_v001", "sx", "crossReferenceSimexMonotonicity_v002", "tm", "tempoDecomposition_v002", "d10", "checkD10VGFCoeffCI_v002");
MATF  = struct("d10", "src/checkD10VGFCoeffCI_v002_results.mat");   % result file outside results/ (small cache of the RDS-only vgfLMM, Finding #217)
SRC   = struct("tm", "analyseTempoDecomposition_v002");   % script name where it differs from the result file name
VARS  = struct("sx", {{"res"}}, "tm", {{"Mt"}}, "d10", {{"maxHW"}});   % load only these variables where the mat is large (the SIMEX join table is not needed)
OUT   = fullfile(ROOT, "results", "checkManuscriptClaims_v001");
RDS_ONLY = "vgfLMM_v001";   % vgfLMM_v001_latest.mat, 15.9 GB, RDS only (Finding #217); declared so it is visible, not silently skipped

%% Load documents (whitespace normalised) and results
dk = string(fieldnames(DOCS))';  txt = struct();
for k = dk
    f = fullfile(ROOT, "docs", DOCS.(k));
    if ~isfile(f), error("claims:doc", "%s", "Missing document: " + f); end
    txt.(k) = regexprep(string(fileread(f)), "\s+", " ");
end
mk = string(fieldnames(MATS))';  M = struct();
for k = mk
    fm = fullfile(ROOT, "results", MATS.(k) + ".mat");  sn = MATS.(k);  if isfield(SRC, k), sn = SRC.(k); end
    if isfield(MATF, k), fm = fullfile(ROOT, MATF.(k)); end
    fs = fullfile(ROOT, "src", sn + ".m");
    if ~isfile(fm), error("claims:mat", "%s", "Missing result file: " + fm); end
    if ~isfile(fs), error("claims:script", "%s", "Missing script: " + fs); end
    a = dir(fm);  b = dir(fs);
    if a.datenum < b.datenum
        warning("claims:stale", "%s", MATS.(k) + ".mat is OLDER than its script: re-run before trusting this check.");
    end
    if isfield(VARS, k), M.(k) = load(fm, VARS.(k){:}); else, M.(k) = load(fm); end
end

%% Part B: claims ledger
C = claims_local();  nC = size(C, 1);
status = strings(nC, 1);  lit = strings(nC, 1);
for i = 1:nC
    id = C{i, 1};  dkey = C{i, 2};  anchor = C{i, 3};
    t = txt.(dkey);  p = strfind(t, anchor);
    lit(i) = regexprep(string(C{i, 4}(M)), "\s+", " ");
    if isempty(p)
        status(i) = "ANCHOR LOST";                                  % prose around the claim was reworded or removed
    else
        t = extractAfter(t, p(1) - 1);
        if contains(t, lit(i)), status(i) = "OK"; else, status(i) = "MISMATCH"; end
    end
end
L = table(string(C(:, 1)), string(C(:, 2)), lit, status, 'VariableNames', ["id" "doc" "literal" "status"]);

%% Part A: static audit of versioned names in the manuscript
ms = string(fileread(fullfile(ROOT, "docs", DOCS.ms)));
nm = unique(string(regexp(ms, '(?<![A-Za-z0-9_])[A-Za-z][A-Za-z0-9_]*_v\d{3}(?![0-9])', 'match')));  nm = nm(:);
gi = splitlines(string(fileread(fullfile(ROOT, ".gitignore"))));
kind = strings(numel(nm), 1);  allow = false(numel(nm), 1);
for i = 1:numel(nm)
    n = nm(i);
    if isfile(fullfile(ROOT, "src", n + ".m")),                kind(i) = "script";   allow(i) = any(gi == "!src/" + n + ".m");
    elseif isfile(fullfile(ROOT, "src", "functions", n + ".m")), kind(i) = "function";
    elseif isfile(fullfile(ROOT, n + ".m")),                   kind(i) = "root script";
    elseif ~isempty(dir(fullfile(ROOT, "docs", n + ".*"))) || ~isempty(dir(fullfile(ROOT, "docs", "obsolete", n + ".*"))) || ~isempty(dir(fullfile(ROOT, n + ".*"))), kind(i) = "document";
    elseif ~isempty(dir(fullfile(ROOT, "figures", "*" + n + "*"))),  kind(i) = "figure";
    elseif ~isempty(dir(fullfile(ROOT, "src", n + "*"))) || ~isempty(dir(fullfile(ROOT, "results", n + "*"))), kind(i) = "output";
    elseif ismember(n, RDS_ONLY),                              kind(i) = "RDS-only (declared)";
    else,                                                      kind(i) = "UNRESOLVED";
    end
end
A = table(nm, kind, allow, 'VariableNames', ["name" "kind" "allowlistedByName"]);

%% Report and save
fprintf("Part B: %d claims, %d OK, %d MISMATCH, %d ANCHOR LOST\n", nC, sum(status == "OK"), sum(status == "MISMATCH"), sum(status == "ANCHOR LOST"));
disp(L(:, ["id" "status" "literal"]))
[uk, ~, jk] = unique(A.kind);  nk = accumarray(jk, 1);
fprintf("Part A: %d versioned names in the manuscript: %s\n", height(A), strjoin(compose("%d %s", nk, uk)', ", "));
fprintf("UNRESOLVED (no script, function, document or output found):\n");
disp(A.name(A.kind == "UNRESOLVED"))
fprintf("Scripts in src/ not allowlisted by name (informational; a directory rule may cover them):\n");
disp(A.name(A.kind == "script" & ~A.allowlistedByName))
writetable(L, OUT + ".csv");
save(OUT + ".mat", "L", "A", "MATS", "DOCS", "-v7.3");
fprintf("Saved: %s.mat and .csv\n", OUT);
if any(status ~= "OK")
    error("claims:mismatch", "%s", sprintf("%d claim(s) do not match the saved results: %s", sum(status ~= "OK"), strjoin(L.id(status ~= "OK"), ", ")));
end

%% =========================================================================
function C = claims_local()
% {id, doc key, anchor that precedes the claim, handle M -> literal string as printed in the prose}
mn = @mn_local;  pn = @pnl_local;
dsN = @(r, s) r.estimator == "constellation" & r.contrast == s;
bootRow = @(M, c) M.unc.rowsD(dsN(M.unc.rowsD, c), :);
bootLit = @(lab, c) @(M) sprintf("%s %+.3f [%+.3f, %+.3f]", lab, bootRow(M, c).diff, bootRow(M, c).ciLo, bootRow(M, c).ciHi);
sub2 = char(8322);
A252 = "Finding #252";
C = {
 "T1 Cook tempo all trials",  "ms", A252, @(M) sprintf("uses all %d trials", M.sd2.rowsD.nTrials(1));
 "T2 Cook tempo gated all",   "ms", A252, @(M) sprintf("%.2f [%s, %.2f] (p = %s, %d trials, %d subjects)", M.sd2.rowsD.b_log2(2), mn(M.sd2.rowsD.ciLo(2), 2), M.sd2.rowsD.ciHi(2), pn(M.sd2.rowsD.p(2), 2), M.sd2.rowsD.nTrials(2), M.sd2.rowsD.nSubjects(2));
 "T3 Cook tempo gated SG-IRLS", "ms", A252, @(M) sprintf("%.2f [%.4f, %.2f] (p = %s, %d trials, %d subjects)", M.sd2.rowsD.b_log2(3), M.sd2.rowsD.ciLo(3), M.sd2.rowsD.ciHi(3), pn(M.sd2.rowsD.p(3), 4), M.sd2.rowsD.nTrials(3), M.sd2.rowsD.nSubjects(3));
 "L1 thinning max shift",     "ms", A252, @(M) sprintf("largest shift %.4f", max(abs(M.sd2.rowsB.shiftFromFull)));
 "L2 Fraser - Cook CTRL",     "ms", A252, bootLit("Cook CTRL", "Fraser - Cook_CTRL");
 "L3 Fraser - Cook ASD",      "ms", A252, bootLit("Cook ASD", "Fraser - Cook_ASD");
 "L4 Fraser - Hickman PLAC",  "ms", A252, bootLit("Hickman PLAC", "Fraser - Hickman_PLAC");
 "L5 Fraser - Hickman HALO",  "ms", A252, bootLit("Hickman HALO", "Fraser - Hickman_HALO");
 "L6 gate-f0 correlation", "ms", A252, @(M) sprintf("log%s f0 %+.2f to %+.2f", sub2, min(rho_local(M, "log2 f0")), max(rho_local(M, "log2 f0")));
 "L7 gate shift at most",     "ms", A252, @(M) sprintf("at most %.2f", max(abs(shift_local(M))));
 "L8 P(PASS) Cook IRLS",      "ms", A252, @(M) sprintf("P(PASS) %.2f to %.2f", min(pass_local(M)), max(pass_local(M)));
 "L9 asymmetric pairs",       "ms", A252, @(M) sprintf("%d of %d single-pipeline pair differences", sum(M.sd2.rowsC.exceedsMDC277 & M.sd2.rowsC.gate ~= "all-pipeline"), sum(M.sd2.rowsC.gate ~= "all-pipeline"));
 "F1 coverage Fraser",        "fr", "## Finding #252", @(M) covLit_local(M, 1, "Fraser");
 "F2 coverage Cook CTRL",     "fr", "## Finding #252", @(M) covLit_local(M, 2, "Cook CTRL");
 "F3 coverage Cook ASD",      "fr", "## Finding #252", @(M) covLit_local(M, 3, "Cook ASD");
 "F4 coverage Hickman PLAC",  "fr", "## Finding #252", @(M) covLit_local(M, 4, "Hickman PLAC");
 "F5 coverage Hickman HALO",  "fr", "## Finding #252", @(M) covLit_local(M, 5, "Hickman HALO");
 "V1 verdict counts N200",     "ms", "N_REPS = 20 (8/2", @(M) sprintf("(%s at N_REPS = 200)", cnt_local(M.pb.V, "v015", false));
 "V2 in-domain counts",       "ms", "Within the validity domain (30", @(M) sprintf("%s (%s at N_REPS = 200", cnt_local(M.pb.V, "v012", true), cnt_local(M.pb.V, "v015", true));
 "V3 in-zone share",          "ms", "With pure monotonicity", @(M) sprintf("%s%% of them do in the five in-domain datasets (%.1f%% in Dhieb)", rngp_local(sz_local(M, true).pctInvInZone), sz_local(M, false).pctInvInZone);
 "V4 clearance",              "ms", "With pure monotonicity", @(M) sprintf("median %.2f-%.2f, at least the MDC in %.0f-%.0f%% of trials (Dhieb %.2f and %.1f%%)", min(sz_local(M, true).clearMed), max(sz_local(M, true).clearMed), min(sz_local(M, true).pctClearMDC), max(sz_local(M, true).pctClearMDC), sz_local(M, false).clearMed, sz_local(M, false).pctClearMDC);
 "V5 rising-branch share",    "ms", "(a) at each trial's own tempo", @(M) sprintf("%s%% of the trials the local check calls invertible", rngp_local(sz_local(M, true).pctInvInZone));
 "V6 clearance (summary)",    "ms", "(a) at each trial's own tempo", @(M) sprintf("median %.2f-%.2f (at least one MDC in %.0f-%.0f%%)", min(sz_local(M, true).clearMed), max(sz_local(M, true).clearMed), min(sz_local(M, true).pctClearMDC), max(sz_local(M, true).pctClearMDC));
 "V7 peak at top, maps only", "ms", "(a) at each trial's own tempo", @(M) sprintf("%.0f-%.0f%% of in-domain trials with a map", min(pm_local(M)), max(pm_local(M)));
 "V9 peak at top (Part 2)",   "ms", "top of the sweep (β_gen 0.75)", @(M) sprintf("in %.0f-%.0f%% of in-domain trials", min(pm_local(M)), max(pm_local(M)));
 "V10 peak at top (Part 1B)", "ms", "tempo, where", @(M) sprintf("%.0f-%.0f%% of in-domain per-trial maps", min(pm_local(M)), max(pm_local(M)));
 "V8 peak at top (#239)",     "ms", "#239's", @(M) sprintf("#239's %.0f-%.0f%% kept the 84", min(sz_local(M, true).pctPeakAtTop), max(sz_local(M, true).pctPeakAtTop));
 "E1 32 of 42 fail",           "ms", "Under that local standard", @(M) sprintf("%d of %d dataset-pipeline combinations (%.1f%%; %d of the %d in the validity domain)", sum(M.pb.V.verdict(M.pb.V.basis == "v015") == "FAIL"), sum(M.pb.V.basis == "v015"), 100 * mean(M.pb.V.verdict(M.pb.V.basis == "v015") == "FAIL"), sum(M.pb.V.verdict(M.pb.V.basis == "v015" & M.pb.V.inDomain) == "FAIL"), sum(M.pb.V.basis == "v015" & M.pb.V.inDomain));
 "E2 84 trials without a map", "ms", "two-orbit window the template", @(M) sprintf("There are %d trials too short to model in the domain (%.1f%%: Fraser %d, Hickman HALO %d, PLAC %d)", sum(nm_local(M, [])), 100 * sum(nm_local(M, [])) / sum(M.em.S.nTrials(M.em.S.inDomain)), nm_local(M, "Fraser"), nm_local(M, "Hickman_HALO"), nm_local(M, "Hickman_PLAC"));
 "E3 HALO IRLS set aside",    "ms", "two-orbit window the template", @(M) sprintf("(%.1f%% and %.1f%%; %s, and %s in the domain)", 100 * em_irls_local(M).cov(1), 100 * em_irls_local(M).cov(2), cnt1_local(M.em.V, "v015", 2, false), cnt1_local(M.em.V, "v015", 2, true));
 "E4 Fraser OLS at N20",      "ms", "two-orbit window the template", @(M) sprintf("also lift Fraser's two OLS cells to %s (%s)", em_fraserols_local(M), cnt1_local(M.em.V, "v012", 2, false));
 "E5 Cook low tempo",         "ms", "two-orbit window the template", @(M) sprintf("%.0f-%.0f%% of them below 0.35 Hz, median %.2f-%.2f Hz against %.2f-%.2f Hz", min(ck_local(M).pctFailSlow), max(ck_local(M).pctFailSlow), min(ck_local(M).f0MedFail), max(ck_local(M).f0MedFail), min(ck_local(M).f0MedRise), max(ck_local(M).f0MedRise));
 "E6 fold branch share",      "ms", "two-orbit window the template", @(M) sprintf("%.0f-%.0f%% of failures in Fraser's and Hickman's LMLS cells, %.0f-%.0f%% in every other FAIL cell", min(fold_local(M, true)), max(fold_local(M, true)), min(fold_local(M, false)), max(fold_local(M, false)));
 "E7 HALO IRLS tempo",        "ms", "two-orbit window the template", @(M) sprintf("failing trials median %.2f Hz against %.2f-%.2f Hz", max(hi_local(M).f0MedFail), min(hi_local(M).f0MedRise), max(hi_local(M).f0MedRise));
 "S1 precise and accurate",    "ms", "The registered dual assessment", @(M) sprintf("%.1f%% of evaluable coordinates are precise and accurate, %.1f%% precise but biased", M.sq.R.pctPreciseAccurate(1), M.sq.R.pctPreciseBiased(1));
 "S2 imprecise shares",       "ms", "The registered dual assessment", @(M) sprintf("%.1f%% biased and imprecise, and %.1f%% accurate but imprecise", M.sq.R.pctImpreciseBiased(1), M.sq.R.pctImpreciseAccurate(1));
 "S3 biased given precise",   "ms", "The registered dual assessment", @(M) sprintf("%.1f%% are biased by more than the MDC, in every pipeline (%.1f-%.1f%%)", M.sq.R.pctBiasedGivenPrecise(1), min(M.sq.P.pctBiasedGivenPrecise), max(M.sq.P.pctBiasedGivenPrecise));
 "S4 registered noise range", "ms", "The registered dual assessment", @(M) sprintf("(%s %s 3) the share is %.1f%%, with %.1f%% of coordinates precise and accurate", char(945), char(8804), M.sq.R.pctBiasedGivenPrecise(2), M.sq.R.pctPreciseAccurate(2));
 "S5 evaluable coordinates",  "ms", "The registered dual assessment", @(M) sprintf("(%s evaluable coordinates)", thou_local(M.sq.R.n(1)));
 "A1 abstract range N20",     "ms", "candidate cut if space requires", @(M) sprintf("gives %.2f-%.2f", min(M.ab.T.medN20), max(M.ab.T.medN20));
 "A2 abstract range N200",    "ms", "candidate cut if space requires", @(M) sprintf("The %.2f-%.2f range", min(M.ab.T.medN200), max(M.ab.T.medN200));
 "G1 below the grid floor",   "ms", "at each trial's own", @(M) sprintf("%.0f%% of Fraser, %.0f%% of Cook CTRL and %.0f%% of Cook ASD trials are slower", tv_local(M, "Fraser", "pctBelowPrim"), tv_local(M, "Cook_CTRL", "pctBelowPrim"), tv_local(M, "Cook_ASD", "pctBelowPrim"));
 "G2 Hickman both ends",      "ms", "at each trial's own", @(M) sprintf("(PLAC %.0f%% below and %.0f%% above; HALO %.0f%% and %.0f%%)", tv_local(M, "Hickman_PLAC", "pctBelowPrim"), tv_local(M, "Hickman_PLAC", "pctAbovePrim"), tv_local(M, "Hickman_HALO", "pctBelowPrim"), tv_local(M, "Hickman_HALO", "pctAbovePrim"));
 "G3 all-trial shares",       "ms", "at each trial's own", @(M) sprintf("the shares below are %.0f%%, %.0f%% and %.0f%%", tv_local(M, "Fraser", "pctBelowRef"), tv_local(M, "Cook_CTRL", "pctBelowRef"), tv_local(M, "Cook_ASD", "pctBelowRef"));
 "G4 grid VGF range",         "ms", "L_0.28 = 245.8", @(M) sprintf("(%.1f-%.1f) is 0.49-1.79 Hz", M.tv.constants.VGF_LO, M.tv.constants.VGF_HI);
 "G5 L_1/3",                  "ms", "L_0.28 = 245.8", @(M) sprintf("L_1/3 = %.1f", M.tv.constants.K13);
 "K1 Cook tempo contrast",     "ms", "The Cook group contrast and tempo", @(M) sprintf("+%.2f log%s units, about %.1f times the tempo; SE %.2f, t(%d) = %.2f, p = %s, 95%% CI [%s, %.2f]", tc_local(M, "Cook").b, char(8322), 2 ^ tc_local(M, "Cook").b, tc_local(M, "Cook").SE, tc_local(M, "Cook").df, tc_local(M, "Cook").t, pn(tc_local(M, "Cook").p, 2), mn(tc_local(M, "Cook").CIlo, 2), tc_local(M, "Cook").CIhi);
 "K2 raw exponent difference", "ms", "The Cook group contrast and tempo", @(M) sprintf("+%.3f SG-IRLS, +%.3f BWFD-OLS", rh_local(ct_local(M, "SG-IRLS", "betaObs").b_only, 3), rh_local(ct_local(M, "BWFD-OLS", "betaObs").b_only, 3));
 "K3 tempo alone",            "ms", "The Cook group contrast and tempo", @(M) sprintf("reduces the difference to +%.3f (+%.3f)", ct_local(M, "SG-IRLS", "betaObs").b_tempo, ct_local(M, "BWFD-OLS", "betaObs").b_tempo);
 "K4 Shapley split",          "ms", "The Cook group contrast and tempo", @(M) sprintf("attributes %.0f%% of the reduction to tempo and %.0f%% to noise, in both pipelines", shap_local(M, false, "tempo"), shap_local(M, false, "noise"));
 "K5 with tempo x group",     "ms", "The Cook group contrast and tempo", @(M) sprintf("tempo still takes %.0f-%.0f%% of the reduction in both pipelines", min(shap_local(M, true, "tempo", true)), max(shap_local(M, true, "tempo", true)));
 "K6 interaction terms",      "ms", "The Cook group contrast and tempo", @(M) sprintf("b = %s, p = %s, SG-IRLS; %s, p = %s, BWFD-OLS", mn(it_local(M, "SG-IRLS", "betaObs").b, 3), pn(it_local(M, "SG-IRLS", "betaObs").p, 2), mn(it_local(M, "BWFD-OLS", "betaObs").b, 3), pn(it_local(M, "BWFD-OLS", "betaObs").p, 2));
 "K7 noise colour x tempo",   "ms", "The Cook group contrast and tempo", @(M) sprintf("(b = %.2f per doubling, p = %s)", it_local(M, "SG-IRLS", "alpha (noise colour)").b, pn(it_local(M, "SG-IRLS", "alpha (noise colour)").p, 2));
 "K8 haloperidol tempo",      "ms", "Haloperidol does not change tempo", @(M) sprintf("(+%.2f log%s units, SE %.3f, t(%d) = %.2f, p = %s, 95%% CI [%s, %.2f]", tc_local(M, "Hickman").b, char(8322), tc_local(M, "Hickman").SE, tc_local(M, "Hickman").df, tc_local(M, "Hickman").t, pn(tc_local(M, "Hickman").p, 2), mn(tc_local(M, "Hickman").CIlo, 2), tc_local(M, "Hickman").CIhi);
 "W1 Fraser coverage v012 to v015", "ms", "Fraser's own coverage status, updated", @(M) fraserShift_local(M);
 "N1 grid ellipse",           "ms", "in real mm to the grid's", @(M) sprintf("semi-axes %.1f %s %.1f harness-mm", M.el.grid.aFitMM, char(215), M.el.grid.bFitMM);
 "N2 empirical a/b",          "ms", "in real mm to the grid's", @(M) sprintf("(a/b %.1f-%.1f in domain, against %.2f)", min(el_local(M, true).ab), max(el_local(M, true).ab), M.el.grid.aFitPx / M.el.grid.bFitPx);
 "N3 Cook and Hickman size",  "ms", "in real mm to the grid's", @(M) sprintf("semi-axes are about %.0f %s %.0f mm", mean(el_local(M, "ch").aMM), char(215), mean(el_local(M, "ch").bMM));
 "N4 Fraser size",            "ms", "in real mm to the grid's", @(M) sprintf("Fraser's are %.1f %s %.1f mm, so it understates Fraser's by a factor of %.1f (%s/a %.3f, against %.3f at the matched grid %s)", el_local(M, "Fraser").aMM, char(215), el_local(M, "Fraser").bMM, el_local(M, "Fraser").gridAoverEmpA, char(963), el_local(M, "Fraser").sigmaOverA, el_local(M, "Fraser").gridSigmaOverA, char(963));
 "N5 adequacy after scaling",  "ms", "in real mm to the grid's", @(M) sprintf("%d of %d cells stay adequate (%d of %d in domain)", sum(M.sl.out.adequateMatched), height(M.sl.out), sum(M.sl.out.adequateMatched(~ismember(M.sl.out.dataset, ["Dhieb" "Zarandi"]))), sum(~ismember(M.sl.out.dataset, ["Dhieb" "Zarandi"])));
 "N6 Fraser SEM rise",         "ms", "in real mm to the grid's", @(M) sprintf("Fraser's SEM rises from %.4f-%.4f to %.4f-%.4f", min(M.sl.out.semCanonical(M.sl.out.dataset == "Fraser")), max(M.sl.out.semCanonical(M.sl.out.dataset == "Fraser")), min(M.sl.out.semMatched(M.sl.out.dataset == "Fraser")), max(M.sl.out.semMatched(M.sl.out.dataset == "Fraser")));
 "D1 within-session slopes",  "ms", "within-session slope is", @(M) sprintf("+%.4f for SG-IRLS %s_obs (t(%d) = %.1f) against +%.4f for %s_gen* (t(%d) = %.1f), a %.1f%% reduction; for BWFD-OLS, +%.4f against +%.4f, a %.1f%% reduction", td_local(M, 2, "SG-IRLS", "beta_obs").b, char(946), td_local(M, 2, "SG-IRLS", "beta_obs").df, td_local(M, 2, "SG-IRLS", "beta_obs").t, td_local(M, 2, "SG-IRLS", "beta_gen*").b, char(946), td_local(M, 2, "SG-IRLS", "beta_gen*").df, td_local(M, 2, "SG-IRLS", "beta_gen*").t, red_local(M, 2, "SG-IRLS"), td_local(M, 2, "BWFD-OLS", "beta_obs").b, td_local(M, 2, "BWFD-OLS", "beta_gen*").b, red_local(M, 2, "BWFD-OLS"));
 "D2 with Dhieb in",          "ms", "within-session slope is", @(M) sprintf("with Dhieb in, as first reported, %.1f%% and %.1f%%", red_local(M, 1, "SG-IRLS"), red_local(M, 1, "BWFD-OLS"));
 "C1 VGF bias range",         "ms", "Bias is large and uniformly signed", @(M) sprintf("overestimates true tempo by %.1f-%.1f%%", min(vc_local(M, 2, "biasPct")), max(vc_local(M, 2, "biasPct")));
 "C2 Zarandi",                "ms", "Bias is large and uniformly signed", @(M) sprintf("(bias %.1f%%, SEM %.1f%%)", vc_local(M, 2, "biasPct", "Zarandi"), vc_local(M, 2, "semPct", "Zarandi"));
 "C3 VGF SEM range",          "ms", "Bias is large and uniformly signed", @(M) sprintf("(%.1f-%.1f%% of each dataset's own tempo)", min(vc_local(M, 2, "semPct")), max(vc_local(M, 2, "semPct")));
 "C4 mean vs median f0",      "ms", "Bias is large and uniformly signed", @(M) sprintf("Cook ASD's bias is %.1f%% rather than %.1f%%", vc_local(M, 3, "biasPct", "Cook ASD"), vc_local(M, 2, "biasPct", "Cook ASD"));
 "C5 retired approximation",  "ms", "The previously published", @(M) sprintf("%.1f-%.1f%%", min(vc_local(M, 0, "biasPct")), max(vc_local(M, 0, "biasPct")));
 "R1 grid tempo range",       "ms", "Tempo rises with", @(M) sprintf("%.2f-%.2f Hz at %s_gen = 1/3 and %.2f-%.2f Hz at 2/3", gt_local(M, 1/3, "f0Min_realised"), gt_local(M, 1/3, "f0Max_realised"), char(946), gt_local(M, 2/3, "f0Min_realised"), gt_local(M, 2/3, "f0Max_realised"));
 "X1 SIMEX cell counts",      "ms", "Joined against the", @(M) sprintf("(%s SIMEX cells; %s joined to a placed trial, %d rows from trials with no %s_gen* excluded)", thou_local(sx_local(M, "V2 PRIMARY").nSimex), thou_local(sx_local(M, "V2 PRIMARY").nInv + sx_local(M, "V2 PRIMARY").nNon), sx_local(M, "V2 PRIMARY").nUnplaced, char(946));
 "X2 SIMEX widths and tests", "ms", "Joined against the", @(M) sprintf("(mean %.3f vs %.3f), a moderate effect in the predicted direction (Cohen's d = %.2f, p < %s Welch, p < %s Wilcoxon rank-sum)", sx_local(M, "V2 PRIMARY").meanNon, sx_local(M, "V2 PRIMARY").meanInv, sx_local(M, "V2 PRIMARY").d, tenpow_local(sx_local(M, "V2 PRIMARY").pWelch), tenpow_local(sx_local(M, "V2 PRIMARY").pWilcoxon));
 "X3 SIMEX overlap",          "ms", "Joined against the", @(M) sprintf("(%.0f%% of invertible cells exceed the non-invertible group's median width; %.0f%% of non-invertible cells exceed the invertible group's median)", sx_local(M, "V2 PRIMARY").pctInvOverNonMed, sx_local(M, "V2 PRIMARY").pctNonOverInvMed);
 "X4 SIMEX retired d",        "ms", "Joined against the", @(M) sprintf("d = %.2f (mean width %.3f vs %.3f; overlap %.0f%% and %.0f%% on re-run", sx_local(M, "V0 RETIRED").d, sx_local(M, "V0 RETIRED").meanNon, sx_local(M, "V0 RETIRED").meanInv, sx_local(M, "V0 RETIRED").pctInvOverNonMed, sx_local(M, "V0 RETIRED").pctNonOverInvMed);
 "X5 SIMEX constant d",       "ms", "Joined against the", @(M) sprintf("the exact L_1/3 alone gives d = %.2f", sx_local(M, "V1 CONSTANT").d);
 "P1 off-map split",          "ms", "What fails is mostly the map", @(M) sprintf("%.0f%% of failing trials observe a value their map cannot produce (%.1f%% above its top, %.1f%% below its floor), and failing maps are flatter (median chord slope %.2f, against %.2f", M.pb.pooled.pctOffMap, M.pb.pooled.pctOffMapAbove, M.pb.pooled.pctOffMapBelow, M.pb.pooled.slopeFail, M.pb.pooled.slopeRise);
 "P2 no map and within",      "ms", "What fails is mostly the map", @(M) sprintf("%.0f%% come from trials too short to model at all", M.pb.pooled.pctNoMap);
 "P3 within range no run",    "ms", "What fails is mostly the map", @(M) sprintf("and %.0f%% fall within the map's range on no monotonic run", M.pb.pooled.pctWithin);
 "P4 lookup overstatement",   "ms", "Cook's and Hickman's semi-axes are", @(M) sprintf("so the lookup overstates their relative noise by about %.0f%%", 100 * (mean(el_local(M, "ch").gridSigmaOverA ./ el_local(M, "ch").sigmaOverA) - 1));
 "P5 per-dataset VGF",        "ms", "Placed this way, Cook CTRL (VGF", @(M) sprintf("Cook CTRL (VGF %.1f) and Cook ASD (%.1f) sit below the grid's floor; Fraser (%.1f), Hickman PLAC (%.1f) and HALO (%.1f) sit inside", vc_local(M, 2, "vgfOwn", "Cook CTRL"), vc_local(M, 2, "vgfOwn", "Cook ASD"), vc_local(M, 2, "vgfOwn", "Fraser"), vc_local(M, 2, "vgfOwn", "Hickman PLAC"), vc_local(M, 2, "vgfOwn", "Hickman HALO"));
 "Q1 matched-tempo contrasts", "ms", "At matched tempo (coarsened exact matching", @(M) sprintf("%s_obs differs by +%.3f (+%.3f) and %s_gen* by %s (%s)", char(946), rh_local(cem_local(M, "SG-IRLS", "betaObs").b_matched, 3), rh_local(cem_local(M, "BWFD-OLS", "betaObs").b_matched, 3), char(946), mn(rh_local(cem_local(M, "SG-IRLS", "betaGen").b_matched, 3), 3), mn(rh_local(cem_local(M, "BWFD-OLS", "betaGen").b_matched, 3), 3));
 "D3 vgfLMM half-width",      "ms", "largest coefficient confidence half-width", @(M) sprintf("is %.2f VGF units", M.d10.maxHW);
 "D4 vgfLMM half-width (#217)", "fr", "Max CI half-width = ", @(M) sprintf("**%.4f mm/s**", M.d10.maxHW);
 "Q2 matched-tempo subjects", "ms", "At matched tempo (coarsened exact matching", @(M) sprintf("%d-%d CTRL and %d-%d ASD subjects", min(M.tm.Mt.nSubjCTRL), max(M.tm.Mt.nSubjCTRL), min(M.tm.Mt.nSubjASD), max(M.tm.Mt.nSubjASD));
 };
end

function s = fraserShift_local(M)
% Fraser's coverage under N_REPS = 20 (v012) and 200 (v015): OLS and IRLS pairs. Errors if a pair rounds to two values
% where the prose prints one.
V = M.pb.V;  F = V(V.dataset == "Fraser", :);
g = @(b, f) 100 * F.coverage(F.basis == b & contains(F.pipeline, f));
one = @(x) sprintf("%.1f", unique(round(x, 1)));
if numel(unique(round(g("v012", "OLS"), 1))) ~= 1 || numel(unique(round(g("v015", "OLS"), 1))) ~= 1 || numel(unique(round(g("v012", "IRLS"), 1))) ~= 1
    error("claims:fraser", "Fraser's pipelines no longer round to one value where the prose prints one.");
end
x15 = g("v015", "IRLS");
s = sprintf("OLS from %s%% to %s%% and IRLS from %s%% to %.1f-%.1f%%", one(g("v012", "OLS")), one(g("v015", "OLS")), one(g("v012", "IRLS")), min(x15), max(x15));
end

function T = el_local(M, which)
% Rows of the ellipse-size table: in-domain datasets (true), Cook and Hickman ("ch"), or one named dataset.
T = M.el.T;
if islogical(which), T = T(T.inDomain, :); elseif which == "ch", T = T(startsWith(T.dataset, "Cook") | startsWith(T.dataset, "Hickman"), :); else, T = T(T.dataset == which, :); end
end

function r = td_local(M, k, pl, outc)
% Pooled within-session slope row (term lfw) from the tempo decomposition: k = 1 as saved (Dhieb in), 2 = validity domain.
T = M.td.D1{k};  r = T(T.pipeline == pl & T.outcome == outc & T.term == "lfw", :);
end

function x = red_local(M, k, pl)
% Percentage reduction of the within-session slope from beta_obs to beta_gen*.
x = 100 * (1 - td_local(M, k, pl, "beta_gen*").b / td_local(M, k, pl, "beta_obs").b);
end

function x = vc_local(M, variant, col, ds)
% VGF-recovery %-of-tempo: six-pipeline mean per dataset (adequate cells), then the in-domain datasets as a vector
% or one named dataset. variant 2 = exact L_1/3 at median f0, 3 = mean f0, 0 = the retired approximation.
C = M.vc.cellsAll;  C = C(C.variant == variant & C.adequate, :);
nm = ["Fraser" "Cook CTRL" "Cook ASD" "Hickman PLAC" "Hickman HALO"];  if nargin > 3, nm = ds; end
x = arrayfun(@(d) mean(C.(col)(C.dataset == d)), nm)';
end

function x = gt_local(M, beta, col)
% Realised grid tempo (all sampling rates) at one generating exponent.
R = M.gt.rows;  x = R.(col)(abs(R.beta_gen - beta) < 1e-3 & R.fs == "all fs");
end

function r = sx_local(M, variant)
% Hickman (pre-fix, cited) SIMEX-versus-monotonicity row: the corpus the manuscript quotes.
R = M.sx.res;  r = R(R.group == "Hickman (pre-fix, cited)" & R.variant == variant, :);
end

function s = tenpow_local(p)
% "10" with a superscript minus exponent: the largest k with p < 10^-k, as the prose prints it.
k = floor(-log10(p));  sup = char("⁰¹²³⁴⁵⁶⁷⁸⁹");  digs = sup(double(char(string(k))) - 47);
s = "10" + string("⁻") + string(digs);
end

function y = rh_local(x, d)
% Round half up. The saved result is stored at 4 decimals (checkCookTempoInteraction_v001 L77), so +0.0455 sits exactly
% on the third-decimal boundary and %.3f would print +0.045; the prose rounds half up (+0.046). The ledger can confirm
% consistency with the stored value, not the unrounded value.
y = sign(x) .* floor(abs(x) * 10 ^ d + 0.5 + 1e-9) / 10 ^ d;
end

function r = tc_local(M, grp)
% First SG-IRLS tempo-contrast row for "Cook" or "Hickman" (the tempo model is pipeline-independent).
T = M.ct.tc;  r = T(startsWith(T.contrast, grp) & T.pipeline == "SG-IRLS", :);
end

function r = ct_local(M, pl, outcome)
% Shapley decomposition row for one pipeline and outcome.
T = M.ct.R;  r = T(T.pipeline == pl & T.outcome == outcome, :);
end

function r = it_local(M, pl, outcome)
% Group-by-tempo interaction row for one pipeline and outcome.
T = M.ct.I;  r = T(T.pipeline == pl & T.outcome == outcome, :);
end

function v = shap_local(M, withInt, part, asVector)
% Tempo or noise percentage of the Shapley reduction, beta_obs, both pipelines (stored as amounts, not fractions).
% Without asVector the two pipelines must round to the same whole percent, else error.
T = M.ct.R(M.ct.R.outcome == "betaObs", :);
if withInt, a = T.tempoShareInt; b = T.noiseShareInt; else, a = T.tempoShare; b = T.noiseShare; end
if part == "tempo", v = 100 * a ./ (a + b); else, v = 100 * b ./ (a + b); end
if nargin < 4 || ~asVector
    if numel(unique(round(v))) ~= 1, error("claims:shapley", "The two pipelines round to different percentages: %s", mat2str(v)); end
    v = v(1);
end
end

function x = tv_local(M, ds, col)
% One field of #240's per-dataset tempo-versus-grid table (checkEmpiricalTempoVsGridVGF_v002).
R = M.tv.results;  x = R.(col)(R.dataset == ds);
end

function s = thou_local(n)
% Integer with thousands separators, as the prose prints it.
s = regexprep(sprintf("%d", round(n)), '(\d)(?=(\d{3})+$)', '$1,');   % a zero-width match would be skipped by regexprep
end

function n = nm_local(M, ds)
% Trials without a forward map (in-domain); a named dataset, or all in-domain datasets when ds is empty.
S = M.em.S;  if isempty(ds), n = S.nNoMap(S.inDomain); else, n = S.nNoMap(S.dataset == ds); end
end

function s = cnt1_local(V, basis, col, domOnly)
% "p/c/f" verdict counts for column col (1 = production, 2 = unevaluable trials set aside).
v = V.verdict(V.basis == basis & (~domOnly | V.inDomain), col);
s = sprintf("%d/%d/%d", sum(v == "PASS"), sum(v == "CONDITIONAL"), sum(v == "FAIL"));
end

function r = em_irls_local(M)
V = M.em.V;  k = V.dataset == "Hickman_HALO" & V.basis == "v015" & contains(V.pipeline, "IRLS");
r = struct("cov", V.cov(k, 2));
if any(V.verdict(k, 2) ~= "CONDITIONAL"), error("claims:em", "Hickman HALO's IRLS cells are not both CONDITIONAL once unevaluable trials are set aside."); end
end

function s = em_fraserols_local(M)
V = M.em.V;  v = V.verdict(V.dataset == "Fraser" & V.basis == "v012" & contains(V.pipeline, "OLS"), 2);
s = unique(v);  if numel(s) ~= 1, error("claims:em", "Fraser's two OLS cells do not share one verdict."); end
end

function r = cem_local(M, pl, outc)
% Cook ASD-CTRL tempo-matched (CEM-weighted) row from tempoDecomposition_v002 (Mt, stored at 4 d.p.).
T = M.tm.Mt;  r = T(T.pipeline == pl & T.outcome == outc, :);
if height(r) ~= 1, error("claims:cem", "%s", "tempo-matched row not unique: " + pl + " " + outc); end
end

function T = ck_local(M)
c = M.pb.comp;  T = c(startsWith(c.dataset, "Cook"), :);
end

function f = fold_local(M, lmls)
c = M.pb.comp;  k = contains(c.pipeline, "LMLS") & (c.dataset == "Fraser" | startsWith(c.dataset, "Hickman"));
if lmls, f = c.pctFailFold(k); else, f = c.pctFailFold(~k); end
end

function T = hi_local(M)
c = M.pb.comp;  T = c(c.dataset == "Hickman_HALO" & contains(c.pipeline, "IRLS"), :);
end

function s = cnt_local(V, basis, domOnly)
% PASS/CONDITIONAL/FAIL counts as "p/c/f" over the dataset-pipeline cells of one N_REPS basis.
v = V.verdict(V.basis == basis & (~domOnly | V.inDomain));
s = sprintf("%d/%d/%d", sum(v == "PASS"), sum(v == "CONDITIONAL"), sum(v == "FAIL"));
end

function T = sz_local(M, domOnly)
% Floor-0 rows of #239's per-dataset summary: the five in-domain datasets, or Dhieb alone.
S = M.pb.S0;  S = S(S.floor == 0, :);
if domOnly, T = S(S.inDomain, :); else, T = S(S.dataset == "Dhieb", :); end
end

function p = pm_local(M)
% Peak-at-top share over trials WITH a map (Finding #248's basis; pctPeakAtTop1), in-domain datasets.
P = M.pm.P;  p = P.pctPeakAtTop1(P.inDomain);
end

function s = rngp_local(x)
% "lo-hi" with lo to one decimal and hi as printed (100 prints as 100).
s = sprintf("%.1f-%g", min(x), round(max(x), 1));
end

function s = covLit_local(M, r, lab)
% Dataset coverage across (saved as fractions, printed as percent) the six pipelines, as printed in Finding #252.
s = lab + " " + strjoin(compose("%.2f", 100 * M.sd2.cov(r, :)), ", ");
end

function r = rho_local(M, g)
A = M.sd2.rowsA;  r = A.rho(A.gate == "all-pipeline" & A.G == g);
end

function r = shift_local(M)
A = M.sd2.rowsA;  r = A.shift(A.gate == "all-pipeline" & A.G == "betaObs");
end

function p = pass_local(M)
B = M.sd2.covB;  p = B.pPass(startsWith(B.target, "Cook") & contains(B.pipeline, "IRLS"));
end

function s = mn_local(x, d)
% Fixed-point string with a true minus sign (U+2212), as the manuscript prints it.
s = strrep(sprintf("%.*f", d, x), "-", char(8722));
end

function s = pnl_local(p, d)
% p-value without a leading zero.
s = regexprep(sprintf("%.*f", d, p), "^0", "");
end
