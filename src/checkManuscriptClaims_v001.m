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
% Updated 2026-10-23 (TODO_PreSubmission A1, A3): D1 split into D1a-D1f on the reworded tempo bullet, D5-D6 for the
% random-slope refit (results/tempoSim/fitTempoRandomSlope_v001.mat), P3 rows for Part 3's cccMed table. The person-level
% bootstrap intervals (bootstrapTempoReduction_v001.mat) are RDS-only and not checked here.
% Updated 2026-10-24: SS1-SS4, Schaal & Sternad (2001) placed on the full-range collapse, now Fig 1B-S1 (plotFilterCutoffCollapseSchaal_v001.mat);
% Z1-Z6, Fig 1B-2 now the empirical tempo range (filterCutoffZoom_v004.mat).
% IP1-IP4, D-T6 inverse-probability weighting for the gate (results/ipwGateTempoSlope_v001.mat).
% RD1-RD7, Rerun D16 stage 1 (results/obsBetaRateCutoff_v002.mat, clean subsampling; _v001.mat, anti-aliased): §2.3 fs-sensitivity, Part 1B cutoff sweep.
% SR1-SR5, §2.8 cross-machine reproducibility and SG curvature residue (results/checkSGCurvatureResidue_v001.mat; Finding #259).
% PG1-PG4, D-T7 sensitivity of the Cook null (results/powerGroupContrast_v001.mat).
% DG1-DG4, §2.5a causal structure (results/tempoDAG_v002.mat); its simulation figures (#257, #258) are not checked here (#258 is RDS-only).
% SB1-SB17, Supplement B (which part of the noise produces the tempo slope; results/tempoSim/analyseR29*, analyseR30, checkR30AlphaVsOwn).
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
    "vc", "checkVGFSEMCentroid_v005", "gt", "checkGridTempoRange_v001", "sx", "crossReferenceSimexMonotonicity_v002", "tm", "tempoDecomposition_v002", "d10", "checkD10VGFCoeffCI_v002", "rs", "fitTempoRandomSlope_v001", "p3", "checkPart3PrimaryMetric_v001_v015", "ss", "plotFilterCutoffCollapseSchaal_v001", "zm", "filterCutoffZoom_v004", "r29", "analyseR29_v001", "r29c", "analyseR29Cropped_v001", "r29s", "analyseR29SeamFree_v001", "r30", "analyseR30_v001", "r30c", "checkR30AlphaVsOwn_v001", "dag", "tempoDAG_v002", "pgc", "powerGroupContrast_v001", "ipw", "ipwGateTempoSlope_v001", "d16", "obsBetaRateCutoff_v002", "d16a", "obsBetaRateCutoff_v001", "sgr", "checkSGCurvatureResidue_v001");
MATF  = struct("d10", "src/checkD10VGFCoeffCI_v002_results.mat", "rs", "results/tempoSim/fitTempoRandomSlope_v001.mat", "r29", "results/tempoSim/analyseR29_v001.mat", "r29c", "results/tempoSim/analyseR29Cropped_v001.mat", "r29s", "results/tempoSim/analyseR29SeamFree_v001.mat", "r30", "results/tempoSim/analyseR30_v001.mat", "r30c", "results/tempoSim/checkR30AlphaVsOwn_v001.mat");   % result file outside results/ (small cache of the RDS-only vgfLMM, Finding #217)
SRC   = struct("tm", "analyseTempoDecomposition_v002", "p3", "checkPart3PrimaryMetric_v001", "zm", "checkFilterCutoffCollapse_v004");   % script name where it differs from the result file name
VARS  = struct("sx", {{"res"}}, "tm", {{"Mt"}}, "d10", {{"maxHW"}}, "zm", {{"Tq", "O", "G", "flt", "M", "F0"}}, "d16", {{"r1", "Q1", "Q3", "SET", "frames"}}, "d16a", {{"r1"}}, "sgr", {{"M", "P", "W"}});   % load only these variables where the mat is large (the SIMEX join table is not needed)
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
    elseif ~isempty(dir(fullfile(ROOT, "src", n + "*"))) || ~isempty(dir(fullfile(ROOT, "results", n + "*"))) || ~isempty(dir(fullfile(ROOT, "results", "tempoSim", n + "*"))), kind(i) = "output";
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
 "V1 verdict counts N200",     "ms", "Across all 42 cells, 36 are SEM-adequate", @(M) sprintf("%s at N_REPS = 200 (v015; %s at", regexprep(cnt_local(M.pb.V, "v015", false), "(\d+)/(\d+)/(\d+)", "$1 PASS / $2 CONDITIONAL / $3 FAIL"), cnt_local(M.pb.V, "v012", false));
 "V2 in-domain counts",       "ms", "Within the validity domain (30", @(M) sprintf("the classification is %s (%s at N_REPS = 20)", cnt_local(M.pb.V, "v015", true), cnt_local(M.pb.V, "v012", true));
 "V3 in-zone share",          "ms", "With pure monotonicity", @(M) sprintf("%s%% of them do in the five in-domain datasets (%.1f%% in Dhieb)", rngp_local(sz_local(M, true).pctInvInZone), sz_local(M, false).pctInvInZone);
 "V4 clearance",              "ms", "With pure monotonicity", @(M) sprintf("median %.2f-%.2f, at least the MDC in %.0f-%.0f%% of trials (Dhieb %.2f and %.1f%%)", min(sz_local(M, true).clearMed), max(sz_local(M, true).clearMed), min(sz_local(M, true).pctClearMDC), max(sz_local(M, true).pctClearMDC), sz_local(M, false).clearMed, sz_local(M, false).pctClearMDC);
 "V5 rising-branch share",    "ms", "(a) at each trial's own tempo", @(M) sprintf("%s%% of the trials the local check calls invertible", rngp_local(sz_local(M, true).pctInvInZone));
 "V6 clearance (summary)",    "ms", "(a) at each trial's own tempo", @(M) sprintf("median %.2f-%.2f (at least one MDC in %.0f-%.0f%%)", min(sz_local(M, true).clearMed), max(sz_local(M, true).clearMed), min(sz_local(M, true).pctClearMDC), max(sz_local(M, true).pctClearMDC));
 "V7 peak at top, maps only", "ms", "(a) at each trial's own tempo", @(M) sprintf("%.0f-%.0f%% of in-domain trials with a map", min(pm_local(M)), max(pm_local(M)));
 "V9 peak at top (Part 2)",   "ms", "top of the sweep (β_gen 0.75)", @(M) sprintf("in %.0f-%.0f%% of in-domain trials", min(pm_local(M)), max(pm_local(M)));
 "V10 peak at top (Part 1B)", "ms", "tempo, where", @(M) sprintf("%.0f-%.0f%% of in-domain per-trial maps", min(pm_local(M)), max(pm_local(M)));
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
 "A1 abstract range N20",     "ms", "Provenance:", @(M) sprintf("gave %.2f-%.2f", min(M.ab.T.medN20), max(M.ab.T.medN20));
 "A2 abstract range N200",    "ms", "Provenance:", @(M) sprintf("The %.2f-%.2f range", min(M.ab.T.medN200), max(M.ab.T.medN200));
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
 "D1a within-session slope obs",  "ms", "within-session slope in the vetted profile (SG-IRLS)", @(M) sprintf("is +%.3f for %s_obs", td_local(M, 2, "SG-IRLS", "beta_obs").b, char(946));
 "D1b within-session slope gen",  "ms", "within-session slope in the vetted profile (SG-IRLS)", @(M) sprintf("against +%.4f for %s_gen*", td_local(M, 2, "SG-IRLS", "beta_gen*").b, char(946));
 "D1c reduction vetted",          "ms", "within-session slope in the vetted profile (SG-IRLS)", @(M) sprintf("a %.0f%% reduction (CI", red_local(M, 2, "SG-IRLS"));
 "D1d legacy slopes",             "ms", "within-session slope in the vetted profile (SG-IRLS)", @(M) sprintf("it is +%.3f (CI", td_local(M, 2, "BWFD-OLS", "beta_obs").b);
 "D1e legacy gen slope",          "ms", "within-session slope in the vetted profile (SG-IRLS)", @(M) sprintf("against +%.4f (CI", td_local(M, 2, "BWFD-OLS", "beta_gen*").b);
 "D1f reduction legacy",          "ms", "within-session slope in the vetted profile (SG-IRLS)", @(M) sprintf("a %.0f%% reduction", red_local(M, 2, "BWFD-OLS"));
 "D2 with Dhieb in",              "ms", "within-session slope in the vetted profile (SG-IRLS)", @(M) sprintf("with Dhieb in, %.1f%% and %.1f%%", red_local(M, 1, "SG-IRLS"), red_local(M, 1, "BWFD-OLS"));
 "D5 random-slope vetted",        "ms", "Refitted with a random tempo slope per person", @(M) sprintf("the vetted slopes are +%.3f for %s_obs (SE %.3f, t(%.0f) = %.1f, p < .001) and +%.3f for %s_gen* (SE %.3f, t(%.0f) = %.2f, p = %s)", rs_local(M, "SG-IRLS", "betaObs").b, char(946), rs_local(M, "SG-IRLS", "betaObs").SE, rs_local(M, "SG-IRLS", "betaObs").df, rs_local(M, "SG-IRLS", "betaObs").t, rs_local(M, "SG-IRLS", "betaGen").b, char(946), rs_local(M, "SG-IRLS", "betaGen").SE, rs_local(M, "SG-IRLS", "betaGen").df, rs_local(M, "SG-IRLS", "betaGen").t, pn(rs_local(M, "SG-IRLS", "betaGen").p, 3));
 "D6 random-slope legacy",        "ms", "Refitted with a random tempo slope per person", @(M) sprintf("the legacy slopes are +%.3f for %s_obs and +%.3f for %s_gen* (p = %s, 95%% CI %s to +%.3f)", rs_local(M, "BWFD-OLS", "betaObs").b, char(946), rs_local(M, "BWFD-OLS", "betaGen").b, char(946), pn(rs_local(M, "BWFD-OLS", "betaGen").p, 2), mn(rs_local(M, "BWFD-OLS", "betaGen").CIlo, 3), rs_local(M, "BWFD-OLS", "betaGen").CIhi);
 "P3 cccMed Fraser",  "ms", "Primary result (cccMed; v015 corpus", @(M) p3row_local(M, "Fraser", "Fraser");
 "P3 cccMed Cook CTRL",  "ms", "Primary result (cccMed; v015 corpus", @(M) p3row_local(M, "Cook CTRL", "Cook CTRL");
 "P3 cccMed Cook ASD",  "ms", "Primary result (cccMed; v015 corpus", @(M) p3row_local(M, "Cook ASD", "Cook ASD");
 "P3 cccMed Hickman HALO",  "ms", "Primary result (cccMed; v015 corpus", @(M) p3row_local(M, "Hickman HALO", "Hickman HALO");
 "P3 cccMed Hickman PLAC",  "ms", "Primary result (cccMed; v015 corpus", @(M) p3row_local(M, "Hickman PLAC", "Hickman PLAC");
 "P3 cccMed Dhieb",  "ms", "Primary result (cccMed; v015 corpus", @(M) p3row_local(M, "Dhieb", "Dhieb (outside validity domain)");
 "SS1 Schaal-Sternad point",   "ms", "demonstration sits on the same curve", @(M) sprintf("read %s = %.3f, an effective gain of %.3f at f0/f_c = %.2f", char(946), M.ss.SS.betaRec, M.ss.gSS, M.ss.xSS);
 "SS2 production crossings",   "ms", "demonstration sits on the same curve", @(M) sprintf("(%.2f, %.2f and %.2f Hz at the production cutoff)", 10 * M.ss.T.x_gain09(M.ss.T.fc == 10), 10 * M.ss.T.x_gain08(M.ss.T.fc == 10), 10 * M.ss.T.x_gain05(M.ss.T.fc == 10));
 "SS3 third-harmonic |H|",     "ms", "demonstration sits on the same curve", @(M) sprintf("third harmonic (%.3f) is close to a single-pass fourth-order filter's (%.3f)", 1 / (1 + (3 * M.ss.xSS)^4), 1 / sqrt(1 + (3 * M.ss.xSS)^8));
 "SS4 our curves at their ratio", "ms", "demonstration sits on the same curve", @(M) sprintf("gain near one-third of %.2f to %.2f and read a true %s of 0 as %.2f to %.2f", min(M.ss.T.gainNear13_atSS), max(M.ss.T.gainNear13_atSS), char(946), min(M.ss.T.betaRec_at_beta0_atSS), max(M.ss.T.betaRec_at_beta0_atSS));
 "Z1 empirical tempo range",    "ms", "Over the tempi at which people draw", @(M) sprintf("(Fig 1B-2; %.2f-%.2f Hz, the lowest 1st", floor(min(M.zm.Tq.p1) * 20) / 20, ceil(max(M.zm.Tq.p99) * 4) / 4);
 "Z2 trials above 2 Hz",        "ms", "Over the tempi at which people draw", @(M) zab_local(M, "long");
 "Z3 gain at top and slow end", "ms", "Over the tempi at which people draw", @(M) zend_local(M);
 "Z4 caption onsets and masks", "ms", "Fig 1B-2. The pull to one-third over the tempi", @(M) zcap_local(M);
 "Z5 caption negative gain",    "ms", "Fig 1B-2. The pull to one-third over the tempi", @(M) sprintf("Above about %.1f Hz the 0.338 s window's gain falls below zero", zneg_local(M));
 "Z6 law bullet",               "ms", "Filtering manufactures one-third only when", @(M) zab_local(M, "short");
 "SB1 seam-free Fraser vetted",  "ms", "Noise size or noise structure?", @(M) sb_local(M, "sfFraserVetted");
 "SB2 cropped-only Fraser",      "ms", "Noise size or noise structure?", @(M) sb_local(M, "cropFraser");
 "SB3 seam-free Fraser legacy",  "ms", "Noise size or noise structure?", @(M) sb_local(M, "sfFraserLegacy");
 "SB4 PLAC three analyses",      "ms", "Noise size or noise structure?", @(M) sb_local(M, "plac");
 "SB5 wrapped share",            "ms", "Noise size or noise structure?", @(M) sb_local(M, "wrapped");
 "SB6 seam-free colour kept",    "ms", "Noise size or noise structure?", @(M) sb_local(M, "sfColour");
 "SB7 table own-noise slope",    "ms", "Which part of the structure?", @(M) sb_local(M, "rowOwn");
 "SB8 table white noise",        "ms", "Which part of the structure?", @(M) sb_local(M, "rowXn");
 "SB9 table colour only",        "ms", "Which part of the structure?", @(M) sb_local(M, "rowXc");
 "SB10 table colour and axes",   "ms", "Which part of the structure?", @(M) sb_local(M, "rowXa");
 "SB11 table axes broken",       "ms", "Which part of the structure?", @(M) sb_local(M, "rowXs");
 "SB12 without covariates",      "ms", "Which part of the structure?", @(M) sb_local(M, "v3");
 "SB13 colour check",            "ms", "Colour check.", @(M) sb_local(M, "colour");
 "SB14 white-noise slopes",      "ms", "What can and cannot be said.", @(M) sb_local(M, "white");
 "SB15 PLAC colour sufficient",  "ms", "What can and cannot be said.", @(M) sb_local(M, "placColour");
 "SB16 limits counts",           "ms", "One simulation seed per version; null truth only", @(M) sb_local(M, "limits");
 "SB16b limits correlation",     "ms", "One simulation seed per version; null truth only", @(M) sb_local(M, "limitsCorr");
 "SB17 Part 4 pointer",          "ms", "Which part of the noise produces the pipeline's slope is examined in Supplement B", @(M) sb_local(M, "pointer");
 "DG1 estimated f0 does not block", "ms", "Adjusting for the estimated f0 is not the same as adjusting for tempo.", @(M) dag_local(M, "route");
 "DG2 identification by reading",   "ms", "What each analysis estimates.", @(M) dag_local(M, "ident");
 "DG3 mediator-outcome confounding", "ms", "What each analysis estimates.", @(M) dag_local(M, "mediator");
 "DG4 gate as collider",            "ms", "The gate selects on tempo.", @(M) dag_local(M, "collider");
 "PG1 Cook power, beta_gen*",      "ms", "The design is insensitive at that scale", @(M) pgc_local(M, "gen");
 "PG2 Cook power, beta_obs",       "ms", "The design is insensitive at that scale", @(M) pgc_local(M, "obs");
 "PG3 Cook matched half-widths",   "ms", "β_gen* shows no group", @(M) pgc_local(M, "hw");
 "PG4 Discussion summary",         "ms", "went mostly with tempo, and β_gen* showed no group difference", @(M) pgc_local(M, "disc");
 "IP1 weighted tempo slopes",      "ms", "Weighting each analysed trial by its inverse probability of passing", @(M) ipw_local(M, "slopes");
 "IP2 weighted Cook contrast",     "ms", "Weighting each analysed trial by its inverse probability of passing", @(M) ipw_local(M, "cook");
 "IP3 weights",                    "ms", "Weighting each analysed trial by its inverse probability of passing", @(M) ipw_local(M, "weights");
 "IP4 Part 4 pointer",             "ms", "the remaining +0.0265 in β_gen* is not attributed", @(M) ipw_local(M, "part4");
 "C1 VGF bias range",         "ms", "Bias is large and uniformly signed", @(M) sprintf("overestimates true tempo by %.1f-%.1f%%", min(vc_local(M, 2, "biasPct")), max(vc_local(M, 2, "biasPct")));
 "C2 Zarandi",                "ms", "Bias is large and uniformly signed", @(M) sprintf("(bias %.1f%%, SEM %.1f%%)", vc_local(M, 2, "biasPct", "Zarandi"), vc_local(M, 2, "semPct", "Zarandi"));
 "C3 VGF SEM range",          "ms", "Bias is large and uniformly signed", @(M) sprintf("(%.1f-%.1f%% of each dataset's own tempo)", min(vc_local(M, 2, "semPct")), max(vc_local(M, 2, "semPct")));
 "C4 mean vs median f0",      "ms", "Bias is large and uniformly signed", @(M) sprintf("Cook ASD's bias is %.1f%% rather than %.1f%%", vc_local(M, 3, "biasPct", "Cook ASD"), vc_local(M, 2, "biasPct", "Cook ASD"));
 "R1 grid tempo range",       "ms", "Tempo rises with", @(M) sprintf("%.2f-%.2f Hz at %s_gen = 1/3 and %.2f-%.2f Hz at 2/3", gt_local(M, 1/3, "f0Min_realised"), gt_local(M, 1/3, "f0Max_realised"), char(946), gt_local(M, 2/3, "f0Min_realised"), gt_local(M, 2/3, "f0Max_realised"));
 "X1 SIMEX cell counts",      "ms", "Joined against the", @(M) sprintf("(%s SIMEX cells; %s joined to a placed trial, %d rows from trials with no %s_gen* excluded)", thou_local(sx_local(M, "V2 PRIMARY").nSimex), thou_local(sx_local(M, "V2 PRIMARY").nInv + sx_local(M, "V2 PRIMARY").nNon), sx_local(M, "V2 PRIMARY").nUnplaced, char(946));
 "X2 SIMEX widths and tests", "ms", "Joined against the", @(M) sprintf("(mean %.3f vs %.3f), a moderate effect in the predicted direction (Cohen's d = %.2f, p < %s Welch, p < %s Wilcoxon rank-sum)", sx_local(M, "V2 PRIMARY").meanNon, sx_local(M, "V2 PRIMARY").meanInv, sx_local(M, "V2 PRIMARY").d, tenpow_local(sx_local(M, "V2 PRIMARY").pWelch), tenpow_local(sx_local(M, "V2 PRIMARY").pWilcoxon));
 "X3 SIMEX overlap",          "ms", "Joined against the", @(M) sprintf("(%.0f%% of invertible cells exceed the non-invertible group's median width; %.0f%% of non-invertible cells exceed the invertible group's median)", sx_local(M, "V2 PRIMARY").pctInvOverNonMed, sx_local(M, "V2 PRIMARY").pctNonOverInvMed);
 "X4 SIMEX retired d",        "ms", "Script: `crossReferenceSimexMonotonicity_v002.m` (Finding #246)", @(M) sprintf("reproduces d = %.2f as its anchor", sx_local(M, "V0 RETIRED").d);
 "P1 off-map split",          "ms", "What fails is mostly the map", @(M) sprintf("%.0f%% of failing trials observe a value their map cannot produce (%.1f%% above its top, %.1f%% below its floor), and failing maps are flatter (median chord slope %.2f, against %.2f", M.pb.pooled.pctOffMap, M.pb.pooled.pctOffMapAbove, M.pb.pooled.pctOffMapBelow, M.pb.pooled.slopeFail, M.pb.pooled.slopeRise);
 "P2 no map and within",      "ms", "What fails is mostly the map", @(M) sprintf("%.0f%% come from trials too short to model at all", M.pb.pooled.pctNoMap);
 "P3 within range no run",    "ms", "What fails is mostly the map", @(M) sprintf("and %.0f%% fall within the map's range on no monotonic run", M.pb.pooled.pctWithin);
 "P4 lookup overstatement",   "ms", "Cook's and Hickman's semi-axes are", @(M) sprintf("so the lookup overstates their relative noise by about %.0f%%", 100 * (mean(el_local(M, "ch").gridSigmaOverA ./ el_local(M, "ch").sigmaOverA) - 1));
 "P5 per-dataset VGF",        "ms", "Placed this way, Cook CTRL (VGF", @(M) sprintf("Cook CTRL (VGF %.1f) and Cook ASD (%.1f) sit below the grid's floor; Fraser (%.1f), Hickman PLAC (%.1f) and HALO (%.1f) sit inside", vc_local(M, 2, "vgfOwn", "Cook CTRL"), vc_local(M, 2, "vgfOwn", "Cook ASD"), vc_local(M, 2, "vgfOwn", "Fraser"), vc_local(M, 2, "vgfOwn", "Hickman PLAC"), vc_local(M, 2, "vgfOwn", "Hickman HALO"));
 "Q1 matched-tempo contrasts", "ms", "At matched tempo (coarsened exact matching", @(M) sprintf("%s_obs differs by +%.3f (+%.3f) and %s_gen* by %s (%s)", char(946), rh_local(cem_local(M, "SG-IRLS", "betaObs").b_matched, 3), rh_local(cem_local(M, "BWFD-OLS", "betaObs").b_matched, 3), char(946), mn(rh_local(cem_local(M, "SG-IRLS", "betaGen").b_matched, 3), 3), mn(rh_local(cem_local(M, "BWFD-OLS", "betaGen").b_matched, 3), 3));
 "D3 vgfLMM half-width",      "ms", "largest coefficient confidence half-width", @(M) sprintf("is %.2f VGF units", M.d10.maxHW);
 "D4 vgfLMM half-width (#217)", "fr", "Max CI half-width = ", @(M) sprintf("**%.4f mm/s**", M.d10.maxHW);
 "Q2 matched-tempo subjects", "ms", "At matched tempo (coarsened exact matching", @(M) sprintf("%d-%d CTRL and %d-%d ASD subjects", min(M.tm.Mt.nSubjCTRL), max(M.tm.Mt.nSubjCTRL), min(M.tm.Mt.nSubjASD), max(M.tm.Mt.nSubjASD));
 "RD1 fs-sensitivity slope",        "ms", "Subsampling the Fraser trials cleanly from 240 Hz", @(M) d16_local(M, "slope");
 "RD2 alpha medians by rate",       "ms", "Subsampling the Fraser trials cleanly from 240 Hz", @(M) d16_local(M, "medians");
 "RD3 anti-aliased slope",          "ms", "With an anti-alias filter applied before the rate change", @(M) d16_local(M, "aa");
 "RD4 Fraser at 120 Hz vs Zarandi", "ms", "Rate accounts for little of the difference", @(M) d16_local(M, "zarandi");
 "RD5 cutoff sweep, observed",      "ms", "Changing the cutoff on the same real trials", @(M) d16_local(M, "sweep");
 "RD6 cutoff sweep, SG-IRLS",       "ms", "Changing the cutoff on the same real trials", @(M) d16_local(M, "sweepSG");
 "RD7 cutoff sweep, zero noise",    "ms", "Changing the cutoff on the same real trials", @(M) d16_local(M, "pred");
 "SR1 machine reproducibility",     "ms", "Reproducibility across machines, and curvature residue", @(M) sgr_local(M, "machine");
 "SR2 residue, machine differences", "ms", "Reproducibility across machines, and curvature residue", @(M) sgr_local(M, "residueDiff");
 "SR3 residue samples removed",     "ms", "Reproducibility across machines, and curvature residue", @(M) sgr_local(M, "drop");
 "SR4 residue by dataset",          "ms", "Reproducibility across machines, and curvature residue", @(M) sgr_local(M, "datasets");
 "SR5 narrow windows",              "ms", "Reproducibility across machines, and curvature residue", @(M) sgr_local(M, "narrow");
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

function r = rs_local(M, pl, oc)
% Random-slope (RS) row, variant V1, from fitTempoRandomSlope_v001 (T2 spec 10.16); oc "betaObs" or "betaGen".
R = M.rs.out.R;  r = R(R.pipeline == pl & R.variant == "V1" & R.outcome == oc & R.model == "RS", :);
if height(r) ~= 1, error("claims:rs", "%s", "Expected one RS row for " + pl + " " + oc); end
end

function s = p3row_local(M, ds, lab)
% One row of Part 3's primary-result table (cccMed, v015), as printed (A1, 2026-10-23).
T = M.p3.T;  r = T(T.dataset == ds, :);
if height(r) ~= 1, error("claims:p3", "%s", "Expected one Part 3 row for " + ds); end
c = @(n) regexprep(sprintf("%d", n), "(\d)(?=(\d{3})+$)", "$1,");
s = sprintf("| %s | %s | %s | %s | %s | %s | %s | %s | %.4f |", lab, c(r.Ntrials), c(r.Ninvertible), c(r.Nvalid), ...
    pnl_local(r.cccMed, 3), r.gradeRegistered, r.gradeMcBride, pnl_local(r.pearsonMed, 3), r.maeMed);
end

function s = mn_local(x, d)
% Fixed-point string with a true minus sign (U+2212), as the manuscript prints it.
s = strrep(sprintf("%.*f", d, x), "-", char(8722));
end

function s = pnl_local(p, d)
% p-value without a leading zero.
s = regexprep(sprintf("%.*f", d, p), "^0", "");
end

function s = zab_local(M, form)
% Share of trials faster than 2 Hz (Fig 1B-2, v015). Fraser and Cook ASD are quoted as one figure,
% Cook CTRL as "none": the literal breaks (MISMATCH) if either stops being true.
T = M.zm.Tq;  v = @(d) T.pctAboveBand(T.dataset == d);
if round(v("Fraser")) ~= round(v("Cook_ASD")), s = "Fraser and Cook ASD no longer share one figure"; return, end
if v("Cook_CTRL") ~= 0, s = "Cook CTRL is no longer none"; return, end
if form == "long"
    s = sprintf("%.0f%% of Fraser and Cook ASD trials and none of Cook CTRL's are faster than 2 Hz, against %.0f%% of Hickman PLAC's and %.0f%% of HALO's", v("Fraser"), v("Hickman_PLAC"), v("Hickman_HALO"));
else
    s = sprintf("%.0f%% of Fraser and Cook ASD trials and %.0f-%.0f%% of Hickman's are faster than 2 Hz", v("Fraser"), min(v("Hickman_HALO"), v("Hickman_PLAC")), max(v("Hickman_HALO"), v("Hickman_PLAC")));
end
end

function [gBW, gSG] = zprod_local(M)
% Gain rows of the two production profiles: BWFD-OLS (f_c 10 Hz) and SG-IRLS (T_w 41/240 s).
f = M.zm.flt;  G = M.zm.G;
iBW = find(f.family == "BW" & f.fc == 10);  iSG = find(f.family == "SG" & abs(f.Tw - 41/240) < 1e-12);
if numel(iBW) ~= 1 || numel(iSG) ~= 1, error("claims:zm", "%s", "Production filters not found in filterCutoffZoom_v004"); end
gBW = sortrows(G(G.filt == iBW, :), "f0");  gSG = sortrows(G(G.filt == iSG, :), "f0");
end

function s = zend_local(M)
% Gain at the top simulated tempo and the slow-end minimum, production profiles (Part 1B (b)).
[gBW, gSG] = zprod_local(M);  T = M.zm.Tq;  pb = T.pctBelowBand(ismember(T.dataset, ["Fraser" "Cook_CTRL" "Cook_ASD"]));
lowMin = min([gBW.gainOLS(gBW.f0 < 0.5); gSG.gainIRLS(gSG.f0 < 0.5)]);
s = sprintf("At the fastest tempo simulated (%.2f Hz) the production settings keep a gain of %.2f (BWFD-OLS) and %.2f (SG-IRLS). At zero noise nothing is lost at the slow end (gain at least %.3f below 0.5 Hz at the production settings), where %.0f-%.0f%% of Fraser and Cook trials sit", ...
    gBW.f0(end), gBW.gainOLS(end), gSG.gainIRLS(end), floor(lowMin * 1000) / 1000, min(pb), max(pb));
end

function s = zcap_local(M)
% Fig 1B-2 caption: top tempo, production onsets, masked IRLS tempi.
O = M.zm.O;  K = M.zm.M;
on = @(fam, set, reg) O.onsetHz(O.family == fam & O.setting == set & O.regression == reg);
if any(K.maskedOLS ~= 0), s = "OLS gains are masked too"; return, end
s = sprintf("Tempi were simulated to %.2f Hz. Shaded: 0.5-2 Hz. Bottom: per-trial f0 by dataset (median, interquartile range, 10th to 90th percentile). The dashed lines mark the production onsets (gain 0.9 at %.2f Hz, Butterworth-OLS at 10 Hz, and %.2f Hz, Savitzky-Golay-OLS at 0.171 s). Gain is left undefined at a tempo where any fit in the %s window failed (IRLS only, %d to %d of %d tempi per setting)", ...
    max(M.zm.F0), on("BW", "f_c 10 Hz", "OLS"), on("SG", "T_w 0.171 s", "OLS"), char(946) + "_gen", min(K.maskedIRLS), max(K.maskedIRLS), numel(M.zm.F0));
end

function f = zneg_local(M)
% First tempo (linear interpolation) at which the 0.338 s Savitzky-Golay window's OLS or IRLS gain falls below zero.
fl = M.zm.flt;  i = find(fl.family == "SG" & abs(fl.Tw - 81/240) < 1e-12);
if numel(i) ~= 1, error("claims:zm", "%s", "0.338 s window not found"); end
g = sortrows(M.zm.G(M.zm.G.filt == i, :), "f0");  f = Inf;
for c = ["gainOLS" "gainIRLS"]
    y = g.(c);  k = find(y < 0, 1);
    if ~isempty(k) && k > 1, f = min(f, g.f0(k-1) - y(k-1) * (g.f0(k) - g.f0(k-1)) / (y(k) - y(k-1))); end
end
if ~isfinite(f), error("claims:zm", "%s", "0.338 s gain never falls below zero"); end
end

function s = sb_local(M, which)
% Supplement B literals (2026-10-24). Signed values print as the manuscript does: "+0.004", "−0.008" (U+2212).
sg = @(x, d) strrep(sprintf("%+.*f", d, x), "-", char(8722));
r29 = @(T, ds, p) sbpick_local(T, ds, p);
r30 = @(ds, p, v) M.r30.R30a(M.r30.R30a.dataset == ds & M.r30.R30a.pipeline == p & M.r30.R30a.variant == v, :);
al  = @(ds, arm) M.r30c.A(M.r30c.A.dataset == ds & M.r30c.A.arm == arm, :);
cell3 = @(f) [r30("Fraser", "SG-IRLS", "V1").(f) r30("Hickman_PLAC", "SG-IRLS", "V1").(f) r30("Fraser", "BWFD-OLS", "V1").(f)];
switch which
    case "sfFraserVetted"
        r = r29(M.r29s.rows, "Fraser", "SG-IRLS");
        s = sprintf("swapping the whole noise raised the slope by %s per doubling (z %.1f); swapping the structure alone by %s (z %.1f); swapping the size alone by %s (z %.1f)", ...
            sg(r.Delta, 3), r.zDelta, sg(r.T, 3), r.zT, sg(r.L, 3), r.zL);
    case "cropFraser"
        r = r29(M.r29c.rows, "Fraser", "SG-IRLS");
        s = sprintf("as in the cropped-only refit (%s, z %.1f) and in the first run. The size contribution is smaller and not stable across the three analyses (%s in the cropped-only refit)", sg(r.T, 3), r.zT, sg(r.L, 3));
    case "sfFraserLegacy"
        r = r29(M.r29s.rows, "Fraser", "BWFD-OLS");
        s = sprintf("Legacy profile: structure (%s, z %.1f; whole noise %s, z %.1f)", sg(r.T, 3), r.zT, sg(r.Delta, 3), r.zDelta);
    case "plac"
        o = r29(M.r29.R29a, "Hickman_PLAC", "SG-IRLS");  c = r29(M.r29c.rows, "Hickman_PLAC", "SG-IRLS");  f = r29(M.r29s.rows, "Hickman_PLAC", "SG-IRLS");
        k = M.r29c.kept.nKeptRows(M.r29c.kept.dataset == "Hickman_PLAC");
        s = sprintf("(whole noise %s, z %.1f), but neither robustness analysis did (cropped-only %s, z %s, on %d trials; seam-free %s, z %s)", ...
            sg(o.Delta, 3), o.zDelta, sg(c.Delta, 3), sg(c.zDelta, 1), k, sg(f.Delta, 3), sprintf("%.1f", f.zDelta));
    case "wrapped"
        w = 100 * (1 - M.r29c.kept.shareCropped);
        s = sprintf("In the first run, %.0f to %.0f%% of the borrowed segments were shorter", min(w), max(w));
    case "sfColour"
        d = M.r30c.A.dAlpha(ismember(M.r30c.A.arm, ["AXdF" "AXoF" "AXlF"]));
        s = sprintf("(median spectral exponent within %.2f)", ceil(max(abs(d)) * 100) / 100);
    case "rowOwn"
        x = cell3("sXr");  s = sprintf("| Slope with the trial's own noise | %s | %s | %s |", sg(x(1), 3), sg(x(2), 3), sg(x(3), 3));
    case {"rowXn" "rowXc" "rowXa" "rowXs"}
        lab = dictionary(["rowXn" "rowXc" "rowXa" "rowXs"], ["White noise" "Colour only" "Colour and relation between axes" "Relation between axes broken"]);
        a = extractAfter(which, "row");  d = cell3("d" + a);  z = cell3("z" + a);
        s = sprintf("| %s | %s (z %s) | %s (z %s) | %s (z %s) |", lab(which), sg(d(1), 3), sg(z(1), 1), sg(d(2), 3), sg(z(2), 1), sg(d(3), 3), sg(z(3), 1));
        s = strrep(s, "(z +", "(z ");
    case "v3"
        s = sprintf("Fraser's colour-only difference falls to %s (z %.1f) and PLAC's axis-swap difference reaches z %.1f", ...
            sg(r30("Fraser", "SG-IRLS", "V3").dXc, 3), r30("Fraser", "SG-IRLS", "V3").zXc, r30("Hickman_PLAC", "SG-IRLS", "V3").zXs);
    case "colour"
        fc = al("Fraser", "AXcF");  fa = al("Fraser", "AXaF");  pc = al("Hickman_PLAC", "AXcF");  pa = al("Hickman_PLAC", "AXaF");
        fs = al("Fraser", "AXsF");  ps = al("Hickman_PLAC", "AXsF");
        s = sprintf("On Fraser their median spectral exponent fell %.2f and %.2f below the trial's own (%.2f and %.2f, against %.2f), outside the 0.05 criterion. On PLAC colour only kept it (%s) and colour and axes did not (%s). The axis-swap replacement kept it in both datasets (%s and %s)", ...
            -fc.dAlpha, -fa.dAlpha, fc.medianAlpha, fa.medianAlpha, fc.medAXr, sg(pc.dAlpha, 2), sg(pa.dAlpha, 2), sg(fs.dAlpha, 2), sg(ps.dAlpha, 2));
        if ~(~fc.colourPreserved && ~fa.colourPreserved && pc.colourPreserved && ~pa.colourPreserved && fs.colourPreserved && ps.colourPreserved), s = "colour-check verdicts changed"; end
    case "white"
        n = cell3("sXn");  o = cell3("sXr");
        s = sprintf("to about zero in Hickman PLAC (%s) and in Fraser's legacy profile (%s), and to less than half in Fraser's vetted profile (%s against %s)", sg(n(2), 3), sg(n(3), 3), sg(n(1), 3), sg(o(1), 3));
        if ~(abs(n(2)) < 0.0108 && abs(n(3)) < 0.0108 && n(1) / o(1) < 0.5), s = "white-noise reading no longer holds"; end
    case "placColour"
        r = r30("Hickman_PLAC", "SG-IRLS", "V1");
        s = "the colour-only difference is inside the SEM threshold of 0.0108 per doubling and not distinguishable from zero";
        if ~(abs(r.dXc) < 0.0108 && abs(r.zXc) < 2), s = "PLAC colour-only no longer within threshold"; end
    case "limits"
        k = M.r29s.corr;  n = @(ds) unique(k.nSkipped(k.dataset == ds));
        s = sprintf("(%d in Fraser, %d in PLAC)", n("Fraser"), n("Hickman_PLAC"));
    case "limitsCorr"
        c = M.r29s.corr.corrLfwAlpha(M.r29s.corr.dataset == "Hickman_PLAC" & ismember(M.r29s.corr.arm, ["AXdF" "AXoF"]));
        s = sprintf("(%s to %s), which is unexplained", sg(max(c), 2), sg(min(c), 2));
    case "pointer"
        n = cell3("sXn");  o = cell3("sXr");
        s = "with white noise in place of each trial's real noise the slope falls in every case, to about zero in Hickman PLAC and to less than half in Fraser's vetted profile";
        if ~(all(n < o) && abs(n(2)) < 0.0108 && n(1) / o(1) < 0.5), s = "pointer reading no longer holds"; end
    otherwise
        error("claims:sb", "%s", "Unknown Supplement B literal: " + which);
end
end

function r = sbpick_local(T, ds, p)
% One V1 row of a saved R29 table (the first-run table has no variant column).
k = T.dataset == ds & T.pipeline == p;
if ismember("variant", T.Properties.VariableNames), k = k & T.variant == "V1"; end
r = T(k, :);
if height(r) ~= 1, error("claims:sb", "%s", sprintf("Expected one row for %s %s, found %d", ds, p, height(r))); end
end

function s = dag_local(M, which)
% §2.5a literals (2026-10-24): each returns the manuscript's wording only if tempoDAG_v002's saved result supports it.
r = M.dag.res;  none = "{} (no adjustment needed)";
switch which
    case "route"
        s = "in the graph the route is blocked given true tempo, but not given the estimated f0, nor given the estimated f0 and the estimated noise together";
        ok = isequal(double(r.b.routeBlocked), [1 0 0]);
    case "ident"
        s = "Under reading (b) it needs no adjustment. Under reading (a) no observed variable identifies it: it needs the subject traits that drive both group and tempo, and if tempo itself sets the condition, it needs true tempo";
        ok = r.b.adj_Bg_obs == none && r.a.adj_Bg_obs == "NOT IDENTIFIED" && r.a.adj_Bg_obs_T == "NOT IDENTIFIED" && r.a.adj_Bg_obs_T_U == "{U}" && r.a2.adj_Bg_obs_T == "{T}";
    case "mediator"
        s = "it is itself biased if subject traits act on both tempo and " + char(946) + "_gen, because conditioning on tempo then links the group to those traits; measuring the traits removes that bias";
        ok = isequal(double(r.b.mediatorConfounding), [0 0 1]);
    case "collider"
        s = "tempo and " + char(946) + "_gen are independent given the group and subject traits, but become dependent once trials are selected by the gate, or once the analysis conditions on " + char(946) + "_obs, the estimated f0, " + char(946) + "_gen* or the record";
        ok = all(arrayfun(@(k) isequal(double(r.(k).collider), [1 0 0 0 0 0]), ["a" "a2" "b"]));
    otherwise
        error("claims:dag", "%s", "Unknown §2.5a literal: " + which);
end
if ~ok, s = "tempoDAG_v002 result no longer supports: " + which; end
end

function s = pgc_local(M, which)
% D-T7 literals (2026-10-24), from powerGroupContrast_v001.mat (R).
R = M.pgc.R;  v = @(p, o, a, f) R.(f)(R.pipeline == p & R.outcome == o & R.analysis == a);
g = R(R.outcome == "betaGen", :);  m = g(g.analysis ~= "all", :);
switch which
    case "gen"
        nS = v("SG-IRLS", "betaGen", "all", "nSubjCTRL") + v("SG-IRLS", "betaGen", "all", "nSubjASD");
        L = g(g.pipeline == "BWFD-OLS", :);
        s = sprintf("with %d subjects it had 80%% power only for a %s_gen* difference of at least %.3f over all trials (%.3f and %.3f at half-octave and 0.1-doubling matching) in the vetted profile and %.3f to %.3f in the legacy one", ...
            nS, char(946), v("SG-IRLS", "betaGen", "all", "MDE80"), v("SG-IRLS", "betaGen", "matched0.5", "MDE80"), v("SG-IRLS", "betaGen", "matched0.1", "MDE80"), min(L.MDE80), max(L.MDE80));
    case "obs"
        o = R(R.outcome == "betaObs", :);
        s = sprintf("and for a %s_obs difference of %.3f to %.3f. Equivalence within %s0.03 is not shown for any of these contrasts", char(946), min(o.MDE80), max(o.MDE80), char(177));
        if any(R.equivalentAtMDC), s = "an equivalence call changed"; end
    case "hw"
        hw = (m.CIhi_pub - m.CIlo_pub) / 2;
        s = sprintf("though the matched intervals (95%% half-widths %.2f to %.2f) are too wide", min(hw), max(hw));
    case "disc"
        hw = (m.CIhi_pub - m.CIlo_pub) / 2;
        s = sprintf("95%% intervals %s%.2f to %s%.2f, too wide for equivalence, and 80%% power only for differences of about %.2f to %.2f", char(177), min(hw), char(177), max(hw), min(g.MDE80), max(g.MDE80));
    otherwise
        error("claims:pgc", "%s", "Unknown D-T7 literal: " + which);
end
end

function s = ipw_local(M, which)
% D-T6 literals (2026-10-24), from ipwGateTempoSlope_v001.mat (S slopes, C Cook contrasts, W weights).
S = M.ipw.S;  C = M.ipw.C;  W = M.ipw.W;
sl = @(p) S.lfw(S.pipeline == p & S.propensity ~= "none");  s0 = @(p) S.lfw(S.pipeline == p & S.propensity == "none");
switch which
    case "slopes"
        s = sprintf("it moves from %+.4f to %+.3f to %+.3f in the vetted profile and from %+.4f to %+.3f to %+.3f in the legacy one", ...
            s0("SG-IRLS"), min(sl("SG-IRLS")), max(sl("SG-IRLS")), s0("BWFD-OLS"), min(sl("BWFD-OLS")), max(sl("BWFD-OLS")));
        if any(sl("SG-IRLS") < s0("SG-IRLS")) || any(sl("BWFD-OLS") < s0("BWFD-OLS")), s = "a weighted slope fell below the unweighted one"; end
    case "cook"
        d = [];
        for p = ["SG-IRLS" "BWFD-OLS"], for an = ["all" "matched0.5"]
            k = C.pipeline == p & C.analysis == an;  d = [d; abs(C.b(k & C.propensity ~= "none") - C.b(k & C.propensity == "none"))]; %#ok<AGROW>
        end, end
        s = sprintf("and the Cook %s_gen* contrast changes by at most %.3f", char(946), ceil(max(d) * 1000) / 1000);
    case "weights"
        k = W.pipeline == "SG-IRLS";
        s = sprintf("(largest weight %.1f, effective sample size about %s of %s trials)", max(W.maxW), ...
            regexprep(sprintf("%d", round(min(W.ESS(k)), -2)), "(\d)(?=(\d{3})+$)", "$1,"), regexprep(sprintf("%d", W.nAnalysed(find(k, 1))), "(\d)(?=(\d{3})+$)", "$1,"));
    case "part4"
        x = sl("SG-IRLS");
        s = sprintf("weighting for the gate's selection on tempo and noise does not reduce it (%+.3f to %+.3f", min(x), max(x));
        if any(x < s0("SG-IRLS")), s = "a weighted slope fell below the unweighted one"; end
    otherwise
        error("claims:ipw", "%s", "Unknown D-T6 literal: " + which);
end
end

function s = d16_local(M, which)
% Rerun D16 stage 1 literals (2026-10-24): obsBetaRateCutoff_v002.mat (clean subsampling), _v001.mat (anti-aliased).
r = M.d16.r1;  Q1 = M.d16.Q1;  Q3 = M.d16.Q3;  sg3 = @(x) strrep(sprintf("%+.3f", x), "-", char(8722));
oBW = Q3.medianObsChange(Q3.pipeline == "BWFD-OLS");  oSG = Q3.medianObsChange(Q3.pipeline == "SG-IRLS");
switch which
    case "slope"
        s = sprintf("lowers %s by %.3f per halving of the rate (95%% person-bootstrap interval %.3f to %.3f; %s trials, %d people)", ...
            char(945), r.mean, r.lo, r.hi, thou_local(r.nTrials), r.nPersons);
        if r.negligible, s = "R1 says negligible: the correction is wrong"; end
    case "medians"
        a = @(f) Q1.alphaMedian(Q1.fs == f);
        s = sprintf("the median %s is %.2f, %.2f and %.2f", char(945), a(240), a(120), a(60));
    case "aa"
        r0 = M.d16a.r1;
        s = sprintf("the slope is %s (%s to %s)", mn_local(r0.mean, 3), mn_local(r0.lo, 3), sg3(r0.hi));
    case "zarandi"
        s = sprintf("Fraser still reads %.2f against", Q1.alphaMedian(Q1.fs == 120));
    case "sweep"
        if ~isequal(string(Q3.setting(Q3.pipeline == "BWFD-OLS"))', ["S240fc05" "S240fc07" "S240fc14" "S240fc20"]), s = "sweep settings out of order"; return, end
        s = sprintf("changes %s_obs by a median of %s, %s, %s and %s per trial", char(946), sg3(oBW(1)), sg3(oBW(2)), sg3(oBW(3)), sg3(oBW(4)));
    case "sweepSG"
        s = "(Butterworth-OLS; Savitzky-Golay-IRLS within 0.0005 of each)";
        if max(abs(oSG - oBW)) >= 0.0005, s = sprintf("SG-IRLS differs from BWFD-OLS by up to %.4f", max(abs(oSG - oBW))); end
    case "pred"
        s = sprintf("predicts median changes below 4 %s 10%s%s", char(215), char(8315), char(8309));
        if max(abs(Q3.medianPredChange)) >= 4e-5, s = sprintf("largest median predicted change %.2g", max(abs(Q3.medianPredChange))); end
    otherwise
        error("claims:d16", "%s", "Unknown D16 literal: " + which);
end
end

function s = sgr_local(M, which)
% §2.8 literals (2026-10-24), from checkSGCurvatureResidue_v001.mat (M machine table, P per dataset, W Fraser sweep).
T = M.sgr.M;  P = M.sgr.P;  W = M.sgr.W;  c = T(~T.residue, :);  r = T(T.residue, :);  g = @(d) P(P.dataset == d, :);
switch which
    case "machine"
        s = sprintf("to within 4 %s 10%s%s in every pipeline on %s of the %s Fraser trials. The other %d expose", char(215), char(8315), char(8312), ...
            thou_local(c.nTrials), thou_local(sum(T.nTrials)), r.nTrials);
        if max(c.maxDiff) >= 4e-8, s = sprintf("clean-trial machine difference %.2g", max(c.maxDiff)); end
    case "residueDiff"
        s = sprintf("by up to %.3f (OLS), %.3f (LMLS) and %.4f (IRLS)", r.maxDiff(2), r.maxDiff(4), r.maxDiff(6));
        if max(r.maxDiff([1 3 5])) >= 1e-6, s = "BWFD differs between machines on residue trials"; end
    case "drop"
        f = g("Fraser");
        s = sprintf("Savitzky-Golay-OLS by a median of %.3f on those trials and Savitzky-Golay-IRLS by a median of %.4f, with %d trials moving by at least the SEM threshold", ...
            f.medianAbsChangeIfDropped(1), f.medianAbsChangeIfDropped(3), f.nSGIRLSchangeAtLeastSEM);
    case "datasets"
        s = sprintf("occurs in %d of Dhieb's %d trials, %d of Hickman PLAC's %d and %d of Cook ASD's %d, and in none of Cook CTRL's, Hickman HALO's or Zarandi's", ...
            g("Dhieb").nResidue, g("Dhieb").nTrials, g("Hickman_PLAC").nResidue, g("Hickman_PLAC").nTrials, g("Cook_ASD").nResidue, g("Cook_ASD").nTrials);
        if any([g("Cook_CTRL").nResidue g("Hickman_HALO").nResidue g("Zarandi").nResidue] > 0), s = "a dataset said to have none has residue trials"; end
    case "narrow"
        k = [M.d16.SET.name] == "S240fc20";
        s = sprintf("%d of the Fraser trials with a %d-sample window", W.nResidue(W.setting == "S240fc20"), M.d16.frames(k));
    otherwise
        error("claims:sgr", "%s", "Unknown residue literal: " + which);
end
end
