%% checkClaimsV006_v002.m
% Claims ledger for the numbers added in skeleton v006 (Part 3 addenda 1-3, Supplement C), same mechanism as checkManuscriptClaims_v001: each row rebuilds a literal
% from a saved .mat and requires it in the manuscript after an anchor (whitespace normalised; U+2212 read as '-'). MISMATCH = prose and file disagree;
% ANCHOR LOST = the anchor phrase was reworded. Ends with an error if any claim fails, after reporting every claim. Items printed only in logs are listed as MANUAL, not checked.
% Run from src/. Reads docs/ (highest PaperDraft_IdentifiabilitySkeleton_v*.md) and results/. Writes results/checkClaimsV006_v002.mat.   Fraser, D.S. (2026)  v002 (2026-11-04: rows HS.1-HS.26, harmonic stripping)
ROOT = fileparts(fileparts(mfilename("fullpath")));  R = fullfile(ROOT, "results");
ml = dir(fullfile(ROOT, "docs", "PaperDraft_IdentifiabilitySkeleton_v*.md"));  nm = sort(string({ml.name}));  msName = nm(end);
txt = regexprep(strrep(string(fileread(fullfile(ROOT, "docs", msName))), char(8722), '-'), "\s+", " ");  fprintf("Manuscript checked: %s\n", msName);
H1 = load(fullfile(R, "analyseH1_v001.mat")).T;  S2 = load(fullfile(R, "D16_S2_v001.mat")).rows;  PC = load(fullfile(R, "D16_PairedCeiling_v001.mat")).T;
PP = load(fullfile(R, "summariseI2PerProfile_v001.mat")).T;  CT = load(fullfile(R, "summariseI2Contrast_v001.mat")).T;  TL = load(fullfile(R, "templateTailStratify_v002.mat")).T;
XC = load(fullfile(R, "analyseCrossAxisCoherence_v002.mat")).T;  GF = load(fullfile(R, "analyseGeneratorFidelity_v001.mat")).T;  SP = load(fullfile(R, "analyseA60pSupport_v001.mat")).T;  PT = load(fullfile(R, "analyseA60pPattern_v001.mat")).T;  DC = load(fullfile(R, "analyseD16Checks_v001.mat")).T;
SR = load(fullfile(R, "checkSGResidueSensitivity_v001.mat"));  SC = load(fullfile(R, "checkSGResidueCCC_v001_v015.mat")).T;  EA = load(fullfile(R, "templateEndAnchorTest_v001.mat"));  TO = load(fullfile(R, "checkTemplateTailOmission_v001.mat"));
NV = load(fullfile(R, "constellationMetrics_v005.mat"));  NB = load(fullfile(R, "checkCCCNullBaseline_v001.mat"));  V42 = load(fullfile(R, "compare42CellVerdict_v002_counts.mat")).V;
T6 = load(fullfile(R, "buildTable6aV015_v001_v015.mat"));
t6c = @(k, j) sprintf("%.3f (%.3f±%.3f)", T6.T6a.med((j-1)*6+k), T6.T6a.mean((j-1)*6+k), T6.T6a.sd((j-1)*6+k));
t6row = @(k) "| " + T6.info.pipelines(k) + " | " + join(arrayfun(@(j) string(t6c(k, j)), 1:5), " | ") + " |";   % no bold (the bold SG-LMLS row was removed 2026-10-31)
ncrow = @(p) "| " + p + " | " + join(arrayfun(@(j) string(sprintf("%.1f%%", T6.TNC{char(p), j})), 1:5), " | ") + " |";
PH = load(fullfile(R, "checkPoolCCCHeterogeneity_v004.mat")).T;
ph = @(ds, v) PH.(v)(PH.dataset == ds);
B6 = load(fullfile(R, "buildTable6bV015_v001.mat")).T;  B6 = B6(B6.basis == "N200", :);
b6row = @(ds, nm) "| " + nm + " | " + sprintf("%.4f (%.4f) [%.3f, %.3f]", B6.cluster(B6.dataset == ds), B6.flat(B6.dataset == ds), B6.ciLo(B6.dataset == ds), B6.ciHi(B6.dataset == ds)) + " | " + sprintf("%.3f ± %.3f", B6.mean(B6.dataset == ds), B6.sd(B6.dataset == ds)) + " | " + addcomma(B6.Nfin(B6.dataset == ds)) + " |";
EAH = load(fullfile(R, "analyseEndAnchorHeadline_v001.mat")).out;  EAm = load(fullfile(R, "compareMatchedTrialsEndAnchor_v001.mat")).out.T;   % end-anchored rerun (Addenda 74-77)
ee = @(p, n, v) EAH.E.(v)(EAH.E.pipeline == p & EAH.E.estimand == n);
ec = @(ds, v) EAH.C.(v)(EAH.C.dataset == ds);
c1row = @(ds, nm) "| " + nm + " | " + sprintf("%d / %d", EAm.nA(EAm.dataset == ds), EAm.nB(EAm.dataset == ds)) + " | " + dot_local(EAm.allA(EAm.dataset == ds)) + " → " + dot_local(EAm.allB(EAm.dataset == ds)) + " | " + dot_local(EAm.matchedA(EAm.dataset == ds)) + " → " + dot_local(EAm.matchedB(EAm.dataset == ds)) + " | " + sd_local(EAm.medPairedD(EAm.dataset == ds)) + " [" + sd_local(EAm.pdLo(EAm.dataset == ds)) + ", " + sd_local(EAm.pdHi(EAm.dataset == ds)) + "] | " + dot_local(EAm.onlyB_ccc(EAm.dataset == ds)) + " |";
NMd = load(fullfile(R, "printDatasetNoiseMedians_v001.mat")).T;  nmv = @(ds, b, v) NMd.(v)(NMd.dataset == ds & NMd.basis == b);
AP = load(fullfile(R, "checkAlphaIRAProvenance_v001.mat")).T;
RS = load(fullfile(R, "constellationRSA_HPC_v009.mat")).results;  RT = struct2table(RS.perDataset);  RT.dataset = string(RT.dataset);
sgn = @(x) strrep(strrep(sprintf("%+.3f", x), "0.", "."), "-", char(8722));   % signed R, three decimals, no leading zero
pf = @(p) strrep(sprintf("%.4f", max(p, 1e-4)), "0.", ".");
rsrow = @(ds, nm) "| " + nm + " | " + sprintf("%d", RT.nTrials(RT.dataset == ds)) + " | " + sgn(RT.mantelR(RT.dataset == ds)) + " | " + pf(RT.pValue(RT.dataset == ds)) + " |";
v42 = @(a) sprintf("%d/%d/%d", a(1), a(2), a(3));
nv = @(ds, v) NV.T.(v)(NV.T.dataset == ds & NV.T.noiseModel == "shaped_xu");   % shaped_xu row of constellationMetrics_v005
nvm = @(ds, m, v) NV.T.(v)(NV.T.dataset == ds & NV.T.noiseModel == m);
pt = @(sc, v) NV.PT.(v)(NV.PT.scope == sc & NV.PT.noiseModel == "shaped_xu");
sk = @(sc, v) NB.S.(v)(NB.S.scope == sc & NB.S.noiseModel == "shaped_xu");
SA = load(fullfile(ROOT, "src", "semAgreement_v003.mat")).results.ratio;
sar = @(ds, v) SA.(v)(SA.dataset == ds);
PB = load(fullfile(R, "plotSecondaryMetricsV015_v001.mat")).BA;
pb = @(ds, v) PB.(v)(PB.dataset == ds);
sp2 = @(x) strrep(sprintf("%+.2f", x), "0.", ".");   % signed, two decimals, no leading zero (as the prose)
sg = @(x, d) string(sprintf("%+.*f", d, x));   % signed, ASCII minus
f2 = @(x, d) string(sprintf("%.*f", d, x));
h1 = @(ds) H1(strcmp(strrep(string({H1.dataset}), "_", " "), ds));  xc = @(ds) XC(strcmp(string({XC.dataset}), ds));
gf = @(ds, m) GF(strcmp(string({GF.dataset}), ds) & strcmp(string({GF.model}), m));
rw = @(st, pr) S2(string({S2.setting}) == st & string({S2.profile}) == pr & string({S2.g}) == "own" & string({S2.subset}) == "all");
pc = @(st, pr) PC(string({PC.setting}) == st & string({PC.profile}) == pr);
pp = @(file) PP(contains(string({PP.file}), file));   % e.g. "Fraser_v005", "nonEllipse", "PLAC", "Zarandi"
col = @(t, nme) find(t.names == nme);
dsl = ["Cook CTRL" "Cook ASD" "Hickman PLAC" "Hickman HALO" "Dhieb" "Fraser"];
bootD = arrayfun(@(d) h1(d).dMed(2), dsl);  xuD = arrayfun(@(d) h1(d).dMed(1), ["Cook CTRL" "Cook ASD" "Dhieb" "Fraser"]);
z = h1("Zarandi");  zc = sort(z.dCI(2, :));
A1 = "Part 3 addendum 1";  A2 = "Part 3 addendum 2";  A3 = "Part 3 addendum 3";
pl = pp("PLAC");  zr = pp("Zarandi");  fb = pp("Fraser_v005");  ne = pp("nonEllipse");  np = pp("nonPowerLaw");  sv = pp("shapeVar");  sd = pp("signalDep");  c1 = CT(contains(string({CT.file}), "Fraser_v005"));  c6 = CT(contains(string({CT.file}), "PLAC"));
cc = @(t) t.cost(1:6)';  hrow = @(t, a) t.hit(end, a);
irl = col(pl, "SG-IRLS");  lml = col(pl, "SG-LMLS");  ols = col(pl, "BWFD-OLS");
ceil = @(st) sprintf("%s and %s", f2(pc(st, "SG-IRLS").cccMaxPair, 2), f2(pc(st, "BWFD-OLS").cccMaxPair, 2));
allc = [PC.cccObs];  allm = [PC.cccMaxPair];  cis = vertcat(PC.cccCI);
tl = @(ds) TL(strcmp(string({TL.dataset}), ds));
U15 = load(fullfile(R, "checkDatasetMedianUncertainty_v002_v015.mat"));  U19 = load(fullfile(R, "checkDatasetMedianUncertainty_v002_v019.mat"));   % Table C2 (GD12)
K15 = load(fullfile(ROOT, "src", "loopClosureVarDecomp_v013_gated_all6_v015.mat")).synthMed_k3;  K15 = K15(1);
K19 = load(fullfile(ROOT, "src", "loopClosureVarDecomp_v013_gated_all6_v019.mat")).synthMed_k3;  K19 = K19(1);
nl0 = @(x, d) regexprep(string(sprintf("%.*f", d, x)), "^0\.", ".");   % no leading zero
f3 = @(x, varargin) string(sprintf("%.3f", x));   % three decimals, leading zero (prose style)
c2m = @(U, ds) nl0(U.rowsM.median(U.rowsM.dataset == ds & U.rowsM.estimator == "constellation"), 3);
c2d = @(U, nm, est) c2cell(U.rowsD(U.rowsD.contrast == nm & U.rowsD.estimator == est, :));
c2p = @(K, ci, m) nl0(K.(m), 3) + " [" + nl0(K.(ci)(1), 3) + ", " + nl0(K.(ci)(2), 3) + "]";
c2h = @(K) "[" + nl0(K.ciHKSJ(1), 3) + ", " + nl0(K.ciHKSJ(2), 3) + "]";
V15 = load(fullfile(R, "checkDatasetMedianUncertainty_v003_cluster_v015.mat"));  V19 = load(fullfile(R, "checkDatasetMedianUncertainty_v003_cluster_v019.mat"));   % cluster basis (GD5)
c3d = @(U, nm, est) c3cell(U.rowsD(U.rowsD.contrast == nm & U.rowsD.estimator == est, :));
VC = load(fullfile(R, "checkVGFInterceptCoupling_v001.mat"));
PB = load(fullfile(R, "tempoSim", "checkSuppBPairedSE_v001.mat")).rows;  PB.dataset = string(PB.dataset);  PB.arm = string(PB.arm);  PB.pipeline = string(PB.pipeline);  PB.spec = string(PB.spec);
pbr = @(ds, arm, p) PB(PB.dataset == ds & PB.arm == arm & PB.pipeline == p & PB.spec == "V1", :);
pbc = @(ds, arm, p) c3cell(pbr(ds, arm, p));
pbrow = @(nm, arm) "| " + nm + " | " + pbc("Fraser", arm, "SG-IRLS") + " | " + pbc("Hickman_PLAC", arm, "SG-IRLS") + " | " + pbc("Fraser", arm, "BWFD-OLS") + " |";
PB1 = PB(PB.spec == "V1", :);  PBn = PB1(~(PB1.arm == "AXnF" & PB1.dataset == "Fraser"), :);  PBw = PB1(PB1.arm == "AXnF" & PB1.dataset == "Fraser", :);
PBf = PB1(PB1.dataset == "Fraser", :);  PBp = PB1(PB1.dataset == "Hickman_PLAC", :);
PS = load(fullfile(R, "checkPipelineSettings_v001.mat"));
GM = load(fullfile(R, "checkGenModelSpec_v001.mat"));
BS = load(fullfile(R, "checkBootstrapSpec_v001.mat"));
SF = load(fullfile(R, "checkSensorFloor_v001.mat")).out;
SS = load(fullfile(R, "checkSpecSetup_v001.mat"));
LC = load(fullfile(R, "checkLmmSpec_v001.mat"));
SJ = load(fullfile(R, "checkSimexJoinSpec_v001.mat")).out;
TB = load(fullfile(R, "checkTempoBinSpec_v001.mat"));
CC = load(fullfile(R, "checkCookContrastSpec_v001.mat"));
RP = load(fullfile(R, "checkReplayDiffSpec_v001.mat"));
PF = load(fullfile(R, "plotConstellationMetrics_v003.mat"));  DN = load(fullfile(R, "checkDhiebSnapSpec_v001.mat"));  NI = load(fullfile(R, "checkNodeInterpolationSpec_v001.mat"));  MS = load(fullfile(R, "checkMedianShareV015Spec_v001.mat"));  msRow = @(d) find(MS.dataset == d);
CF = load(fullfile(R, "checkV015ConfigSpec_v001.mat"));  cfFs = @(n) CF.fsDS(CF.datasets == n);
HS = load(fullfile(R, "checkHalfSampleSpec_v001.mat"));
NP = load(fullfile(R, "checkNoisePhaseSpec_v001.mat"));
DCs = load(fullfile(R, "checkDatasetCountsSpec_v001.mat"));  numWords = ["zero" "one" "two" "three" "four" "five" "six" "seven"];
rsaF = dir(fullfile(R, "constellationRSA_HPC_v009.mat"));  rsaName = rsaF(1).name;
GMd = @(d) GM.Dstar(GM.Dstar.dataset == d, :);
mcF = dir(fullfile(R, "positiveControlMatchedByDataset_v002_*.mat"));  [~, mcI] = max([mcF.datenum]);  MC = load(fullfile(mcF(mcI).folder, mcF(mcI).name)).out;
MCd = @(d) MC(MC.dataset == d, :);
m3F = dir(fullfile(R, "positiveControlMatchedByDataset_v003_*.mat"));  [~, m3I] = max([m3F.datenum]);  M3 = load(fullfile(m3F(m3I).folder, m3F(m3I).name)).out;
M3d = @(d) M3(M3.dataset == d, :);
m4F = dir(fullfile(R, "positiveControlMatchedByDataset_v004_*.mat"));  [~, m4I] = max([m4F.datenum]);  M4 = load(fullfile(m4F(m4I).folder, m4F(m4I).name)).out;
M4d = @(d) M4(M4.dataset == d, :);
in5m = ismember(M4.dataset, ["Fraser","Cook CTRL","Cook ASD","Hickman PLAC","Hickman HALO"]);
clin4 = ismember(M4.dataset, ["Cook CTRL","Cook ASD","Hickman PLAC","Hickman HALO"]);
in5 = ["Fraser","Cook CTRL","Cook ASD","Hickman PLAC","Hickman HALO"];
DG = load(fullfile(R, "checkMatchedControlDiagnosis_v001.mat"));  DGm = @(c) median(DG.TS.lambda(DG.TS.dataset == c), "omitnan");
inDomRaw = max(MC.subjMed_lambda(~ismember(MC.dataset, ["Zarandi","Dhieb"])));
% Critique 2026-11-01 (Abstract coverage of one-third, v019 range, Supplement C template-sensitivity paragraph)
cvB = @(ds) B6.ciLo(B6.dataset == ds) <= 1/3 & B6.ciHi(B6.dataset == ds) >= 1/3;   % Table 6b (v015, inversion draw included)
cvV = @(V) V.rowsM.ciLo(V.rowsM.estimator == "constellation") <= 1/3 & V.rowsM.ciHi(V.rowsM.estimator == "constellation") >= 1/3;
vci = @(V, ds, f) V.rowsM.(f)(V.rowsM.dataset == ds & V.rowsM.estimator == "constellation");
zin = @(V, nm, est) V.rowsD.ciLo(V.rowsD.contrast == nm & V.rowsD.estimator == est) <= 0 & V.rowsD.ciHi(V.rowsD.contrast == nm & V.rowsD.estimator == est) >= 0;
CV = load(fullfile(R, "compareTable6bCoverage_v001.mat")).T;  tv = @(b, ds, f) CV.(f)(CV.template == b & CV.dataset == ds);   % Table 6b intervals, v015 and v019
ci3 = @(b, ds) string(sprintf("[%.3f, %.3f]", tv(b, ds, "ciLo"), tv(b, ds, "ciHi")));  n19 = nnz(CV.includesOneThird(CV.template == "v019"));
% Sampling-rate bullets of section 2.3 (formerly [GENERATED, UNVERIFIED]; ledgered 2026-11-01)
SRL = load(fullfile(R, "reportSamplingRateLmm_v001.mat"));  SRE = load(fullfile(R, "reportSamplingRateEffects_v001.mat"));  SRA = load(fullfile(R, "reportSyntheticAliasing_v001.mat")).R;
OB1 = load(fullfile(R, "obsBetaRateCutoff_v001.mat"), "r1").r1;  OB2 = load(fullfile(R, "obsBetaRateCutoff_v002.mat"), "r1").r1;  EC = load(fullfile(R, "checkEdgeClipShare_v001.mat"));
L9n = load(fullfile(R, "L9ModelInspect_v002.mat")).out.nObs;  ali = load(fullfile(R, "reportSyntheticAliasing_v001.mat"), "aObs", "lo");
qv = @(Tb, nm) Tb.estimate(Tb.quantity == nm);  qc = @(Tb, nm, f) Tb.(f)(Tb.quantity == nm);
lt = @(nm, f) SRL.lmmTerms.(f)(SRL.lmmTerms.term == nm);  im = @(pp) SRL.impl.impliedChange60to240(SRL.impl.pipeline == pp);  bf = @(fs) SRL.byFs.shareAdequate(SRL.byFs.fs == fs);
BWp = ["BWFD_OLS" "BWFD_LMLS" "BWFD_IRLS"];  SGp = ["SG_OLS" "SG_LMLS" "SG_IRLS"];  ALLp = [BWp SGp];
c120 = arrayfun(@(pp) qv(SRE.Chg, "dbetaObs_" + pp + "_120"), ALLp);  s120 = arrayfun(@(pp) qv(SRE.Shr, "shareAtTol_dbetaObs_" + pp + "_120"), ALLp);
cBW60 = arrayfun(@(pp) qv(SRE.Chg, "dbetaObs_" + pp + "_60"), BWp);  cSG60 = arrayfun(@(pp) qv(SRE.Chg, "dbetaObs_" + pp + "_60"), SGp);  sSG60 = arrayfun(@(pp) qv(SRE.Shr, "shareAtTol_dbetaObs_" + pp + "_60"), SGp);
[srMin, iMin] = min(SRL.rows.shareAdequate);  [srMax, iMax] = max(SRL.rows.shareAdequate);
[imMin, jMin] = min(SRL.impl.impliedChange60to240);  [imMax, jMax] = max(SRL.impl.impliedChange60to240);
ecM = @(fs, v) EC.T.(v)(EC.T.fs == fs);
% Harmonic stripping, simulation side (Part 1B (a)-(d)) and the sliding target (Part 4); ledgered 2026-11-04 (Findings #225, #226, #233, #234, #249)
HA = load(fullfile(R, "tempoFoldAnalysis_v001.mat"), "At", "Ecc", "Top");   % analyseTempoFoldSweep_v001 (Block A zero-noise maps; map top by tempo)
haA = @(col, f) HA.At.(col)(HA.At.f0 == f);  haS = @(col, fs) strjoin(arrayfun(@(f) string(sprintf("%.3f", haA(col, f))), fs), ", ");
haTop = @(ds, f) sprintf("%.2f-%.2f", HA.Top.min_top(HA.Top.dataset == ds & HA.Top.f0 == f), HA.Top.max_top(HA.Top.dataset == ds & HA.Top.f0 == f));
haRow = @(f, lab) "| " + lab + " | " + haTop("Fraser", f) + " | " + haTop("Cook_CTRL", f) + " | " + haTop("Hickman_PLAC", f) + " |";
HG1 = load(fullfile(R, "predictHarmonicGain_v001.mat"), "TW", "T_rms");   % harmonic weights, theory/emulation RMS
HG3 = load(fullfile(R, "predictHarmonicGain_v003.mat"), "T_233", "T_fine", "S233");   % #249 fine-grid, edge-safe onsets
hgF = HG3.T_fine;  hgT = HG3.T_233;  hgS = HG3.S233;
hgOn = @(fam, flt, ab, reg, v) hgF.(v)(hgF.family == fam & hgF.filt == flt & abs(hgF.ab - ab) < 1e-9 & hgF.regression == reg);
bwO = hgT.scaledFine(hgT.family == "BW" & hgT.regression == "OLS");  sgO = hgT.scaledFineSafe(hgT.family == "SG" & hgT.regression == "OLS");
w3O = HG1.TW.w3(HG1.TW.regression == "OLS");  rmsE = HG1.T_rms.rmsEmulVsMeasured;  relO = hgF.relDiff(hgF.regression == "OLS");
GT = load(fullfile(R, "gainTarget_v002.mat"), "P", "S3", "B3");   % extractGainTarget_v002 (gain and target, simulated and per-trial)
gtA = GT.P.target(GT.P.source == "A" & GT.P.reliable);  gtS = GT.S3;  gtB = GT.B3;
gtIn = ~ismember(gtS.dataset, ["Zarandi" "Dhieb"]);   % the five in-domain datasets (Finding #242)
gtb = @(ds, b, v) gtB.(v)(gtB.dataset == ds & gtB.bin == b);
gtBins = @(ds) strjoin(arrayfun(@(b) string(sprintf("%.3f", gtb(ds, b, "median_target"))), ["<0.35" "0.35-0.7" "0.7-1.4" "1.4-2.8"]), ", ");
GL = load(fullfile(R, "gainTargetLMMDomain_v001.mat"), "t");  gl = GL.t(GL.t.basis == "domain5_F242" & GL.t.term == "lf", :);   % checkGainTargetLMMDomain_v001
C = {
 "H1.1 Zarandi medians",  A1, "Zarandi was the exception", "bootstrap gave " + f2(z.med(3), 3) + " against " + f2(z.med(1), 3) + ", a difference of " + sg(z.dMed(2), 3) + " (subject-clustered 95% CI " + sg(zc(1), 3) + " to " + sg(zc(2), 3) + "; " + z.n + " trials, " + z.nSubj + " subjects)"
 "H1.2 six datasets",     A1, "within about 0.02 (largest +0.021) in six of seven datasets:", "Cook CTRL " + sg(bootD(1), 3) + ", Cook ASD " + sg(bootD(2), 3) + ", Hickman PLAC " + sg(bootD(3), 3) + ", Hickman HALO " + sg(bootD(4), 3) + ", Dhieb " + sg(bootD(5), 3) + ", Fraser " + sg(bootD(6), 3)
 "H1.3 xu range",         A1, "Straight xu lowered", "by " + f2(min(abs(xuD)), 2) + " to " + f2(max(abs(xuD)), 2) + " in four datasets"
 "H1.4 anisotropy",       A1, "extremely anisotropic", "ratio " + f2(gf("Zarandi", "bootstrap").sdRatioReal, 2) + " against " + f2(min([gf("Hickman PLAC", "bootstrap").sdRatioReal gf("Fraser", "bootstrap").sdRatioReal]), 1) + " to " + f2(max([gf("Hickman PLAC", "bootstrap").sdRatioReal gf("Fraser", "bootstrap").sdRatioReal]), 1) + " elsewhere; tablet-frame x and y residuals correlated at " + f2(xc("Zarandi").rhoXY, 2)
 "H1.5 generated ratio",  A1, "extremely anisotropic", "(" + f2(gf("Zarandi", "bootstrap").sdRatioReal, 2) + " real, " + f2(gf("Zarandi", "bootstrap").sdRatioGen, 2) + " generated)"
 "H1.6 anis Spearman",    A1, "extremely anisotropic", "(Spearman " + sg(xc("Zarandi").spAnis, 2) + ")"
 "H1.7 cross-axis",       A1, "Cross-axis coherence", "(median |" + char(961) + "| " + f2(xc("Zarandi").rhoMajMin, 2) + " against " + f2(xc("Fraser").rhoMajMin, 2) + " in Fraser), coherence not higher (" + f2(xc("Zarandi").coh, 2) + " against " + f2(min([xc("Hickman PLAC").coh xc("Fraser").coh]), 2) + " to " + f2(max([xc("Hickman PLAC").coh xc("Fraser").coh]), 2) + ")"
 "H1.8 cross-axis Sp",    A1, "Cross-axis coherence", "(Spearman " + sg(xc("Zarandi").spRho, 2) + " and " + sg(xc("Zarandi").spCoh, 2) + ")"
 "H1.9 alpha deficit",    A1, "The surrogate's spectral slope", "an IRASA exponent " + f2(min(-[gf("Zarandi", "bootstrap").dAlpha gf("Hickman PLAC", "bootstrap").dAlpha gf("Fraser", "bootstrap").dAlpha]), 1) + " to " + f2(max(-[gf("Zarandi", "bootstrap").dAlpha gf("Hickman PLAC", "bootstrap").dAlpha gf("Fraser", "bootstrap").dAlpha]), 1) + " below"
 "R2.1 pass count",       A2, "Result: the registered rule is not met", sprintf("%d of %d setting-by-profile cells pass", sum([S2(string({S2.g}) == "own" & string({S2.subset}) == "all").pass]), nnz(string({S2.g}) == "own" & string({S2.subset}) == "all"))
 "R2.2 60 Hz CCC",        A2, "Result: the registered rule is not met", "(CCC " + sg(rw("A60p", "SG-IRLS").ccc, 3) + " [" + sg(rw("A60p", "SG-IRLS").cccLo, 3) + ", " + sg(rw("A60p", "SG-IRLS").cccHi, 3) + "] at 60 Hz and " + sg(rw("A120p", "SG-IRLS").ccc, 3) + " at 120 Hz for SG-IRLS)"
 "R2.3 cutoff CCC range", A2, "Result: the registered rule is not met", "(CCC " + f2(min(arrayfun(@(s) min([rw(s, "SG-IRLS").ccc rw(s, "BWFD-OLS").ccc]), ["S240fc05" "S240fc07" "S240fc14" "S240fc20"])), 2) + " to " + f2(max(arrayfun(@(s) max([rw(s, "SG-IRLS").ccc rw(s, "BWFD-OLS").ccc]), ["S240fc05" "S240fc07" "S240fc14" "S240fc20"])), 2) + ", intervals excluding 0)"
 "R2.4 ceiling 5 Hz",     A2, "What the cutoff failure means", "SG-IRLS and BWFD-OLS: " + f2(pc("S240fc05", "SG-IRLS").cccMaxPair, 2) + " and " + f2(pc("S240fc05", "BWFD-OLS").cccMaxPair, 2) + " at 5 Hz, " + f2(pc("S240fc07", "SG-IRLS").cccMaxPair, 2) + " and " + f2(pc("S240fc07", "BWFD-OLS").cccMaxPair, 2) + " at 7 Hz, " + f2(pc("S240fc14", "SG-IRLS").cccMaxPair, 2) + " and " + f2(pc("S240fc14", "BWFD-OLS").cccMaxPair, 2) + " at 14 Hz, " + f2(pc("S240fc20", "SG-IRLS").cccMaxPair, 2) + " and " + f2(pc("S240fc20", "BWFD-OLS").cccMaxPair, 2) + " at 20 Hz"
 "R2.5 ceiling range",    A2, "What the cutoff failure means", "was " + f2(min(allm), 2) + " to " + f2(max(allm), 2) + " ("
 "R2.6 observed range",   A2, "What the cutoff failure means", "against observed " + f2(min(allc), 2) + " to " + f2(max(allc), 2) + " (subset CIs up to " + f2(max(cis(:, 2)), 2) + ")"
 "R2.7 cells under 0.70", A2, "What the cutoff failure means", sprintf("out of reach even for a perfect model in %d of %d cutoff cells", nnz(allm < 0.70), numel(allm))
 "R2.8 60Hz usable",      A2, "Coverage limit at 60 Hz", "rests on " + sprintf("%s (SG-IRLS) and %s (BWFD-OLS) trials", addcomma(rw("A60p", "SG-IRLS").n), addcomma(rw("A60p", "BWFD-OLS").n))
 "R2.9 60Hz whole-curve", A2, "Coverage limit at 60 Hz", sprintf("for %s of %s trials", addcomma(PT(1).anyNaN), addcomma(SP(1).n))
 "R2.10 betaObs terciles",A2, "Coverage limit at 60 Hz", sprintf("(failure rate %.0f%% in the lowest tercile, %.0f%% in the highest)", PT(1).nanByObsTercile(1), PT(1).nanByObsTercile(3))
 "I2.1 PLAC costs",       A3, "Hickman PLAC: every profile loses hits", "+" + f2(pl.cost(ols), 1) + " pp (legacy BWFD-OLS), +" + f2(pl.cost(irl), 1) + " (SG-IRLS) and +" + f2(pl.cost(lml), 1) + " (SG-LMLS)"
 "I2.2 PLAC pooled",      A3, "Hickman PLAC: every profile loses hits", "the pooled median loses more (+" + f2(pl.cost(7), 1) + ") and has the lowest estimated hit rate (" + f2(pl.hit(7, 1), 1) + "% against " + f2(min(pl.hit(1:6, 1)), 0) + " to " + f2(max(pl.hit(1:6, 1)), 0) + "%)"
 "I2.3 Zarandi",          A3, "Hickman PLAC: every profile loses hits", "legacy " + sg(zr.cost(ols), 1) + ", SG-IRLS " + sg(zr.cost(irl), 1)
 "I2.4 nonEllipse",       A3, "Hickman PLAC: every profile loses hits", "costs " + sg(min(cc(ne)), 0) + " to " + sg(max(cc(ne)), 0) + " pp in every profile"
 "I2.5 nonPowerLaw",      A3, "Hickman PLAC: every profile loses hits", "(pooled " + sg(np.cost(7), 1) + ": estimating beats knowing)"
 "I2.6 shapeVar",         A3, "Hickman PLAC: every profile loses hits", "shape variation " + sg(min(cc(sv)), 0) + " to " + sg(max(cc(sv)), 0)
 "I2.7 signalDep",        A3, "Hickman PLAC: every profile loses hits", "signal-dependent noise " + sg(min(sd.cost([1 2 3 4])), 0) + " to " + sg(max(sd.cost([1 2 3 4])), 0) + " for the OLS and LMLS profiles but about " + sg(mean(sd.cost([5 6])), 0) + " for both IRLS profiles"
 "I2.8 LMLS contrast",    A3, "Legacy against the vetted family", "(Fraser base " + sg(c1.dHit(1, 1), 1) + " [" + sg(c1.dHitCI(1, 1, 1), 1) + ", " + sg(c1.dHitCI(2, 1, 1), 1) + "])"
 "I2.9 IRLS contrast",    A3, "Legacy against the vetted family", "(" + sg(c1.dHit(2, 1), 1) + " [" + sg(c1.dHitCI(1, 2, 1), 1) + ", " + sg(c1.dHitCI(2, 2, 1), 1) + "] in the base world)"
 "I2.10 PLAC attribution",A3, "Attribution", "est " + f2(pl.hit(7, 1), 1) + ", aor " + f2(pl.hit(7, 3), 1) + ", gor " + f2(pl.hit(7, 4), 1) + ", nor " + f2(pl.hit(7, 5), 1) + ", orc " + f2(pl.hit(7, 2), 1)
 "I2.11 geometry world",  A3, "Attribution", "wrong-geometry world " + f2(ne.hit(7, 1), 1) + ", " + f2(ne.hit(7, 3), 1) + ", " + f2(ne.hit(7, 4), 1) + ", " + f2(ne.hit(7, 5), 1) + ", " + f2(ne.hit(7, 2), 1)
 "T1 pooled slopes",      "Supplement C", "Pooled over settings", "each cycle of added zero stretch cost " + f2(-100 * tl("Hickman PLAC").fit.b(2), 0) + " pp (PLAC), " + f2(-100 * tl("Zarandi").fit.b(2), 0) + " pp (Zarandi) and " + f2(-100 * tl("Fraser").fit.b(2), 0) + " pp (Fraser)"
 "T2 PLAC [8,8]",         "Supplement C", "Wherever a setting added", "Hickman PLAC [8,8] " + tab_local(tl("Hickman PLAC"), "est[8,8]")
 "T3 Zarandi [8,5]",      "Supplement C", "Wherever a setting added", "Zarandi [8,5] " + tab_local(tl("Zarandi"), "est[8,5]")
 "T4 Fraser [8,2]",       "Supplement C", "Wherever a setting added", "Fraser [8,2] " + tab_local(tl("Fraser"), "est[8,2]")
 "R2.11 pairing check",   A2, "What the cutoff failure means", sprintf("(pairing check: median replicate-wise correlation %.2f)", DC.pairing.median)
 "R2.12 rate-arm spreads",A2, "The rate arms fail", "predicted standard deviation about " + f2(DC.rate(1).sdPred, 3) + " (60 Hz) and " + f2(DC.rate(3).sdPred, 3) + " (120 Hz) against observed " + f2(DC.rate(1).sdObs, 3) + " and " + f2(DC.rate(3).sdObs, 3)
 "T5 baseline tails",     "Supplement C", "The baseline itself carries", sprintf("such a tail in %.0f%% (Zarandi), %.0f%% (PLAC) and %.0f%% (Fraser) of trials", tl("Zarandi").baseTailShare, tl("Hickman PLAC").baseTailShare, tl("Fraser").baseTailShare)
 "V1 residue slope",      "Residue", "Leaving the residue trials out of the analyses", "from " + f2(red_local(SR.SUM, "V1 Dhieb out", "SG-IRLS"), 1) + "% to " + f2(red_local(SR.SUM, "V2 Dhieb out, residue trials out", "SG-IRLS"), 1) + "% (legacy profile " + f2(red_local(SR.SUM, "V1 Dhieb out", "BWFD-OLS"), 1) + "% to " + f2(red_local(SR.SUM, "V2 Dhieb out, residue trials out", "BWFD-OLS"), 1) + "%)"
 "V2 residue outputs",    "Residue", "Leaving the residue trials out of the analyses", "no dataset's median β_gen* by more than " + f2(max(abs(SR.DS.dMedGen)), 3) + ", the invertible share by at most " + f2(max(abs(SR.DS.dPctInv)), 1) + " pp and the primary cccMed by at most " + f2(max(abs(SC.dCCC)), 3)
 "E1 anchoring effect",   "Supplement C", "The causal test is a variant", "raised the baseline hit rate by " + ea_local(EA.R1, "Hickman PLAC", "anchoring effect [4,4]") + " pp (PLAC), " + ea_local(EA.R1, "Zarandi", "anchoring effect [4,4]") + " (Zarandi) and " + ea_local(EA.R1, "Fraser", "anchoring effect [4,4]") + " (Fraser)"
 "E2 window effects",     "Supplement C", "The causal test is a variant", "Hickman PLAC [8,8] " + ea_local(EA.R1, "Hickman PLAC", "window effect, production: [8,8] vs [4,4]") + " without anchoring, " + ea_local(EA.R1, "Hickman PLAC", "window effect, anchored: [8,8]+EA vs [4,4]+EA") + " with; Zarandi [8,5] " + ea_local(EA.R1, "Zarandi", "window effect, production: [8,5] vs [4,4]") + " and " + ea_local(EA.R1, "Zarandi", "window effect, anchored: [8,5]+EA vs [4,4]+EA")
 "E3 Fraser residual",    "Supplement C", "The causal test is a variant", "Fraser [8,2] (" + ea_local(EA.R1, "Fraser", "window effect, production: [8,2] vs [4,4]") + " and " + ea_local(EA.R1, "Fraser", "window effect, anchored: [8,2]+EA vs [4,4]+EA") + ")"
 "E4 non-tail control",   "Supplement C", "On real trials (20 repetitions", "trials without a tail changed by a median of at most " + f2(max(arrayfun(@(r) abs(r.tab.medDeltaBG(r.tab.group == "non-tail")), EA.R2)), 4)
 "E5 tail deltas",        "Supplement C", "On real trials (20 repetitions", "lowered β_gen* by a median of " + f2(-r2_local(EA.R2, "Fraser", "tail"), 3) + " (Fraser), " + f2(-r2_local(EA.R2, "Cook_CTRL", "tail"), 3) + " (Cook CTRL), " + f2(-r2_local(EA.R2, "Cook_ASD", "tail"), 3) + " (Cook ASD), " + f2(-r2_local(EA.R2, "Hickman_PLAC", "tail"), 3) + " (Hickman PLAC) and " + f2(-r2_local(EA.R2, "Hickman_HALO", "tail"), 3) + " (Hickman HALO)"
 "E6 dataset medians",    "Supplement C", "On real trials (20 repetitions", "no dataset's median changed by more than " + f2(max(arrayfun(@(r) abs(r.tab.medBGea(r.tab.group == "all") - r.tab.medBGprod(r.tab.group == "all")), EA.R2)), 3)
 "E7 omission contrast",  "Supplement C", "Leaving the tail trials out instead", "moved the vetted tempo-slope reduction from " + f2(red_local(TO.SUM, "domain, all trials", "SG-IRLS"), 1) + "% to " + f2(red_local(TO.SUM, "domain, tail trials out", "SG-IRLS"), 1) + "% and Cook CTRL's cccMed from " + dot_local(TO.A.cccMedAll(TO.A.dataset == "Cook_CTRL")) + " to " + dot_local(TO.A.cccMedDrop(TO.A.dataset == "Cook_CTRL"))
 "N1 comparator Fraser",  "Part 3", "Post hoc check: a trial-agnostic comparator", "median per-trial CCC of " + dot_local(nv("Fraser", "nullCccMed")) + " in Fraser (model " + dot_local(nv("Fraser", "cccMed")) + ")"
 "N2 comparator pooled",  "Part 3", "Post hoc check: a trial-agnostic comparator", "the model's median is " + dot_local(pt("inDomain", "cccMed")) + " [" + dot_local(pt("inDomain", "cccMedLo")) + ", " + dot_local(pt("inDomain", "cccMedHi")) + "] against " + dot_local(pt("inDomain", "nullCccMed")) + " for the comparator (model above it in " + sprintf("%.0f", 100 * NB.P.fracBeatNull(NB.P.scope == "inDomain" & NB.P.noiseModel == "shaped_xu")) + "% of trials)"
 "N3 comparator exFraser",  "Part 3", "Post hoc check: a trial-agnostic comparator", "it is " + dot_local(pt("exFraser", "cccMed")) + " against " + dot_local(pt("exFraser", "nullCccMed"))
 "N4 Fraser share",        "Part 3", "Post hoc check: a trial-agnostic comparator", sprintf("which supplies %.0f%% of the valid trials", 100 * nv("Fraser", "Nvalid") / pt("inDomain", "Nvalid"))
 "N5 main table Fraser",   "Supplement A", "Main constellation table (all seven, v015)", "| Fraser | " + dot_local(nv("Fraser", "loopCCC")) + " | " + dot_local(nv("Fraser", "cccMed")) + " [" + dot_local(nv("Fraser", "cccMedLo")) + ", " + dot_local(nv("Fraser", "cccMedHi")) + "] | " + sprintf("%.4f", nv("Fraser", "maeMed")) + " |"
 "N6 main table PLAC",     "Supplement A", "Main constellation table (all seven, v015)", "| Hickman PLAC | " + dot_local(nv("Hickman_PLAC", "loopCCC")) + " | " + dot_local(nv("Hickman_PLAC", "cccMed")) + " [" + dot_local(nv("Hickman_PLAC", "cccMedLo")) + ", " + dot_local(nv("Hickman_PLAC", "cccMedHi")) + "] | " + sprintf("%.4f", nv("Hickman_PLAC", "maeMed")) + " |"
 "N7 comparator table PLAC","Supplement A", "Trial-agnostic comparator for the pattern metrics", "| Hickman PLAC | " + dot_local(nv("Hickman_PLAC", "cccMed")) + " | " + dot_local(nv("Hickman_PLAC", "nullCccMed")) + " | " + sprintf("%.0f", 100 * nv("Hickman_PLAC", "fracBeatNull")) + "% | " + dot_local(nv("Hickman_PLAC", "rhoMed")) + " | " + dot_local(nv("Hickman_PLAC", "nullRhoMed")) + " | " + sprintf("%.1f", nv("Hickman_PLAC", "cpMDC")) + "% | " + sprintf("%.1f", nv("Hickman_PLAC", "nullCpMDC")) + "% | " + sprintf("%.3f", nv("Hickman_PLAC", "tdiEmp")) + " | " + sprintf("%.3f", nv("Hickman_PLAC", "nullTdiEmp")) + " |"
 "N8 stacked decomposition","Supplement A", "The stacked agreement over all trial-by-contrast pairs", "stacked CCC is " + dot_local(sk("inDomain", "stackedCCC")) + " (Pearson r " + dot_local(sk("inDomain", "pearsonR")) + ", SD ratio predicted/observed " + strrep(sprintf("%.2f", sk("inDomain", "sdRatio")), "0.", ".") + ", mean shift " + sg(sk("inDomain", "meanShift"), 3) + ")"
 "N9 three-model Fraser",  "Supplement A", "Three-model comparison (Finding #92; model selection", "| Fraser | " + dot_local(nvm("Fraser", "xu", "cccMed")) + " | " + dot_local(nv("Fraser", "cccMed")) + " | **" + dot_local(nvm("Fraser", "bootstrap", "cccMed")) + "** | " + dot_local(nvm("Fraser", "xu", "loopCCC")) + " | " + dot_local(nv("Fraser", "loopCCC")) + " | **" + dot_local(nvm("Fraser", "bootstrap", "loopCCC")) + "** |"
 "N10 Bland-Altman PLAC",  "Supplement A", "Bland-Altman limits of agreement (Finding #126; v015)", "| Hickman PLAC | " + addcomma(nv("Hickman_PLAC", "nPairs")) + " | " + sg(nv("Hickman_PLAC", "baBias"), 3) + " | [" + sg(nv("Hickman_PLAC", "loaLo"), 3) + ", " + sg(nv("Hickman_PLAC", "loaHi"), 3) + "] |"
 "N11 TDI multiples",      "Supplement A", "TDI and Coverage Probability (Finding #129; v015)", sprintf("exceeds the MDC by a wide margin: %.1fx even for Fraser, up to %.1fx for Hickman PLAC", nv("Fraser", "tdiEmp") / 0.03, nv("Hickman_PLAC", "tdiEmp") / 0.03)
 "S1 SEM table Fraser",    "Supplement A", "SEM agreement (Finding #128; v015)", "| Fraser | " + f2(sar("Fraser", "obsPredMed"), 1) + "x (" + f2(sar("Fraser", "obsPredMin"), 1) + "-" + f2(sar("Fraser", "obsPredMax"), 1) + "x) |"
 "S2 SEM table PLAC",      "Supplement A", "SEM agreement (Finding #128; v015)", "| Hickman PLAC | " + f2(sar("Hickman PLAC", "obsPredMed"), 1) + "x (" + f2(sar("Hickman PLAC", "obsPredMin"), 1) + "-" + f2(sar("Hickman PLAC", "obsPredMax"), 1) + "x) |"
 "S3 SEM closing ratios",  "Supplement A", "Fails everywhere it could meaningfully pass", "closes most of the gap for five of seven datasets (" + f2(min(SA.lcPredMed(~ismember(SA.dataset, ["Fraser" "Dhieb"]))), 1) + "-" + f2(max(SA.lcPredMed(~ismember(SA.dataset, ["Fraser" "Dhieb"]))), 1) + "x); **Fraser does not join that band " + char(8212) + " " + f2(sar("Fraser", "lcPredMed"), 2) + "x"
 "S4 SEM Part 3",          "Part 3", "SEM agreement (Finding #128) fails everywhere", "(" + f2(min(SA.lcPredMed(~ismember(SA.dataset, ["Fraser" "Dhieb"]))), 1) + "-" + f2(max(SA.lcPredMed(~ismember(SA.dataset, ["Fraser" "Dhieb"]))), 1) + "x), except Fraser at " + f2(sar("Fraser", "lcPredMed"), 2) + "x"
 "B1 BA slopes Hickman",   "Supplement A", "Fig 10 (v015) shows the difference", "(OLS slope " + sg(pb("Hickman_PLAC", "slope"), 2) + " for PLAC and " + sg(pb("Hickman_HALO", "slope"), 2) + " for HALO, Spearman " + sp2(pb("Hickman_PLAC", "spearman")) + " for both)"
 "B2 BA slopes Cook",      "Supplement A", "Fig 10 (v015) shows the difference", "in Cook ASD (" + sg(pb("Cook_ASD", "slope"), 2) + "; " + sp2(pb("Cook_ASD", "spearman")) + ") and, more moderately, in Cook CTRL (" + sg(pb("Cook_CTRL", "slope"), 2) + "; " + sp2(pb("Cook_CTRL", "spearman")) + ")"
 "B3 BA slopes flat",      "Supplement A", "Fig 10 (v015) shows the difference", "not in Fraser (" + sg(pb("Fraser", "slope"), 2) + "; " + sp2(pb("Fraser", "spearman")) + ") or Dhieb (" + sg(pb("Dhieb", "slope"), 2) + "; " + sp2(pb("Dhieb", "spearman")) + ")"
 "V1 42-cell N200",        "Part 1", "Across all 42 cells, 36 are SEM-adequate", "coverage classification is " + V42.all200(1) + " PASS / " + V42.all200(2) + " CONDITIONAL / " + V42.all200(3) + " FAIL at N_REPS = 200 (v015; " + v42(V42.all20) + " at the first-reported N_REPS = 20)"
 "V2 42-cell domain",      "Part 1", "Across all 42 cells, 36 are SEM-adequate", "the classification is " + v42(V42.dom200) + " (" + strjoin(string(V42.dom20), "/") + " at N_REPS = 20)"
 "V3 42-cell Part 2",      "Part 2", "Coverage classification result (current", "32 of 42 empirical dataset-by-pipeline combinations FAIL outright; " + V42.all200(2) + " are CONDITIONAL; " + V42.all200(1) + " PASS"
 "V4 42-cell Part 3",      "Part 3", "all 30 are SEM-adequate; coverage is", "coverage is " + V42.dom200(1) + " PASS, " + V42.dom200(2) + " CONDITIONAL, " + V42.dom200(3) + " FAIL at N_REPS = 200 (" + strjoin(string(V42.dom20), ", ") + " at N_REPS = 20)"
 "T6a-1 beta_obs BWFD-OLS", "Part 4", "Table 6a — β_obs, all six pipelines", t6row(1)
 "T6a-2 beta_obs SG-OLS", "Part 4", "Table 6a — β_obs, all six pipelines", t6row(2)
 "T6a-3 beta_obs BWFD-IRLS", "Part 4", "Table 6a — β_obs, all six pipelines", t6row(3)
 "T6a-4 beta_obs SG-IRLS", "Part 4", "Table 6a — β_obs, all six pipelines", t6row(4)
 "T6a-5 beta_obs BWFD-LMLS", "Part 4", "Table 6a — β_obs, all six pipelines", t6row(5)
 "T6a-6 beta_obs SG-LMLS", "Part 4", "Table 6a — β_obs, all six pipelines", t6row(6)
 "NC BWFD-OLS", "Part 4", "Table 6a's companion", ncrow("BWFD-OLS")
 "NC SG-OLS", "Part 4", "Table 6a's companion", ncrow("SG-OLS")
 "NC BWFD-LMLS", "Part 4", "Table 6a's companion", ncrow("BWFD-LMLS")
 "NC SG-LMLS", "Part 4", "Table 6a's companion", ncrow("SG-LMLS")
 "NC BWFD-IRLS", "Part 4", "Table 6a's companion", ncrow("BWFD-IRLS")
 "NC SG-IRLS", "Part 4", "Table 6a's companion", ncrow("SG-IRLS")
 "H1 PLAC per-pipeline",  "Part 3", "Per-pipeline heterogeneity (Finding #113)", "OLS both Poor at " + d2_local(ph("Hickman_PLAC", "BWFD_OLS")) + "/" + d2_local(ph("Hickman_PLAC", "SG_OLS")) + "; LMLS and IRLS " + d2_local(min([ph("Hickman_PLAC", "BWFD_LMLS") ph("Hickman_PLAC", "SG_LMLS") ph("Hickman_PLAC", "BWFD_IRLS") ph("Hickman_PLAC", "SG_IRLS")])) + "-" + d2_local(max([ph("Hickman_PLAC", "BWFD_LMLS") ph("Hickman_PLAC", "SG_LMLS") ph("Hickman_PLAC", "BWFD_IRLS") ph("Hickman_PLAC", "SG_IRLS")]))
 "H2 OLS worst others",   "Part 3", "Per-pipeline heterogeneity (Finding #113)", "(Cook CTRL " + d2_local(ph("Cook_CTRL", "BWFD_OLS")) + "/" + d2_local(ph("Cook_CTRL", "SG_OLS")) + ", Cook ASD " + d2_local(ph("Cook_ASD", "BWFD_OLS")) + "/" + d2_local(ph("Cook_ASD", "SG_OLS")) + ", Hickman HALO " + d2_local(ph("Hickman_HALO", "BWFD_OLS")) + "/" + d2_local(ph("Hickman_HALO", "SG_OLS"))
 "H3 Fraser reversal",    "Part 3", "Per-pipeline heterogeneity (Finding #113)", "OLS is best (" + d2_local(ph("Fraser", "BWFD_OLS")) + "/" + d2_local(ph("Fraser", "SG_OLS")) + ") and LMLS worst (" + d2_local(ph("Fraser", "BWFD_LMLS")) + "/" + d2_local(ph("Fraser", "SG_LMLS")) + ")"
 "T6b Fraser", "Part 4", "Table 6b — validity gradient", b6row("Fraser", "Fraser")
 "T6b Cook ASD", "Part 4", "Table 6b — validity gradient", b6row("Cook ASD", "Cook ASD")
 "T6b Cook CTRL", "Part 4", "Table 6b — validity gradient", b6row("Cook CTRL", "Cook CTRL")
 "T6b Hickman HALO", "Part 4", "Table 6b — validity gradient", b6row("Hickman HALO", "Hickman HALO")
 "T6b Hickman PLAC", "Part 4", "Table 6b — validity gradient", b6row("Hickman PLAC", "Hickman PLAC")
 "RSA Fraser", "Supplement A", "Cross-trial RSA (Finding #125) — the literal", rsrow("Fraser", "Fraser")
 "RSA Zarandi", "Supplement A", "Cross-trial RSA (Finding #125) — the literal", rsrow("Zarandi", "Zarandi")
 "RSA Cook CTRL", "Supplement A", "Cross-trial RSA (Finding #125) — the literal", rsrow("Cook CTRL", "Cook CTRL")
 "RSA Cook ASD", "Supplement A", "Cross-trial RSA (Finding #125) — the literal", rsrow("Cook ASD", "Cook ASD")
 "RSA Hickman PLAC", "Supplement A", "Cross-trial RSA (Finding #125) — the literal", rsrow("Hickman PLAC", "Hickman PLAC")
 "RSA Hickman HALO", "Supplement A", "Cross-trial RSA (Finding #125) — the literal", rsrow("Hickman HALO", "Hickman HALO")
 "RSA pooled", "Supplement A", "Cross-trial RSA (Finding #125) — the literal", "| Pooled (6, Dhieb excl.) | " + sprintf("%d", RS.pooled.nTrials) + " | " + sgn(RS.pooled.mantelR) + " | " + pf(RS.pooled.pValue) + " |"
 "N12 pattern ranks",      "Supplement A", "Pattern preservation (Finding #127; v015)", "in both Hickman arms the median ρ is " + dot_local(nv("Hickman_PLAC", "rhoMed")) + " against " + dot_local(nv("Hickman_PLAC", "nullRhoMed"))

 "ER1 tempo vetted", "Supplement C", "*Tempo.* The within-session", sprintf("is %.1f%% (95%% CI %.1f to %.1f%%) against %.1f%% (%.1f to %.1f%%) with the original template", 100*ee("SG-IRLS","redV1","est19"), 100*ee("SG-IRLS","redV1","lo19"), 100*ee("SG-IRLS","redV1","hi19"), 100*ee("SG-IRLS","redV1","est15"), 100*ee("SG-IRLS","redV1","lo15"), 100*ee("SG-IRLS","redV1","hi15"))
 "ER2 tempo slope", "Supplement C", "*Tempo.* The within-session", sprintf("the β_gen* slope is %+.4f against %+.4f (CI %+.3f to %+.3f", ee("SG-IRLS","bGenV1","est19"), ee("SG-IRLS","bGenV1","est15"), ee("SG-IRLS","bGenV1","lo19"), ee("SG-IRLS","bGenV1","hi19"))
 "ER3 tempo legacy", "Supplement C", "*Tempo.* The within-session", sprintf("gives %.1f%% (%.1f to %.1f%%) against %.1f%% (%.1f to %.1f%%)", 100*ee("BWFD-OLS","redV1","est19"), 100*ee("BWFD-OLS","redV1","lo19"), 100*ee("BWFD-OLS","redV1","hi19"), 100*ee("BWFD-OLS","redV1","est15"), 100*ee("BWFD-OLS","redV1","lo15"), 100*ee("BWFD-OLS","redV1","hi15"))
 "ER4 legacy slope", "Supplement C", "*Tempo.* The within-session", sprintf("(%+.4f, CI %+.3f to %+.3f)", ee("BWFD-OLS","bGenV1","est19"), ee("BWFD-OLS","bGenV1","lo19"), ee("BWFD-OLS","bGenV1","hi19"))
 "ER5 loopCCC", "Supplement C", "*Validation.* Per-dataset cccMed", "(Fraser " + dot_local(ec("Fraser","loopCCC15")) + " to " + dot_local(ec("Fraser","loopCCC19")) + ", Cook CTRL " + dot_local(ec("Cook_CTRL","loopCCC15")) + " to " + dot_local(ec("Cook_CTRL","loopCCC19")) + ", Cook ASD " + dot_local(ec("Cook_ASD","loopCCC15")) + " to " + dot_local(ec("Cook_ASD","loopCCC19")) + ", Hickman PLAC " + dot_local(ec("Hickman_PLAC","loopCCC15")) + " to " + dot_local(ec("Hickman_PLAC","loopCCC19")) + ", Hickman HALO " + dot_local(ec("Hickman_HALO","loopCCC15")) + " to " + dot_local(ec("Hickman_HALO","loopCCC19")) + ")"
 "ER6 null range", "Supplement C", "*Validation.* Per-dataset cccMed", strrep(sprintf("(cccMed %.2f to %.2f against a comparator of %.2f to %.2f)", ec("Hickman_PLAC","ccc19"), ec("Hickman_HALO","ccc19"), ec("Hickman_PLAC","null19"), ec("Hickman_HALO","null19")), "0.", ".")
 "ER7 C1 Cook CTRL", "Supplement C", "Table C1.", c1row("Cook_CTRL", "Cook CTRL")
 "ER8 C1 Cook ASD", "Supplement C", "Table C1.", c1row("Cook_ASD", "Cook ASD")
 "ER9 C1 Hickman PLAC", "Supplement C", "Table C1.", c1row("Hickman_PLAC", "Hickman PLAC")
 "ER10 C1 Hickman HALO", "Supplement C", "Table C1.", c1row("Hickman_HALO", "Hickman HALO")
 "ER11 C1 Fraser", "Supplement C", "Table C1.", c1row("Fraser", "Fraser")
 "ER12 Part 3 tempo sentence", "Part 3", "Rerunning the simulations with the end-anchored template", sprintf("%.0f%% (CI %.0f to %.0f%%) in the vetted profile and %.0f%% (CI %.0f to %.0f%%) in the legacy one", 100*ee("SG-IRLS","redV1","est19"), 100*ee("SG-IRLS","redV1","lo19"), 100*ee("SG-IRLS","redV1","hi19"), 100*ee("BWFD-OLS","redV1","est19"), 100*ee("BWFD-OLS","redV1","lo19"), 100*ee("BWFD-OLS","redV1","hi19"))
 "ER13 compound criterion", "Part 4", "With the end-anchored template of Supplement C the same two pass:", dot_local(ec("Fraser","ccc19")) + " and " + dot_local(ec("Cook_CTRL","ccc19"))
 "NM1 Fraser alpha sigma v015", "Part 2", "whose dataset medians differ slightly", sprintf("(Fraser α %.2f, σ %.2f mm in the v015 corpus)", nmv("Fraser","v015","alphaMed"), nmv("Fraser","v015","sigmaMed"))
 "NM2 alpha by basis", "Part 2", "Script: `printDatasetNoiseMedians_v001.m`", sprintf("for example Cook CTRL %.2f against %.2f", nmv("Cook_CTRL","N20","alphaMed"), nmv("Cook_CTRL","v015","alphaMed"))
 "NM3 alpha fields identical", "Part 2", "Script: `printDatasetNoiseMedians_v001.m`", iff_local(all(AP.maxDiff_alphaMaj == 0 & AP.maxDiff_alphaMin == 0 & AP.maxDiff_v015_vs_char == 0), "identical in both runs", "NOT IDENTICAL")
 "C2.1 medians v015 to v019", "Part 3", "*Table C3. As Table C2", "| Median β_gen*: Fraser | " + c2m(U15, "Fraser") + " | " + c2m(U19, "Fraser") + " |"
 "C2.1 medians CTRL", "Part 3", "*Table C3. As Table C2", "| Median β_gen*: Cook CTRL | " + c2m(U15, "Cook_CTRL") + " | " + c2m(U19, "Cook_CTRL") + " |"
 "C2.1 medians ASD", "Part 3", "*Table C3. As Table C2", "| Median β_gen*: Cook ASD | " + c2m(U15, "Cook_ASD") + " | " + c2m(U19, "Cook_ASD") + " |"
 "C2.1 medians PLAC", "Part 3", "*Table C3. As Table C2", "| Median β_gen*: Hickman PLAC | " + c2m(U15, "Hickman_PLAC") + " | " + c2m(U19, "Hickman_PLAC") + " |"
 "C2.1 medians HALO", "Part 3", "*Table C3. As Table C2", "| Median β_gen*: Hickman HALO | " + c2m(U15, "Hickman_HALO") + " | " + c2m(U19, "Hickman_HALO") + " |"
 "C2.2 Fraser minus Cook CTRL", "Part 3", "*Table C3. As Table C2", "| Fraser minus Cook CTRL | " + c2d(U15, "Fraser - Cook_CTRL", "constellation") + " | " + c2d(U19, "Fraser - Cook_CTRL", "constellation") + " |"
 "C2.3 Fraser minus Cook ASD", "Part 3", "*Table C3. As Table C2", "| Fraser minus Cook ASD | " + c2d(U15, "Fraser - Cook_ASD", "constellation") + " | " + c2d(U19, "Fraser - Cook_ASD", "constellation") + " |"
 "C2.4 Fraser minus PLAC", "Part 3", "*Table C3. As Table C2", "| Fraser minus Hickman PLAC | " + c2d(U15, "Fraser - Hickman_PLAC", "constellation") + " | " + c2d(U19, "Fraser - Hickman_PLAC", "constellation") + " |"
 "C2.5 Fraser minus HALO", "Part 3", "*Table C3. As Table C2", "| Fraser minus Hickman HALO | " + c2d(U15, "Fraser - Hickman_HALO", "constellation") + " | " + c2d(U19, "Fraser - Hickman_HALO", "constellation") + " |"
 "C2.6 SG-IRLS PLAC", "Part 3", "*Table C3. As Table C2", "| SG-IRLS only: Fraser minus Hickman PLAC | " + c2d(U15, "Fraser - Hickman_PLAC", "SG-IRLS") + " | " + c2d(U19, "Fraser - Hickman_PLAC", "SG-IRLS") + " |"
 "C2.7 SG-IRLS HALO", "Part 3", "*Table C3. As Table C2", "| SG-IRLS only: Fraser minus Hickman HALO | " + c2d(U15, "Fraser - Hickman_HALO", "SG-IRLS") + " | " + c2d(U19, "Fraser - Hickman_HALO", "SG-IRLS") + " |"
 "C2.8 pool RE", "Part 3", "*Table C2. Dataset medians, Fraser separation", "| Pool (k = 3), random effects | " + c2p(K15, "ciRE", "mPoolRE") + " | " + c2p(K19, "ciRE", "mPoolRE") + " |"
 "C2.9 pool HKSJ", "Part 3", "*Table C2. Dataset medians, Fraser separation", "| Pool (k = 3), HKSJ interval | " + c2h(K15) + " | " + c2h(K19) + " |"
 "C2c.med Fraser", "Part 3", "*Table C2. Dataset medians, Fraser separation", "| Median of subject medians, β_gen*: Fraser | " + c2m(V15, "Fraser") + " | " + c2m(V19, "Fraser") + " |"
 "C2c.med CTRL", "Part 3", "*Table C2. Dataset medians, Fraser separation", "| Median of subject medians, β_gen*: Cook CTRL | " + c2m(V15, "Cook_CTRL") + " | " + c2m(V19, "Cook_CTRL") + " |"
 "C2c.med ASD", "Part 3", "*Table C2. Dataset medians, Fraser separation", "| Median of subject medians, β_gen*: Cook ASD | " + c2m(V15, "Cook_ASD") + " | " + c2m(V19, "Cook_ASD") + " |"
 "C2c.med PLAC", "Part 3", "*Table C2. Dataset medians, Fraser separation", "| Median of subject medians, β_gen*: Hickman PLAC | " + c2m(V15, "Hickman_PLAC") + " | " + c2m(V19, "Hickman_PLAC") + " |"
 "C2c.med HALO", "Part 3", "*Table C2. Dataset medians, Fraser separation", "| Median of subject medians, β_gen*: Hickman HALO | " + c2m(V15, "Hickman_HALO") + " | " + c2m(V19, "Hickman_HALO") + " |"
 "C2c.cons CTRL", "Part 3", "*Table C2. Dataset medians, Fraser separation", "| Fraser minus Cook CTRL | " + c2d(V15, "Fraser - Cook_CTRL", "constellation") + " | " + c2d(V19, "Fraser - Cook_CTRL", "constellation") + " |"
 "C2c.SG-I CTRL", "Part 3", "*Table C2. Dataset medians, Fraser separation", "| SG-IRLS only: Fraser minus Cook CTRL | " + c2d(V15, "Fraser - Cook_CTRL", "SG-IRLS") + " | " + c2d(V19, "Fraser - Cook_CTRL", "SG-IRLS") + " |"
 "C2c.cons ASD", "Part 3", "*Table C2. Dataset medians, Fraser separation", "| Fraser minus Cook ASD | " + c2d(V15, "Fraser - Cook_ASD", "constellation") + " | " + c2d(V19, "Fraser - Cook_ASD", "constellation") + " |"
 "C2c.SG-I ASD", "Part 3", "*Table C2. Dataset medians, Fraser separation", "| SG-IRLS only: Fraser minus Cook ASD | " + c2d(V15, "Fraser - Cook_ASD", "SG-IRLS") + " | " + c2d(V19, "Fraser - Cook_ASD", "SG-IRLS") + " |"
 "C2c.cons PLAC", "Part 3", "*Table C2. Dataset medians, Fraser separation", "| Fraser minus Hickman PLAC | " + c2d(V15, "Fraser - Hickman_PLAC", "constellation") + " | " + c2d(V19, "Fraser - Hickman_PLAC", "constellation") + " |"
 "C2c.SG-I PLAC", "Part 3", "*Table C2. Dataset medians, Fraser separation", "| SG-IRLS only: Fraser minus Hickman PLAC | " + c2d(V15, "Fraser - Hickman_PLAC", "SG-IRLS") + " | " + c2d(V19, "Fraser - Hickman_PLAC", "SG-IRLS") + " |"
 "C2c.cons HALO", "Part 3", "*Table C2. Dataset medians, Fraser separation", "| Fraser minus Hickman HALO | " + c2d(V15, "Fraser - Hickman_HALO", "constellation") + " | " + c2d(V19, "Fraser - Hickman_HALO", "constellation") + " |"
 "C2c.SG-I HALO", "Part 3", "*Table C2. Dataset medians, Fraser separation", "| SG-IRLS only: Fraser minus Hickman HALO | " + c2d(V15, "Fraser - Hickman_HALO", "SG-IRLS") + " | " + c2d(V19, "Fraser - Hickman_HALO", "SG-IRLS") + " |"
 "P4.ii CTRL", "Part 4", "(ii) On each dataset's own subject-cluster bootstrap", "Cook CTRL " + c3d(V15, "Fraser - Cook_CTRL", "constellation")
 "P4.ii ASD", "Part 4", "(ii) On each dataset's own subject-cluster bootstrap", "Cook ASD " + c3d(V15, "Fraser - Cook_ASD", "constellation")
 "P4.ii PLAC", "Part 4", "(ii) On each dataset's own subject-cluster bootstrap", "Hickman PLAC " + c3d(V15, "Fraser - Hickman_PLAC", "constellation")
 "P4.ii HALO", "Part 4", "(ii) On each dataset's own subject-cluster bootstrap", "Hickman HALO " + c3d(V15, "Fraser - Hickman_HALO", "constellation")
 "P4.ii SGIRLS CTRL", "Part 4", "(ii) On each dataset's own subject-cluster bootstrap", "Cook CTRL " + c3d(V15, "Fraser - Cook_CTRL", "SG-IRLS")
 "P4.ii SGIRLS ASD", "Part 4", "(ii) On each dataset's own subject-cluster bootstrap", "Cook ASD " + c3d(V15, "Fraser - Cook_ASD", "SG-IRLS")
 "P4.ii SGIRLS HALO", "Part 4", "(ii) On each dataset's own subject-cluster bootstrap", "Hickman HALO " + c3d(V15, "Fraser - Hickman_HALO", "SG-IRLS")
 "P4.ii SGIRLS PLAC", "Part 4", "(ii) On each dataset's own subject-cluster bootstrap", "(" + c3d(V15, "Fraser - Hickman_PLAC", "SG-IRLS") + ")"
 "ABS.1 Fraser vs others range", "Abstract", "Fraser's corpus sat above the others", "(median of subject medians " + f3(V15.rowsM.median(V15.rowsM.dataset == "Fraser" & V15.rowsM.estimator == "constellation"), 3) + " against " + f3(min(V15.rowsM.median(V15.rowsM.dataset ~= "Fraser" & V15.rowsM.estimator == "constellation")), 3) + "-" + f3(max(V15.rowsM.median(V15.rowsM.dataset ~= "Fraser" & V15.rowsM.estimator == "constellation")), 3) + ")"
 "ABS.2 Fraser vs others Part 4", "Part 4", "Fraser's estimate is the highest", "(median of subject medians " + f3(V15.rowsM.median(V15.rowsM.dataset == "Fraser" & V15.rowsM.estimator == "constellation"), 3) + " against " + f3(min(V15.rowsM.median(V15.rowsM.dataset ~= "Fraser" & V15.rowsM.estimator == "constellation")), 3) + "-" + f3(max(V15.rowsM.median(V15.rowsM.dataset ~= "Fraser" & V15.rowsM.estimator == "constellation")), 3) + ")"
 "ABS.3 provenance range", "Abstract", "*Provenance:", f3(min(V15.cMed)) + "-" + f3(max(V15.cMed)) + " as medians of subject medians"
 "ABS.4 Cook cohorts", "Discussion", "The two cohorts' close individual agreement", f3(V15.rowsM.median(V15.rowsM.dataset == "Cook_CTRL" & V15.rowsM.estimator == "constellation"), 3) + " vs " + f3(V15.rowsM.median(V15.rowsM.dataset == "Cook_ASD" & V15.rowsM.estimator == "constellation"), 3) + " as medians of subject medians on v015"
 "ABS.5 comparator beats model outside Fraser", "Abstract", "scored higher than the model on it in the four datasets other than Fraser", iff_local(all(arrayfun(@(d) nv(d, "nullCccMed") > nv(d, "cccMed"), ["Cook_CTRL" "Cook_ASD" "Hickman_PLAC" "Hickman_HALO"])) && nv("Fraser", "cccMed") > nv("Fraser", "nullCccMed"), "scored higher than the model on it in the four datasets other than Fraser", "COMPARATOR CLAIM NOT SUPPORTED")
 "GD6.1 VGF slope", "Part 1", "Two coupled mechanisms can produce this", sprintf("(d ln L / dβ = −%.2f", -VC.slope)
 "GD6.2 VGF inflation", "Part 1", "Two coupled mechanisms can produce this", sprintf("a deficit of 0.05 in β alone inflates the recovered VGF by %.0f%% and a deficit of 0.10 by %.0f%%", 100 * VC.inflFirstOrder(1), 100 * VC.inflFirstOrder(2))
 "B2 AXdF", "Supplement B", "*Table B2. Paired difference", pbrow("Whole noise replaced", "AXdF")
 "B2 AXoF", "Supplement B", "*Table B2. Paired difference", pbrow("Structure replaced, size kept", "AXoF")
 "B2 AXlF", "Supplement B", "*Table B2. Paired difference", pbrow("Size replaced, structure kept", "AXlF")
 "B2 AXnF", "Supplement B", "*Table B2. Paired difference", pbrow("White noise", "AXnF")
 "B2 AXcF", "Supplement B", "*Table B2. Paired difference", pbrow("Colour only", "AXcF")
 "B2 AXaF", "Supplement B", "*Table B2. Paired difference", pbrow("Colour and relation between axes", "AXaF")
 "B2 AXsF", "Supplement B", "*Table B2. Paired difference", pbrow("Relation between axes broken", "AXsF")
 "B2 SE ratio range", "Supplement B", "*Paired intervals (post hoc).*", sprintf("The paired standard errors are %.2f to %.1f times the model-based ones, except for the white-noise replacement in Fraser, where they are %.1f (vetted) and %.1f (legacy) times as large", min(PBn.seRatio), max(PBn.seRatio), PBw.seRatio(PBw.pipeline == "SG-IRLS"), PBw.seRatio(PBw.pipeline == "BWFD-OLS"))
 "B2 agreement", "Supplement B", "*Paired intervals (post hoc).*", sprintf("In Fraser its differences reproduce the contrasts above to within %.3f (in Hickman PLAC, where the paired samples are %d to %d trials, to within %.3f)", -floor(-1000 * max(abs(PBf.diff - PBf.pubDiff))) / 1000, min(PBp.nTrials), max(PBp.nTrials), -floor(-1000 * max(abs(PBp.diff - PBp.pubDiff))) / 1000)
 "B2 white vetted slopes", "Supplement B", "*Paired intervals (post hoc).*", sprintf("(the slope halves, from %+.3f to %+.3f in this estimator)", PBw.slopeRef(PBw.pipeline == "SG-IRLS"), PBw.slopeArm(PBw.pipeline == "SG-IRLS"))
 "B2 reference slopes", "Supplement B", "*Paired intervals (post hoc).*", sprintf("%+.3f to %+.3f (Fraser, vetted), %+.3f to %+.3f (Hickman PLAC, vetted) and %+.3f to %+.3f (Fraser, legacy)", min(PBf.slopeRef(PBf.pipeline == "SG-IRLS")), max(PBf.slopeRef(PBf.pipeline == "SG-IRLS")), min(PBp.slopeRef), max(PBp.slopeRef), min(PBf.slopeRef(PBf.pipeline == "BWFD-OLS")), max(PBf.slopeRef(PBf.pipeline == "BWFD-OLS")))
 "D.1 SG frames", "Supplement D", "Savitzky-Golay differentiation (SG).", sprintf("The resulting frame lengths are %d, %d, %d, %d and %d samples at 60, 100, 120, 133 and 240 Hz (%.3f, %.3f, %.3f, %.3f and %.3f s)", PS.sg.frameLen, PS.sg.frameSec)
 "D.2 BW response", "Supplement D", "Butterworth differentiation (BWFD).", sprintf("the forward-backward pass has 1/%d there (%.2f dB at 10 Hz", round(1 / PS.bw.zeroLagPower(1)), PS.bw.zeroLagDB(1))
 "D.3 floor", "Supplement D", "Edge clip and length floor.", sprintf("shorter than 4 × (%d + 50) = %d samples", PS.edgeClip, PS.floorSamples)
 "D.4 clip", "Supplement D", "Edge clip and length floor.", sprintf("the clip is %d samples at each end of the differentiated signal", PS.edgeClip)
 "D.5 SG order", "Supplement D", "Savitzky-Golay differentiation (SG).", sprintf("Polynomial order %d. The reference frame is %d samples at 100 Hz", PS.sgOrder, PS.sgRef)
 "D.6 BW order", "Supplement D", "Butterworth differentiation (BWFD).", sprintf("A second-order Butterworth low-pass filter (`butter(%d, fc/(fs/2))`, fc = %d Hz)", PS.bwOrder, PS.bwFc)
 "D.7 seeds", "Supplement D", "Regressions.", "seeded with the OLS estimate from the same trial and filter",
 "D.8 grid", "Supplement D", "Inversion of the observed exponent.", sprintf("25 equally spaced generative exponents on [0, 0.75] (spacing %.5f)", PS.grid.spacing),
 "D.9 smoothing", "Supplement D", "Inversion of the observed exponent.", sprintf("after a %d-point moving mean (a window of 0.0625", PS.grid.windowPts),
 "D.10 min run", "Supplement D", "Inversion of the observed exponent.", sprintf("at least 0.0625 (%d grid spacings)", PS.grid.minSegSpacings),
 "D.11 slope trim", "Supplement D", "Inversion of the observed exponent.", "mean absolute slope above 0.05, and it ends at the first step that falls against its direction by at least 1.0 times the mean noise scale",
 "D.12 ambiguous", "Supplement D", "Inversion of the observed exponent.", "it lies inside more than one run (rising or falling), the trial is classed as ambiguous",
 "D.13 grid separate", "Supplement D", "Inversion of the observed exponent.", "separate from the 22-value sweep of §2.1 (β_gen up to 0.70)",
 "D.14 D5 grid", "Deviations", "| D5 |", "25 nodes over β_gen ∈ [0, 0.75] (the registered grid stops at 2/3",
 "D.15 GM cal", "Supplement D", "Generative-model check (λ, ρ, λ*; D20).", sprintf("(medA %.3f)", GM.medA),
 "D.16 GM windows", "Supplement D", "Generative-model check (λ, ρ, λ*; D20).", sprintf("the enacted end read %.3f with sensor noise, 0.006 above its window, and %.3f without", GM.enactedSensor, GM.enactedClean),
 "D.17 GM windows set", "Supplement D", "Generative-model check (λ, ρ, λ*; D20).", "|medA| ≤ 0.30 and 0.70 ≤ medR ≤ 1.30",
 "D.18 GM NSurr", "Supplement D", "Generative-model check (λ, ρ, λ*; D20).", sprintf("%d surrogates each in the full-corpus run (seed 1729)", GM.nSurr),
 "D.19 GM template", "Supplement D", "Generative-model check (λ, ρ, λ*; D20).", "the first four harmonics of the trial's fundamental frequency f0, fitted by Hann-weighted least squares in windows of four cycles with 50% overlap",
 "D.20 GM Dhieb", "Part 3", "Dhieb's raw λ", sprintf("(%.2f) was inconclusive", GMd("Dhieb").subjMed_lambda),
 "D.21 GM Dhieb star", "Part 3", "On the calibrated scale below its", sprintf("λ* (%.2f) exceeds the 0.5 boundary", GMd("Dhieb").subjMed_lambdaStar_APPROX),
 "D.22 MC met", "Part 3", "Controls matched to each dataset", sprintf("met the same pre-set windows for %s and %s", strjoin(MC.dataset(find(MC.calOK, 2, "first")), ", "), MC.dataset(find(MC.calOK, 1, "last"))),
 "D.22b MC CTRL", "Part 3", "Cook CTRL, Hickman PLAC and Hickman HALO", sprintf("(additive %.2f, %.2f, %.2f; enacted %.2f, %.2f, %.2f)", MCd("Cook CTRL").medA, MCd("Hickman PLAC").medA, MCd("Hickman HALO").medA, MCd("Cook CTRL").medR, MCd("Hickman PLAC").medR, MCd("Hickman HALO").medR),
 "D.23 MC near", "Part 3", "narrowly missed", sprintf("them for Cook ASD (enacted %.2f) and Dhieb (%.2f), overshot for Fraser (%.2f) and failed for Zarandi (−%.2f)", MCd("Cook ASD").medR, MCd("Dhieb").medR, MCd("Fraser").medR, -MCd("Zarandi").medR),
 "D.24 MC lamstar", "Part 3", "On each dataset's own control", sprintf("λ* is %.2f in Fraser, −%.2f to −%.2f in the four clinical arms and %.2f in Dhieb", MCd("Fraser").lambdaStar_own, -min(MC.lambdaStar_own(ismember(MC.dataset, ["Cook CTRL","Cook ASD","Hickman PLAC","Hickman HALO"]))), -max(MC.lambdaStar_own(ismember(MC.dataset, ["Cook CTRL","Cook ASD","Hickman PLAC","Hickman HALO"]))), MCd("Dhieb").lambdaStar_own),
 "D.25 MC diag", "Part 3", "a diagnosis on 25 Zarandi", sprintf("trials read −%.2f with the matched parameters, %.2f with only the colour", -DGm("matched"), DGm("alpha3.18")),
 "D.26 MC diag all", "Part 3", "the minor) and", sprintf("%.2f with all four original", DGm("allOriginal")),
 "D.27 raw split", "Part 3", "the raw median λ is at", sprintf("most %.2f in the five in-domain datasets and %.2f in Dhieb and %.2f in", inDomRaw, MCd("Dhieb").subjMed_lambda, MCd("Zarandi").subjMed_lambda),
 "D.28 MC colour", "Part 3", "fitted colour is", sprintf("%.1f on the major axis and %.1f on the minor", median(DG.P(:,5)), median(DG.P(:,7))),
 "D.29 floor plateau", "Supplement D", "The sensor floor is bracketed, not measured", sprintf("with a plateau only in %s", strjoin(SF.dataset(SF.plateau), " and ")),
 "D.35 boot SD", "Part 4 methods", "The inversion-uncertainty draw is Gaussian", sprintf("divided by 2 × %.3f", BS.Z90),
 "D.36 boot inner", "Part 4 methods", "the SD of the median of", sprintf("%d inner draws", BS.NInner),
 "D.37 boot floor", "Part 4 methods", "floored at", "1.2533 times the between-pipeline SD",
 "D.38 boot n", "Part 4 methods", "Each stage uses", sprintf("%d resamples (seed %d)", BS.NBoot, BS.RngSeed),
 "D.39 boot subj-only", "Part 4 methods", "separate bootstrap that resamples subjects only", sprintf("(%d resamples, seed %d)", BS.datasetMedianB, BS.datasetMedianSeed),
 "D.40 boot median", "Part 4 methods", "across-subject median"" is", "the median of per-subject medians",
 "D.41 length floor", "Supplement D", "Edge clip and length floor", sprintf("4 × (%d + 50) = %d samples", SS.edgeClip, SS.mFloor),
 "D.42 vgf grid", "Supplement D", "velocity-gain-factor grid", sprintf("%.1f to %.1f in steps of %.1f in the natural logarithm", SS.vgfLo, SS.vgfHi, SS.vgfLnStep),
 "D.43 template cols", "Supplement D", "Harmonic template and noise residual", sprintf("%d columns", SS.nTemplateCols),
 "D.44 ira band", "Supplement D", "The residual's exponent α is estimated by IRASA", sprintf("%g Hz to min(%g Hz, fs/2 − 1) and dilation factors %.1f to %.1f in steps of %.2f", SS.iraFLow, SS.iraFHighMax, SS.iraH(1), SS.iraH(end), SS.iraH(2)-SS.iraH(1)),
 "D.45 harness clip", "Supplement D", "Edge clip and length floor", sprintf("the edge clip is %d samples at every sampling rate", SS.harnessEdgeClip),
 "D.46 harness guard", "Supplement D", "Edge clip and length floor", sprintf("fewer than %d valid (finite, positive-curvature) points", SS.harnessMinValid),
 "D.47 harness fail N", "Supplement D", "Edge clip and length floor", sprintf("N ≤ %d samples therefore fails", SS.harnessFailNmax),
 "D.48 lmm count", "Supplement D", "LMM coding (Stage 1)", sprintf("gives %d coefficients: the intercept, %d main effects and %d interactions", LC.nCoef, LC.nMain, LC.nInter),
 "D.49 lmm groups", "Supplement D", "LMM coding (Stage 1)", sprintf("The random intercept has %s groups", regexprep(sprintf("%d", LC.nGroups), "\d(?=(\d{3})+$)", "$0,")),
 "D.50 lmm signif", "Supplement D", "LMM coding (Stage 1)", sprintf("for all %d main effects and %d of the %d interactions", LC.mainSig, LC.interSig, LC.nInter),
 "D.51 lmm scaled", "Supplement D", "LMM coding (Stage 1)", "sampling rate was not scaled, so its coefficients are per Hz",
 "D.52 lmm 2.4", "Part 2 methods", "× σ × α + (1|paramComboID).", "σ and α were centred and scaled to unit SD; sampling rate was left in Hz",
 "D.53 simex Hickman", "Part 3", "(4,033 SIMEX cells; 3,829 joined to a placed trial;", sprintf("the other %d are %d trials × 6 pipelines", SJ.nUnplaced(1), SJ.nUnplacedTrials(1)),
 "D.54 simex Cook", "Part 3", "Cook CTRL+ASD (n=1,174 SIMEX cells,", sprintf("%d of them on placed trials; the %d left out are %d trials × 6 pipelines", SJ.nJoined(4), SJ.nUnplaced(4), SJ.nUnplacedTrials(4)),
 "D.55 bin edges", "Part 4", "SG-IRLS medians of the reliable targets", sprintf("(edges %g, %g, %g and %g Hz", TB.edges(2:5)),
 "D.56 Cook CTRL bins", "Part 4", "SG-IRLS medians of the reliable targets", sprintf("Cook CTRL %.3f, %.3f, %.3f (%d, %d and %d trials)", TB.CTRL, TB.CTRLn),
 "D.57 high bin", "Part 4", "A fifth bin above 2.8 Hz", sprintf("%d (PLAC, %.3f), %d (HALO, %.3f) and %d (Dhieb, %.3f)", TB.HI.hickmanPLAC(2), TB.HI.hickmanPLAC(1), TB.HI.hickmanHALO(2), TB.HI.hickmanHALO(1), TB.HI.dhieb(2), TB.HI.dhieb(1)),
 "D.58 Cook raw", "Part 4", "had a higher raw exponent (β_obs", sprintf("%+.4f SG-IRLS, %+.4f BWFD-OLS", CC.rawObs),
 "D.59 Cook tempo", "Part 4", "Adding tempo alone reduces the difference to", sprintf("%+.3f (%+.3f)", CC.tempoObs),
 "D.60 Cook half-octave", "Part 4", "At matched tempo (coarsened exact matching", sprintf("β_obs differs by %+.3f (%+.3f) and β_gen* by %+.3f (%+.3f)", CC.halfOctObs, CC.halfOctGen),
 "D.61 Cook tenth", "Part 4", "At matched tempo (coarsened exact matching", sprintf("%+.3f (%+.3f) and %+.3f (%+.3f)", CC.tenthObs, CC.tenthGen),
 "D.62 GD SEs", "General Discussion", "variances (SE", sprintf("%.4f and %.4f; the bootstrap SE of their difference is %.4f)", CC.gdSeCtrl, CC.gdSeAsd, CC.gdSeDiff),
 "D.63 GD sign", "General Discussion", "variances (SE", sprintf("(%.4f, SE %.4f; ASD minus CTRL is %+.4f", CC.gdCtrlMinusAsd, CC.gdSeDiff, CC.gdAsdMinusCtrl),
 "D.64 replay noise", "Part 1B", "Across those 1,320 nodes", sprintf("by a median of %.4f and at most %.3f, with %d nodes above 0.01", RP.noiseMedian, RP.noiseMax, RP.noiseOver01),
 "D.65 P4 B2", "Part 1B", "P4 (replay within 0.01) held in substance", sprintf("a median of %.4f (at most %.4f in any pipeline) and by up to %.3f (%s), with %d of %d nodes above 0.01", RP.b2Median, RP.b2PipeMedianMax, RP.b2Max, RP.b2MaxPipeline, RP.b2Over01, RP.nB2),
 "D.66 replay zero", "Part 1B", "P4 (replay within 0.01) held in substance", sprintf("matches the grid to 5 × 10⁻⁶ at all %d σ = 0 nodes", RP.nZero),
 "D.67 cell-mean gate", "Part 1B", "OLS on cell means; applied to the registered VGF formula", sprintf("all %d fitted fixed effects to within %.3f standard errors, the largest absolute difference being %.1f × 10⁻⁶", RP.gateN, RP.gateMaxSE, RP.gateMaxAbs*1e6),
 "D.30 emp enacted", "Part 3", "unreadable: it read", sprintf("%.2f to %.2f in the five in-domain datasets, %.2f in Zarandi and −%.2f in", min(M3.medR(ismember(M3.dataset, in5))), max(M3.medR(ismember(M3.dataset, in5))), M3d("Zarandi").medR, -M3d("Dhieb").medR),
 "D.31 hybrid", "Part 3", "adding the rest after the retiming as sensor noise, gave", sprintf("%.2f to %.2f in the five in-domain datasets, %.2f in Zarandi and %.2f in Dhieb", min(M4.medR20(in5m)), max(M4.medR20(in5m)), M4d("Zarandi").medR20, M4d("Dhieb").medR20),
 "D.32 hybrid sens", "Part 3", "(splits at 10 and 35 Hz:", sprintf("%.2f to %.2f in the five)", min([M4.medR10(in5m); M4.medR35(in5m)]), max([M4.medR10(in5m); M4.medR35(in5m)])),
 "D.33 below additive", "Part 3", "the clinical arms (", sprintf("−%.2f to −%.2f) lie below their additive controls (%.2f to %.2f)", -min(M4.subjMed_lambda(clin4)), -max(M4.subjMed_lambda(clin4)), min(M4.medA(clin4)), max(M4.medA(clin4))),
 "D.34 Fraser lam", "Part 3", "Fraser's λ of", sprintf("%.2f is what partial enactment", MC.subjMed_lambda(MC.dataset == "Fraser")),
 "D.68 cfg filters", "Supplement D", "Settings the v015 result files do not record", sprintf("Butterworth order 2 with cutoff %d Hz, Savitzky-Golay order 4 with a %d-sample reference frame", CF.bwCutoff, CF.sgRef),
 "D.69 cfg clip", "Supplement D", "Settings the v015 result files do not record", sprintf("the option `EdgeClip` with default %d", CF.edgeClip),
 "D.70 cfg fs", "Supplement D", "Settings the v015 result files do not record", sprintf("Fraser %d Hz, Cook (both groups) %d Hz, Hickman (both arms) %d Hz, Zarandi %d Hz, Dhieb %d Hz", cfFs("Fraser"), cfFs("Cook_CTRL"), cfFs("Hickman_PLAC"), cfFs("Zarandi"), cfFs("Dhieb")),
 "D.71 half spread", "Part 1B", "The collapse of the Butterworth curves across", sprintf("spread over %.4f in the pipeline", HS.spreadMeasured),
 "D.72 half emul", "Part 1B", "The collapse of the Butterworth curves across", sprintf("the spread is %.4f, and with that delay removed it falls to %.4f (theory %.4f)", HS.spreadEmul, HS.spreadEmulZP, HS.spreadTheory),
 "D.73 half effect", "Part 1B", "The collapse of the Butterworth curves across", sprintf("by %.2f%%, %.2f%% and %.2f%% at cutoffs of 5, 10 and 20 Hz", HS.delayEffectPct),
 "D.74 phase script", "Section 2.5", "Stages 1 and 5 therefore differ only in Fourier phase", sprintf("`%s`", NP.generator),
 "D.75 phase formula", "Section 2.5", "Stages 1 and 5 therefore differ only in Fourier phase", "F = |FFT(resid)| · exp(i · angle(FFT(Xu_base)))",
 "D.76 flip count", "Section 3b", "the flip to the expected direction is", sprintf("It appears at %s of the %s datasets", numWords(DCs.nInDomain + 1), numWords(DCs.nDatasets + 1)),
 "D.77 in-domain count", "Section 3b", "The sign holds at every empirical centroid.", sprintf("At the %s in-domain datasets", numWords(DCs.nInDomain + 1)),
 "D.78 alpha count", "Section 3b", "That the four Cook and Hickman datasets", sprintf("%s of them above 5", numWords(DCs.nAlphaGt5 + 1)),
 "D.80 share Fraser", "Part 3", "Checked directly on the v015 corpus", sprintf("| Fraser | **%.1f%%** | %.1f%% | %.1f%% |", 100*MS.olsShare(msRow("Fraser")), 100*MS.lmlsShare(msRow("Fraser")), 100*MS.irlsShare(msRow("Fraser"))),
 "D.81 share CTRL", "Part 3", "Checked directly on the v015 corpus", sprintf("| Cook CTRL | %.1f%% | **%.1f%%** | %.1f%% |", 100*MS.olsShare(msRow("Cook CTRL")), 100*MS.lmlsShare(msRow("Cook CTRL")), 100*MS.irlsShare(msRow("Cook CTRL"))),
 "D.82 share ASD", "Part 3", "Checked directly on the v015 corpus", sprintf("| Cook ASD | %.1f%% | %.1f%% | **%.1f%%** |", 100*MS.olsShare(msRow("Cook ASD")), 100*MS.lmlsShare(msRow("Cook ASD")), 100*MS.irlsShare(msRow("Cook ASD"))),
 "D.83 share HALO", "Part 3", "Checked directly on the v015 corpus", sprintf("| Hickman HALO | %.1f%% | %.1f%% | **%.1f%%** |", 100*MS.olsShare(msRow("Hickman HALO")), 100*MS.lmlsShare(msRow("Hickman HALO")), 100*MS.irlsShare(msRow("Hickman HALO"))),
 "D.84 share PLAC", "Part 3", "Checked directly on the v015 corpus", sprintf("| Hickman PLAC | %.1f%% | %.1f%% | **%.1f%%** |", 100*MS.olsShare(msRow("Hickman PLAC")), 100*MS.lmlsShare(msRow("Hickman PLAC")), 100*MS.irlsShare(msRow("Hickman PLAC"))),
 "D.85 share Zarandi", "Part 3", "Checked directly on the v015 corpus", sprintf("| Zarandi (outside validity domain) | **%.1f%%** | %.1f%% | %.1f%% |", 100*MS.olsShare(msRow("Zarandi")), 100*MS.lmlsShare(msRow("Zarandi")), 100*MS.irlsShare(msRow("Zarandi"))),
 "D.86 fig7 outside", "Part 3", "Panel 7c has", sprintf("%s of its %s points outside the axes", regexprep(sprintf("%d", PF.nOutside7c), "(\d)(?=(\d{3})+$)", "$1,"), regexprep(sprintf("%d", PF.nPoints7c), "(\d)(?=(\d{3})+$)", "$1,")),
 "D.87 fig7 files", "Part 3", "Fig 7, on the v015 basis", "`src/plotConstellationMetrics_v003.m` reads the numbers from",
 "D.88 Dhieb node", "Section 3c", "Dhieb (0.0189; 0.0171;", sprintf("α = %.2f, σ = %.2f mm; the centroid itself is α = %.2f, σ = %.2f mm", DN.nodeAlpha, DN.nodeSigmaMM, DN.centroidAlpha, DN.centroidSigmaMM),
 "D.89 interp median", "Part 1B", "Snapping to the nearest node does not drive the result", sprintf("by a median of %.4f to %.4f in every dataset and by at most %.3f (Fraser), never by as much as the MDC (0.03)", min(NI.medAbs), max(NI.medAbs), NI.maxAbs(1)),
 "D.90 interp tol", "Part 1B", "Snapping to the nearest node does not drive the result", sprintf("%.1f%% of Fraser's and %.1f%% of Hickman HALO's trial-by-contrast cells move by more than the SEM threshold (0.0108), and none elsewhere", 100*NI.fracOverTol(1), 100*NI.fracOverTol(5)),
 "D.91 interp ccc", "Part 1B", "Snapping to the nearest node does not drive the result", sprintf("%+.4f (Fraser), %+.4f (Cook CTRL), %+.4f (Cook ASD), %+.4f (Hickman PLAC), %+.4f (Hickman HALO) and %+.4f (Dhieb)", NI.cccMedInterp - NI.cccMedNearest),
 "D.79 RSA basis", "Supplement A", "The cross-trial RSA table is also on v015", sprintf("saved in `results/%s`", rsaName),
 "CR.1 v019 range", "Abstract", "A sensitivity rerun with an end-anchored harmonic template", "lowered the estimates to " + f3(min(V19.cMed)) + "-" + f3(max(V19.cMed))
 "CR.2 below 1/3", "Abstract", "and gave model-implied estimates of", iff_local(all(V15.cMed < 1/3), "all below one-third", "BELOW-1/3 CLAIM NOT SUPPORTED")
 "CR.3 6b coverage", "Abstract", "with inversion uncertainty included", iff_local(cvB("Fraser") && cvB("Cook CTRL") && cvB("Cook ASD") && ~cvB("Hickman PLAC") && ~cvB("Hickman HALO"), "the intervals of three datasets (Fraser, Cook CTRL and Cook ASD) included one-third and those of both haloperidol-study arms did not", "TABLE 6B COVERAGE CLAIM NOT SUPPORTED")
 "CR.4 v019 separation", "Abstract", "in the end-anchored rerun not against Hickman PLAC", iff_local(zin(V19, "Fraser - Hickman_PLAC", "constellation") && zin(V19, "Fraser - Hickman_PLAC", "SG-IRLS") && zin(V19, "Fraser - Hickman_HALO", "SG-IRLS") && ~zin(V19, "Fraser - Hickman_HALO", "constellation"), "for the constellation median, nor against either Hickman arm for the vetted profile", "V019 SEPARATION CLAIM NOT SUPPORTED")
 "CR.5 subject-only coverage", "Supplement C", "On this subject-only bootstrap, which carries no inversion draw", iff_local(nnz(cvV(V15)) == 1 && vci(V15, "Fraser", "ciLo") <= 1/3 && vci(V15, "Fraser", "ciHi") >= 1/3 && nnz(cvV(V19)) == 0, "only Fraser's interval includes one-third with the original template ([" + dot_local(vci(V15, "Fraser", "ciLo")) + ", " + dot_local(vci(V15, "Fraser", "ciHi")) + "]), and none does with the anchored one (Fraser [" + dot_local(vci(V19, "Fraser", "ciLo")) + ", " + dot_local(vci(V19, "Fraser", "ciHi")) + "])", "SUBJECT-ONLY COVERAGE CLAIM NOT SUPPORTED")
 "CR.6 pool vs 1/3", "Supplement C", "On this subject-only bootstrap, which carries no inversion draw", iff_local((K19.ciRE(1) > 1/3 || K19.ciRE(2) < 1/3) && K19.ciHKSJ(1) <= 1/3 && K19.ciHKSJ(2) >= 1/3, "The pool's random-effects interval excludes one-third with the anchored template ([" + nl0(K19.ciRE(1), 3) + ", " + nl0(K19.ciRE(2), 3) + "]); its HKSJ interval does not.", "POOL CLAIM NOT SUPPORTED")
 "CR.7 6b v019 coverage", "Abstract", "lowered the estimates to", iff_local(n19 == 0 && tv("v019", "Cook ASD", "oneThirdMinusCiHi") == min(CV.oneThirdMinusCiHi(CV.template == "v019")), "and no dataset's interval then included one-third (Cook ASD's upper limit fell " + f3(tv("v019", "Cook ASD", "oneThirdMinusCiHi")) + " short)", "V019 6B COVERAGE CLAIM NOT SUPPORTED")
 "CR.8 6b v019 CIs", "Part 2", "**Per-dataset intervals and 1/3.**", iff_local(n19 == 0, "none includes 1/3: Fraser " + ci3("v019", "Fraser") + ", Cook CTRL " + ci3("v019", "Cook CTRL") + ", Cook ASD " + ci3("v019", "Cook ASD") + ", Hickman PLAC " + ci3("v019", "Hickman PLAC") + " and Hickman HALO " + ci3("v019", "Hickman HALO"), "V019 6B COVERAGE CLAIM NOT SUPPORTED")
 "CR.9 6b v019 nearest", "Supplement C", "On this subject-only bootstrap, which carries", iff_local(n19 == 0, "for no dataset with the anchored one (nearest, Cook ASD [" + dot_local(tv("v019", "Cook ASD", "ciLo")) + ", " + dot_local(tv("v019", "Cook ASD", "ciHi")) + "] and Fraser [" + dot_local(tv("v019", "Fraser", "ciLo")) + ", " + dot_local(tv("v019", "Fraser", "ciHi")) + "]", "V019 6B COVERAGE CLAIM NOT SUPPORTED")
 "SR.1 alpha per halving", "Section 2.3", "Subsampling the Fraser trials cleanly from 240 Hz", "lowers α by " + f3(OB2.mean) + " per halving of the rate (95% person-bootstrap interval " + f3(OB2.lo) + " to " + f3(OB2.hi) + "; " + addcomma(OB2.nTrials) + " trials, " + OB2.nPersons + " people)"
 "SR.2 alpha by rate", "Section 2.3", "Subsampling the Fraser trials cleanly from 240 Hz", "the median α is " + f2(qv(SRE.Lev, "alpha_240"), 2) + ", " + f2(qv(SRE.Lev, "alpha_120"), 2) + " and " + f2(qv(SRE.Lev, "alpha_60"), 2)
 "SR.3 anti-aliased", "Section 2.3", "With an anti-alias filter applied before the rate change", sprintf("the slope is %.3f (%.3f to %+.3f)", OB1.mean, OB1.lo, OB1.hi)
 "SR.4 Fraser vs Zarandi", "Section 2.3", "Rate accounts for little of the difference between", "subsampled to 120 Hz, Fraser still reads " + f2(qv(SRE.Lev, "alpha_120"), 2) + " against Zarandi's " + f2(nmv("Zarandi", "v015", "alphaMed"), 2) + " at 100 Hz"
 "SR.5 sigma by rate", "Section 2.3", "Subsampling leaves the noise amplitude essentially unchanged", iff_local(round(qv(SRE.Lev, "sigmaMM_240"), 2) == round(qv(SRE.Lev, "sigmaMM_120"), 2), sprintf("median σ is %.2f mm at 240 and 120 Hz and %.2f mm at 60 Hz, with paired per-trial ratios to 240 Hz of %.4f (%.4f to %.4f) and %.4f (%.4f to %.4f)", qv(SRE.Lev, "sigmaMM_240"), qv(SRE.Lev, "sigmaMM_60"), qv(SRE.Chg, "sigmaRatio_120"), qc(SRE.Chg, "sigmaRatio_120", "lo95"), qc(SRE.Chg, "sigmaRatio_120", "hi95"), qv(SRE.Chg, "sigmaRatio_60"), qc(SRE.Chg, "sigmaRatio_60", "lo95"), qc(SRE.Chg, "sigmaRatio_60", "hi95")), "SIGMA 240/120 DIFFER AT 2 DP")
 "SR.6 synthetic folding", "Section 2.3", "Folding of a pure power law does not by itself reproduce", sprintf("interpolated to the observed α of %.2f between generating exponents %g and %g) changes by %.3f per doubling of the rate (%.3f to %.3f), against %+.3f per doubling in the real trials", ali.aObs, SRA.A.alphaGen(ali.lo), SRA.A.alphaGen(ali.lo + 1), SRA.A.slopeAtObs, SRA.A.slopeAtObsCI(1), SRA.A.slopeAtObsCI(2), OB2.mean)
 "SR.7 independent noise", "Section 2.3", "Folding of a pure power law does not by itself reproduce", sprintf("noise generated independently at each rate by %+.3f per doubling (%.3f to %.3f)", SRA.B.slopeAtObs, SRA.B.slopeAtObsCI(1), SRA.B.slopeAtObsCI(2))
 "SR.8 120 Hz change", "Section 2.3", "The observed exponent is far more stable than α", iff_local(max(abs(c120)) < 0.0011, "At 120 Hz the median paired change in β_obs is below 0.0011 in all six pipelines", "120 HZ CHANGE NOT BELOW 0.0011") + sprintf(", inside the 0.0108 tolerance of D27, although %.0f%% to %.0f%% of individual trials", 100*min(s120), 100*max(s120))
 "SR.9 60 Hz change", "Section 2.3", "The observed exponent is far more stable than α", sprintf("At 60 Hz the Butterworth-differentiation pipelines rise by at most %.3f, whereas the Savitzky-Golay pipelines fall by %.4f to %.4f (SG-IRLS %.3f to %.3f), and %.0f%% to %.0f%% of their trials exceed the tolerance", max(cBW60), -max(cSG60), -min(cSG60), qv(SRE.Lev, "betaObs_SG_IRLS_240"), qv(SRE.Lev, "betaObs_SG_IRLS_60"), 100*min(sSG60), 100*max(sSG60)) + iff_local(all(cBW60 > 0) && all(cSG60 < 0), "", " SIGN PATTERN NOT SUPPORTED")
 "SR.10 residue excluded", "Section 2.3", "The observed exponent is far more stable than α", sprintf("excluding the %d trials with Savitzky-Golay curvature residue changes this little (SG-IRLS %.3f to %.3f)", SRE.nT - SRE.nKeep, qv(SRE.LevC, "betaObs_SG_IRLS_240"), qv(SRE.LevC, "betaObs_SG_IRLS_60"))
 "SR.11 LMM main", "Section 2.3", "In the fitted Stage 1 model (N =", "(N = " + addcomma(L9n) + sprintf(" simulated trajectories) the sampling-rate main effect on δβ is %+.2f × 10⁻⁵ per Hz (SE %.1f × 10⁻⁶;", 1e5*lt("samplingRate", "estimatePerHz"), 1e6*lt("samplingRate", "SE"))
 "SR.12 LMM span and SG", "Section 2.3", "In the fitted Stage 1 model (N =", sprintf("or %+.4f over the 180 Hz from 60 to 240 Hz. The derivation interaction is %.2f × 10⁻⁵ per Hz (SE %.1f × 10⁻⁶), so the Savitzky-Golay pipelines rise less with rate than the Butterworth ones (%+.4f against %+.4f under OLS); the regression interactions are positive (LMLS %+.2f × 10⁻⁵, IRLS %+.2f × 10⁻⁵ per Hz)", lt("samplingRate", "span60to240"), 1e5*lt("samplingRate:filterType_6", "estimatePerHz"), 1e6*lt("samplingRate:filterType_6", "SE"), im("SG-OLS"), im("BWFD-OLS"), 1e5*lt("samplingRate:regressionType_4", "estimatePerHz"), 1e5*lt("samplingRate:regressionType_5", "estimatePerHz"))
 "SR.13 LMM pipeline range", "Section 2.3", "In the fitted Stage 1 model (N =", iff_local(SRL.impl.pipeline(jMin) == "SG-OLS" && SRL.impl.pipeline(jMax) == "BWFD-LMLS" && imMin > 0, sprintf("falls further below the generated one as rate rises, by %.4f (SG-OLS) to %.4f (BWFD-LMLS) from 60 to 240 Hz", imMin, imMax), "PIPELINE RANGE NOT SUPPORTED")
 "SR.14 SEM shares", "Section 2.3", "In the fitted Stage 1 model (N =", iff_local(SRL.rows.pipeline(iMin) == "SG-IRLS" && SRL.rows.fs(iMin) == 60 && SRL.rows.pipeline(iMax) == "BWFD-LMLS" && SRL.rows.fs(iMax) == 240, sprintf("%.1f%%, %.1f%% and %.1f%% at 60, 120 and 240 Hz (by pipeline %.1f%% to %.1f%%, lowest for SG-IRLS at 60 Hz and highest for BWFD-LMLS at 240 Hz)", 100*bf(60), 100*bf(120), 100*bf(240), 100*srMin, 100*srMax), "SEM SHARE EXTREMES NOT SUPPORTED")
 "SR.15 edge clip", "Section 2.3", "All of this is confounded with the fixed", sprintf("%d-sample edge clip, which at β_gen = 1/3 and the log-midpoint VGF removes %.1f%%, %.1f%% and %.1f%% of samples at 60, 120 and 240 Hz (across the grid's VGF values %.1f%% to %.1f%% at 60 Hz and %.1f%% to %.1f%% at 240 Hz", EC.edgeClip, ecM(60, "pctDiscarded_mid"), ecM(120, "pctDiscarded_mid"), ecM(240, "pctDiscarded_mid"), ecM(60, "pctDiscarded_gridMin"), ecM(60, "pctDiscarded_gridMax"), ecM(240, "pctDiscarded_gridMin"), ecM(240, "pctDiscarded_gridMax")) + iff_local(all(abs(EC.ratio - 4) < 0.05), ", always about four times as much at 60 Hz", " RATIO NOT ABOUT FOUR")
 "HS.1 one-third fixed", "Part 1B", "Every other exponent moves towards", iff_local(all(round(HA.At.("b0.333_SGIRLS"), 3) == 0.333), "β_gen = 1/3 reads 0.333 at every tempo from 0.25 to 6 Hz", "1/3 DOES NOT READ 0.333 AT EVERY TEMPO")
 "HS.2 two-thirds by tempo", "Part 1B", "Every other exponent moves towards", "β_gen = 2/3 reads " + regexprep(haS("b0.667_SGIRLS", [0.25 1 1.5 2.5 4 6]), ", (\S+)$", " and $1") + " at 0.25, 1, 1.5, 2.5, 4 and 6 Hz"
 "HS.3 zero by tempo", "Part 1B", "Every other exponent moves towards", "β_gen = 0 (SG-OLS) reads " + regexprep(haS("b0.000_SGOLS", [0.25 1 1.5 2.5 4 6]), ", (\S+)$", " and $1")
 "HS.4 eccentricity", "Part 1B", "Eccentricity matters less than tempo", sprintf("(β_gen 2/3 at 2.5 Hz: %.3f at a/b 1.25, %.3f at 4.1)", HA.Ecc.SG_IRLS(HA.Ecc.ab == 1.25), HA.Ecc.SG_IRLS(HA.Ecc.ab == 4.1))
 "HS.5 weight at 3f0", "Part 1B", "Harmonic traversal of an ellipse obeys", sprintf("(weight %.2f to %.1f depending on eccentricity", min(w3O), max(w3O))
 "HS.6 theory RMS", "Part 1B", "This first-order", sprintf("to RMS ≤ 0.003 in %d of %d filter-regression cells", nnz(rmsE <= 0.003), numel(rmsE))
 "HS.7 fine-grid onsets", "Part 1B", "This first-order", sprintf("the theory's onsets lie within %.2f%% of the emulated ones (OLS)", 100 * max(abs(relO)))
 "HS.8 geometry onsets", "Part 1B", "This first-order", sprintf("(Butterworth-OLS onsets %.2f, %.2f and %.2f Hz at a/b 1.25, 2.16 and 4.1, against %.2f, %.2f and %.2f Hz predicted)", hgOn("BW", 2, 1.25, "OLS", "pipe"), hgOn("BW", 2, 2.16, "OLS", "pipe"), hgOn("BW", 2, 4.1, "OLS", "pipe"), hgOn("BW", 2, 1.25, "OLS", "theory"), hgOn("BW", 2, 2.16, "OLS", "theory"), hgOn("BW", 2, 4.1, "OLS", "theory"))
 "HS.9 Butterworth onset", "Part 1B", "Gain against tempo therefore collapses", sprintf("Gain falls to 0.9 at f0/f_c = %.3f-%.3f for Butterworth-OLS, a cutoff multiple of %.2f-%.2f (Butterworth-IRLS %.3f)", min(bwO), max(bwO), 1 / max(bwO), 1 / min(bwO), hgS.mean_scaledFine(hgS.family == "BW" & hgS.regression == "IRLS"))
 "HS.10 SG onset", "Part 1B", "Gain against tempo therefore collapses", sprintf("at f0·T_w = %.3f-%.3f for Savitzky-Golay-OLS in all three windows (SG-IRLS %.3f)", min(sgO), max(sgO), hgS.mean_scaledFineSafe(hgS.family == "SG" & hgS.regression == "IRLS"))
 "HS.11 production onsets", "Part 1B", "The production settings lose gain from", sprintf("lose gain from about %.1f Hz (BWFD, 10 Hz) and %.1f Hz (SG-IRLS, 0.17 s)", hgOn("BW", 2, 2.16, "OLS", "pipe"), hgOn("SG", 5, 2.16, "IRLS", "pipe"))
 "HS.12 map top 0.25 Hz", "Part 1B", "flattened on both sides of about 1 Hz", haRow(0.25, "0.25 Hz")
 "HS.13 map top 1 Hz", "Part 1B", "flattened on both sides of about 1 Hz", haRow(1, "1 Hz")
 "HS.14 map top 4 Hz", "Part 1B", "flattened on both sides of about 1 Hz", haRow(4, "4 Hz")
 "HS.15 curves analysed", "Part 1B", "A target is", addcomma(height(GT.P)) + " curves were analysed (" + addcomma(nnz(GT.P.source ~= "v015")) + " simulated, " + addcomma(nnz(GT.P.source == "v015")) + " empirical per-trial maps); " + addcomma(nnz(GT.P.reliable)) + " were reliable"
 "HS.16 zero-noise target", "Part 1B", "At σ = 0 and fixed", sprintf("(%d reliable targets: median %.4f, IQR %.3f-%.3f)", numel(gtA), median(gtA), prctile(gtA, 25), prctile(gtA, 75))
 "HS.17 per-trial gain", "Part 4", "Tempo in the per-trial maps: the target slides", sprintf("median gain runs from %.2f (Dhieb) to %.2f (Hickman PLAC, SG-IRLS)", min(gtS.median_gain(gtS.dataset == "Dhieb")), gtS.median_gain(gtS.dataset == "Hickman_PLAC" & gtS.pipeline == "SG-IRLS"))
 "HS.18 reliable share", "Part 4", "Tempo in the per-trial maps: the target slides", sprintf("reliable targets cover %.0f-%.0f%% of trials in domain", 100 * min(gtS.fun1_reliable(gtIn)), 100 * max(gtS.fun1_reliable(gtIn)))
 "HS.19 bins Fraser", "Part 4", "Tempo in the per-trial maps: the target slides", "Fraser " + gtBins("Fraser") + "; Cook ASD " + gtBins("Cook_ASD")
 "HS.20 bins Cook CTRL", "Part 4", "Tempo in the per-trial maps: the target slides", sprintf("Cook CTRL %.3f, %.3f, %.3f (%d, %d and %d trials)", gtb("Cook_CTRL", "<0.35", "median_target"), gtb("Cook_CTRL", "0.35-0.7", "median_target"), gtb("Cook_CTRL", "0.7-1.4", "median_target"), gtb("Cook_CTRL", "<0.35", "fun1_target"), gtb("Cook_CTRL", "0.35-0.7", "fun1_target"), gtb("Cook_CTRL", "0.7-1.4", "fun1_target"))
 "HS.21 bins Hickman", "Part 4", "Tempo in the per-trial maps: the target slides", "Hickman PLAC " + gtBins("Hickman_PLAC") + "; Hickman HALO " + gtBins("Hickman_HALO")
 "HS.22 fifth bin", "Part 4", "Tempo in the per-trial maps: the target slides", sprintf("A fifth bin above 2.8 Hz holds only %d (PLAC, %.3f), %d (HALO, %.3f) and %d (Dhieb, %.3f) reliable targets", gtb("Hickman_PLAC", ">2.8", "fun1_target"), gtb("Hickman_PLAC", ">2.8", "median_target"), gtb("Hickman_HALO", ">2.8", "fun1_target"), gtb("Hickman_HALO", ">2.8", "median_target"), gtb("Dhieb", ">2.8", "fun1_target"), gtb("Dhieb", ">2.8", "median_target"))
 "HS.23 LMM sample", "Part 4", "In an LMM on reliable in-domain SG-IRLS targets", "(trial as unit, " + addcomma(gl.nTrials) + " trials, " + gl.nSubjects + " subjects"
 "HS.24 LMM log f0", "Part 4", "In an LMM on reliable in-domain SG-IRLS targets", sprintf("log f0 has b = %.3f (SE %.4f, t(%d) = %.1f, 95%% CI %.3f-%.3f)", gl.b, gl.SE, gl.df, gl.t, gl.CIlo, gl.CIhi)
 "HS.25 fast bin", "Part 4", "In an LMM on reliable in-domain SG-IRLS targets", sprintf("Cook ASD on %d trials", gtb("Cook_ASD", "1.4-2.8", "fun1_target"))
 "HS.26 Dhieb fast", "Part 4", "In an LMM on reliable in-domain SG-IRLS targets", sprintf("(%.3f at 1.4-2.8 Hz)", gtb("Dhieb", "1.4-2.8", "median_target"))
};
MANUAL = strings(0, 1);   % none: every v006 number is rebuilt from a saved result (a generating script saves it)
nC = size(C, 1);  st = strings(nC, 1);  lit = strings(nC, 1);
for k = 1:nC
    lit(k) = strrep(string(C{k, 4}), char(8722), '-');  a = strfind(txt, C{k, 3});
    if isempty(a), st(k) = "ANCHOR LOST"; elseif contains(extractBetween(txt, a(1), min(strlength(txt), a(1) + 2500)), lit(k)), st(k) = "ok"; else, st(k) = "MISMATCH"; end
    fprintf("%-22s %-12s %s\n", C{k, 1}, st(k), lit(k));
end
fprintf("\nMANUAL (not checked):\n");  fprintf("  %s\n", MANUAL);
nBad = nnz(st ~= "ok");  save(fullfile(R, "checkClaimsV006_v002.mat"), "C", "st", "lit", "MANUAL");  fprintf("\n%d claims, %d not ok\n", nC, nBad);
if nBad > 0, error("claimsV006:fail", "%d claim(s) not ok (see report above)", nBad); end
fprintf("=== CLAIMS V006 OK ===\n");

function s = addcomma(n)
    s = string(regexprep(sprintf("%d", n), "(\d)(?=(\d{3})+$)", "$1,"));   % every thousands group (the one-comma version broke above 999,999)
end

function s = tab_local(t, lab)
% "-30.5 [-41.0, -22.0]" : the dz>0 group D and CI of the setting in the tail table
    k = find(strcmp(t.tab(:, 1), lab));  ci = t.tab{k, 6};  s = sprintf("%+.1f [%+.1f, %+.1f]", t.tab{k, 5}, ci(1), ci(2));
end

function x = red_local(S, variant, pipe)
% reduction (%) of the within-session tempo slope for a variant and pipeline
    x = S.reductionPct(S.variant == variant & S.pipeline == pipe);
end

function s = ea_local(R1, ds, nm)
% "-18.3 [-25.8, -10.9]" : paired contrast (pp) with its bootstrap CI, signed
    r = R1(string({R1.dataset}) == ds);  c = r.contrast(r.contrast.contrast == nm, :);
    s = string(sprintf("%+.1f [%+.1f, %+.1f]", c.diffPP, c.ciLo, c.ciHi));
end

function x = r2_local(R2, ds, grp)
% median paired change in beta_gen* (anchored - production) for a group
    r = R2(string({R2.dataset}) == ds);  x = r.tab.medDeltaBG(r.tab.group == grp);
end

function s = dot_local(x)
% ".776" : three decimals, no leading zero
    s = strrep(string(sprintf("%.3f", x)), "0.", ".");
end

function y = iff_local(c, a, b)
% MATLAB has no ternary operator.
if c, y = a; else, y = b; end
end

function s = d2_local(x)
% ".38" : two decimals, no leading zero (as the prose)
    s = strrep(string(sprintf("%.2f", x)), "0.", ".");
end

function s = sd_local(x)
% signed, no leading zero; four decimals only when |x| < .001 (as Table C1)
    if abs(x) < 0.001, s = sprintf("%+.4f", x); else, s = sprintf("%+.3f", x); end
    s = regexprep(s, "0\.", ".");
end

function s = c2cell(r)
% "+.056 [+.024, +.105]" : bootstrap difference with percentile CI (three decimals; four when |x| < .001), as Table C2
    s = sd_local(r.diff) + " [" + sd_local(r.ciLo) + ", " + sd_local(r.ciHi) + "]";
end

function s = c3cell(r)
% "+0.056 [+0.020, +0.133]" : bootstrap difference with percentile CI, three decimals, as the Part 4 prose
    s = string(sprintf("%+.3f [%+.3f, %+.3f]", r.diff, r.ciLo, r.ciHi));
end
