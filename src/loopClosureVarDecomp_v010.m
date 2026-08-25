function out = loopClosureVarDecomp_v010(opts)
% loopClosureVarDecomp_v010  Adds a k=3 synthesis view -- Fraser, Cook CTRL,
% Hickman PLAC only -- to v009's existing k=4/k=5/k=6 views. Nothing else
% changes: v009's per-dataset computation, cluster bootstrap, and DL/HKSJ/
% mKH synthesis machinery are byte-identical here; k=3 is one more call to
% the same crossDatasetSynth_local already used for the other three views.
%
% Motivation (2026-08-08): working out which datasets belong in a
% defensible headline pool, prompted directly by Dagmar's own question
% ("we cannot justify including Dhieb, Dagenais, ASD, HALO, PD on/off or
% Zarandi -- am I missing any?"). Each exclusion below has its own,
% different justification -- not one rule applied six times:
%   - Dagenais, James: out of methodological scope entirely (non-elliptical
%     paradigm, never eligible for loop closure at all; R7).
%   - Dhieb: measurement contamination, three converging signals (R7 gate
%     fail, lowest non-degenerate Mantel, elevated beta_gen*; Findings
%     #32/#65).
%   - PD ON/OFF (Hickman Study 1): excluded on data-quality grounds, D13.
%   - Zarandi: Scenario A, non-exchangeable with the other datasets,
%     k=6 pooling heterogeneity I2=97.3% (Finding #89).
%   - Hickman HALO: NOT dropped for the same reason as ASD below. PLAC and
%     HALO are the same subjects under two drug conditions (crossover,
%     28/35 shared) with an established null (Finding #30, reconfirmed
%     properly paired in Finding #139: constellation-median estimator
%     shows equivalence within MDC=0.03). Using PLAC alone instead of the
%     k=4 PLAC+HALO merge is a cleaner way to reach the same place -- no
%     shared-subject question even arises if only one arm is ever touched.
%   - Cook ASD: a genuinely different kind of exclusion, not the same logic
%     reapplied. CTRL and ASD are not repeated measures of the same people;
%     there is no null result licensing a merge or a drop the way HALO's
%     exclusion is licensed. Finding #139's own properly-powered comparison
%     (independent-groups cluster bootstrap, not a flat t-test) landed
%     genuinely INCONCLUSIVE, not equivalent -- underpowered even with the
%     correct machinery (M=18 vs M=21 subjects). Finding #114 already
%     flagged that Cook CTRL vs ASD is this literature's own primary
%     clinical contrast; dropping it from the headline pool is therefore
%     not "cleaning up the pool", it is declining to fold a genuinely
%     unresolved case-control question into an average that would partly
%     cancel out the very difference it exists to test. Handled instead as
%     its own dedicated comparison (Finding #139), not silently discarded --
%     the pooled headline answers "is beta_gen* < 1/3"; the CTRL/ASD
%     comparison answers a different question entirely, and pooling them
%     together would answer neither cleanly.
%
% See v008's header for the full cluster-bootstrap rationale and citations
% (Davison & Hinkley 1997; Field & Welsh 2007; Cameron, Gelbach & Miller
% 2008) -- not repeated here; v008 remains the canonical methodology
% reference, v009/v010 are data-source and scope updates only.
%
% Reads loopClosureResults_<tag>_all_shaped_xu_v007.mat for Cook CTRL/ASD,
% Hickman PLAC/HALO, Zarandi; loopClosureResults_Fraser_all_shaped_xu_v008.mat
% for Fraser -- unchanged from v009.
%
% Saves loopClosureVarDecomp_v010_<subsetTag>.mat, same subsetTag convention
% as v006-v009.
%
% Usage:
%   out = loopClosureVarDecomp_v010;
%
% Fraser, D.S. (2026)  v010

    arguments
        opts.NBoot          (1,1) double {mustBeInteger, mustBePositive} = 2000
        opts.NInner         (1,1) double {mustBeInteger, mustBePositive} = 200
        opts.RngSeed        (1,1) double {mustBeInteger}                 = 20260619
        opts.PrimaryPP      (1,1) double {mustBeInRange(opts.PrimaryPP, 1, 6)} = 4
        opts.MDC            (1,1) double {mustBePositive}                = 0.03
        opts.PipelineSubset (1,:) double {mustBeInteger, mustBeInRange(opts.PipelineSubset, 1, 6)} = 1:6
    end

    if numel(unique(opts.PipelineSubset)) ~= numel(opts.PipelineSubset)
        error("loopClosureVarDecomp_v010:duplicateSubsetEntries", "%s", ...
            "PipelineSubset contains duplicate pipeline indices.");
    end
    if numel(opts.PipelineSubset) < 2
        error("loopClosureVarDecomp_v010:subsetTooSmall", "%s", ...
            "PipelineSubset must contain at least 2 pipelines -- a median of one point is not a constellation-median.");
    end

    %% --- Config ----------------------------------------------------------
    PP_NAMES_V007 = ["BWFD-OLS","SG-OLS","BWFD-LMLS","SG-LMLS","BWFD-IRLS","SG-IRLS"];
    PP_PRIMARY    = opts.PrimaryPP;
    N_PIPELINES   = 6;
    Z90           = 1.644853626951472;
    ONE_THIRD     = 1/3;
    MDC           = opts.MDC;
    SUBSET        = opts.PipelineSubset;

    isDefaultSubset = isequal(sort(SUBSET), 1:6);
    subsetTag = "all6";
    if ~isDefaultSubset
        subsetTag = strjoin(lower(erase(PP_NAMES_V007(SUBSET), "-")), "");
    end
    fprintf("PipelineSubset = [%s]  (%s)\n\n", num2str(SUBSET), strjoin(PP_NAMES_V007(SUBSET), ", "));

    DATASETS = { ...
        "Cook_CTRL",    "Cook CTRL",    133; ...
        "Cook_ASD",     "Cook ASD",     133; ...
        "Hickman_PLAC", "Hickman PLAC", 133; ...
        "Hickman_HALO", "Hickman HALO", 133; ...
        "Fraser",       "Fraser",       240; ...
        "Zarandi",      "Zarandi",      100 };
    nDS = size(DATASETS, 1);

    srcDir  = fileparts(mfilename("fullpath"));
    outFile = fullfile(srcDir, sprintf("loopClosureVarDecomp_v010_%s.mat", subsetTag));
    rng(opts.RngSeed, "twister");

    cfg = struct("N_PIPELINES",N_PIPELINES, "PP_PRIMARY",PP_PRIMARY, "Z90",Z90, ...
        "NInner",opts.NInner, "NBoot",opts.NBoot, "ONE_THIRD",ONE_THIRD, "MDC",MDC, ...
        "SUBSET",SUBSET);

    %% --- Per-dataset computation -----------------------------------------
    rAll = cell(nDS, 1);
    for k = 1:nDS
        tag   = DATASETS{k,1};
        label = DATASETS{k,2};
        r = loadResults_local(srcDir, tag);
        rAll{k} = r;
        fprintf("Computing %-14s (N=%d) ...\n", label, numel(r));
        out(k) = computeDatasetEntry_local(r, label, tag, cfg);
    end

    %% --- Hickman PLAC+HALO merged entry (independence-respecting) --------
    tagsCol = string(DATASETS(:,1));
    placIdx = find(tagsCol == "Hickman_PLAC", 1);
    haloIdx = find(tagsCol == "Hickman_HALO", 1);
    rMerged = [rAll{placIdx}(:); rAll{haloIdx}(:)];
    fprintf("Computing %-14s ...\n", "Hickman (combined)");
    hickmanCombined = computeDatasetEntry_local(rMerged, "Hickman (combined)", ...
        "Hickman_PLAC+HALO", cfg);
    hickmanCombined.note = ["PLAC and HALO share n=28 of 35 subjects (crossover " ...
        "design, confirmed by subjectID intersection). The cluster bootstrap treats " ...
        "each real person's PLAC and HALO trials as one subject-cluster. Justified " ...
        "by Finding #30 (PLAC approx HALO, no haloperidol effect)."];

    %% --- Cross-dataset synthesis, three views (UNCHANGED machinery) ------
    zarIdx = find(tagsCol == "Zarandi", 1);
    idxK5  = setdiff(1:nDS, zarIdx);
    idxK6  = 1:nDS;
    outK4  = [out(idxK5(1)), out(idxK5(2)), hickmanCombined, out(idxK5(end))];

    % k=3 (v010 addition): Fraser, Cook CTRL, Hickman PLAC only -- Cook ASD,
    % Hickman HALO, Zarandi all excluded outright, not merged/downweighted.
    % See this file's header for the reasoning behind each exclusion.
    ctrlIdx = find(tagsCol == "Cook_CTRL", 1);
    fraserIdx = find(tagsCol == "Fraser", 1);
    idxK3 = [ctrlIdx, placIdx, fraserIdx];

    synthMed     = crossDatasetSynth_local(out(idxK5), "median",  ONE_THIRD, MDC);
    synthPrim    = crossDatasetSynth_local(out(idxK5), "primary", ONE_THIRD, MDC);
    synthMed_k4  = crossDatasetSynth_local(outK4,      "median",  ONE_THIRD, MDC);
    synthPrim_k4 = crossDatasetSynth_local(outK4,      "primary", ONE_THIRD, MDC);
    synthMed_k6  = crossDatasetSynth_local(out(idxK6), "median",  ONE_THIRD, MDC);
    synthPrim_k6 = crossDatasetSynth_local(out(idxK6), "primary", ONE_THIRD, MDC);
    synthMed_k3  = crossDatasetSynth_local(out(idxK3), "median",  ONE_THIRD, MDC);
    synthPrim_k3 = crossDatasetSynth_local(out(idxK3), "primary", ONE_THIRD, MDC);

    %% --- Report ------------------------------------------------------------
    printBanner_local("E3  CLUSTER BOOTSTRAP: SUBJECT COUNTS AND WITHIN-SUBJECT TRIAL COUNTS", ...
        "M = subject count (the cluster-bootstrap unit). Fraser is now all-trials (v009), no longer 1 trial/subject.");
    printClusterCounts_local(out);

    printBanner_local("E3  FLAT (PSEUDO-REPLICATED) vs CLUSTER (RIGOROUS) SE  --  constellation median", ...
        "flatSE = old v001-v007 method. clusterSE = two-stage subject x trial bootstrap. ratio = clusterSE/flatSE.");
    printFlatVsCluster_local("constellation median", out, "cluster_median", "flat_median");
    printFlatVsCluster_local("SG-LMLS (primary)",    out, "cluster_primary", "flat_primary");

    printBanner_local("E3  3-STAGE ABLATION VARIANCE DECOMPOSITION  --  constellation median", ...
        "stage0 = inversion only | stage1 = +within-subject trial resample | stage2 = +between-subject resample (=primary CI SE).");
    printAblation_local("constellation median", out, "cluster_median");
    printAblation_local("SG-LMLS (primary)",    out, "cluster_primary");

    printBanner_local("E3  SURROGATE-NOISE JUSTIFICATION  (shaped_xu vs fftnoise per dataset)", ...
        "LC-residual alpha (Finding #83) > ~3 means fftnoise breaks (Finding #67); shaped_xu required.");
    printAlphaJustification_local(out);

    printBanner_local("E3  PER-PIPELINE FORWARD-MAP RELIABILITY  (%monotone among resolved trials; %undef of all trials)", "");
    printPipeRel_local(out, PP_NAMES_V007);

    printBanner_local("E3  LEGACY 2-WAY DECOMPOSITION  (trial-level; does NOT correct for subject clustering)", "");
    printDecomp_local("constellation median", out, "dec_median");
    printDecomp_local("SG-LMLS (primary)",   out, "dec_primary");

    printBanner_local(sprintf("E3  MEDIAN CI + FORMAL OVERLAP TEST vs 1/3  (cluster bootstrap, t(M-1)-referenced, NBoot=%d, MDC=%.3f)", opts.NBoot, MDC), "");
    printTest_local("constellation median", out, "test_median");
    printTest_local("SG-LMLS (primary)",   out, "test_primary");
    printTest_local("Hickman combined (median, independence check)", hickmanCombined, "test_median");

    printBanner_local("E3  CROSS-DATASET SYNTHESIS  [k=5, FRASER ALL-TRIALS]  (Cook CTRL/ASD, Hickman PLAC, Hickman HALO, Fraser)", ...
        "Same five datasets as v008 (Fraser already replaced Pilot there); the change here is Fraser's own N (48->2829, per_subject->all).");
    printSynth_local("constellation median (primary)", synthMed);
    printSynth_local("SG-LMLS (sensitivity)",          synthPrim);

    printBanner_local("E3  CROSS-DATASET SYNTHESIS  [k=4, INDEPENDENCE-RESPECTING]  (Hickman PLAC+HALO merged)", "");
    printSynth_local("constellation median (primary)", synthMed_k4);
    printSynth_local("SG-LMLS (sensitivity)",          synthPrim_k4);

    printBanner_local("E3  CROSS-DATASET SYNTHESIS  [k=6, WITH ZARANDI]  (inclusion sensitivity)", "");
    printSynth_local("constellation median (primary)", synthMed_k6);
    printSynth_local("SG-LMLS (sensitivity)",          synthPrim_k6);

    printBanner_local("E3  CROSS-DATASET SYNTHESIS  [k=3, HEADLINE CANDIDATE]  (Fraser, Cook CTRL, Hickman PLAC only)", ...
        "Cook ASD, Hickman HALO, Zarandi excluded outright -- see header for the reasoning behind each. v010 addition, otherwise byte-identical to v009.");
    printSynth_local("constellation median (primary)", synthMed_k3);
    printSynth_local("SG-LMLS (sensitivity)",          synthPrim_k3);

    %% --- Save --------------------------------------------------------------
    meta = struct("NBoot", opts.NBoot, "NInner", opts.NInner, "RngSeed", opts.RngSeed, ...
        "PP_PRIMARY", PP_PRIMARY, "pipelineNames", PP_NAMES_V007, "Z90", Z90, ...
        "MDC", MDC, "oneThird", ONE_THIRD, "createdAt", datetime("now"), ...
        "source", "loopClosureResults_<tag>_all_shaped_xu_v007.mat (Fraser: _all_shaped_xu_v008.mat, CHANGED from v008's _per_subject_)", ...
        "datasetOrder", string(DATASETS(:,2)'), "zarandiAlphaLC_median", out(zarIdx).medAlphaLC, ...
        "pipelineSubset", SUBSET, "pipelineSubsetNames", PP_NAMES_V007(SUBSET), ...
        "isDefaultSubset", isDefaultSubset, ...
        "methodChange", "Data-source update only vs v008: Fraser's all-trials result (N=2829) replaces the per_subject one (N=48). Cluster-bootstrap methodology unchanged (Finding #136).");
    save(outFile, "out", "hickmanCombined", "meta", ...
        "synthMed", "synthPrim", "synthMed_k4", "synthPrim_k4", "synthMed_k6", "synthPrim_k6", ...
        "synthMed_k3", "synthPrim_k3", "-v7.3");
    fprintf("\nSaved: %s\n", outFile);
end

% ======================================================================== %
%  Per-dataset compute                                                     %
% ======================================================================== %

function ds = computeDatasetEntry_local(r, label, tag, cfg)
    N = numel(r);
    N_PIPELINES = cfg.N_PIPELINES;
    SUBSET = cfg.SUBSET;

    bgStarPP  = nan(N, N_PIPELINES);
    ciHiMat   = nan(N, N_PIPELINES);
    ciLoMat   = nan(N, N_PIPELINES);
    alphaLC   = nan(N, 1);
    subjID    = strings(N,1);

    for ti = 1:N
        bgStarPP(ti,:) = getRowVec_local(r(ti), "betaGenStar", N_PIPELINES);
        ciHiMat(ti,:)  = getRowVec_local(r(ti), "ciHi",        N_PIPELINES);
        ciLoMat(ti,:)  = getRowVec_local(r(ti), "ciLo",        N_PIPELINES);
        alphaLC(ti)    = getfieldScalar_local(r(ti), "alphaMaj");
        subjID(ti)     = getSubjID_local(r(ti));
    end
    if any(subjID == "")
        error("loopClosureVarDecomp_v010:missingSubjectID", "%s", ...
            sprintf("%d/%d trials in dataset '%s' have missing subjectID -- cluster bootstrap requires it for every trial.", ...
            sum(subjID==""), N, label));
    end

    ciWidth  = ciLoMat - ciHiMat;
    sdInvPP  = ciWidth / (2*cfg.Z90);

    bgStarMed = median(bgStarPP(:, SUBSET), 2, "omitnan");

    sdNaive = innerMedianSD_local(bgStarPP(:, SUBSET), sdInvPP(:, SUBSET), cfg.NInner);
    [sdInvMed, calib] = floorMedianSD_local(bgStarPP(:, SUBSET), sdNaive);

    clM = clusterEstimate_local(bgStarMed,                sdInvMed,                subjID, cfg);
    clP = clusterEstimate_local(bgStarPP(:,cfg.PP_PRIMARY), sdInvPP(:,cfg.PP_PRIMARY), subjID, cfg);

    flatM = bootstrapBeta_local(bgStarMed,                sdInvMed,                cfg.NBoot, cfg.ONE_THIRD);
    flatP = bootstrapBeta_local(bgStarPP(:,cfg.PP_PRIMARY), sdInvPP(:,cfg.PP_PRIMARY), cfg.NBoot, cfg.ONE_THIRD);

    decP = varDecomp_local(bgStarPP(:,cfg.PP_PRIMARY), sdInvPP(:,cfg.PP_PRIMARY));
    decM = varDecomp_local(bgStarMed,                  sdInvMed);

    pipeRel = pipelineReliability_local(ciWidth, bgStarPP);

    ds = struct();
    ds.label          = label;
    ds.tag            = tag;
    ds.N_total        = N;
    ds.N_subjects     = clM.M;
    ds.medAlphaLC     = median(alphaLC, "omitnan");
    ds.dec_primary    = decP;
    ds.dec_median     = decM;
    ds.cluster_median  = clM;
    ds.cluster_primary = clP;
    ds.flat_median     = flatM;
    ds.flat_primary     = flatP;
    ds.test_median    = overlapTestFromCluster_local(clM, cfg.ONE_THIRD, cfg.MDC);
    ds.test_primary   = overlapTestFromCluster_local(clP, cfg.ONE_THIRD, cfg.MDC);
    ds.boot_median    = struct("N", clM.M);
    ds.boot_primary   = struct("N", clP.M);
    ds.pipeRel        = pipeRel;
    ds.sdMedCalib     = calib;
    ds.note           = "";
end

function r = loadResults_local(srcDir, tag)
% CHANGED from v008: Fraser now loads the all-trials file, not per_subject.
    if tag == "Fraser"
        matFile = fullfile(srcDir, "loopClosureResults_Fraser_all_shaped_xu_v008.mat");
    else
        matFile = fullfile(srcDir, sprintf("loopClosureResults_%s_all_shaped_xu_v007.mat", tag));
    end
    if ~isfile(matFile)
        error("loopClosureVarDecomp_v010:notFound", "%s", ...
            sprintf("Loop-closure mat not found (FAILED PATH): %s", matFile));
    end
    D = load(matFile, "results");
    r = D.results(:);
end

function s = getSubjID_local(r)
    s = "";
    if isfield(r, "subjectID") && ~isempty(r.subjectID)
        s = string(r.subjectID);
    end
end

% ======================================================================== %
%  Two-stage cluster bootstrap (UNCHANGED FROM v008)                       %
% ======================================================================== %

function cl = clusterEstimate_local(point, sdInv, subjID, cfg)
    fin = isfinite(point);
    p = point(fin); s = sdInv(fin); sID = subjID(fin);
    s(~isfinite(s) | s < 0) = 0;

    uSubj = unique(sID);
    M = numel(uSubj);
    subjTrialIdx = cell(M,1);
    subjN = zeros(M,1);
    for i = 1:M
        subjTrialIdx{i} = find(sID == uSubj(i));
        subjN(i) = numel(subjTrialIdx{i});
    end

    subjMed = nan(M,1);
    for i = 1:M
        subjMed(i) = median(p(subjTrialIdx{i}), "omitnan");
    end
    pointEst = median(subjMed, "omitnan");

    NBoot  = cfg.NBoot;
    refVal = cfg.ONE_THIRD;

    stage0     = nan(NBoot,1);
    stage1     = nan(NBoot,1);
    stage2     = nan(NBoot,1);
    trialOnly2 = nan(NBoot,1);

    for bi = 1:NBoot
        pert0 = p + s .* randn(numel(p),1);
        s0med = nan(M,1);
        for i = 1:M
            s0med(i) = median(pert0(subjTrialIdx{i}), "omitnan");
        end
        stage0(bi) = median(s0med, "omitnan");

        s1med = nan(M,1);
        for i = 1:M
            idxI = subjTrialIdx{i};
            nI   = numel(idxI);
            draw = idxI(randi(nI, nI, 1));
            pert1 = p(draw) + s(draw) .* randn(nI,1);
            s1med(i) = median(pert1, "omitnan");
        end
        stage1(bi) = median(s1med, "omitnan");

        subjDraw = randi(M, M, 1);
        s2med   = nan(M,1);
        s2medNo = nan(M,1);
        for i = 1:M
            idxI = subjTrialIdx{subjDraw(i)};
            nI   = numel(idxI);
            draw = idxI(randi(nI, nI, 1));
            pert2 = p(draw) + s(draw) .* randn(nI,1);
            s2med(i)   = median(pert2,  "omitnan");
            s2medNo(i) = median(p(draw), "omitnan");
        end
        stage2(bi)     = median(s2med,   "omitnan");
        trialOnly2(bi) = median(s2medNo, "omitnan");
    end

    cl = struct();
    cl.M = M;
    cl.subjN_min = min(subjN); cl.subjN_med = median(subjN); cl.subjN_max = max(subjN);
    cl.N_trials  = numel(p);
    cl.pointEst  = pointEst;

    cl.seStage0 = std(stage0, 0, "omitnan");
    cl.seStage1 = std(stage1, 0, "omitnan");
    cl.seStage2 = std(stage2, 0, "omitnan");
    cl.seTrialOnly = std(trialOnly2, 0, "omitnan");

    cl.varInversion      = cl.seStage0^2;
    cl.varWithinSubject  = max(0, cl.seStage1^2 - cl.seStage0^2);
    cl.varBetweenSubject = max(0, cl.seStage2^2 - cl.seStage1^2);
    vTot = cl.varInversion + cl.varWithinSubject + cl.varBetweenSubject;
    if vTot > 0
        cl.fracInversion      = cl.varInversion      / vTot;
        cl.fracWithinSubject  = cl.varWithinSubject  / vTot;
        cl.fracBetweenSubject = cl.varBetweenSubject / vTot;
    else
        cl.fracInversion = NaN; cl.fracWithinSubject = NaN; cl.fracBetweenSubject = NaN;
    end

    cl.ciPercentile          = prctile(stage2, [2.5 97.5]);
    cl.ciPercentileTrialOnly = prctile(trialOnly2, [2.5 97.5]);
    if M >= 2
        tcrit = tinv(0.975, M-1);
        cl.dfT = M-1;
        cl.ciT          = [pointEst - tcrit*cl.seStage2,     pointEst + tcrit*cl.seStage2];
        cl.ciTTrialOnly = [pointEst - tcrit*cl.seTrialOnly,  pointEst + tcrit*cl.seTrialOnly];
    else
        cl.dfT = NaN; cl.ciT = [NaN NaN]; cl.ciTTrialOnly = [NaN NaN];
    end

    cl.pEmp = mean(stage2 >= refVal);
    if cl.seStage2 > 0
        cl.pGauss = normcdf((pointEst - refVal) / cl.seStage2);
    else
        cl.pGauss = NaN;
    end
end

function t = overlapTestFromCluster_local(cl, refVal, MDC)
    t = struct("median",cl.pointEst,"ciFM",cl.ciT,"seFM",cl.seStage2, ...
        "pEmp",cl.pEmp,"pGauss",cl.pGauss,"marginBeta",NaN,"marginMDC",NaN, ...
        "ciExcludes",false);
    t.marginBeta = refVal - cl.pointEst;
    t.marginMDC  = t.marginBeta / MDC;
    if all(isfinite(cl.ciT))
        t.ciExcludes = cl.ciT(2) < refVal;
    end
end

% ======================================================================== %
%  Legacy 2-way decomposition + flat bootstrap (UNCHANGED FROM v008)       %
% ======================================================================== %

function d = varDecomp_local(b, s)
    d = struct("N",0,"nDrop",0,"sigTotal",NaN,"sigWithin",NaN,"sigBio",NaN, ...
        "varTotal",NaN,"varWithin",NaN,"varBio",NaN,"fWithin",NaN,"fBio",NaN, ...
        "bioClamped",false);
    mask    = isfinite(b) & isfinite(s) & (s >= 0);
    d.nDrop = sum(isfinite(b)) - sum(mask);
    bb = b(mask); ss = s(mask);
    d.N = numel(bb);
    if d.N < 3, return; end
    d.varTotal  = var(bb, 0);
    d.varWithin = mean(ss.^2);
    vbio = d.varTotal - d.varWithin;
    if vbio < 0, d.bioClamped = true; vbio = 0; end
    d.varBio   = vbio;
    d.sigTotal  = sqrt(d.varTotal);
    d.sigWithin = sqrt(d.varWithin);
    d.sigBio    = sqrt(d.varBio);
    if d.varTotal > 0
        d.fWithin = d.varWithin / d.varTotal;
        d.fBio    = d.varBio    / d.varTotal;
    end
end

function b = bootstrapBeta_local(point, sdInv, NBoot, refVal)
    fin = isfinite(point);
    p   = point(fin);
    s   = sdInv(fin);
    s(~isfinite(s) | s < 0) = 0;
    N = numel(p);

    b = struct("N", N, "median", NaN, "ciFM", [NaN NaN], "seFM", NaN, "pGE", NaN, "refVal", refVal);
    if N < 3, return; end

    b.median = median(p);
    medFM    = nan(NBoot, 1);
    for bi = 1:NBoot
        idx = randi(N, N, 1);
        pr  = p(idx);
        medFM(bi) = median(pr + s(idx) .* randn(N, 1));
    end
    b.ciFM = prctile(medFM, [2.5 97.5]);
    b.seFM = std(medFM, 0);
    b.pGE  = mean(medFM >= refVal);
end

% ======================================================================== %
%  Decomposition / synthesis / calibration / reliability / readers /       %
%  printers -- ALL UNCHANGED FROM v006/v007/v008 BELOW THIS POINT          %
% ======================================================================== %

function S = crossDatasetSynth_local(dsList, which, refVal, MDC)
    nDS = numel(dsList);
    m = nan(nDS,1); se = nan(nDS,1); ub = nan(nDS,1); lab = strings(nDS,1);
    for k = 1:nDS
        lab(k) = dsList(k).label;
        if which == "median"
            tt = dsList(k).test_median;
        else
            tt = dsList(k).test_primary;
        end
        m(k) = tt.median; se(k) = tt.seFM;
        if all(isfinite(tt.ciFM)), ub(k) = tt.ciFM(2); end
    end
    ok = isfinite(m) & isfinite(se) & se > 0;

    S = struct("which",which,"labels",lab,"medians",m,"seFM",se,"ciUpper",ub, ...
        "nUsed",sum(ok),"range",NaN,"rangeMDC",NaN,"Q",NaN,"dfQ",NaN,"pHet",NaN, ...
        "mPool",NaN,"sePool",NaN,"poolMarginBeta",NaN,"poolMarginMDC",NaN, ...
        "poolZ",NaN,"poolP",NaN,"maxCiUpper",NaN,"allCIbelow",false, ...
        "tau2",NaN,"I2",NaN,"mPoolRE",NaN,"sePoolRE",NaN,"ciRE",[NaN NaN], ...
        "poolMarginBetaRE",NaN,"poolMarginMDC_RE",NaN,"poolZ_RE",NaN, ...
        "poolP_RE",NaN,"reExcludes",false, ...
        "dfHKSJ",NaN,"sePoolHKSJ",NaN,"ciHKSJ",[NaN NaN], ...
        "poolMarginBetaHKSJ",NaN,"poolMarginMDC_HKSJ",NaN,"poolT_HKSJ",NaN, ...
        "poolP_HKSJ",NaN,"hksjExcludes",false, ...
        "qRaw",NaN,"qStar",NaN,"sePoolMKH",NaN,"ciMKH",[NaN NaN], ...
        "poolMarginBetaMKH",NaN,"poolMarginMDC_MKH",NaN,"poolT_MKH",NaN, ...
        "poolP_MKH",NaN,"mkhExcludes",false);
    if sum(ok) < 2, return; end

    S.range      = max(m(ok)) - min(m(ok));
    S.rangeMDC   = S.range / MDC;
    S.maxCiUpper = max(ub(ok));
    S.allCIbelow = S.maxCiUpper < refVal;

    w  = 1 ./ se(ok).^2;
    mp = sum(w .* m(ok)) / sum(w);
    S.mPool  = mp;
    S.sePool = sqrt(1 / sum(w));
    S.Q   = sum(w .* (m(ok) - mp).^2);
    S.dfQ = sum(ok) - 1;
    if S.dfQ >= 1, S.pHet = 1 - chi2cdf(S.Q, S.dfQ); end

    S.poolMarginBeta = refVal - mp;
    S.poolMarginMDC  = S.poolMarginBeta / MDC;
    S.poolZ = (refVal - mp) / S.sePool;
    S.poolP = normcdf((mp - refVal) / S.sePool);

    C    = sum(w) - sum(w.^2)/sum(w);
    tau2 = 0;
    if C > 0, tau2 = max(0, (S.Q - S.dfQ) / C); end
    S.tau2 = tau2;
    if S.Q > 0, S.I2 = 100 * max(0, (S.Q - S.dfQ) / S.Q); end
    wRE  = 1 ./ (se(ok).^2 + tau2);
    mpRE = sum(wRE .* m(ok)) / sum(wRE);
    zc   = 1.959963984540054;
    S.mPoolRE  = mpRE;
    S.sePoolRE = sqrt(1 / sum(wRE));
    S.ciRE     = [mpRE - zc*S.sePoolRE, mpRE + zc*S.sePoolRE];
    S.poolMarginBetaRE = refVal - mpRE;
    S.poolMarginMDC_RE = S.poolMarginBetaRE / MDC;
    S.poolZ_RE = (refVal - mpRE) / S.sePoolRE;
    S.poolP_RE = normcdf((mpRE - refVal) / S.sePoolRE);
    S.reExcludes = S.ciRE(2) < refVal;

    kStudies = sum(ok);
    dfH = kStudies - 1;
    S.dfHKSJ = dfH;
    if dfH >= 1
        qStat   = sum(wRE .* (m(ok) - mpRE).^2) / dfH;
        seHKSJ  = sqrt(qStat / sum(wRE));
        tcrit   = tinv(0.975, dfH);
        S.sePoolHKSJ = seHKSJ;
        S.ciHKSJ     = [mpRE - tcrit*seHKSJ, mpRE + tcrit*seHKSJ];
        S.poolMarginBetaHKSJ = refVal - mpRE;
        S.poolMarginMDC_HKSJ = S.poolMarginBetaHKSJ / MDC;
        if seHKSJ > 0
            S.poolT_HKSJ = (refVal - mpRE) / seHKSJ;
            S.poolP_HKSJ = tcdf((mpRE - refVal) / seHKSJ, dfH);
        end
        S.hksjExcludes = S.ciHKSJ(2) < refVal;

        qStar = max(1, qStat);
        S.qRaw  = qStat;
        S.qStar = qStar;
        seMKH   = sqrt(qStar / sum(wRE));
        S.sePoolMKH = seMKH;
        S.ciMKH     = [mpRE - tcrit*seMKH, mpRE + tcrit*seMKH];
        S.poolMarginBetaMKH = refVal - mpRE;
        S.poolMarginMDC_MKH = S.poolMarginBetaMKH / MDC;
        if seMKH > 0
            S.poolT_MKH = (refVal - mpRE) / seMKH;
            S.poolP_MKH = tcdf((mpRE - refVal) / seMKH, dfH);
        end
        S.mkhExcludes = S.ciMKH(2) < refVal;
    end
end

function [sdMed, calib] = floorMedianSD_local(bgStarPPSubset, sdNaive)
    MEDIAN_SE_CONST = 1.2533141373155003;
    nConv  = sum(isfinite(bgStarPPSubset), 2);
    empSD  = std(bgStarPPSubset, 0, 2, "omitnan");
    empSD(nConv < 2) = 0;
    floorSD = empSD ./ sqrt(max(nConv,1)) * MEDIAN_SE_CONST;
    floorSD(nConv < 2) = 0;

    sdMed = max(sdNaive, floorSD);

    calib = struct();
    calib.medNaive       = median(sdNaive,  "omitnan");
    calib.medFloor       = median(floorSD,  "omitnan");
    calib.medCorrected   = median(sdMed,    "omitnan");
    floorBound = floorSD > sdNaive;
    calib.pctFloorBound  = 100 * mean(floorBound);
end

function pr = pipelineReliability_local(ciWidth, bgStarPP)
    nPP = size(ciWidth, 2);
    pr = struct("pctMono",nan(1,nPP), "pctUndef",nan(1,nPP), "medBgStar",nan(1,nPP));
    for pp = 1:nPP
        finite = isfinite(ciWidth(:,pp));
        pr.pctMono(pp)   = 100 * mean(ciWidth(finite,pp) > 0);
        pr.pctUndef(pp)  = 100 * mean(~finite);
        pr.medBgStar(pp) = median(bgStarPP(:,pp), "omitnan");
    end
end

function sdMed = innerMedianSD_local(bgStarPPSubset, sdInvPPSubset, NInner)
    N = size(bgStarPPSubset, 1);
    sdMed = zeros(N, 1);
    for ti = 1:N
        pts  = bgStarPPSubset(ti, :);
        sds  = sdInvPPSubset(ti, :);
        conv = isfinite(pts);
        if ~any(conv), sdMed(ti) = 0; continue; end
        p = pts(conv);
        d = sds(conv);
        d(~isfinite(d) | d < 0) = 0;
        if all(d == 0), sdMed(ti) = 0; continue; end
        draws = p(:)' + d(:)' .* randn(NInner, numel(p));
        meds  = median(draws, 2);
        sdMed(ti) = std(meds, 0);
    end
end

function v = getfieldScalar_local(s, fld)
    v = NaN;
    if isfield(s, fld)
        x = s.(fld);
        if isnumeric(x) && isscalar(x), v = double(x); end
    end
end

function row = getRowVec_local(s, fld, nExpected)
    row = nan(1, nExpected);
    if isfield(s, fld)
        x = s.(fld);
        if isnumeric(x) && numel(x) == nExpected
            row = double(x(:)');
        end
    end
end

function printClusterCounts_local(out)
    fprintf("%-24s %8s %10s %10s %10s %10s\n", "Dataset", "N_trials", "M(subj)", "min/subj", "med/subj", "max/subj");
    fprintf("%s\n", repmat('-', 1, 76));
    for k = 1:numel(out)
        c = out(k).cluster_median;
        fprintf("%-24s %8d %10d %10d %10.1f %10d\n", ...
            out(k).label, out(k).N_total, out(k).N_subjects, c.subjN_min, c.subjN_med, c.subjN_max);
    end
end

function printFlatVsCluster_local(estLabel, out, clFld, flatFld)
    fprintf("\n--- Flat vs cluster SE [%s] ---\n", estLabel);
    fprintf("%-24s %10s %10s %10s   %s\n", "Dataset", "flatSE", "clusterSE", "ratio", "");
    fprintf("%s\n", repmat('-', 1, 70));
    for k = 1:numel(out)
        cl = out(k).(clFld);
        fl = out(k).(flatFld);
        ratio = NaN;
        if isfinite(fl.seFM) && fl.seFM > 0, ratio = cl.seStage2 / fl.seFM; end
        flag = "";
        if isfinite(ratio) && ratio > 1.2, flag = "flat UNDERSTATED uncertainty"; end
        fprintf("%-24s %10.4f %10.4f %10.2fx   %s\n", ...
            out(k).label, fl.seFM, cl.seStage2, ratio, flag);
    end
end

function printAblation_local(estLabel, out, clFld)
    fprintf("\n--- 3-stage ablation [%s] ---\n", estLabel);
    fprintf("%-24s %8s %8s %8s   %6s %6s %6s\n", ...
        "Dataset", "SE_inv", "SE_wSub", "SE_bSub", "inv%", "wSub%", "bSub%");
    fprintf("%s\n", repmat('-', 1, 76));
    for k = 1:numel(out)
        c = out(k).(clFld);
        fprintf("%-24s %8.4f %8.4f %8.4f   %6.1f %6.1f %6.1f\n", ...
            out(k).label, sqrt(c.varInversion), sqrt(c.varWithinSubject), sqrt(c.varBetweenSubject), ...
            100*c.fracInversion, 100*c.fracWithinSubject, 100*c.fracBetweenSubject);
    end
end

function printAlphaJustification_local(out)
    fprintf("%-14s %10s   %s\n", "Dataset", "medAlphaLC", "surrogate");
    fprintf("%s\n", repmat('-', 1, 60));
    for k = 1:numel(out)
        a = out(k).medAlphaLC;
        verdict = "shaped_xu OK (fftnoise also viable, alpha<3)";
        if isfinite(a) && a >= 3, verdict = "shaped_xu REQUIRED (fftnoise breaks, alpha>=3)"; end
        fprintf("%-14s %10.3f   %s\n", out(k).label, a, verdict);
    end
end

function printPipeRel_local(out, ppNames)
    nPP = numel(ppNames);
    fprintf("%-14s", "Dataset");
    for pp = 1:nPP, fprintf(" %19s", ppNames(pp)); end
    fprintf("\n%-14s", "");
    for pp = 1:nPP, fprintf(" %9s %9s", "%mono", "%undef"); end
    fprintf("\n%s\n", repmat('-', 1, 14 + 20*nPP));
    for k = 1:numel(out)
        fprintf("%-14s", out(k).label);
        for pp = 1:nPP
            fprintf(" %8.1f%% %8.1f%%", out(k).pipeRel.pctMono(pp), out(k).pipeRel.pctUndef(pp));
        end
        fprintf("\n");
    end
end

function printDecomp_local(estLabel, out, fld)
    fprintf("\n--- Legacy 2-way decomposition [%s] ---\n", estLabel);
    fprintf("%-14s %6s  %8s %8s %8s   %6s %6s  %5s  %s\n", ...
        "Dataset","N","sig_tot","sig_inv","sig_bio","fBio%","fInv%","nDrp","flag");
    fprintf("%s\n", repmat('-', 1, 78));
    for k = 1:numel(out)
        d = out(k).(fld);
        flag = "";
        if d.bioClamped, flag = "CLAMP(bio=0)"; end
        fprintf("%-14s %6d  %8.4f %8.4f %8.4f   %6.1f %6.1f  %5d  %s\n", ...
            out(k).label, d.N, d.sigTotal, d.sigWithin, d.sigBio, ...
            100*d.fBio, 100*d.fWithin, d.nDrop, flag);
    end
end

function printTest_local(estLabel, out, fld)
    fprintf("\n--- Median CI + overlap test [%s] ---\n", estLabel);
    fprintf("%-14s %6s   %8s  %-20s %8s %7s   %8s %10s   %s\n", ...
        "Dataset","M(subj)","median","t(M-1) 95% CI","margin","/MDC","pEmp","pGauss","CI<1/3");
    fprintf("%s\n", repmat('-', 1, 100));
    for k = 1:numel(out)
        t = out(k).(fld);
        nThis = out(k).boot_median.N;
        if fld == "test_primary", nThis = out(k).boot_primary.N; end
        ciStr = sprintf("[%7.4f,%8.4f]", t.ciFM(1), t.ciFM(2));
        yn = "no ";
        if t.ciExcludes, yn = "YES"; end
        fprintf("%-14s %6d   %8.4f  %-20s %8.4f %7.2f   %8.4f %10.2e   %s\n", ...
            out(k).label, nThis, t.median, ciStr, t.marginBeta, t.marginMDC, ...
            t.pEmp, t.pGauss, yn);
    end
end

function printSynth_local(estLabel, S)
    fprintf("\n--- Cross-dataset synthesis [%s] ---\n", estLabel);
    if S.nUsed < 2
        fprintf("  insufficient datasets with finite SE (nUsed=%d)\n", S.nUsed);
        return
    end
    fprintf("  convergence range across %d datasets : %.4f  (= %.2f x MDC)\n", ...
        S.nUsed, S.range, S.rangeMDC);
    hetMsg = "medians consistent with one shared value";
    if isfinite(S.pHet) && S.pHet < 0.05, hetMsg = "genuine between-dataset differences"; end
    fprintf("  heterogeneity (Cochran Q)           : Q=%.2f, df=%d, p_het=%.4f  [%s]\n", ...
        S.Q, S.dfQ, S.pHet, hetMsg);
    fprintf("  fixed-effect pooled beta_gen*       : %.4f +/- %.4f\n", S.mPool, S.sePool);
    fprintf("  between-dataset heterogeneity       : tau=%.4f, I2=%.1f%%\n", sqrt(S.tau2), S.I2);
    fprintf("  random-effects pooled beta_gen* (DL): %.4f +/- %.4f   95%% CI [%.4f, %.4f]\n", ...
        S.mPoolRE, S.sePoolRE, S.ciRE(1), S.ciRE(2));
    fprintf("  RE distance below 1/3               : %.4f  (= %.2f MDCs, %.1f sigma; P>=1/3 = %.2e)\n", ...
        S.poolMarginBetaRE, S.poolMarginMDC_RE, S.poolZ_RE, S.poolP_RE);
    if isfinite(S.sePoolHKSJ)
        fprintf("  HKSJ (small-k) pooled beta_gen*      : %.4f +/- %.4f   95%% CI [%.4f, %.4f]  (t, df=%d)\n", ...
            S.mPoolRE, S.sePoolHKSJ, S.ciHKSJ(1), S.ciHKSJ(2), S.dfHKSJ);
        fprintf("  HKSJ distance below 1/3              : %.4f  (= %.2f MDCs, t=%.2f; P>=1/3 = %.2e)\n", ...
            S.poolMarginBetaHKSJ, S.poolMarginMDC_HKSJ, S.poolT_HKSJ, S.poolP_HKSJ);
    end
    if isfinite(S.sePoolMKH)
        fprintf("  mKH (Roever et al 2015) q=%.3f->q*=%.3f  : %.4f +/- %.4f   95%% CI [%.4f, %.4f]  (t, df=%d)\n", ...
            S.qRaw, S.qStar, S.mPoolRE, S.sePoolMKH, S.ciMKH(1), S.ciMKH(2), S.dfHKSJ);
        fprintf("  mKH distance below 1/3               : %.4f  (= %.2f MDCs, t=%.2f; P>=1/3 = %.2e)\n", ...
            S.poolMarginBetaMKH, S.poolMarginMDC_MKH, S.poolT_MKH, S.poolP_MKH);
    end
    verdict   = "NO"; if S.allCIbelow,  verdict   = "YES"; end
    reVerdict = "NO"; if S.reExcludes,  reVerdict = "YES"; end
    hkVerdict = "NO"; if S.hksjExcludes, hkVerdict = "YES"; end
    fprintf("  every dataset 95%% CI upper < 1/3 ?   : %s  (max upper bound = %.4f)\n", verdict, S.maxCiUpper);
    fprintf("  RE pooled 95%% CI upper < 1/3 ?       : %s  (upper = %.4f)\n", reVerdict, S.ciRE(2));
    if isfinite(S.sePoolHKSJ)
        fprintf("  HKSJ pooled 95%% CI upper < 1/3 ?     : %s  (upper = %.4f)\n", hkVerdict, S.ciHKSJ(2));
    end
    if isfinite(S.sePoolMKH)
        mkVerdict = "NO"; if S.mkhExcludes, mkVerdict = "YES"; end
        fprintf("  mKH pooled 95%% CI upper < 1/3 ?      : %s  (upper = %.4f)\n", mkVerdict, S.ciMKH(2));
    end
end

function printBanner_local(title, body)
    bar = repmat('=', 1, 78);
    fprintf("\n%s\n%s\n%s\n", bar, title, bar);
    if strlength(body) > 0, fprintf("%s\n", body); end
end
