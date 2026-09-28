function tier3 = tier3DesignPrecision_v001(batchResults, options)
% TIER3DESIGNPRECISION_V001  Tier 3: dataset/design-level precision.
%
% Answers the question Tier 1/Tier 2 don't: given a batch of Tier 1
% results from the researcher's OWN pilot or full dataset, is the
% resulting dataset estimate -- or, for a two-arm study, the resulting
% contrast -- precise enough on the generator scale to support the
% clinical inference being attempted, and if not, what would help.
%
% This is a direct, faithful generalisation of the main project's own
% cluster-bootstrap machinery (Finding #216) to arbitrary user-supplied
% data: the three-stage ablation variance decomposition of
% loopClosureVarDecomp_v012_gated.m's clusterEstimate_local (inversion /
% within-subject / between-subject shares), and, for a two-arm contrast,
% either compareHickmanPlacHalo_v002.m's one-sample paired bootstrap or
% compareCookCtrlAsd_v002.m's two-sample independent-groups bootstrap,
% whichever the supplied design calls for. The bootstrap logic below is a
% line-for-line port of those three functions (not reimplemented from a
% description), generalised only to read from this tool's own
% runToolLoopClosure_v001/v002 batch shape (.subLabel/.sesLabel/.trial)
% instead of the main project's saved dataset .mat files. This is why it
% is a tool-build item, not new analysis: the method, its constants
% (Z90, the MEDIAN_SE_CONST floor, the NBoot/NInner defaults) and its
% validation are already established; only the data source is new.
%
% MODES:
%   "single"       - report this one batch's own precision and variance
%                    decomposition. No contrast.
%   "paired"       - within-subject contrast between two conditions
%                    recorded in each trial's own .sesLabel (e.g.
%                    "PLAC" vs "HALO"). Requires batchResults to come
%                    from runToolLoopClosure_v002 or later (v001's
%                    batch has no .sesLabel field -- this function
%                    checks for the field and fails loud, not silently,
%                    if it is absent).
%   "independent"  - two-sample contrast between two subject groups
%                    supplied via options.GroupMap (subLabel -> group
%                    name), e.g. a diagnosis or genotype label that has
%                    no natural place in a per-trial BIDS filename.
%
% WHAT "REQUIRED SE" MEANS, stated precisely because it is easy to
% over-read: options.MDC and the two-sample/one-sample bootstrap already
% in hand give a required SE for equivalence UNDER THE ASSUMPTION the
% true difference is zero -- req = MDC / tinv(0.95, df) -- exactly the
% calculation the manuscript's own General Discussion states for Cook
% CTRL-vs-ASD ("~0.018 ... even if the true difference were zero"). This
% is not a claim about what SE would be needed to detect a real,
% nonzero difference; the corpus this method was built on has no
% positive control for that (Finding #216's own stated limit), and
% neither does this function.
%
% WHAT THIS DOES NOT DO: decide how many trials or subjects to collect.
% It reports precision and variance shares for the data actually
% supplied, and, when a contrast fails, the SE a design of the SAME
% subject/trial structure would need for equivalence -- a diagnostic
% for pilot data, not a prospective sample-size calculator. Turning that
% into a sample-size recommendation would be a genuinely new statistical
% method, out of scope here (see docs/TODO_ToolBuild_v001.md Part K,
% "not yet built" list).
%
% SYNTAX:
%   tier3 = tier3DesignPrecision_v001(batchResults, Mode="single")
%   tier3 = tier3DesignPrecision_v001(batchResults, Mode="paired", Conditions=["PLAC","HALO"])
%   tier3 = tier3DesignPrecision_v001(batchResults, Mode="independent", GroupMap=gm, Groups=["CTRL","ASD"])
%
% INPUTS:
%   batchResults - Output of runToolLoopClosure_v001/v002: Nx1 struct
%                  array, each with .ok, .subLabel, (.sesLabel if v002+),
%                  and (if ok) .trial, a processTrialLoopClosure_v001
%                  result. Trials with ok==false are excluded, counted
%                  and reported, not silently dropped.
%   Name-value options:
%     Mode       "single" (default) | "paired" | "independent"
%     Pipeline   "median" (default; the six-pipeline constellation
%                median per trial, matching the main project's own
%                primary cluster_median estimator) or one of
%                {"BWFD-OLS","SG-OLS","BWFD-LMLS","SG-LMLS","BWFD-IRLS",
%                "SG-IRLS"} to use a single named pipeline's betaGenStar
%                instead (matching cluster_primary).
%     Conditions 1x2 string, required for Mode="paired". Matched
%                case-sensitively against each trial's own .sesLabel.
%     GroupMap   containers.Map, string subLabel -> string group name.
%                Required for Mode="independent". Every subLabel
%                appearing in batchResults must have an entry, or this
%                fails loud rather than silently excluding that subject.
%     Groups     1x2 string, required for Mode="independent": the two
%                group names to contrast (must both appear as values in
%                GroupMap).
%     MDC        0.03 (default), the clinical threshold.
%     NBoot      2000 (default), matching loopClosureVarDecomp_v012_gated.m
%                and both compare*_v002.m scripts.
%     NInner     200 (default), matching the same three scripts.
%     RngSeed    20260619 (default, matching the same three scripts, so a
%                cross-check against them at identical inputs is exact
%                rather than approximate).
%
% OUTPUT: tier3 struct.
%   Mode="single" fields: label, pipeline, M (subject count), N_trials,
%     subjN_min/med/max, pointEst, seStage0/1/2, fracInversion,
%     fracWithinSubject, fracBetweenSubject, ciT, dfT, MDC, semAdequate
%     (seStage2 < MDC/2.77).
%   Mode="paired"/"independent" additionally: conditions or groups,
%     armA/armB (each a full Mode="single" struct for that arm alone),
%     diffPoint, seDiff, df, ci95_t, ci90_t, diffExcludesZero,
%     equivalentWithinMDC, verdict, and, only when NOT equivalentWithinMDC,
%     requiredSEForEquivalence (see note above) and designNote (a plain-
%     language string naming the between-subject share just computed and
%     stating that a paired/within-subject design of the same size would
%     remove it directly, when armA/armB's own fracBetweenSubject exceeds
%     fracInversion + fracWithinSubject combined -- omitted otherwise,
%     since the recommendation only follows when between-subject
%     variance is genuinely the dominant term).
%
% VALIDATION STATUS: NOT YET DONE. This file has been written and
% checked by hand against the three reference scripts it ports
% (loopClosureVarDecomp_v012_gated.m's clusterEstimate_local,
% compareHickmanPlacHalo_v002.m's pairedOneEstimator_local,
% compareCookCtrlAsd_v002.m's compareIndependent_local -- constants,
% bootstrap structure and stage definitions compared line by line, not
% reimplemented from a description) and passes a bracket/paren balance
% check, but it has NOT been run, NOT been through check_matlab_code, and
% NOT been cross-checked against real output of any kind. Required before
% this function, or any number it produces, is cited anywhere:
%   1. check_matlab_code, both this file and runToolLoopClosure_v002.m.
%   2. Mode="single" against tool/examples/fraser_sub-103 (one subject,
%      so this only exercises the inversion/within-subject stages
%      meaningfully -- M=1 makes the between-subject share and the
%      t(M-1) CI undefined by construction, which this file's own
%      singleDatasetPrecision_local warns about rather than hides).
%   3. Mode="paired" and Mode="independent" cross-checked against the
%      main project's own frozen output -- read
%      loopClosureResults_Hickman_PLAC/HALO_all_shaped_xu_v015.mat and
%      loopClosureResults_Cook_CTRL/ASD_all_shaped_xu_v015.mat directly
%      (bypassing Motion-BIDS ingest, which those datasets were never
%      run through) into this function's own (point,sd,subjectID) shape,
%      and confirm diffPoint/seDiff/verdict reproduce
%      compareHickmanPlacHalo_v002.m's and compareCookCtrlAsd_v002.m's
%      own saved .mat output exactly at the same RngSeed. This is a
%      code-correctness check against already-frozen numbers, not new
%      analysis. A first cut of the adapter exists
%      (checkTier3AgainstMainProject_v001.m, drafted alongside this
%      file) but is complete only for the Hickman half -- the Cook
%      half stops after loading the reference .mat, since
%      compareCookCtrlAsd_v002.m's own output struct field names have
%      not been confirmed against it. Neither half has been run.
% Until all three steps are done and their results recorded here, treat
% every number this function returns as unverified.
%
% See also: runToolLoopClosure_v001, runToolLoopClosure_v002,
%           processTrialLoopClosure_v001, tier2BatchValidation_v001
%
% Fraser, D.S. (2026)

arguments
    batchResults (:,1) struct
    options.Mode (1,1) string {mustBeMember(options.Mode, ["single","paired","independent"])} = "single"
    options.Pipeline (1,1) string = "median"
    options.Conditions (1,2) string = ["",""]
    options.GroupMap = containers.Map('KeyType', 'char', 'ValueType', 'char')
    options.Groups (1,2) string = ["",""]
    options.MDC (1,1) double {mustBePositive} = 0.03
    options.NBoot (1,1) double {mustBeInteger, mustBePositive} = 2000
    options.NInner (1,1) double {mustBeInteger, mustBePositive} = 200
    options.RngSeed (1,1) double {mustBeInteger} = 20260619
end

Z90 = 1.644853626951472;
cfg = struct('MDC', options.MDC, 'NBoot', options.NBoot, 'NInner', options.NInner, 'Z90', Z90);

okMask = [batchResults.ok];
if ~any(okMask)
    error('tier3DesignPrecision_v001:NoSuccessfulTrials', '%s', ...
        'batchResults contains zero trials with ok==true -- nothing to compute.');
end

switch options.Mode
case "single"
    [pt, sub] = extractPointsAndSubjects_local(batchResults(okMask), options.Pipeline);
    rng(options.RngSeed, 'twister');
    tier3 = singleDatasetPrecision_local(pt, sub, cfg, "single-dataset", options.Pipeline);
    tier3.mode = "single";

case "paired"
    if any(options.Conditions == "")
        error('tier3DesignPrecision_v001:MissingConditions', '%s', ...
            'Mode="paired" requires Conditions=[condA,condB] naming two values of .sesLabel.');
    end
    if ~isfield(batchResults, 'sesLabel')
        error('tier3DesignPrecision_v001:NoSesLabel', '%s', ...
            'batchResults has no .sesLabel field -- it must come from runToolLoopClosure_v002 or later, not v001.');
    end
    condA = options.Conditions(1); condB = options.Conditions(2);
    okA = okMask & (string({batchResults.sesLabel}) == condA);
    okB = okMask & (string({batchResults.sesLabel}) == condB);
    if ~any(okA)
        error('tier3DesignPrecision_v001:NoTrialsCondA', '%s', sprintf('No successful trials with sesLabel="%s".', condA));
    end
    if ~any(okB)
        error('tier3DesignPrecision_v001:NoTrialsCondB', '%s', sprintf('No successful trials with sesLabel="%s".', condB));
    end
    [pA, sA] = extractPointsAndSubjects_local(batchResults(okA), options.Pipeline);
    [pB, sB] = extractPointsAndSubjects_local(batchResults(okB), options.Pipeline);

    shared = intersect(unique(sA), unique(sB));
    if isempty(shared)
        error('tier3DesignPrecision_v001:NoSharedSubjects', '%s', sprintf( ...
            'No subject has trials in both "%s" and "%s" -- Mode="paired" needs within-subject data. Use Mode="independent" for separate cohorts.', ...
            condA, condB));
    end
    onlyA = setdiff(unique(sA), shared); onlyB = setdiff(unique(sB), shared);

    rng(options.RngSeed, 'twister');
    armA = singleDatasetPrecision_local(pA, sA, cfg, condA, options.Pipeline);
    armB = singleDatasetPrecision_local(pB, sB, cfg, condB, options.Pipeline);

    rng(options.RngSeed, 'twister');
    contrast = pairedContrast_local(pA, sA, pB, sB, shared, cfg);

    tier3 = struct();
    tier3.mode = "paired";
    tier3.conditions = [condA, condB];
    tier3.M = numel(shared);
    tier3.onlyInA = onlyA; tier3.onlyInB = onlyB;
    tier3.armA = armA; tier3.armB = armB;
    tier3 = mergeStruct_local(tier3, contrast);
    tier3 = addDesignNote_local(tier3, armA, armB, condA, condB);

case "independent"
    if any(options.Groups == "")
        error('tier3DesignPrecision_v001:MissingGroups', '%s', ...
            'Mode="independent" requires Groups=[groupA,groupB].');
    end
    gm = options.GroupMap;
    if ~isa(gm, 'containers.Map') || gm.Count == 0
        error('tier3DesignPrecision_v001:MissingGroupMap', '%s', ...
            'Mode="independent" requires GroupMap, a containers.Map from subLabel (char) to group name (char).');
    end
    groupA = options.Groups(1); groupB = options.Groups(2);

    subsAll = unique(string({batchResults(okMask).subLabel}));
    missing = subsAll(~isKey(gm, cellstr(subsAll)));
    if ~isempty(missing)
        error('tier3DesignPrecision_v001:GroupMapIncomplete', '%s', sprintf( ...
            'GroupMap has no entry for subject(s): %s. Every subject in batchResults must be mapped, or excluded before calling this function -- not silently skipped here.', ...
            strjoin(missing, ', ')));
    end
    groupOf = @(s) string(gm(char(s)));
    subGroups = arrayfun(groupOf, subsAll);

    okA = okMask & ismember(string({batchResults.subLabel}), subsAll(subGroups == groupA));
    okB = okMask & ismember(string({batchResults.subLabel}), subsAll(subGroups == groupB));
    if ~any(okA)
        error('tier3DesignPrecision_v001:NoTrialsGroupA', '%s', sprintf('No successful trials in group "%s".', groupA));
    end
    if ~any(okB)
        error('tier3DesignPrecision_v001:NoTrialsGroupB', '%s', sprintf('No successful trials in group "%s".', groupB));
    end
    [pA, sA] = extractPointsAndSubjects_local(batchResults(okA), options.Pipeline);
    [pB, sB] = extractPointsAndSubjects_local(batchResults(okB), options.Pipeline);

    rng(options.RngSeed, 'twister');
    armA = singleDatasetPrecision_local(pA, sA, cfg, groupA, options.Pipeline);
    armB = singleDatasetPrecision_local(pB, sB, cfg, groupB, options.Pipeline);

    rng(options.RngSeed, 'twister');
    contrast = independentContrast_local(pA, sA, pB, sB, cfg);

    tier3 = struct();
    tier3.mode = "independent";
    tier3.groups = [groupA, groupB];
    tier3.armA = armA; tier3.armB = armB;
    tier3 = mergeStruct_local(tier3, contrast);
    tier3 = addDesignNote_local(tier3, armA, armB, groupA, groupB);
end

printSummary_local(tier3, options.MDC);

end

% ======================================================================== %
%  Extraction: this tool's own batch shape -> (point, subjectID) vectors   %
% ======================================================================== %

function [pt, sub] = extractPointsAndSubjects_local(okBatch, pipelineChoice)
    N = numel(okBatch);
    pt = struct('val', nan(N,1), 'sd', nan(N,1));
    sub = strings(N,1);
    for i = 1:N
        tr = okBatch(i).trial;
        sub(i) = string(okBatch(i).subLabel);
        if pipelineChoice == "median"
            pt.val(i) = median(tr.betaGenStar, 'omitnan');
            ciW = tr.ciHi - tr.ciLo;                       % 1x6
            % Same treatment as loopClosureVarDecomp_v012_gated.m's
            % constellation-median path: a per-trial SD for the median
            % is derived by inner resampling across the six pipelines'
            % own point estimates and CI-derived SDs, not a plain mean.
            sdPP = ciW / (2*1.644853626951472);
            pts6 = tr.betaGenStar; conv = isfinite(pts6);
            if ~any(conv)
                pt.sd(i) = 0;
            else
                p6 = pts6(conv); d6 = sdPP(conv); d6(~isfinite(d6)|d6<0) = 0;
                if all(d6 == 0)
                    pt.sd(i) = 0;
                else
                    draws = p6(:)' + d6(:)' .* randn(200, numel(p6));
                    pt.sd(i) = std(median(draws, 2), 0);
                end
            end
        else
            pp = find(tr.pipelineLabels == pipelineChoice, 1);
            if isempty(pp)
                error('tier3DesignPrecision_v001:UnknownPipeline', '%s', sprintf( ...
                    'Pipeline "%s" not recognised. Use "median" or one of: %s', ...
                    pipelineChoice, strjoin(tr.pipelineLabels, ', ')));
            end
            pt.val(i) = tr.betaGenStar(pp);
            pt.sd(i)  = (tr.ciHi(pp) - tr.ciLo(pp)) / (2*1.644853626951472);
        end
    end
    % No further floor is applied at this level, unlike
    % floorMedianSD_local's empirical-spread-across-pipelines floor:
    % single-value-per-trial input here (this tool reports one point per
    % trial already, unlike the main project's per-pipeline matrix) --
    % the "spread across pipelines" term such a floor would add does not
    % apply here; the inner-resample sd above already reflects that
    % spread for Pipeline="median". Stated rather than silently omitted.
end

% ======================================================================== %
%  Single-dataset three-stage ablation (clusterEstimate_local, ported)     %
% ======================================================================== %

function cl = singleDatasetPrecision_local(pt, sub, cfg, label, pipelineChoice)
    fin = isfinite(pt.val);
    p = pt.val(fin); s = pt.sd(fin); s(~isfinite(s) | s < 0) = 0; sID = sub(fin);
    if isempty(p)
        error('tier3DesignPrecision_v001:NoFiniteEstimates', '%s', sprintf( ...
            'Dataset "%s" has zero trials with a finite betaGenStar for pipeline "%s".', label, pipelineChoice));
    end

    uSubj = unique(sID);
    M = numel(uSubj);
    subjTrialIdx = cell(M,1); subjN = zeros(M,1);
    for i = 1:M
        subjTrialIdx{i} = find(sID == uSubj(i));
        subjN(i) = numel(subjTrialIdx{i});
    end

    subjMed = nan(M,1);
    for i = 1:M
        subjMed(i) = median(p(subjTrialIdx{i}), 'omitnan');
    end
    pointEst = median(subjMed, 'omitnan');

    NBoot = cfg.NBoot;
    stage0 = nan(NBoot,1); stage1 = nan(NBoot,1); stage2 = nan(NBoot,1);

    for bi = 1:NBoot
        pert0 = p + s .* randn(numel(p),1);
        s0med = nan(M,1);
        for i = 1:M, s0med(i) = median(pert0(subjTrialIdx{i}), 'omitnan'); end
        stage0(bi) = median(s0med, 'omitnan');

        s1med = nan(M,1);
        for i = 1:M
            idxI = subjTrialIdx{i}; nI = numel(idxI); draw = idxI(randi(nI,nI,1));
            pert1 = p(draw) + s(draw) .* randn(nI,1);
            s1med(i) = median(pert1, 'omitnan');
        end
        stage1(bi) = median(s1med, 'omitnan');

        subjDraw = randi(M, M, 1); s2med = nan(M,1);
        for i = 1:M
            idxI = subjTrialIdx{subjDraw(i)}; nI = numel(idxI); draw = idxI(randi(nI,nI,1));
            pert2 = p(draw) + s(draw) .* randn(nI,1);
            s2med(i) = median(pert2, 'omitnan');
        end
        stage2(bi) = median(s2med, 'omitnan');
    end

    cl = struct();
    cl.label = label; cl.pipeline = pipelineChoice;
    cl.M = M; cl.subjN_min = min(subjN); cl.subjN_med = median(subjN); cl.subjN_max = max(subjN);
    cl.N_trials = numel(p);
    cl.pointEst = pointEst;
    cl.seStage0 = std(stage0, 0, 'omitnan');
    cl.seStage1 = std(stage1, 0, 'omitnan');
    cl.seStage2 = std(stage2, 0, 'omitnan');

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

    if M >= 2
        tcrit = tinv(0.975, M-1);
        cl.dfT = M-1;
        cl.ciT = [pointEst - tcrit*cl.seStage2, pointEst + tcrit*cl.seStage2];
    else
        cl.dfT = NaN; cl.ciT = [NaN NaN];
        warning('tier3DesignPrecision_v001:SingleSubject', '%s', sprintf( ...
            'Dataset "%s" has only one subject -- between-subject variance and the t(M-1) CI are undefined (df=0). seStage2 still reported; treat it as an inversion+within-subject figure only.', label));
    end

    cl.MDC = cfg.MDC;
    cl.semAdequate = cl.seStage2 < (cfg.MDC / 2.77);
end

% ======================================================================== %
%  Paired one-sample bootstrap (pairedOneEstimator_local, ported)          %
% ======================================================================== %

function res = pairedContrast_local(pA, sA, pB, sB, sharedSubj, cfg)
    finA = isfinite(pA.val); vA = pA.val(finA); dA = pA.sd(finA); dA(~isfinite(dA)|dA<0)=0; sidA = sA(finA);
    finB = isfinite(pB.val); vB = pB.val(finB); dB = pB.sd(finB); dB(~isfinite(dB)|dB<0)=0; sidB = sB(finB);

    M = numel(sharedSubj);
    idxA = cell(M,1); idxB = cell(M,1);
    for i = 1:M
        idxA{i} = find(sidA == sharedSubj(i));
        idxB{i} = find(sidB == sharedSubj(i));
    end
    if any(cellfun(@isempty, idxA)) || any(cellfun(@isempty, idxB))
        error('tier3DesignPrecision_v001:MissingSharedTrials', '%s', ...
            'A shared subject has zero finite trials in one condition -- cannot pair.');
    end

    subjA_pt = nan(M,1); subjB_pt = nan(M,1);
    for i = 1:M
        subjA_pt(i) = median(vA(idxA{i}), 'omitnan');
        subjB_pt(i) = median(vB(idxB{i}), 'omitnan');
    end
    subjDiff = subjA_pt - subjB_pt;
    pointEstDiff = median(subjDiff, 'omitnan');

    NBoot = cfg.NBoot;
    diffBoot = nan(NBoot,1);
    for bi = 1:NBoot
        subjDraw = randi(M, M, 1); dBoot = nan(M,1);
        for i = 1:M
            src = subjDraw(i);
            iiA = idxA{src}; nA = numel(iiA); ddA = iiA(randi(nA,nA,1));
            iiB = idxB{src}; nB = numel(iiB); ddB = iiB(randi(nB,nB,1));
            mA = median(vA(ddA) + dA(ddA).*randn(nA,1), 'omitnan');
            mB = median(vB(ddB) + dB(ddB).*randn(nB,1), 'omitnan');
            dBoot(i) = mA - mB;
        end
        diffBoot(bi) = median(dBoot, 'omitnan');
    end

    res = finishContrast_local(pointEstDiff, diffBoot, M-1, cfg.MDC);
    res.M = M;
end

% ======================================================================== %
%  Independent-groups two-sample bootstrap (compareIndependent_local)      %
% ======================================================================== %

function res = independentContrast_local(pA, sA, pB, sB, cfg)
    finA = isfinite(pA.val); vA = pA.val(finA); dA = pA.sd(finA); dA(~isfinite(dA)|dA<0)=0; sidA = sA(finA);
    finB = isfinite(pB.val); vB = pB.val(finB); dB = pB.sd(finB); dB(~isfinite(dB)|dB<0)=0; sidB = sB(finB);

    uA = unique(sidA); Ma = numel(uA);
    uB = unique(sidB); Mb = numel(uB);
    idxA = cell(Ma,1); for i=1:Ma, idxA{i} = find(sidA==uA(i)); end
    idxB = cell(Mb,1); for i=1:Mb, idxB{i} = find(sidB==uB(i)); end

    subjMedA = nan(Ma,1); for i=1:Ma, subjMedA(i) = median(vA(idxA{i}), 'omitnan'); end
    subjMedB = nan(Mb,1); for i=1:Mb, subjMedB(i) = median(vB(idxB{i}), 'omitnan'); end
    pointEstA = median(subjMedA, 'omitnan'); pointEstB = median(subjMedB, 'omitnan');
    diffPoint = pointEstA - pointEstB;

    NBoot = cfg.NBoot;
    diffBoot = nan(NBoot,1); aBoot = nan(NBoot,1); bBoot = nan(NBoot,1);
    for bi = 1:NBoot
        drawA = randi(Ma, Ma, 1); mA = nan(Ma,1);
        for i = 1:Ma
            ii = idxA{drawA(i)}; nI = numel(ii); dd = ii(randi(nI,nI,1));
            mA(i) = median(vA(dd) + dA(dd).*randn(nI,1), 'omitnan');
        end
        drawB = randi(Mb, Mb, 1); mB = nan(Mb,1);
        for i = 1:Mb
            ii = idxB{drawB(i)}; nI = numel(ii); dd = ii(randi(nI,nI,1));
            mB(i) = median(vB(dd) + dB(dd).*randn(nI,1), 'omitnan');
        end
        aBoot(bi) = median(mA, 'omitnan'); bBoot(bi) = median(mB, 'omitnan');
        diffBoot(bi) = aBoot(bi) - bBoot(bi);
    end

    seA = std(aBoot, 0, 'omitnan'); seB = std(bBoot, 0, 'omitnan');
    dfWelch = (seA^2 + seB^2)^2 / ( (seA^4)/(Ma-1) + (seB^4)/(Mb-1) );

    res = finishContrast_local(diffPoint, diffBoot, dfWelch, cfg.MDC);
    res.Ma = Ma; res.Mb = Mb;
end

function res = finishContrast_local(pointEstDiff, diffBoot, df, MDC)
    seDiff = std(diffBoot, 0, 'omitnan');
    tcrit95 = tinv(0.975, df); tcrit90 = tinv(0.95, df);
    ci95t = [pointEstDiff - tcrit95*seDiff, pointEstDiff + tcrit95*seDiff];
    ci90t = [pointEstDiff - tcrit90*seDiff, pointEstDiff + tcrit90*seDiff];

    diffExcludesZero = (ci95t(1) > 0) || (ci95t(2) < 0);
    equivalentWithinMDC = (ci90t(1) > -MDC) && (ci90t(2) < MDC);

    if diffExcludesZero
        verdict = "REAL DIFFERENCE (95% CI excludes 0)";
    elseif equivalentWithinMDC
        verdict = "EQUIVALENT within MDC (90% TOST CI inside +/-MDC)";
    else
        verdict = "INCONCLUSIVE -- neither excludes 0 nor fits within +/-MDC";
    end

    res = struct('diffPoint', pointEstDiff, 'seDiff', seDiff, 'df', df, ...
        'ci95_t', ci95t, 'ci90_t', ci90t, ...
        'diffExcludesZero', diffExcludesZero, 'equivalentWithinMDC', equivalentWithinMDC, ...
        'verdict', verdict);

    if ~equivalentWithinMDC
        % Required SE for equivalence ASSUMING the true difference is
        % zero -- the same calculation, and the same stated assumption,
        % as the manuscript's own General Discussion gives for Cook
        % CTRL-vs-ASD. Not a claim about detecting a real difference.
        res.requiredSEForEquivalence = MDC / tcrit90;
    end
end

% ======================================================================== %
%  Design-rule note: only offered when the numbers actually support it     %
% ======================================================================== %

function tier3 = addDesignNote_local(tier3, armA, armB, labelA, labelB)
    if isfield(tier3, 'equivalentWithinMDC') && ~tier3.equivalentWithinMDC
        fBSA = armA.fracBetweenSubject; fBSB = armB.fracBetweenSubject;
        if isfinite(fBSA) && isfinite(fBSB) && fBSA > 0.5 && fBSB > 0.5
            tier3.designNote = sprintf([ ...
                'Between-subject variance is %.0f%%-%.0f%% of the total in %s/%s ' ...
                'individually -- the dominant term. A paired or within-subject ' ...
                'design cancels this directly in the contrast; an independent-' ...
                'groups design of the same size can only dilute it by averaging ' ...
                'over participants. If a within-subject structure is feasible for ' ...
                'this comparison, it is the highest-leverage design choice available.'], ...
                100*fBSA, 100*fBSB, labelA, labelB);
        end
    end
end

% ======================================================================== %
%  Small helpers                                                           %
% ======================================================================== %

function a = mergeStruct_local(a, b)
    fn = fieldnames(b);
    for k = 1:numel(fn), a.(fn{k}) = b.(fn{k}); end
end

function printSummary_local(tier3, MDC)
    fprintf('\n=== tier3DesignPrecision_v001 (%s) ===\n', tier3.mode);
    if tier3.mode == "single"
        printArm_local(tier3, MDC);
        return
    end
    if tier3.mode == "paired"
        fprintf('Conditions: %s vs %s (M=%d shared subjects)\n', tier3.conditions(1), tier3.conditions(2), tier3.M);
    else
        fprintf('Groups: %s (Ma=%d) vs %s (Mb=%d)\n', tier3.groups(1), tier3.armA.M, tier3.groups(2), tier3.armB.M);
    end
    fprintf('-- Arm A --\n'); printArm_local(tier3.armA, MDC);
    fprintf('-- Arm B --\n'); printArm_local(tier3.armB, MDC);
    fprintf('Contrast: diff=%.4f  SE=%.4f  df=%.1f\n', tier3.diffPoint, tier3.seDiff, tier3.df);
    fprintf('  95%% CI [%.4f, %.4f]   90%% CI (TOST vs +/-%.3f) [%.4f, %.4f]\n', ...
        tier3.ci95_t(1), tier3.ci95_t(2), MDC, tier3.ci90_t(1), tier3.ci90_t(2));
    fprintf('  VERDICT: %s\n', tier3.verdict);
    if isfield(tier3, 'requiredSEForEquivalence')
        fprintf('  Required SE for equivalence (true diff assumed zero): %.4f\n', tier3.requiredSEForEquivalence);
    end
    if isfield(tier3, 'designNote')
        fprintf('  %s\n', tier3.designNote);
    end
end

function printArm_local(cl, MDC)
    fprintf('  %-20s M=%-4d N_trials=%-5d pointEst=%.4f\n', cl.label, cl.M, cl.N_trials, cl.pointEst);
    fprintf('    seStage2 (cluster SE) = %.4f   %s MDC/2.77=%.4f\n', ...
        cl.seStage2, ternary_local(cl.semAdequate, '<', '>='), MDC/2.77);
    if isfinite(cl.fracInversion)
        fprintf('    Variance shares: inversion=%.0f%%  within-subject=%.0f%%  between-subject=%.0f%%\n', ...
            100*cl.fracInversion, 100*cl.fracWithinSubject, 100*cl.fracBetweenSubject);
    end
end

function out = ternary_local(cond, a, b)
    if cond, out = a; else, out = b; end
end
