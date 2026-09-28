function out = checkTier3AgainstMainProject_v001(opts)
% CHECKTIER3AGAINSTMAINPROJECT_V001  Cross-check tier3DesignPrecision_v001
% against the main project's own frozen contrast scripts, on the SAME
% already-computed data -- a code-correctness check, not new analysis.
%
% Reads loopClosureResults_Hickman_PLAC/HALO_all_shaped_xu_v015.mat and
% loopClosureResults_Cook_CTRL/ASD_all_shaped_xu_v015.mat directly (the
% main project's own saved per-trial results, NOT through Motion-BIDS
% ingest -- these datasets have never been converted to Motion-BIDS and
% this check does not require that), reshapes each into the
% (point,sd,subjectID) form tier3DesignPrecision_v001's own extraction
% step produces, and calls the paired/independent internals directly.
% Compares diffPoint/seDiff/verdict against
% compareHickmanPlacHalo_v002.m's and compareCookCtrlAsd_v002.m's own
% saved .mat output at the same RngSeed. diffPoint is deterministic (no
% RNG involved) and should match exactly; seDiff need not match to the
% last digit -- this file's own inline comment at the Hickman comparison
% explains why a reference script's incidental prior RNG consumption can
% offset the bootstrap stream without indicating an algorithmic error.
% It says nothing about Motion-BIDS ingest, which this check
% deliberately bypasses.
%
% REQUIRES: run from a working directory (or with srcDir set) where the
% main project's src/ .mat files are reachable -- this is a one-off
% diagnostic, not part of the tool itself, and is not called by anything
% in tool/.
%
% Fraser, D.S. (2026)

arguments
    opts.SrcDir (1,1) string = "."
    opts.Pipeline (1,1) string = "SG-LMLS"   % matches compareHickmanPlacHalo_v002.m's PrimaryPP=4
end

fprintf('=== Hickman PLAC vs HALO (paired) ===\n');
[pP, sP] = loadMainProjectResults_local(opts.SrcDir, "Hickman_PLAC", "v015", opts.Pipeline);
[pH, sH] = loadMainProjectResults_local(opts.SrcDir, "Hickman_HALO", "v015", opts.Pipeline);
shared = intersect(unique(sP), unique(sH));
fprintf('  Shared subjects: %d\n', numel(shared));

% Direct calls to tier3DesignPrecision_v001's own internal functions are
% not possible (MATLAB local functions are not exported) -- this check
% instead reuses the same public entry point via a synthetic batchResults
% array shaped exactly like runToolLoopClosure_v002's own output, one
% pseudo-trial per real trial, sesLabel="PLAC"/"HALO" carrying the
% condition. betaGenStar/ciLo/ciHi are copied straight from the main
% project's own saved fields, so no new numbers are computed here beyond
% what the main project already produced. subLabel is the RAW subject ID
% here (PrefixSubLabel=false) -- pairing depends on the SAME subLabel
% appearing in both arms for a given real subject.
batchPH = buildSyntheticBatch_local(pP, sP, "PLAC", pH, sH, "HALO", PrefixSubLabel=false);
t3PH = tier3DesignPrecision_v001(batchPH, Mode="paired", Conditions=["PLAC","HALO"], Pipeline=opts.Pipeline);

ref = load(fullfile(opts.SrcDir, "compareHickmanPlacHalo_v002.mat"), "out"); ref = ref.out;
if opts.Pipeline == "SG-LMLS"
    refBlock = ref.primary.paired;
else
    refBlock = ref.median.paired;
end
% seDiff is not expected to match to the last digit even at the same
% RngSeed: compareHickmanPlacHalo_v002.m's own extractDataset_local
% (via extractDataset_local -> innerMedianSD_local) unconditionally
% consumes RNG draws computing a constellation-median SD for BOTH PLAC
% and HALO before this paired bootstrap ever runs -- even when, as here,
% Pipeline="SG-LMLS" means that median-SD value is never used. This
% function's own extraction skips that irrelevant computation for a
% named pipeline, so its bootstrap loop reads from a different offset
% into the same seeded stream. diffPoint has no such dependency (it is
% a plain median of real data, no RNG) and should still match exactly.
reportDiff_local("Hickman PLAC-HALO (paired)", t3PH.diffPoint, t3PH.seDiff, refBlock.diffPoint, refBlock.seDiff);

fprintf('\n=== Cook CTRL vs ASD (independent) ===\n');
[pC, sC] = loadMainProjectResults_local(opts.SrcDir, "Cook_CTRL", "v015", opts.Pipeline);
[pA, sA] = loadMainProjectResults_local(opts.SrcDir, "Cook_ASD", "v015", opts.Pipeline);
gm = containers.Map('KeyType','char','ValueType','char');
for i = 1:numel(unique(sC)), u = unique(sC); gm(char("CTRL_"+u(i))) = 'CTRL'; end
for i = 1:numel(unique(sA)), u = unique(sA); gm(char("ASD_"+u(i))) = 'ASD'; end
% Cook subLabel IS arm-prefixed here (PrefixSubLabel=true): CTRL/ASD are
% different real people, and their raw subject-ID namespaces have not
% been confirmed disjoint, so GroupMap keys must not collide across arms.
batchCA = buildSyntheticBatch_local(pC, sC, "CTRL", pA, sA, "ASD", PrefixSubLabel=true);
t3CA = tier3DesignPrecision_v001(batchCA, Mode="independent", GroupMap=gm, Groups=["CTRL","ASD"], Pipeline=opts.Pipeline);

refCA = load(fullfile(opts.SrcDir, "compareCookCtrlAsd_v002.mat"), "out"); refCA = refCA.out;
if opts.Pipeline == "SG-LMLS"
    refBlockCA = refCA.primary;
else
    refBlockCA = refCA.median;
end
reportDiff_local("Cook CTRL-ASD (independent)", t3CA.diffPoint, t3CA.seDiff, refBlockCA.diffPoint, refBlockCA.seDiff);

out = struct('hickman', t3PH, 'cook', t3CA);
end

function [pt, sub] = loadMainProjectResults_local(srcDir, tag, version, pipeline)
    f = fullfile(srcDir, sprintf('loopClosureResults_%s_all_shaped_xu_%s.mat', tag, version));
    if ~isfile(f)
        error('checkTier3AgainstMainProject_v001:NotFound', '%s', sprintf('FAILED PATH: %s', f));
    end
    D = load(f, 'results'); r = D.results(:);
    N = numel(r);
    Z90 = 1.644853626951472;
    ppNames = ["BWFD-OLS","SG-OLS","BWFD-LMLS","SG-LMLS","BWFD-IRLS","SG-IRLS"];
    pp = find(ppNames == pipeline, 1);
    pt = struct('val', nan(N,1), 'sd', nan(N,1));
    sub = strings(N,1);
    for i = 1:N
        pt.val(i) = getScalar_local(r(i), 'betaGenStar', pp);
        ciHi = getScalar_local(r(i), 'ciHi', pp); ciLo = getScalar_local(r(i), 'ciLo', pp);
        pt.sd(i) = (ciHi - ciLo) / (2*Z90);
        sub(i) = string(r(i).subjectID);
    end
end

function v = getScalar_local(s, fld, idx)
    v = NaN;
    if isfield(s, fld)
        x = s.(fld);
        if isnumeric(x) && numel(x) >= idx, v = double(x(idx)); end
    end
end

function batch = buildSyntheticBatch_local(pA, sA, labelA, pB, sB, labelB, opts)
% Embeds each arm's (already pipeline-selected) point/sd pair into a
% single fixed column (SG-LMLS) of an otherwise-NaN six-pipeline trial
% struct, matching the shape processTrialLoopClosure_v001 actually
% produces. Only that one named column is populated -- callers must pass
% the SAME pipeline name (not "median") to tier3DesignPrecision_v001, or
% every trial will be treated as having only one finite pipeline, which
% Pipeline="median" mode handles correctly via omitnan but was not the
% comparison this diagnostic is built to run.
%
% PrefixSubLabel controls whether .subLabel is arm-prefixed:
%   false (Hickman/paired) -- .subLabel MUST be the raw subject ID,
%     unprefixed, identical across both arms for the same real subject,
%     since Mode="paired" pairs on matching .subLabel. Getting this
%     wrong makes every subject look arm-unique and Mode="paired" fail
%     loud with "No subject has trials in both..." -- exactly what the
%     first run of this file did, with PrefixSubLabel effectively always
%     true before this fix (Finding: caught by running, not by review).
%   true (Cook/independent) -- .subLabel is arm-prefixed to guarantee
%     disjoint IDs across groups for GroupMap lookup, since Cook's own
%     CTRL/ASD subject-ID namespaces have not been confirmed disjoint.
    arguments
        pA, sA, labelA, pB, sB, labelB
        opts.PrefixSubLabel (1,1) logical
    end
    NA = numel(pA.val); NB = numel(pB.val);
    ppNames = ["BWFD-OLS","SG-OLS","BWFD-LMLS","SG-LMLS","BWFD-IRLS","SG-IRLS"];
    idxCol = find(ppNames == "SG-LMLS", 1);
    Z90 = 1.644853626951472;

    if opts.PrefixSubLabel
        subA = labelA + "_" + sA; subB = labelB + "_" + sB;
    else
        subA = sA; subB = sB;
    end

    batch = repmat(struct('ok', true, 'subLabel', "", 'sesLabel', "", 'runLabel', "", 'trial', []), NA+NB, 1);
    for i = 1:NA
        batch(i) = oneSyntheticTrial_local(subA(i), labelA, i, ppNames, idxCol, pA.val(i), pA.sd(i), Z90);
    end
    for i = 1:NB
        batch(NA+i) = oneSyntheticTrial_local(subB(i), labelB, i, ppNames, idxCol, pB.val(i), pB.sd(i), Z90);
    end
end

function b = oneSyntheticTrial_local(subLabel, sesLabel, runNum, ppNames, idxCol, val, sd, Z90)
    b.ok = true;
    b.subLabel = subLabel;
    b.sesLabel = sesLabel;
    b.runLabel = string(runNum);
    b.trial = struct('pipelineLabels', ppNames, 'betaGenStar', nan(1,6), 'ciLo', nan(1,6), 'ciHi', nan(1,6));
    b.trial.betaGenStar(idxCol) = val;
    b.trial.ciLo(idxCol) = val - sd*Z90;
    b.trial.ciHi(idxCol) = val + sd*Z90;
end

function reportDiff_local(label, myDiff, mySE, refDiff, refSE)
    fprintf('  %s\n', label);
    fprintf('    this function : diff=%.6f  SE=%.6f\n', myDiff, mySE);
    fprintf('    reference     : diff=%.6f  SE=%.6f\n', refDiff, refSE);
    fprintf('    delta         : %.2e / %.2e\n', myDiff-refDiff, mySE-refSE);
end
