function constellationFraser_v002(opts)
% constellationFraser_v002  v001 on the N=48 corpus (2026-11-04). ONLY CHANGES: importDB_fraser_v002 (ID101/102/214 out,
% ID401-403 in, gaps split, direct touches dropped), noiseCharacterisation_fraserN48.mat, SIGMA_TO_MM = 1/ipadScaleM4_v001()
% (UIKit points; console summary only, the constellation itself is scale-free: betas, and VGF in native units),
% output constellationFraser_v002.mat. v001 and its mat are left as they were. v001's header follows.
% constellationFraser_v001  Six-pipeline constellation for Fraser (Analysis_2026).
%
% Brain2Bee MSc confirmatory test-retest study (Analysis_2026), the replacement
% cohort for Pilot. iPad Pro 13" (M4) -- same physical device as Pilot,
% confirmed 2026-08-06. 2 visits (~1 week apart), up to 3 sessions/visit,
% 10 ellipse trials/session (data/fraser/automatic/ is already Ellipse-only).
% 48 raw participant folders, 45 complete (6/6 sessions), 3 incomplete
% (F101 5/6, F102 5/6, F214 3/6) -- included here, not pre-filtered; whether
% to exclude them is a downstream decision (data/fraser/README.md point 2).
%
% Modelled directly on constellationPilot_v001.m with one deliberate
% simplification: this script does NOT recompute noise characterisation
% internally (Pilot's own script runs the now-superseded
% characteriseBiologicalNoise_v003 for this purpose, redundant with the
% separately-maintained noiseCharacterisation_fraser.mat which already
% holds the correct v004 result -- see rerunNoiseCharacterisation_v005.m).
% checkInvertibilityForEmpirical_v2_001.m reads noise from that separate
% file regardless of what a constellation script's own bioResults says, so
% this recompute was never load-bearing downstream. This script therefore
% only ever READS noiseCharacterisation_fraser.mat (for console summary
% context) and never writes to it -- the canonical v004 noise mat is not at
% risk of being clobbered by this script.
%
% USAGE:
%   constellationFraser_v001()                 % both visits, all subjects
%   constellationFraser_v001(Visits=1, MaxSubjects=5)   % smoke test
%
% OUTPUT FILE:
%   src/constellationFraser_v001.mat            constellation results only
%   (noiseCharacterisation_fraser.mat is READ, never written, by this script)
%
% Fraser, D.S. (2026)  v001
% See also: constellationPilot_v001, importDB_fraser_v001,
%           rerunNoiseCharacterisation_v005

    arguments
        opts.Visits      (1,:) double = [1 2]
        opts.MaxSubjects (1,1) double = Inf
    end

    srcDir = fileparts(mfilename("fullpath"));
    addpath(genpath(fullfile(srcDir, "functions")));

    FS          = 240;
    SIGMA_TO_MM = 1 / ipadScaleM4_v001();   % v002: UIKit points per mm (was 1/10.41793, a pixel density)
    %SIGMA_TO_MM_v001 = 1 / 10.41793;   % iPad Pro 13" M4, hardware-derived (data/fraser/README.md point 3)
                                  % NOT Pilot's 9.73 -- coherence-check false start corrected there.

    constellFile = fullfile(srcDir, 'constellationFraser_v002.mat');
    noiseFile    = fullfile(srcDir, 'noiseCharacterisation_fraserN48.mat');

    fprintf('=== CONSTELLATION: Fraser (iPad 240 Hz) -- Visits [%s] ===\n\n', ...
        num2str(opts.Visits));

    %% 1. Import
    dataRoot = fullfile(srcDir, '..', 'data', 'fraser');
    trials   = importDB_fraser_v002(dataRoot, ...
        Visits=opts.Visits, MaxSubjects=opts.MaxSubjects, Verbose=true);   % v002
    nTrials  = numel(trials);
    fprintf('\n');

    %% 2. Noise summary (READ ONLY -- see header note)
    if isfile(noiseFile)
        fprintf('--- Noise summary (read from existing noiseCharacterisation_fraser.mat, v004) ---\n');
        nz = load(noiseFile, 'bioResults');
        bioResults = nz.bioResults;
        fprintf('  IRASA alpha:  %.3f +/- %.3f\n', ...
            mean(bioResults.ira_alphaMean, 'omitnan'), std(bioResults.ira_alphaMean, 'omitnan'));
        fprintf('  sigma (mm):   %.3f +/- %.3f  (x %.6f mm/px)\n', ...
            mean(bioResults.sigmaMean, 'omitnan') * SIGMA_TO_MM, ...
            std(bioResults.sigmaMean,  'omitnan') * SIGMA_TO_MM, SIGMA_TO_MM);
    else
        fprintf('  WARNING: %s not found -- run rerunNoiseCharacterisation_v005(DS="Fraser") first.\n', noiseFile);
    end

    %% 3. Pipeline configuration (identical to Pilot -- same fs=240Hz device)
    derivConfigs = struct( ...
        "label",        {"BWFD",    "SG"}, ...
        "filterType",   {2,         6}, ...
        "filterParams", {[2 10 1],  [4 17]});

    regressConfigs = struct( ...
        "label", {"OLS", "LMLS", "IRLS"}, ...
        "type",  {3,     4,      5});

    canonLabels = ["BWFD-OLS","BWFD-LMLS","BWFD-IRLS","SG-OLS","SG-LMLS","SG-IRLS"];
    edgeClip    = 20;
    nDeriv      = numel(derivConfigs);
    nRegress    = numel(regressConfigs);

    %% 4. Process all trials
    betaAll = NaN(nTrials, nDeriv, nRegress);
    vgfAll  = NaN(nTrials, nDeriv, nRegress);

    fprintf('\nProcessing %d trials x %d derivations x %d regressions...\n', ...
        nTrials, nDeriv, nRegress);
    tic;

    for tIdx = 1:nTrials
        x = trials(tIdx).x;
        y = trials(tIdx).y;

        for dIdx = 1:nDeriv
            [dx, dy] = differentiateKinematicsEBR( ...
                x, y, derivConfigs(dIdx).filterType, derivConfigs(dIdx).filterParams, FS);

            vx = dx(edgeClip:end-edgeClip, 2);
            vy = dy(edgeClip:end-edgeClip, 2);
            ax = dx(edgeClip:end-edgeClip, 3);
            ay = dy(edgeClip:end-edgeClip, 3);

            speed = sqrt(vx.^2 + vy.^2);
            kappa = curvatureKinematicEBR(vx, vy, ax, ay);
            LMseeds = [1 -1/3];

            for rIdx = 1:nRegress
                try
                    [betaEst, vgfEst] = regressDataEBR( ...
                        speed, kappa, regressConfigs(rIdx).type, LMseeds, 0, 0);
                    betaAll(tIdx, dIdx, rIdx) = betaEst;
                    vgfAll(tIdx, dIdx, rIdx)  = vgfEst;
                    if rIdx == 1, LMseeds = [vgfEst betaEst]; end
                catch ME
                    fprintf('  WARN trial %d %s-%s: %s\n', tIdx, ...
                        derivConfigs(dIdx).label, regressConfigs(rIdx).label, ME.message);
                end
            end
        end

        if mod(tIdx, 200) == 0
            fprintf('  %d / %d (%.0f s elapsed)\n', tIdx, nTrials, toc);
        end
    end
    fprintf('Done in %.1f s\n\n', toc);

    %% 5. Flatten to canonical 6-pipeline array
    betaCanon = NaN(nTrials, 6);
    col = 0;
    for dIdx = 1:nDeriv
        for rIdx = 1:nRegress
            col = col + 1;
            betaCanon(:, col) = betaAll(:, dIdx, rIdx);
        end
    end

    %% 6. Pipeline summary (all trials)
    fprintf('=== PIPELINE BETA MEANS (all trials) ===\n');
    fprintf('%-14s  %7s  %7s  %7s  %5s\n', 'Pipeline', 'mean', 'std', 'median', 'N');
    fprintf('%s\n', repmat('-', 1, 50));
    for dIdx = 1:nDeriv
        for rIdx = 1:nRegress
            b = betaAll(:, dIdx, rIdx);
            b = b(~isnan(b));
            fprintf('%-14s  %7.4f  %7.4f  %7.4f  %5d\n', ...
                derivConfigs(dIdx).label + "-" + regressConfigs(rIdx).label, ...
                mean(b), std(b), median(b), numel(b));
        end
    end

    %% 7. Visit 1 vs Visit 2 comparison (paired by subject)
    fprintf('\n=== VISIT 1 vs VISIT 2 COMPARISON ===\n');
    subIDs  = string({trials.subjectID});
    visitNs = [trials.visit];
    uSubs   = unique(subIDs);
    bothVisits = false(numel(uSubs), 1);
    for s = 1:numel(uSubs)
        vv = unique(visitNs(subIDs == uSubs(s)));
        bothVisits(s) = numel(vv) == 2;
    end
    subsBoth = uSubs(bothVisits);
    fprintf('Subjects with both visits: %d / %d\n\n', sum(bothVisits), numel(uSubs));

    fprintf('%-14s  %8s  %8s  %8s  %7s  %7s\n', ...
        'Pipeline', 'V1 mean', 'V2 mean', 'V2-V1', 'V1 std', 'V2 std');
    fprintf('%s\n', repmat('-', 1, 65));

    betaDiffByPipeline = NaN(numel(subsBoth), 6);
    for p = 1:6
        betaV1 = NaN(numel(subsBoth), 1);
        betaV2 = NaN(numel(subsBoth), 1);
        for s = 1:numel(subsBoth)
            sid = subsBoth(s);
            m1  = subIDs == sid & visitNs == 1;
            m2  = subIDs == sid & visitNs == 2;
            betaV1(s) = mean(betaCanon(m1, p), 'omitnan');
            betaV2(s) = mean(betaCanon(m2, p), 'omitnan');
        end
        betaDiffByPipeline(:, p) = betaV2 - betaV1;
        fprintf('%-14s  %8.4f  %8.4f  %8.4f  %7.4f  %7.4f\n', ...
            canonLabels(p), mean(betaV1,'omitnan'), mean(betaV2,'omitnan'), ...
            mean(betaV2-betaV1,'omitnan'), std(betaV1,'omitnan'), std(betaV2,'omitnan'));
    end

    fprintf('\n--- Paired t-test: Visit 2 - Visit 1 (per-subject means) ---\n');
    fprintf('%-14s  %8s  %8s  %6s\n', 'Pipeline', 't', 'p', 'df');
    for p = 1:6
        d = betaDiffByPipeline(:, p);
        d = d(~isnan(d));
        [~, pval, ~, stats] = ttest(d);
        fprintf('%-14s  %8.4f  %8.4f  %6d\n', canonLabels(p), stats.tstat, pval, stats.df);
    end

    %% 8. BWFD vs SG divergence at 240 Hz
    fprintf('\n=== BWFD vs SG DIVERGENCE (fs=240 Hz) ===\n');
    fprintf('%-8s  %10s  %10s  %10s\n', 'Regress', 'max|diff|', 'mean|diff|', 'bit-ident%');
    for rIdx = 1:nRegress
        bw    = betaAll(:, 1, rIdx);
        sg    = betaAll(:, 2, rIdx);
        valid = ~isnan(bw) & ~isnan(sg);
        d     = abs(bw(valid) - sg(valid));
        pctId = 100 * sum(bw(valid) == sg(valid)) / sum(valid);
        fprintf('%-8s  %10.2e  %10.2e  %10.1f%%\n', ...
            regressConfigs(rIdx).label, max(d), mean(d), pctId);
    end

    %% 9. Pairwise constellation
    nPipes     = 6;
    pairLabels = strings(0);
    pairDiffs  = NaN(nTrials, nPipes*(nPipes-1)/2);
    pCol = 0;
    for p = 1:nPipes
        for q = (p+1):nPipes
            pCol = pCol + 1;
            pairDiffs(:, pCol) = betaCanon(:, p) - betaCanon(:, q);
            pairLabels(pCol)   = canonLabels(p) + " - " + canonLabels(q);
        end
    end

    fprintf('\n=== CONSTELLATION: mean pairwise differences ===\n');
    fprintf('%-40s  %8s  %8s\n', 'Pair', 'mean', 'std');
    fprintf('%s\n', repmat('-', 1, 60));
    for pCol = 1:numel(pairLabels)
        d = pairDiffs(~isnan(pairDiffs(:,pCol)), pCol);
        fprintf('%-40s  %8.4f  %8.4f\n', pairLabels(pCol), mean(d), std(d));
    end

    %% 10. Save (constellation file ONLY -- noiseCharacterisation_fraser.mat untouched)
    visitMetadata = struct('visit', {trials.visit}, ...
        'sessionOverall', {trials.sessionOverall}, ...
        'subjectID', {trials.subjectID});

    save(constellFile, ...
        'betaAll', 'vgfAll', 'betaCanon', ...
        'pairDiffs', 'pairLabels', 'canonLabels', ...
        'derivConfigs', 'regressConfigs', ...
        'visitMetadata', 'subsBoth', 'betaDiffByPipeline', ...
        'SIGMA_TO_MM', '-v7.3');

    fprintf('\nSaved: %s\n', constellFile);
    fprintf('(noiseCharacterisation_fraser.mat NOT written by this script -- read only)\n');
    fprintf('=== DONE ===\n');
end
