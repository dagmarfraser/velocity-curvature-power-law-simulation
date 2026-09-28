function results = zarandiLMLSFoldDiagnostic_v001(opts)
% zarandiLMLSFoldDiagnostic_v001  Per-pipeline betaObs/betaGenStar/rise-
% branch diagnostic underlying Finding #206.
%
% Shows Zarandi's BWFD-LMLS/SG-LMLS pipelines infer a generator (beta_gen*)
% at or past the established forward-map fold (beta_gen ~ 0.53), while
% BWFD-OLS/BWFD-IRLS infer a generator well below it from essentially the
% same observed beta -- and that this is specific to Zarandi+LMLS, not a
% general LMLS weakness (Cook CTRL, Fraser, Hickman HALO all infer LMLS
% generators clear of the fold). This is the mechanism behind Finding
% #205's LOPO failure at Zarandi's LMLS pipelines specifically.
%
% METHOD: loads loopClosureResults_<dataset>_all_<NoiseModel>_<Version>.mat
% for each requested dataset, reads per-trial betaObs/betaGenStar/
% invertStatus, and reports the median of each per pipeline, plus the
% fraction of trials landing on the "rise" branch (rather than "desc" or
% "fail") for each.
%
% USAGE:
%   R = zarandiLMLSFoldDiagnostic_v001()
%   R = zarandiLMLSFoldDiagnostic_v001(Datasets=["Zarandi","Fraser"])
%
% Sanity check: warns (does not hard-fail; the table is still returned) if
% Zarandi's BWFD-LMLS median betaGenStar does not reproduce Finding #206's
% documented 0.5551 (tolerance 1e-3).
%
% Fraser, D.S. (2026)
% See also: Finding #205 (the LOPO result this explains), Finding #207
% (the companion coverage diagnostic from the same session)

    arguments
        opts.Datasets (1,:) string = ["Zarandi","Cook_CTRL","Fraser","Hickman_HALO"]
        opts.NoiseModel (1,1) string = "shaped_xu"
        opts.Version (1,1) string = "v015"
        opts.Save (1,1) logical = true
    end

    srcDir = fileparts(mfilename('fullpath'));
    nDS = numel(opts.Datasets);
    rowsOut = {};

    for k = 1:nDS
        ds = opts.Datasets(k);
        f = fullfile(srcDir, sprintf('loopClosureResults_%s_all_%s_%s.mat', ds, opts.NoiseModel, opts.Version));
        if ~isfile(f)
            error('zarandiLMLSFold:MissingFile', '%s', sprintf('%s not found.', f));
        end
        S = load(f, 'results', 'pipelineLabels');
        r = S.results;
        PP = S.pipelineLabels;
        nP = numel(PP);
        n = numel(r);

        betaObs = nan(n, nP); betaGenStar = nan(n, nP); invStat = strings(n, nP);
        for i = 1:n
            betaObs(i,:) = r(i).betaObs;
            betaGenStar(i,:) = r(i).betaGenStar(:).';
            invStat(i,:) = r(i).invertStatus;
        end

        for p = 1:nP
            rowsOut(end+1,:) = {ds, PP(p), median(betaObs(:,p),'omitnan'), ...
                median(betaGenStar(:,p),'omitnan'), ...
                100*mean(invStat(:,p)=="fail"), 100*mean(invStat(:,p)=="rise")}; %#ok<AGROW>
        end
    end

    results = cell2table(rowsOut, 'VariableNames', ...
        {'dataset','pipeline','medObs','medGenStar','failPct','risePct'});

    fprintf('=== zarandiLMLSFoldDiagnostic_v001 (Finding #206) ===\n');
    fprintf('%-14s %-10s %8s %10s %8s %8s\n', 'dataset','pipeline','medObs','medGenStar','fail%','rise%');
    for i = 1:height(results)
        fprintf('%-14s %-10s %8.4f %10.4f %7.1f%% %7.1f%%\n', results.dataset(i), results.pipeline(i), ...
            results.medObs(i), results.medGenStar(i), results.failPct(i), results.risePct(i));
    end

    %% Sanity check against Finding #206's documented figure
    chk = results(results.dataset=="Zarandi" & results.pipeline=="BWFD-LMLS", :);
    if isempty(chk)
        warning('zarandiLMLSFold:SanityCheckSkipped', '%s', ...
            'Zarandi/BWFD-LMLS not in requested Datasets -- sanity check against Finding #206 skipped.');
    elseif abs(chk.medGenStar - 0.5551) > 1e-3
        warning('zarandiLMLSFold:SanityCheckFailed', ...
            'Zarandi BWFD-LMLS medGenStar %.4f does not match Finding #206''s documented 0.5551 (tol 1e-3). Source v015 mat may have changed -- do not cite this run under Finding #206 without checking why.', chk.medGenStar);
    else
        fprintf('\nSanity check PASSED: Zarandi BWFD-LMLS medGenStar %.4f matches Finding #206 (0.5551) within tolerance.\n', chk.medGenStar);
    end

    %% Save
    if opts.Save
        resDir = fullfile(srcDir, 'results');
        if ~exist(resDir, 'dir'), mkdir(resDir); end
        matOut = fullfile(resDir, 'zarandiLMLSFoldDiagnostic_v001.mat');
        save(matOut, 'results', 'opts', '-v7.3');
        fprintf('\nSaved: %s\n', matOut);
    end
end
