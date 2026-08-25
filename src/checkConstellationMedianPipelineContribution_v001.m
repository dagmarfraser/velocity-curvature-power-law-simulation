function report = checkConstellationMedianPipelineContribution_v001()
% CHECKCONSTELLATIONMEDIANPIPELINECONTRIBUTION_V001
%
% Tests directly whether particular pipelines (OLS in particular, given
% its known compression bias) disproportionately determine the
% constellation median beta_gen* reported in Table 6b, rather than
% inferring this from resolution-rate/coverage numbers alone.
%
% For each trial: betaGenStar has up to 6 values, some NaN (that
% pipeline's local neighbourhood wasn't invertible for this trial, per
% runLoopClosureFftnoise_v012.m's own fail-loud NaN convention -- not a
% missing-data problem, a correct exclusion). Of the finite ones, the
% constellation median is the median of whichever subset resolved. This
% script identifies, per trial, WHICH pipeline(s) that finite subset's
% median value actually comes from (the middle one for an odd count, the
% two averaged ones for an even count), and tallies that across all
% trials -- "resolution rate" (how often a pipeline is finite at all) vs
% "median-determining rate" (how often its value is one the median is
% actually built from) are reported separately, since they need not
% match.
%
% Fraser, D.S. (2026)

    srcDir = fileparts(mfilename('fullpath'));
    cd(srcDir);

    PP_NAMES = ["BWFD-OLS","SG-OLS","BWFD-LMLS","SG-LMLS","BWFD-IRLS","SG-IRLS"];
    nPP = 6;

    registry = struct( ...
        'label', {'Fraser','Cook ASD','Cook CTRL','Hickman HALO','Hickman PLAC','Zarandi'}, ...
        'file',  {'loopClosureResults_Fraser_all_shaped_xu_v008.mat', ...
                   'loopClosureResults_Cook_ASD_all_shaped_xu_v007.mat', ...
                   'loopClosureResults_Cook_CTRL_all_shaped_xu_v007.mat', ...
                   'loopClosureResults_Hickman_HALO_all_shaped_xu_v007.mat', ...
                   'loopClosureResults_Hickman_PLAC_all_shaped_xu_v007.mat', ...
                   'loopClosureResults_Zarandi_all_shaped_xu_v007.mat'});

    nDS = numel(registry);
    dsRows = cell(nDS, 1);

    fprintf('%s\n', repmat('=', 1, 118));
    fprintf('PER-PIPELINE: resolution rate (finite at all) vs median-determining rate (contributes to the reported median)\n');
    fprintf('%s\n', repmat('=', 1, 118));

    for d = 1:nDS
        matFile = fullfile(srcDir, registry(d).file);
        if ~isfile(matFile)
            error('checkConstellationMedianPipelineContribution:notFound', 'FAILED PATH: %s', matFile);
        end
        S = load(matFile, 'results');
        R = S.results(:);
        N = numel(R);

        bAll = nan(N, nPP);
        for i = 1:N
            v = R(i).betaGenStar;
            if isnumeric(v) && numel(v) == nPP
                bAll(i,:) = double(v(:)');
            end
        end

        resolvedCount = sum(isfinite(bAll), 2);   % pipelines resolving per trial
        resolveRate   = mean(isfinite(bAll), 1);  % per-pipeline resolution rate

        % Median-determining tally
        medianCredit = zeros(1, nPP);
        nContributingTrials = 0;
        for i = 1:N
            row = bAll(i, :);
            finIdx = find(isfinite(row));
            k = numel(finIdx);
            if k == 0
                continue
            end
            nContributingTrials = nContributingTrials + 1;
            [sortedVals, ord] = sort(row(finIdx));
            sortedPipes = finIdx(ord);
            if mod(k, 2) == 1
                midPipe = sortedPipes((k+1)/2);
                medianCredit(midPipe) = medianCredit(midPipe) + 1;
            else
                p1 = sortedPipes(k/2);
                p2 = sortedPipes(k/2 + 1);
                medianCredit(p1) = medianCredit(p1) + 0.5;
                medianCredit(p2) = medianCredit(p2) + 0.5;
            end
        end
        medianDetermRate = medianCredit / max(nContributingTrials, 1);

        fprintf('\n--- %s (N=%d trials, %d contributing at least one pipeline) ---\n', ...
            registry(d).label, N, nContributingTrials);
        fprintf('%-12s', 'Pipeline');
        for p = 1:nPP, fprintf(' %10s', PP_NAMES(p)); end
        fprintf('\n%-12s', 'Resolve %%');
        for p = 1:nPP, fprintf(' %9.1f%%', 100*resolveRate(p)); end
        fprintf('\n%-12s', 'MedianDet %%');
        for p = 1:nPP, fprintf(' %9.1f%%', 100*medianDetermRate(p)); end
        fprintf('\n');

        % Flag disproportion: median-determining rate vs resolve rate, and
        % OLS-family vs LMLS/IRLS-family aggregate share
        olsShare   = sum(medianCredit([1 2])) / max(nContributingTrials, 1);
        lmlsShare  = sum(medianCredit([3 4])) / max(nContributingTrials, 1);
        irlsShare  = sum(medianCredit([5 6])) / max(nContributingTrials, 1);
        fprintf('  Aggregate median-determining share: OLS family %.1f%%, LMLS family %.1f%%, IRLS family %.1f%%\n', ...
            100*olsShare, 100*lmlsShare, 100*irlsShare);

        dsRows{d} = struct('dataset', string(registry(d).label), 'N', N, ...
            'nContributing', nContributingTrials, ...
            'resolveRate', resolveRate, 'medianDetermRate', medianDetermRate, ...
            'olsShare', olsShare, 'lmlsShare', lmlsShare, 'irlsShare', irlsShare);
    end

    report = struct2table([dsRows{:}]);

    fprintf('\n%s\n', repmat('=', 1, 118));
    fprintf('SUMMARY: aggregate median-determining share by pipeline family, all datasets\n');
    fprintf('%s\n', repmat('=', 1, 118));
    fprintf('%-16s %10s %10s %10s\n', 'Dataset', 'OLS %%', 'LMLS %%', 'IRLS %%');
    for d = 1:nDS
        fprintf('%-16s %9.1f%% %9.1f%% %9.1f%%\n', ...
            dsRows{d}.dataset, 100*dsRows{d}.olsShare, 100*dsRows{d}.lmlsShare, 100*dsRows{d}.irlsShare);
    end
    fprintf('%s\n', repmat('=', 1, 118));

end
