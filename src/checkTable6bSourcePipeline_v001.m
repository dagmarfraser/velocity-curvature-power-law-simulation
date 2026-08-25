function report = checkTable6bSourcePipeline_v001()
% CHECKTABLE6BSOURCEPIPELINE_V001  Which pipeline does Table 6b actually report?
%
% Prompted directly by reportLoopClosureVarDecompNvsN_v002's own SG-LMLS
% column not matching Table 6b's published numbers at all (deltas of
% 0.01-0.05), while the constellation-median column matched closely
% (deltas of 0.0004-0.017) -- despite Sec 7's own prose stating SG-LMLS
% is Table 6b's source pipeline. plotBetaGenStarAllDatasets_v002.m's own
% header states plainly "beta_gen* (SG-IRLS)" and hardcodes PIPE=6 (SG-
% IRLS), directly contradicting that prose.
%
% This script computes Table 6b's own stated quantities (median, mean +/-
% SD) directly from the v007 shaped_xu corpus (the same files
% plotBetaGenStarAllDatasets_v002.m reads) for THREE candidate pipelines
% -- SG-IRLS, SG-LMLS, and the constellation median across all six -- and
% prints all three against the published Table 6b numbers side by side,
% so which one Table 6b actually is can be read off directly rather than
% inferred.
%
% Fraser, D.S. (2026)

    srcDir = fileparts(mfilename('fullpath'));
    cd(srcDir);

    PP_NAMES = ["BWFD-OLS","SG-OLS","BWFD-LMLS","SG-LMLS","BWFD-IRLS","SG-IRLS"];
    SG_IRLS_IDX = 6;
    SG_LMLS_IDX = 4;

    % Table 6b's five core datasets, published values (Sec 7)
    published = struct( ...
        'label',  {'Fraser','Cook ASD','Cook CTRL','Hickman HALO','Hickman PLAC'}, ...
        'median', {0.3292,  0.294,     0.284,       0.299,          0.298}, ...
        'meanSD', {[0.3424 0.1074], [0.276 0.130], [0.283 0.130], [0.293 0.114], [0.303 0.086]});

    % v007 registry, matching plotBetaGenStarAllDatasets_v002.m exactly
    % (Fraser uses v008 there per that era's data source; using the same
    % v008 file here for Fraser specifically, v007 for the rest)
    registry = struct( ...
        'label', {'Fraser','Cook ASD','Cook CTRL','Hickman HALO','Hickman PLAC'}, ...
        'file',  {'loopClosureResults_Fraser_all_shaped_xu_v008.mat', ...
                   'loopClosureResults_Cook_ASD_all_shaped_xu_v007.mat', ...
                   'loopClosureResults_Cook_CTRL_all_shaped_xu_v007.mat', ...
                   'loopClosureResults_Hickman_HALO_all_shaped_xu_v007.mat', ...
                   'loopClosureResults_Hickman_PLAC_all_shaped_xu_v007.mat'});

    nDS = numel(registry);
    rows = cell(nDS, 1);

    fprintf('%-14s %10s | %-22s | %-22s | %-22s\n', ...
        'Dataset', 'Published', 'SG-IRLS (median|mean+-SD)', 'SG-LMLS (median|mean+-SD)', 'ConstMed (median|mean+-SD)');
    fprintf('%s\n', repmat('-', 1, 100));

    for d = 1:nDS
        matFile = fullfile(srcDir, registry(d).file);
        if ~isfile(matFile)
            error('checkTable6bSourcePipeline:notFound', 'FAILED PATH: %s', matFile);
        end
        S = load(matFile, 'results');
        R = S.results(:);
        N = numel(R);

        bAll = nan(N, 6);
        for i = 1:N
            v = R(i).betaGenStar;
            if isnumeric(v) && numel(v) == 6
                bAll(i,:) = double(v(:)');
            end
        end

        bIRLS = bAll(:, SG_IRLS_IDX);
        bLMLS = bAll(:, SG_LMLS_IDX);
        bMed  = median(bAll, 2, 'omitnan');

        medIRLS = median(bIRLS, 'omitnan'); meanIRLS = mean(bIRLS, 'omitnan'); sdIRLS = std(bIRLS, 'omitnan');
        medLMLS = median(bLMLS, 'omitnan'); meanLMLS = mean(bLMLS, 'omitnan'); sdLMLS = std(bLMLS, 'omitnan');
        medCM   = median(bMed,  'omitnan'); meanCM   = mean(bMed,  'omitnan'); sdCM   = std(bMed,  'omitnan');

        pub = published(strcmp({published.label}, registry(d).label));

        fprintf('%-14s %10.4f | %6.4f|%6.4f+-%5.4f | %6.4f|%6.4f+-%5.4f | %6.4f|%6.4f+-%5.4f\n', ...
            registry(d).label, pub.median, ...
            medIRLS, meanIRLS, sdIRLS, medLMLS, meanLMLS, sdLMLS, medCM, meanCM, sdCM);

        rows{d} = struct('dataset', string(registry(d).label), 'N', N, ...
            'publishedMedian', pub.median, 'publishedMean', pub.meanSD(1), 'publishedSD', pub.meanSD(2), ...
            'sgIrlsMedian', medIRLS, 'sgIrlsMean', meanIRLS, 'sgIrlsSD', sdIRLS, ...
            'sgLmlsMedian', medLMLS, 'sgLmlsMean', meanLMLS, 'sgLmlsSD', sdLMLS, ...
            'constMedMedian', medCM, 'constMedMean', meanCM, 'constMedSD', sdCM, ...
            'deltaMedianIRLS', abs(medIRLS - pub.median), ...
            'deltaMedianLMLS', abs(medLMLS - pub.median), ...
            'deltaMedianConstMed', abs(medCM - pub.median));
    end

    report = struct2table([rows{:}]);

    fprintf('\n%s\n', repmat('=', 1, 60));
    fprintf('MEAN ABSOLUTE DEVIATION FROM PUBLISHED MEDIAN, per candidate:\n');
    fprintf('  SG-IRLS      : %.4f\n', mean(report.deltaMedianIRLS));
    fprintf('  SG-LMLS      : %.4f\n', mean(report.deltaMedianLMLS));
    fprintf('  ConstMedian  : %.4f\n', mean(report.deltaMedianConstMed));
    fprintf('%s\n', repmat('=', 1, 60));

end
