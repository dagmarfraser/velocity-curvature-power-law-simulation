function out = queryMonotonicSlopeAtEmpirical_v005(opts)
% queryMonotonicSlopeAtEmpirical_v005  Successor to v004 and v003 (Finding #128 Part B diagnostic).
% v005 (2026-09-30) changes two things against v004: (1) Zarandi and Dhieb use their own trial
% betaGenStarMed too (flagged out of domain, #242), not 1/3; (2) the placed-trial invRiseFracPlaced is
% added (unplaced no-exponent trials counted separately as noExp) because v004's invRiseFrac counted
% them as non-invertible, an artefact. invRiseFrac (v003's denominator) is kept for the V0 anchor.
%
% v003 placed each trial on the grid's VGF axis as f0 x K_conv, K_conv = 175.2636 (beta = 1/3).
% That constant is 5.28% below the grid's true tempo constant K(1/3) = 185.02 and K falls
% steeply with beta (testGridKConv_v001). v004 runs three variants on identical inputs so that
% each change is attributable:
%   V0 RETIRED   f0 x 175.2636         regression anchor: must reproduce v003's saved dsSummary.
%   V1 CONSTANT  f0 x K(1/3)           the constant fixed only.
%   V2 PRIMARY   f0_i x K(beta_i)      per trial. In-domain datasets (Finding #242): beta_i is the
%                                      trial's betaGenStarMed; a trial without one is NOT imputed,
%                                      it is not placed (counted as noExp). Zarandi and Dhieb (outside
%                                      the validity domain) follow the same rule, flagged.
% Everything else is v003 unchanged (registry, alpha/sigma, checkTrialInvertibility_v001,
% aggregation, Part B ratios reused from Finding #128 and disclosed as reused).
%
% Alignment: v003's constellation/bioResults order is matched to the v015 corpus (source of
% betaGenStarMed) by comparing f0 element by element; any mismatch errors.
%
% Reads: src/constellation*_v001.mat, src/noiseCharacterisation_*.mat, importers,
%        src/loopClosureResults_<name>_all_shaped_xu_v015.mat, src/monotonicSegments_v2_002.mat,
%        src/results/queryMonotonicSlopeAtEmpirical_v003.mat (anchor).
% Writes: results/queryMonotonicSlopeAtEmpirical_v005.mat (project-root results/).
% USAGE:  from the project root: R = queryMonotonicSlopeAtEmpirical_v005()   (about minutes; run as
%         a background job, the importers exceed the 60 s device limit). Save=false to skip writing.
% Fraser, D.S. (2026)  v005. NOT preregistered.
    arguments
        opts.Save (1,1) logical = true
    end
    ROOT   = fileparts(fileparts(mfilename('fullpath')));
    srcDir = fullfile(ROOT, 'src');
    addpath(ROOT); addpath(srcDir); addpath(genpath(fullfile(srcDir, 'functions')));

    KCONV_RETIRED = 175.2636;  K13_PUB = 185.0242;  BETA_REF = 1/3;
    ANCHOR_TOL = 1e-9;
    K13 = gridKConv_v001(BETA_REF);
    if abs(K13 - K13_PUB) > 1e-3
        error('queryMonoSlope4:K13', '%s', sprintf('K(1/3) %.4f does not reproduce %.4f', K13, K13_PUB));
    end
    v003Mat = fullfile(srcDir, 'results', 'queryMonotonicSlopeAtEmpirical_v003.mat');
    if ~isfile(v003Mat), error('queryMonoSlope4:NoAnchor', '%s', 'FAILED PATH (anchor): ' + string(v003Mat)); end

    monoFile = fullfile(srcDir, 'monotonicSegments_v2_002.mat');
    if ~isfile(monoFile), error('queryMonoSlope4:MissingBundle', '%s', 'monotonicSegments_v2_002.mat not found.'); end
    mono = load(monoFile, 'params', 'pipeOrder', 'alphaGrid', 'sigmaGrid', 'fsGrid', 'betaGenGrid', 'VGFGrid', 'betaSeg');

    AlphaSource = "ira_alphaMean";
    registry = { ...
        struct('name', "Fraser",       'con', "constellationFraser_v001.mat",       'noise', "noiseCharacterisation_fraser.mat",      'imp', @importDB_fraser_v001, 'args', {{}}, 'sigmaToMM', 1/10.41793); ...
        struct('name', "Zarandi",      'con', "constellationZarandi_v001.mat",      'noise', "noiseCharacterisation_zarandi.mat",     'imp', @importDB_zarandi_v001, 'args', {{ "Verbose", false }}, 'sigmaToMM', 10.0); ...
        struct('name', "Cook_CTRL",    'con', "constellationCook_v001.mat",         'noise', "noiseCharacterisation_cook.mat",        'imp', @importDB_cook_v002, 'args', {{ "Group", "CTRL", "Tasks", 7, "Verbose", false }}, 'sigmaToMM', 0.248); ...
        struct('name', "Cook_ASD",     'con', "constellationCookASD_v001.mat",      'noise', "noiseCharacterisation_cookASD.mat",     'imp', @importDB_cook_v002, 'args', {{ "Group", "ASD", "Tasks", 7, "Verbose", false }}, 'sigmaToMM', 0.248); ...
        struct('name', "Dhieb",        'con', "constellationDhieb_v001.mat",        'noise', "noiseCharacterisation_dhieb.mat",       'imp', @importDB_dhieb_v001, 'args', {{}}, 'sigmaToMM', 0.1478); ...
        struct('name', "Hickman_PLAC", 'con', "constellationHickmanPLAC_v001.mat",  'noise', "noiseCharacterisation_hickmanPLAC.mat", 'imp', @importDB_hickman_v002, 'args', {{ "Study", 2, "Group", "PLAC", "Shapes", 3, "Verbose", false }}, 'sigmaToMM', 0.248); ...
        struct('name', "Hickman_HALO", 'con', "constellationHickmanHALO_v001.mat",  'noise', "noiseCharacterisation_hickmanHALO.mat", 'imp', @importDB_hickman_v002, 'args', {{ "Study", 2, "Group", "HALO", "Shapes", 3, "Verbose", false }}, 'sigmaToMM', 0.248); ...
        };
    nSets = numel(registry);

    %% Gather per-dataset inputs once (identical across variants)
    DS = cell(nSets, 1);
    for k = 1:nSets
        sp = registry{k};
        cf = fullfile(srcDir, sp.con);  nf = fullfile(srcDir, sp.noise);
        if ~isfile(cf) || ~isfile(nf), error('queryMonoSlope4:MissingFile', '%s', 'FAILED PATH: ' + string(cf) + ' | ' + string(nf)); end
        emp = load(cf, 'betaCanon');
        bio = load(nf, 'bioResults').bioResults;
        nTri = size(emp.betaCanon, 1);
        if ~ismember('f0', bio.Properties.VariableNames) || ~ismember(AlphaSource, bio.Properties.VariableNames)
            error('queryMonoSlope4:NoColumn', '%s', sp.name + ': bioResults lacks f0 or ' + AlphaSource);
        end
        if height(bio) ~= nTri
            error('queryMonoSlope4:RowMismatch', '%s', sprintf('%s: bioResults %d rows, constellation %d trials', sp.name, height(bio), nTri));
        end
        trials = sp.imp(sp.args{:});
        if numel(trials) ~= nTri
            error('queryMonoSlope4:LengthMismatch', '%s', sprintf('%s: importer %d trials, constellation %d', sp.name, numel(trials), nTri));
        end
        C = load(fullfile(srcDir, "loopClosureResults_" + sp.name + "_all_shaped_xu_v015.mat"), "results").results;
        if numel(C) ~= nTri
            error('queryMonoSlope4:CorpusCount', '%s', sprintf('%s: v015 corpus %d trials, constellation %d', sp.name, numel(C), nTri));
        end
        f0Bio = double(bio.f0);  f0C = arrayfun(@(r) double(r.f0), C(:));
        if ~isequaln(isfinite(f0Bio), isfinite(f0C)) || any(abs(f0Bio(isfinite(f0Bio)) ./ f0C(isfinite(f0C)) - 1) > 1e-6)
            error('queryMonoSlope4:Align', '%s', sp.name + ': f0 differs between bioResults and the v015 corpus; trial order not aligned');
        end
        Dx = datasetOwnTempoExponent_v002(sp.name);
        DS{k} = struct('name', string(sp.name), 'inDomain', Dx.inDomain, 'nTri', nTri, ...
            'alpha', double(bio.(AlphaSource)), 'sigmaMM', double(bio.sigmaMean) * sp.sigmaToMM, ...
            'fs', arrayfun(@(t) double(t.fs), trials).', 'f0', f0Bio, ...
            'beta', arrayfun(@(r) double(r.betaGenStarMed), C(:)));
        fprintf('%-13s n=%4d inDomain=%d  no exponent: %d\n', sp.name, nTri, Dx.inDomain, ...
            nnz(isfinite(f0Bio) & ~(isfinite(DS{k}.beta) & DS{k}.beta >= 0)));
    end

    %% Variants: VGF per trial
    vname = ["V0 RETIRED (175.2636)", "V1 CONSTANT (K(1/3))", "V2 PRIMARY (K(beta_i))"];
    nV = numel(vname);
    for k = 1:nSets
        d = DS{k};  has = isfinite(d.f0) & d.f0 > 0;
        v0 = nan(d.nTri, 1);  v0(has) = d.f0(has) * KCONV_RETIRED;
        v1 = nan(d.nTri, 1);  v1(has) = d.f0(has) * K13;
        v2 = nan(d.nTri, 1);                               % all datasets: own exponent, no imputation
        ok = has & isfinite(d.beta) & d.beta >= 0;
        v2(ok) = d.f0(ok) .* gridKConv_v001(d.beta(ok));
        DS{k}.vgf = {v0, v1, v2};
    end

    %% Run each variant
    ratioTable = table(["Fraser";"Zarandi";"Cook_CTRL";"Cook_ASD";"Hickman_PLAC";"Hickman_HALO";"Dhieb"], ...
        [7.43; 3.13; 2.98; 3.11; 2.89; 3.96; 0.58], 'VariableNames', {'dataset','semPartBRatio'});   % Finding #128, reused
    perTrial = cell(1, nV);  perDataPipe = cell(1, nV);  dsSummary = cell(1, nV);  nOutGrid = zeros(nSets, nV);
    for v = 1:nV
        fprintf('\n##### %s #####\n', vname(v));
        L = cell(nSets, 1);
        for k = 1:nSets
            d = DS{k};
            L{k} = buildLong(d, d.vgf{v}, mono);
            vg = d.vgf{v};  nOutGrid(k, v) = nnz(vg < min(mono.VGFGrid) | vg > max(mono.VGFGrid));
        end
        perTrial{v} = vertcat(L{:});
        [perDataPipe{v}, dsSummary{v}] = aggregateVariant(perTrial{v}, ratioTable);
    end

    %% Regression anchor: V0 must reproduce v003's saved dsSummary
    A = load(v003Mat, 'dsSummary').dsSummary;  B = dsSummary{1};
    if ~isequal(sort(A.dataset), sort(B.dataset)), error('queryMonoSlope4:AnchorKeys', '%s', 'V0 datasets differ from v003'); end
    B = B(arrayfun(@(x) find(B.dataset == x), A.dataset), :);
    for col = ["meanGRiseAcrossPipe" "meanGDescAcrossPipe" "invRiseFrac" "medianAlpha" "medianSigmaMM"]
        a = A.(col);  b = B.(col);
        if ~isequal(isnan(a), isnan(b)) || max(abs(a - b), [], 'omitnan') > ANCHOR_TOL
            error('queryMonoSlope4:Anchor', '%s', sprintf('V0 %s does not reproduce v003 (max diff %.3g)', col, max(abs(a - b), [], 'omitnan')));
        end
    end
    fprintf('\nREGRESSION ANCHOR passed: V0 reproduces v003 dsSummary (tolerance %.0e).\n', ANCHOR_TOL);

    %% Report
    dsNames = string(cellfun(@(x) x.name, DS));
    for v = 1:nV
        T = dsSummary{v};
        fprintf('\n=== %s: dataset summary x Finding #128 Part B ratio ===\n', vname(v));
        fprintf('%-13s %6s %6s %8s %8s %9s %7s %6s %8s\n', 'dataset', 'alpha', 'sigMM', 'gRise', 'gDesc', 'invRiseFr', 'ratioB', 'nOut', 'unplaced');
        for r = 1:height(T)
            k = find(dsNames == T.dataset(r), 1);
            fprintf('%-13s %6.3f %6.3f %8.4f %8.4f %9.3f %7.2f %6d %8d\n', T.dataset(r), T.medianAlpha(r), T.medianSigmaMM(r), ...
                T.meanGRiseAcrossPipe(r), T.meanGDescAcrossPipe(r), T.invRiseFracPlaced(r), T.semPartBRatio(r), nOutGrid(k, v), round(T.nUnplaced(r)));
        end
        [rS, pS] = corr(T.meanGRiseAcrossPipe, T.semPartBRatio, 'Type', 'Spearman', 'Rows', 'complete');
        [rA, pA] = corr(T.medianAlpha, T.semPartBRatio, 'Type', 'Spearman', 'Rows', 'complete');
        [rG, pG] = corr(T.medianSigmaMM, T.semPartBRatio, 'Type', 'Spearman', 'Rows', 'complete');
        [rI, pI] = corr(T.invRiseFracPlaced, T.semPartBRatio, 'Type', 'Spearman', 'Rows', 'complete');
        fprintf('Spearman vs Part B ratio (n=%d): gRise r=%+.3f p=%.4f | alpha r=%+.3f p=%.4f | sigma r=%+.3f p=%.4f | invRise r=%+.3f p=%.4f\n', ...
            height(T), rS, pS, rA, pA, rG, pG, rI, pI);
    end
    fprintf('\nCHANGE AGAINST V0 (gRise, invRiseFrac):\n');
    for v = 2:nV
        fprintf('  %s\n', vname(v));
        T0 = dsSummary{1};  T1 = dsSummary{v};
        for r = 1:height(T0)
            j = find(T1.dataset == T0.dataset(r));
            fprintf('    %-13s gRise %8.4f -> %8.4f   invRiseFr %6.3f -> %6.3f\n', T0.dataset(r), ...
                T0.meanGRiseAcrossPipe(r), T1.meanGRiseAcrossPipe(j), T0.invRiseFracPlaced(r), T1.invRiseFracPlaced(j));
        end
    end

    out = struct('perTrial', {perTrial}, 'perDataPipe', {perDataPipe}, 'dsSummary', {dsSummary}, ...
        'vname', vname, 'nOutGrid', nOutGrid, 'datasets', dsNames);
    if opts.Save
        resDir = fullfile(ROOT, 'results');
        if ~exist(resDir, 'dir'), error('queryMonoSlope4:NoResultsDir', '%s', 'results/ not found at the project root'); end
        constants = struct('K13', K13, 'KCONV_RETIRED', KCONV_RETIRED, 'runDate', string(datetime('now', 'Format', 'yyyy-MM-dd')));
        matOut = fullfile(resDir, 'queryMonotonicSlopeAtEmpirical_v005.mat');
        save(matOut, 'perTrial', 'perDataPipe', 'dsSummary', 'ratioTable', 'vname', 'nOutGrid', 'constants', 'monoFile', '-v7.3');
        fprintf('\nSaved: %s\n', matOut);
    end
end

function T = buildLong(d, vgf, mono)
% Long-format rows (trial x pipeline), as v003. A trial with no VGF is recorded, not dropped.
    nP = numel(mono.pipeOrder);  dummyBeta = NaN(1, nP);
    rows = cell(d.nTri * nP, 1);  ri = 0;
    for t = 1:d.nTri
        res = checkTrialInvertibility_v001(d.alpha(t), d.sigmaMM(t), d.fs(t), vgf(t), dummyBeta, mono);
        for p = 1:nP
            ri = ri + 1;
            base = struct('dataset', d.name, 'trialIdx', t, 'pipeline', mono.pipeOrder(p), ...
                'alphaTrial', d.alpha(t), 'sigmaMMTrial', d.sigmaMM(t));
            if isnan(res.aIdx)
                rows{ri} = catstruct2(base, struct('snapWarn', "coord NaN, not snapped", 'gRise', NaN, 'gDesc', NaN, 'invRise', false, 'invDesc', false, 'placed', false));
                continue;
            end
            a = res.aIdx; s = res.sIdx; f = res.fIdx; v = res.vIdx;
            rows{ri} = catstruct2(base, struct('snapWarn', res.snapWarnings, ...
                'gRise', mono.betaSeg.rise.meanSlope(p, a, s, f, v), 'gDesc', mono.betaSeg.desc.meanSlope(p, a, s, f, v), ...
                'invRise', mono.betaSeg.rise.invertible(p, a, s, f, v), 'invDesc', mono.betaSeg.desc.invertible(p, a, s, f, v), 'placed', true));
        end
    end
    T = struct2table([rows{:}]);
end

function s = catstruct2(a, b)
    s = a;  fn = fieldnames(b);
    for i = 1:numel(fn), s.(fn{i}) = b.(fn{i}); end
end

function [perDataPipe, dsSummary] = aggregateVariant(perTrial, ratioTable)
% v003's aggregation, unchanged.
    [G, gKeys] = findgroups(perTrial(:, {'dataset','pipeline'}));
    perDataPipe = gKeys;
    perDataPipe.nTrials       = splitapply(@numel, perTrial.trialIdx, G);
    perDataPipe.nInvRise      = splitapply(@(x) sum(x), perTrial.invRise, G);
    perDataPipe.nPlaced       = splitapply(@(x) sum(x), perTrial.placed, G);
    perDataPipe.medianAlpha   = splitapply(@(x) median(x, 'omitnan'), perTrial.alphaTrial, G);
    perDataPipe.medianSigmaMM = splitapply(@(x) median(x, 'omitnan'), perTrial.sigmaMMTrial, G);
    perDataPipe.medianGRise   = splitapply(@(x) median(x, 'omitnan'), perTrial.gRise, G);
    perDataPipe.medianGDesc   = splitapply(@(x) median(x, 'omitnan'), perTrial.gDesc, G);
    [Gd, dKeys] = findgroups(perDataPipe.dataset);
    dsSummary = table(dKeys, 'VariableNames', {'dataset'});
    dsSummary.meanGRiseAcrossPipe = splitapply(@(x) mean(x, 'omitnan'), perDataPipe.medianGRise, Gd);
    dsSummary.meanGDescAcrossPipe = splitapply(@(x) mean(x, 'omitnan'), perDataPipe.medianGDesc, Gd);
    dsSummary.medianAlpha   = splitapply(@(x) mean(x, 'omitnan'), perDataPipe.medianAlpha, Gd);
    dsSummary.medianSigmaMM = splitapply(@(x) mean(x, 'omitnan'), perDataPipe.medianSigmaMM, Gd);
    dsSummary.invRiseFrac   = splitapply(@(x) mean(x, 'omitnan'), perDataPipe.nInvRise ./ perDataPipe.nTrials, Gd);
    dsSummary.invRiseFracPlaced = splitapply(@(x) mean(x, 'omitnan'), perDataPipe.nInvRise ./ perDataPipe.nPlaced, Gd);
    dsSummary.nUnplaced = splitapply(@(x) mean(x), perDataPipe.nTrials - perDataPipe.nPlaced, Gd);
    dsSummary = innerjoin(dsSummary, ratioTable, 'Keys', 'dataset');
end
