function out = investigateSEMPartBTrialLevel_v002(opts)
% investigateSEMPartBTrialLevel_v002  Trial-level follow-up to the Part B ecological correlation
% (Finding #128), Fraser alone, rebuilt on the grid's true tempo constant.
%
% v001 placed every trial on the VGF axis as f0 x K_conv, K_conv = 175.2636 (beta = 1/3), 5.28%
% below K(1/3) = 185.02 (testGridKConv_v001), with K itself falling steeply in beta. v002 runs
% three variants on identical inputs so each change is attributable:
%   V0 RETIRED   f0 x 175.2636      regression anchor: must reproduce v001's saved table T.
%   V1 CONSTANT  f0 x K(1/3)        the constant fixed only.
%   V2 PRIMARY   f0 x K(beta_i)     beta_i = the trial's own betaGenStarMed (from the same file).
% Every retained row already has a finite betaGenStarMed (v001 required it for semTrial), so
% nothing is imputed and no trial is dropped by the exponent rule.
% Deliberately unchanged from v001: the source file loopClosureResults_Fraser_all_shaped_xu_v008
% (20 replicates per trial, the basis of the published diagnostic; the current v015 corpus holds 200
% and would change semTrial itself, a separate question), fs = 240 (importDB_fraser uniform clock),
% the semTrial definition, the predictors and the analysis. NOT preregistered.
%
% Reads:  src/loopClosureResults_Fraser_all_shaped_xu_v008.mat, src/monotonicSegments_v2_002.mat,
%         src/results/investigateSEMPartBTrialLevel_v001.mat (anchor).
% Writes: results/investigateSEMPartBTrialLevel_v002.mat (project-root results/).
% USAGE:  from the project root: R = investigateSEMPartBTrialLevel_v002()   (Save=false to skip)
% Fraser, D.S. (2026)  v002
    arguments
        opts.Save (1,1) logical = true
    end
    ROOT = fileparts(fileparts(mfilename('fullpath')));
    srcDir = fullfile(ROOT, 'src');
    addpath(ROOT); addpath(srcDir); addpath(genpath(fullfile(srcDir, 'functions')));
    KCONV_RETIRED = 175.2636;  K13_PUB = 185.0242;  FS_FIXED = 240;  ANCHOR_TOL = 1e-9;
    K13 = gridKConv_v001(1/3);
    if abs(K13 - K13_PUB) > 1e-3
        error('investSEMB2:K13', '%s', sprintf('K(1/3) %.4f does not reproduce %.4f', K13, K13_PUB));
    end
    v001Mat = fullfile(srcDir, 'results', 'investigateSEMPartBTrialLevel_v001.mat');
    lcPath  = fullfile(srcDir, 'loopClosureResults_Fraser_all_shaped_xu_v008.mat');
    monoFile = fullfile(srcDir, 'monotonicSegments_v2_002.mat');
    for f = [string(v001Mat) string(lcPath) string(monoFile)]
        if ~isfile(f), error('investSEMB2:NoFile', '%s', 'FAILED PATH: ' + f); end
    end
    mono = load(monoFile, 'params', 'pipeOrder', 'alphaGrid', 'sigmaGrid', 'fsGrid', 'betaGenGrid', 'VGFGrid', 'betaSeg');
    lc = load(lcPath, 'results').results;

    vname = ["V0 RETIRED (175.2636)", "V1 CONSTANT (K(1/3))", "V2 PRIMARY (K(beta_i))"];
    nV = numel(vname);  T = cell(1, nV);  models = cell(1, nV);
    for v = 1:nV
        T{v} = buildTable(lc, mono, FS_FIXED, v, KCONV_RETIRED, K13);
        models{v} = fitPerPipeline(T{v});
    end

    %% Regression anchor: V0 must reproduce v001's saved T
    A = load(v001Mat, 'T').T;  B = T{1};
    if height(A) ~= height(B), error('investSEMB2:AnchorRows', '%s', sprintf('V0 has %d rows, v001 %d', height(B), height(A))); end
    if ~isequal(string(A.trialID), string(B.trialID)) || ~isequal(string(A.pipeline), string(B.pipeline))
        error('investSEMB2:AnchorKeys', '%s', 'V0 trialID/pipeline order differs from v001');
    end
    for col = ["semTrial" "gRise" "gDesc" "alpha" "sigmaMM" "f0"]
        d = max(abs(A.(col) - B.(col)), [], 'omitnan');
        if d > ANCHOR_TOL, error('investSEMB2:Anchor', '%s', sprintf('V0 %s differs from v001 by %.3g', col, d)); end
    end
    if ~isequal(A.invRise, B.invRise) || ~isequal(A.inGridVGF, B.inGridVGF)
        error('investSEMB2:AnchorLogical', '%s', 'V0 invRise/inGridVGF differ from v001');
    end
    fprintf('REGRESSION ANCHOR passed: V0 reproduces v001 T (%d rows, tolerance %.0e).\n', height(A), ANCHOR_TOL);

    %% Report
    predictors = ["gRise","gDesc","invRise","betaGenStarMed","alpha","sigmaMM"];
    for v = 1:nV
        Tv = T{v};  Tin = Tv(Tv.inGridVGF, :);
        fprintf('\n=== %s: rows %d, in-grid VGF %d (%.1f%%) ===\n', vname(v), height(Tv), height(Tin), 100 * height(Tin) / height(Tv));
        fprintf('%-16s %10s %10s | %12s\n', 'predictor', 'Pearson r', 'Spearman r', 'in-grid r');
        for pr = predictors
            rP = corr(double(Tv.(pr)), Tv.semTrial, 'Rows', 'complete');
            rS = corr(double(Tv.(pr)), Tv.semTrial, 'Rows', 'complete', 'Type', 'Spearman');
            rI = corr(double(Tin.(pr)), Tin.semTrial, 'Rows', 'complete');
            fprintf('%-16s %10.3f %10.3f | %12.3f\n', pr, rP, rS, rI);
        end
        fprintf('fitlm per pipeline (semTrial ~ gRise + gDesc + invRise + betaGenStarMed): R2 and gRise coefficient (p)\n');
        for i = 1:numel(models{v})
            m = models{v}(i);
            fprintf('  %-8s n=%d R2=%.3f  gRise %+.4f (p=%.3g)\n', m.pipeline, m.n, m.R2, m.gRiseB, m.gRiseP);
        end
    end
    fprintf('\nCHANGE AGAINST V0 (pooled Pearson r, semTrial vs predictor):\n%-16s %8s %8s %8s\n', 'predictor', 'V0', 'V1', 'V2');
    for pr = predictors
        r = arrayfun(@(v) corr(double(T{v}.(pr)), T{v}.semTrial, 'Rows', 'complete'), 1:nV);
        fprintf('%-16s %8.3f %8.3f %8.3f\n', pr, r);
    end

    out = struct('T', {T}, 'models', {models}, 'vname', vname);
    if opts.Save
        resDir = fullfile(ROOT, 'results');
        if ~exist(resDir, 'dir'), error('investSEMB2:NoResultsDir', '%s', 'results/ not found at the project root'); end
        constants = struct('K13', K13, 'KCONV_RETIRED', KCONV_RETIRED, 'runDate', string(datetime('now', 'Format', 'yyyy-MM-dd')));
        matOut = fullfile(resDir, 'investigateSEMPartBTrialLevel_v002.mat');
        save(matOut, 'T', 'models', 'vname', 'constants', 'monoFile', 'lcPath', '-v7.3');
        fprintf('\nSaved: %s\n', matOut);
    end
end

function T = buildTable(lc, mono, fs, variant, kRetired, k13)
% v001's row construction; only the VGF placement depends on the variant.
    nP = numel(mono.pipeOrder);  dummyBeta = NaN(1, nP);
    rows = cell(numel(lc) * nP, 1);  ri = 0;
    for t = 1:numel(lc)
        alphaT = lc(t).alphaIRA;  sigmaT = lc(t).sigmaMM;  f0T = lc(t).f0;  bG = lc(t).betaGenStarMed;
        if ~isfinite(alphaT) || ~isfinite(sigmaT) || ~isfinite(f0T), continue; end
        switch variant
            case 1, vgfT = f0T * kRetired;
            case 2, vgfT = f0T * k13;
            case 3
                if ~isfinite(bG) || bG < 0, continue; end      % no exponent: not placed (v001 also dropped these rows)
                vgfT = f0T * gridKConv_v001(bG);
        end
        res = checkTrialInvertibility_v001(alphaT, sigmaT, fs, vgfT, dummyBeta, mono);
        if isnan(res.aIdx), continue; end
        a = res.aIdx; s = res.sIdx; f = res.fIdx; v = res.vIdx;
        for p = 1:nP
            slice = lc(t).betaRecSlice(:, p);
            if ~(isfinite(bG) && sum(isfinite(slice)) >= 2), continue; end
            ri = ri + 1;
            rows{ri} = struct('trialID', string(lc(t).trialID), 'pipeline', mono.pipeOrder(p), ...
                'alpha', alphaT, 'sigmaMM', sigmaT, 'f0', f0T, 'betaGenStarMed', bG, ...
                'semTrial', std(slice, 'omitnan'), ...
                'gRise', mono.betaSeg.rise.meanSlope(p, a, s, f, v), 'gDesc', mono.betaSeg.desc.meanSlope(p, a, s, f, v), ...
                'invRise', mono.betaSeg.rise.invertible(p, a, s, f, v), 'invDesc', mono.betaSeg.desc.invertible(p, a, s, f, v), ...
                'inGridVGF', ~contains(res.snapWarnings, "VGF="));
        end
    end
    T = struct2table([rows{1:ri}]);
    T = T(isfinite(T.semTrial), :);
end

function M = fitPerPipeline(T)
    pipes = unique(T.pipeline, 'stable');
    M = repmat(struct('pipeline', "", 'n', 0, 'R2', NaN, 'gRiseB', NaN, 'gRiseP', NaN, 'mdl', []), numel(pipes), 1);
    for i = 1:numel(pipes)
        R = T(T.pipeline == pipes(i), :);
        tbl = table(R.semTrial, R.gRise, R.gDesc, double(R.invRise), R.betaGenStarMed, ...
            'VariableNames', {'semTrial','gRise','gDesc','invRise','betaGenStarMed'});
        mdl = fitlm(tbl, 'semTrial ~ gRise + gDesc + invRise + betaGenStarMed');
        M(i) = struct('pipeline', pipes(i), 'n', height(tbl), 'R2', mdl.Rsquared.Ordinary, ...
            'gRiseB', mdl.Coefficients{'gRise', 'Estimate'}, 'gRiseP', mdl.Coefficients{'gRise', 'pValue'}, 'mdl', mdl);
    end
end
