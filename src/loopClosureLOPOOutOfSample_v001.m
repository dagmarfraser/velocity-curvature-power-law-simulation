function out = loopClosureLOPOOutOfSample_v001(opts)
%LOOPCLOSURELOPOOUTOFSAMPLE_V001  Strong-form leave-one-pipeline-out.
%
%   THE TEST
%   For each trial and each held-out pipeline p:
%     1. infer beta_gen* from the OTHER five pipelines (median of their finite
%        betaGenStar, requiring >= MinOthers of them)
%     2. push that value through pipeline p's OWN forward map to predict what
%        p should have observed:  predicted = interp1(betaGenVec,
%        betaRecCurveMed(p,:), betaGenStar_minus_p)
%     3. compare against what p actually observed: resid = predicted - betaObs(p)
%
%   Neither p's own betaGenStar nor p's own betaObs enters the prediction, so
%   this is genuinely out of sample for p. That is what distinguishes it from
%   the scoped form (loopClosureLOPO_v001.m, Finding #199), which compares
%   inverted values against each other and therefore cannot detect a bias
%   shared by all six pipelines.
%
%   WHY IT NEEDS v015
%   v013 discarded the forward map, persisting only betaRecSlice -- a
%   [N_REPS x 6] replicate cloud at ONE grid point, not a curve. v015 persists
%   betaRecCurveMed / betaRecCurveSD, the exact two reductions
%   invertBeta_local itself consumes, so the prediction runs against the
%   identical curve the inversion used rather than a reconstruction of it.
%   Running this against a v013 corpus is an error, not a fallback.
%
%   WHAT A RESULT WOULD MEAN
%   This is the residue of GPT SOL's circularity critique (section 2) and the
%   manuscript's own most substantive open item. Small residuals mean the
%   framework predicts a pipeline's behaviour it was not shown; large ones mean
%   the per-trial forward maps do not generalise across pipelines, which would
%   bound how far pointwise inversion can be trusted.
%
%   BENCHMARKS, both computed here, because a residual with nothing to compare
%   against is uninterpretable:
%     - IN-SAMPLE floor: predict p from p's OWN betaGenStar. This is what the
%       inversion fits, so it is the best the curve can do. Out-of-sample
%       residuals should be read as a multiple of this, not against zero.
%     - SPREAD reference: the scoped-form LOPO residual (Finding #199) was
%       median 0.0184 with 66.5% inside MDC.
%
%   Fail Loud, Never Fake: a corpus without betaRecCurveMed is an error. A
%   predicted value falling outside the curve's beta range is NaN and COUNTED,
%   never clamped to the nearest endpoint -- clamping would manufacture
%   agreement exactly where the map has stopped covering the data.
%
%   Provenance: written 2026-09-16 while the v015 re-run was in flight, against
%   the field contract makeRunLoopClosureV015_v001.m defines. Not yet executed.
%
%   See also LOOPCLOSURELOPO_V001, MAKERUNLOOPCLOSUREV015_V001.

arguments
    opts.Datasets   (1,:) string = ["Zarandi","Cook_CTRL","Cook_ASD", ...
                                    "Hickman_PLAC","Hickman_HALO","Dhieb","Fraser"]
    opts.NoiseModel (1,1) string  = "shaped_xu"
    opts.Version    (1,1) string  = "v015"
    opts.MinOthers  (1,1) double  = 3
    opts.MDC        (1,1) double  = 0.03
    opts.Save       (1,1) logical = true
end

PP = ["BWFD-OLS","SG-OLS","BWFD-LMLS","SG-LMLS","BWFD-IRLS","SG-IRLS"];
srcDir = fileparts(mfilename("fullpath"));

rows = [];   % [residOut, residIn, pipelineIdx, datasetIdx]
nOutOfRange = 0; nAttempt = 0;
perDS = struct("dataset",{},"n",{},"medOut",{},"pctMDC",{},"medIn",{},"ratio",{},"oor",{});

fprintf("STRONG-FORM LOPO: predict held-out pipeline's betaObs through its own forward map\n");
fprintf("%-14s %7s %9s %8s %9s %7s %6s\n", ...
    "dataset","nCells","medOut","%<MDC","medIn","ratio","OOR%");

for k = 1:numel(opts.Datasets)
    tag = strrep(opts.Datasets(k), " ", "_");
    f = fullfile(srcDir, sprintf("loopClosureResults_%s_all_%s_%s.mat", ...
        tag, opts.NoiseModel, opts.Version));
    if ~isfile(f)
        error("lopoOOS:missingCorpus", "Not found:\n  %s", f);
    end
    S = load(f, "results", "betaGenVec");
    r = S.results;
    if ~isfield(r, "betaRecCurveMed")
        error("lopoOOS:noCurve", ...
            "%s has no betaRecCurveMed. This is a pre-v015 corpus; the " + ...
            "strong-form test is not possible against it.", f);
    end
    bg = S.betaGenVec(:).';

    dsRows = []; dsOOR = 0;
    for ti = 1:numel(r)
        bgs   = r(ti).betaGenStar(:).';        % 1 x 6
        bobs  = r(ti).betaObs(:).';            % 1 x 6
        curve = r(ti).betaRecCurveMed;         % 6 x N_BETA
        if size(curve,2) ~= numel(bg)
            error("lopoOOS:curveShape", "%s trial %d: curve is %s, betaGenVec is %d long.", ...
                tag, ti, mat2str(size(curve)), numel(bg));
        end

        for p = 1:6
            if ~isfinite(bobs(p)), continue; end
            others = bgs(setdiff(1:6, p));
            others = others(isfinite(others));
            if numel(others) < opts.MinOthers, continue; end

            nAttempt = nAttempt + 1;
            predOut = interpCurve_local(bg, curve(p,:), median(others));
            if isnan(predOut)
                nOutOfRange = nOutOfRange + 1; dsOOR = dsOOR + 1;
                continue
            end
            % in-sample floor: p's own inverted value through p's own curve
            predIn = NaN;
            if isfinite(bgs(p))
                predIn = interpCurve_local(bg, curve(p,:), bgs(p));
            end
            dsRows(end+1,:) = [predOut - bobs(p), predIn - bobs(p), p, k]; %#ok<AGROW>
        end
    end

    if isempty(dsRows)
        fprintf("%-14s %7d   -- no evaluable cells --\n", opts.Datasets(k), 0);
        continue
    end
    ao = abs(dsRows(:,1));
    ai = abs(dsRows(:,2)); ai = ai(isfinite(ai));
    medIn = median(ai);
    ratio = median(ao) / max(medIn, eps);
    oorPct = 100*dsOOR / max(dsOOR + size(dsRows,1), 1);
    fprintf("%-14s %7d %9.4f %7.1f%% %9.4f %7.1f %5.1f%%\n", ...
        opts.Datasets(k), size(dsRows,1), median(ao), 100*mean(ao < opts.MDC), ...
        medIn, ratio, oorPct);

    perDS(end+1) = struct("dataset",opts.Datasets(k), "n",size(dsRows,1), ...
        "medOut",median(ao), "pctMDC",100*mean(ao < opts.MDC), "medIn",medIn, ...
        "ratio",ratio, "oor",oorPct); %#ok<AGROW>
    rows = [rows; dsRows]; %#ok<AGROW>
end

if isempty(rows)
    error("lopoOOS:noCells", "No evaluable cells in any dataset.");
end

fprintf("\nPOOLED per held-out pipeline:\n%-12s %8s %9s %8s %9s\n", ...
    "pipeline","n","medOut","%<MDC","medIn");
perPP = struct("pipeline",{},"n",{},"medOut",{},"pctMDC",{},"medIn",{});
for p = 1:6
    m = rows(:,3) == p;
    if ~any(m), continue; end
    ao = abs(rows(m,1)); ai = abs(rows(m,2)); ai = ai(isfinite(ai));
    fprintf("%-12s %8d %9.4f %7.1f%% %9.4f\n", PP(p), nnz(m), median(ao), ...
        100*mean(ao < opts.MDC), median(ai));
    perPP(end+1) = struct("pipeline",PP(p), "n",nnz(m), "medOut",median(ao), ...
        "pctMDC",100*mean(ao < opts.MDC), "medIn",median(ai)); %#ok<AGROW>
end

ao = abs(rows(:,1)); ai = abs(rows(:,2)); ai = ai(isfinite(ai));
fprintf("\nALL: n=%d  medOut=%.4f  %%<MDC=%.1f%%  medIn=%.4f  ratio=%.1f\n", ...
    size(rows,1), median(ao), 100*mean(ao < opts.MDC), median(ai), median(ao)/max(median(ai),eps));
fprintf("out-of-range predictions excluded: %d of %d attempts (%.1f%%)\n", ...
    nOutOfRange, nAttempt, 100*nOutOfRange/max(nAttempt,1));
fprintf("\nRead medOut against medIn, not against zero. Scoped-form comparator\n");
fprintf("(Finding #199): median 0.0184, 66.5%% inside MDC.\n");

out = struct("perDataset",perDS, "perPipeline",perPP, "rows",rows, ...
    "pipelineNames",PP, "opts",opts, "nOutOfRange",nOutOfRange, ...
    "nAttempt",nAttempt, "generated",string(datetime("now")));

if opts.Save
    o = fullfile(srcDir, "loopClosureLOPOOutOfSample_v001.mat");
    save(o, "-struct", "out");
    fprintf("\nSaved: %s\n", o);
end
end

% =========================================================================
function y = interpCurve_local(bg, c, x)
% Linear interpolation on the forward map, NaN outside its support.
% No clamping: a beta_gen* outside the curve's range means the map does not
% cover that trial, which is a result, not something to round away.
m = isfinite(c) & isfinite(bg);
if nnz(m) < 2 || ~isfinite(x), y = NaN; return; end
bgm = bg(m); cm = c(m);
if x < min(bgm) || x > max(bgm), y = NaN; return; end
[bgu, iu] = unique(bgm);              % interp1 needs strictly increasing x
y = interp1(bgu, cm(iu), x, "linear");
end
