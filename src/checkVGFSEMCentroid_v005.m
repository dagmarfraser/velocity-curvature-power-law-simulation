%% checkVGFSEMCentroid_v005.m
% v005 (2026-09-30) supersedes v004 only in the out-of-domain rule: Zarandi and Dhieb use their own
% median betaGenStarMed (flagged, not a generator estimate) instead of the 1/3 reference. In-domain unchanged.
% VGF bias and SEM at each dataset's own tempo (Part 1 Results (c); Findings #194, D7),
% rebuilt on the grid's true tempo constant. v003 placed each dataset on the VGF axis with
% K_conv = 175.2636 (perimeter / mean(kappa^-1/3)) at beta = 1/3 and mean f0. That constant
% is 5.28% below K(1/3) = 185.02 (testGridKConv_v001), and K falls steeply with beta.
% v004 runs four variants on identical inputs so each change is attributable:
%   V0  RETIRED   mean f0 (noiseCharacterisation) x 175.2636          regression anchor: must
%                 reproduce the published 26.6-67.9% (bias) and 1.7-9.2% (SEM) exactly.
%   V1  CONSTANT  mean f0 (noiseCharacterisation) x K(1/3)            the constant fixed only.
%   V2  PRIMARY   median f0, ALL trials (v015 corpus) x K(beta_ds)    decided 2026-09-29:
%                 beta_ds = the dataset's median betaGenStarMed (in-domain datasets);
%                 outside the validity domain (Zarandi, Dhieb; #242) the same rule, flagged
%                 (their beta_gen* is a model response, not a generator estimate).
%   V3  MEAN      mean f0, all trials (v015) x K(beta_ds)             continuity with the
%                 manuscript's "mean f0" wording.
% In every variant the Hz conversions of bias and SEM use that variant's own K, so
% bias%f0 = bias_VGF / VGF_own.
% Unchanged from v003: alpha/sigma centroids (noiseCharacterisation means), fs, nearest-node
% snapping in alpha, sigma, fs and VGF, the beta SEM gate (semAdequacyThreshold_v001), and
% INCLUDED cells only in the six-pipeline mean.
% Reporting rule (Dagmar, 2026-09-29): the headline range is in-domain only; Zarandi (outside the
% validity domain, #242) is reported beside it, not folded in, because its difference is itself informative.
% The all-datasets range is kept for the V0 regression anchor.
% Self-checks (Fail Loud): K(1/3) = 185.0242; V0 reproduces the published ranges; in-domain
%   trial counts equal Table 6a; v015 and noiseCharacterisation mean f0 agree within 2% (warn).
% Reads:  src/noiseCharacterisation_*.mat, src/loopClosureResults_*_all_shaped_xu_v015.mat,
%         src/perCoordinateSEM_v2_001.mat, src/perCoordinateSEM_VGF_v001.mat.
% Writes: results/checkVGFSEMCentroid_v005.mat
% USAGE:  from the project root: checkVGFSEMCentroid_v005
% Fraser, D.S. (2026)  v005

%% CONFIG
ROOT   = fileparts(fileparts(mfilename("fullpath")));
SRC    = fullfile(ROOT, "src");
addpath(ROOT); addpath(SRC); addpath(genpath(fullfile(SRC, "functions")));
% label, noiseCharacterisation file, fs, native-sigma-to-mm factor, v015 name, in validity domain
CONFIG = {
    "Fraser",       "noiseCharacterisation_fraser.mat",       240, 0.095988, "Fraser",       true;
    "Cook CTRL",    "noiseCharacterisation_cook.mat",         133, 0.248,    "Cook_CTRL",    true;
    "Cook ASD",     "noiseCharacterisation_cookASD.mat",      133, 0.248,    "Cook_ASD",     true;
    "Hickman PLAC", "noiseCharacterisation_hickmanPLAC.mat",  133, 0.248,    "Hickman_PLAC", true;
    "Hickman HALO", "noiseCharacterisation_hickmanHALO.mat",  133, 0.248,    "Hickman_HALO", true;
    "Zarandi",      "noiseCharacterisation_zarandi.mat",      100, 10.0,     "Zarandi",      false;
    "Dhieb",        "noiseCharacterisation_dhieb.mat",        100, 0.1478,   "Dhieb",        false;
};
N_TABLE6A = [2829 94 102 359 338 NaN NaN];
BETA_REF  = 1/3;
KCONV_RETIRED = 175.2636;
K13_PUB = 185.0242;  K13_TOL = 1e-3;
PUB_BIAS = [26.6 67.9];  PUB_SEM = [1.7 9.2];        % v004 Part 1 Results (c), from checkVGFSEMCentroid_v003
F0_AGREE_TOL = 0.02;
OUT_MAT = fullfile(ROOT, "results", "checkVGFSEMCentroid_v005.mat");
nD = size(CONFIG, 1);

K13 = gridKConv_v001(BETA_REF);
if abs(K13 - K13_PUB) > K13_TOL
    error("vgfCentroid4:K13", "%s", sprintf("K(1/3) %.4f does not reproduce %.4f", K13, K13_PUB));
end

%% Centroids from noiseCharacterisation, and f0 / exponent from the v015 corpora
label = strings(nD, 1); alpha = nan(nD, 1); sigmaMM = nan(nD, 1); fs = nan(nD, 1);
f0MeanBio = nan(nD, 1); f0Mean15 = nan(nD, 1); f0Med15 = nan(nD, 1);
betaDs = nan(nD, 1); betaGenMed = nan(nD, 1); nTr = nan(nD, 1); nExp = nan(nD, 1); inDom = false(nD, 1);
for i = 1:nD
    ncFile = fullfile(SRC, CONFIG{i, 2});
    if ~isfile(ncFile), error("vgfCentroid4:NotFound", "%s", "FAILED PATH: " + ncFile); end
    T = load(ncFile, "bioResults").bioResults;
    label(i)   = CONFIG{i, 1};   fs(i) = CONFIG{i, 3};   inDom(i) = CONFIG{i, 6};
    alpha(i)   = mean(T.ira_alphaMean, "omitnan");
    sigmaMM(i) = mean(T.sigmaMean, "omitnan") * CONFIG{i, 4};
    f0MeanBio(i) = mean(T.f0, "omitnan");

    cf = fullfile(SRC, "loopClosureResults_" + CONFIG{i, 5} + "_all_shaped_xu_v015.mat");
    if ~isfile(cf), error("vgfCentroid4:NoCorpus", "%s", "FAILED PATH: " + cf); end
    R = load(cf, "results").results;
    nTr(i) = numel(R);
    if ~isnan(N_TABLE6A(i)) && nTr(i) ~= N_TABLE6A(i)
        error("vgfCentroid4:Count", "%s", sprintf("%s: %d trials, Table 6a says %d", label(i), nTr(i), N_TABLE6A(i)));
    end
    f0 = arrayfun(@(r) double(r.f0), R(:));
    g  = arrayfun(@(r) double(r.betaGenStarMed), R(:));
    okF = isfinite(f0) & f0 > 0;
    if any(~okF)
        warning("vgfCentroid4:NoF0", "%s", sprintf("%s: %d trials without a usable f0 excluded", label(i), nnz(~okF)));
    end
    f0Mean15(i) = mean(f0(okF));  f0Med15(i) = median(f0(okF));
    betaGenMed(i) = median(g(okF & isfinite(g)));
    nExp(i) = nnz(okF & isfinite(g));
    if isnan(betaGenMed(i)) || betaGenMed(i) < 0
        error("vgfCentroid5:Beta", "%s", sprintf("%s: no usable dataset exponent", label(i)));
    end
    betaDs(i) = betaGenMed(i);
    rel = f0Mean15(i) / f0MeanBio(i) - 1;
    if abs(rel) > F0_AGREE_TOL
        warning("vgfCentroid4:F0Source", "%s", sprintf("%s: mean f0 differs by %+.1f%% between v015 (%.4f) and noiseCharacterisation (%.4f)", ...
            label(i), 100 * rel, f0Mean15(i), f0MeanBio(i)));
    end
end
centroids = table(label, inDom, nTr, nExp, alpha, sigmaMM, fs, f0MeanBio, f0Mean15, f0Med15, betaGenMed, betaDs);
disp(centroids);

%% Grid tables
ctx.Tbeta = load(fullfile(SRC, "perCoordinateSEM_v2_001.mat"), "coordTable").coordTable;
vf = fullfile(SRC, "perCoordinateSEM_VGF_v001.mat");
if ~isfile(vf), error("vgfCentroid4:NotFound", "%s", "FAILED PATH: " + vf + ". Run computePerCoordinateSEM_VGF_v001() first."); end
ctx.Tvgf = load(vf, "coordTable").coordTable;
ctx.SEM_ADEQUATE = semAdequacyThreshold_v001();
ctx.allAlpha = sort(unique(ctx.Tbeta.alpha));  ctx.allSigma = sort(unique(ctx.Tbeta.sigma));
ctx.allFs = sort(unique(ctx.Tbeta.fs));        ctx.allVGF = sort(unique(ctx.Tvgf.VGF));
ctx.VGF_LO = exp(4.5);  ctx.VGF_HI = exp(5.8);

%% Variants
V = struct("name", {"V0 RETIRED  (mean f0 x 175.2636)", "V1 CONSTANT (mean f0 x K(1/3))", ...
                    "V2 PRIMARY  (median f0, all trials x K(beta_ds))", "V3 MEAN     (mean f0, all trials x K(beta_ds))"}, ...
           "f0", {f0MeanBio, f0MeanBio, f0Med15, f0Mean15}, ...
           "K",  {repmat(KCONV_RETIRED, nD, 1), repmat(K13, nD, 1), gridKConv_v001(betaDs), gridKConv_v001(betaDs)});
cellsAll = table();  summ = cell(numel(V), 1);
for k = 1:numel(V)
    [C, S] = evalVariant(V(k).f0, V(k).K, label, centroids, inDom, ctx);
    C.variant = repmat(k - 1, height(C), 1);  cellsAll = [cellsAll; C]; %#ok<AGROW>
    summ{k} = S;
end

%% Regression anchor: V0 must reproduce the published ranges
S0 = summ{1};
got = [round(S0.biasRange, 1); round(S0.semRange, 1)];
if any(abs(got(1, :) - PUB_BIAS) > 0) || any(abs(got(2, :) - PUB_SEM) > 0)
    error("vgfCentroid4:Anchor", "%s", sprintf("V0 gives bias %.1f-%.1f, SEM %.1f-%.1f; published %.1f-%.1f, %.1f-%.1f", ...
        S0.biasRange, S0.semRange, PUB_BIAS, PUB_SEM));
end
fprintf("\nREGRESSION ANCHOR passed: V0 reproduces bias %.1f-%.1f%% and SEM %.1f-%.1f%%.\n", S0.biasRange, S0.semRange);

%% Report
for k = 1:numel(V)
    S = summ{k};
    fprintf("\n=== %s ===\n", V(k).name);
    fprintf("%-13s %6s %8s %9s %6s | %9s %9s %s\n", "dataset", "f0(Hz)", "K", "VGF_own", "node", "bias%f0", "sem%f0", "");
    for d = 1:nD
        fprintf("%-13s %6.3f %8.2f %9.1f %6.1f | %9s %9s %s\n", label(d), V(k).f0(d), V(k).K(d), S.vgfOwn(d), S.node(d), ...
            fmt(S.meanBias(d)), fmt(S.meanSem(d)), extrapTag(S.extrap(d), inDom(d)));
    end
    fprintf("HEADLINE range (in-domain, INCLUDED cells): bias %.1f-%.1f%%   SEM %.1f-%.1f%%   (datasets: %d)\n", ...
        S.biasRangeDom, S.semRangeDom, S.nDatasetsDom);
    fprintf("With Zarandi, as v003 published:            bias %.1f-%.1f%%   SEM %.1f-%.1f%%   (datasets: %d)\n", ...
        S.biasRange, S.semRange, S.nDatasets);
end

%% Save
constants = struct("K13", K13, "KCONV_RETIRED", KCONV_RETIRED, "BETA_REF", BETA_REF, ...
    "runDate", string(datetime("now", "Format", "yyyy-MM-dd")));
save(OUT_MAT, "centroids", "cellsAll", "summ", "constants");
fprintf("\nSaved %s\n", OUT_MAT);

%% Local functions
function s = fmt(x)
    if isnan(x), s = "--"; else, s = sprintf("%.1f", x); end
end

function t = extrapTag(e, dom)
    t = "";
    if e, t = t + "[VGF EXTRAPOLATED] "; end
    if ~dom, t = t + "[outside domain: own beta is a model response]"; end
end

function [cells, S] = evalVariant(f0v, Kv, label, cen, dom, ctx)
    nD = numel(label);
    rows = cell(0, 1);
    vgfOwn = nan(nD, 1); node = nan(nD, 1); extrap = false(nD, 1);
    meanBias = nan(nD, 1); meanSem = nan(nD, 1);
    for di = 1:nD
        tV = f0v(di) * Kv(di);
        [~, ai] = min(abs(ctx.allAlpha - cen.alpha(di)));
        [~, si] = min(abs(ctx.allSigma - cen.sigmaMM(di)));
        [~, fi] = min(abs(ctx.allFs - cen.fs(di)));
        [~, vi] = min(abs(ctx.allVGF - tV));
        sA = ctx.allAlpha(ai); sS = ctx.allSigma(si); sF = ctx.allFs(fi); sV = ctx.allVGF(vi);
        vgfOwn(di) = tV;  node(di) = sV;  extrap(di) = tV < ctx.VGF_LO || tV > ctx.VGF_HI;
        subBeta = ctx.Tbeta(ctx.Tbeta.alpha == sA & ctx.Tbeta.sigma == sS & ctx.Tbeta.fs == sF, :);
        subVGF  = ctx.Tvgf(ctx.Tvgf.alpha == sA & ctx.Tvgf.sigma == sS & ctx.Tvgf.fs == sF & ctx.Tvgf.VGF == sV, :);
        if isempty(subBeta)
            error("vgfCentroid4:NoRows", "%s", sprintf("%s: no beta rows at the snapped coordinate", label(di)));
        end
        [Gb, pipes] = findgroups(subBeta.pipeline);
        betaSEM = splitapply(@(x) mean(x, "omitnan"), subBeta.sem, Gb);
        bp = [];  sp = [];
        for p = 1:numel(pipes)
            adequate = betaSEM(p) < ctx.SEM_ADEQUATE;
            vr = subVGF(subVGF.pipeline == pipes(p), :);
            if isempty(vr), biasHz = NaN; semHz = NaN;
            else
                biasHz = mean(vr.meanBias, "omitnan") / Kv(di);
                semHz  = mean(vr.sem, "omitnan") / Kv(di);
            end
            biasPct = 100 * abs(biasHz) / f0v(di);
            semPct  = 100 * semHz / f0v(di);
            if adequate, bp(end+1) = biasPct; sp(end+1) = semPct; end %#ok<AGROW>
            rows{end+1, 1} = struct("dataset", label(di), "pipeline", string(pipes(p)), "adequate", adequate, ...
                "vgfOwn", tV, "node", sV, "extrap", extrap(di), "biasPct", biasPct, "semPct", semPct); %#ok<AGROW>
        end
        if ~isempty(bp), meanBias(di) = mean(bp); meanSem(di) = mean(sp); end
    end
    cells = struct2table(vertcat(rows{:}));
    okD = ~isnan(meanBias);
    okIn = okD & dom;                                      % headline range: validity domain only
    S = struct("vgfOwn", vgfOwn, "node", node, "extrap", extrap, "meanBias", meanBias, "meanSem", meanSem, ...
        "biasRange", [min(meanBias(okD)) max(meanBias(okD))], "semRange", [min(meanSem(okD)) max(meanSem(okD))], ...
        "nDatasets", nnz(okD), ...
        "biasRangeDom", [min(meanBias(okIn)) max(meanBias(okIn))], "semRangeDom", [min(meanSem(okIn)) max(meanSem(okIn))], ...
        "nDatasetsDom", nnz(okIn));
end
