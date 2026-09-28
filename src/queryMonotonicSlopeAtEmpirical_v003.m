function queryMonotonicSlopeAtEmpirical_v003()
% queryMonotonicSlopeAtEmpirical_v003.m
% Extends v002 (f0-based VGF matching) two ways, both in service of testing
% Dagmar's SEM-agreement-Part-B hypothesis (Finding #128 update note,
% 2026-08-08 session): "higher noise moves out onto the flats" -- refined
% in-session to "lower alpha means sigma has more energy at higher
% frequencies", i.e. the forward map's local steepness is jointly governed
% by (alpha, sigma), not sigma alone.
%
%   1. Fraser added to the registry (constellationFraser_v001.mat,
%      noiseCharacterisation_fraser.mat, importDB_fraser_v001,
%      sigmaToMM = 1/10.41793, Findings #134/#135). Not present in v001/v002.
%   2. Per-trial alpha (ira_alphaMean) and sigma (mm) are now carried into
%      the long-format table and aggregated alongside slope -- v002 computed
%      both internally but never surfaced them. A dataset-level summary
%      (mean slope across pipelines, median alpha, median sigma, rise-
%      invertible fraction) is cross-tabulated against Finding #128's own
%      Part B ratio (loop-closure-internal SEM / grid SEM_predicted,
%      hardcoded from that finding's update note, disclosed as a reuse, not
%      recomputed here) and Spearman-correlated against it.
%
% This does NOT touch the loop-closure replicate data itself (betaRecSlice)
% -- it stays on the static monotonicSegments_v2_002.mat grid, same as
% v001/v002. Whether the mechanism is "local slope magnitude" or "proximity
% to the rise/desc segment boundary" (branch-switching across N_REPS
% replicates) is NOT distinguished by this script -- flagged as a follow-up,
% not answered here.
%
% Fraser (2026) -- SEM-agreement-Part-B diagnostic. NOT preregistered.

clearvars -except ans
clc;

srcDir = fileparts(mfilename('fullpath'));
resDir = fullfile(srcDir, 'results');
if ~exist(resDir, 'dir'), mkdir(resDir); end
addpath(genpath(fullfile(srcDir, 'functions')));

fprintf('==============================================\n');
fprintf('  queryMonotonicSlopeAtEmpirical_v003\n');
fprintf('  (Fraser added; alpha/sigma reported alongside slope;\n');
fprintf('   cross-tabulated against Finding #128 Part B ratios)\n');
fprintf('==============================================\n\n');

%% K_conv -- identical derivation to v002, self-contained
BETA_REF = 1/3;
shapeFile = fullfile(srcDir, 'functions', 'baselineShp6_120Hz.mat');
if ~isfile(shapeFile)
    shapeFile = fullfile(srcDir, 'functions', 'baselineShp6_60Hz.mat');
end
if ~isfile(shapeFile)
    error('queryMonoSlope3:NoShapeFile', '%s', ...
        'baselineShp6_*Hz.mat not found in src/functions/.');
end
S = load(shapeFile, 'pathXYresample', 'K');
perimeter_px = sum(sqrt(sum(diff(S.pathXYresample, 1, 1).^2, 2)));
kappa  = S.K(:);
kappa  = kappa(kappa > 0 & isfinite(kappa));
K_conv = perimeter_px / mean(kappa .^ (-BETA_REF));
fprintf('K_conv = %.4f (VGF / f0[Hz])\n\n', K_conv);

%% Load the monotonic segments bundle
monoFile = fullfile(srcDir, 'monotonicSegments_v2_002.mat');
if ~isfile(monoFile)
    error('queryMonoSlope3:MissingBundle', '%s', ...
        'monotonicSegments_v2_002.mat not found.');
end
monoBundle = load(monoFile, 'params', 'pipeOrder', ...
    'alphaGrid', 'sigmaGrid', 'fsGrid', 'betaGenGrid', 'VGFGrid', 'betaSeg');
nP = numel(monoBundle.pipeOrder);
f0GridLo = min(monoBundle.VGFGrid) / K_conv;
f0GridHi = max(monoBundle.VGFGrid) / K_conv;
fprintf('VGF grid [%.2f, %.2f] -> f0 grid [%.4f, %.4f] Hz\n\n', ...
    min(monoBundle.VGFGrid), max(monoBundle.VGFGrid), f0GridLo, f0GridHi);

%% Registry -- v002's six, plus Fraser (Findings #134/#135)
opts.VGFStrategy   = "conditional_median";
opts.BetaStrategy  = "conditional_median";
opts.VGFConversion = @(b, s) s^(1 - b);   % OLD method, kept for parity/diagnostics only
opts.AlphaSource   = "ira_alphaMean";

registry = { ...
    struct('name', "Fraser", ...
        'constellationMat', fullfile(srcDir, "constellationFraser_v001.mat"), ...
        'noiseMat',         fullfile(srcDir, "noiseCharacterisation_fraser.mat"), ...
        'importerHandle',   @importDB_fraser_v001, ...
        'importerArgs',     {{}}, ...
        'sigmaToMM',        1/10.41793); ...
    struct('name', "Zarandi", ...
        'constellationMat', fullfile(srcDir, "constellationZarandi_v001.mat"), ...
        'noiseMat',         fullfile(srcDir, "noiseCharacterisation_zarandi.mat"), ...
        'importerHandle',   @importDB_zarandi_v001, ...
        'importerArgs',     {{ "Verbose", false }}, ...
        'sigmaToMM',        10.0); ...
    struct('name', "Cook_CTRL", ...
        'constellationMat', fullfile(srcDir, "constellationCook_v001.mat"), ...
        'noiseMat',         fullfile(srcDir, "noiseCharacterisation_cook.mat"), ...
        'importerHandle',   @importDB_cook_v002, ...
        'importerArgs',     {{ "Group", "CTRL", "Tasks", 7, "Verbose", false }}, ...
        'sigmaToMM',        0.248); ...
    struct('name', "Cook_ASD", ...
        'constellationMat', fullfile(srcDir, "constellationCookASD_v001.mat"), ...
        'noiseMat',         fullfile(srcDir, "noiseCharacterisation_cookASD.mat"), ...
        'importerHandle',   @importDB_cook_v002, ...
        'importerArgs',     {{ "Group", "ASD", "Tasks", 7, "Verbose", false }}, ...
        'sigmaToMM',        0.248); ...
    struct('name', "Dhieb", ...
        'constellationMat', fullfile(srcDir, "constellationDhieb_v001.mat"), ...
        'noiseMat',         fullfile(srcDir, "noiseCharacterisation_dhieb.mat"), ...
        'importerHandle',   @importDB_dhieb_v001, ...
        'importerArgs',     {{}}, ...
        'sigmaToMM',        0.1478); ...
    struct('name', "Hickman_PLAC", ...
        'constellationMat', fullfile(srcDir, "constellationHickmanPLAC_v001.mat"), ...
        'noiseMat',         fullfile(srcDir, "noiseCharacterisation_hickmanPLAC.mat"), ...
        'importerHandle',   @importDB_hickman_v002, ...
        'importerArgs',     {{ "Study", 2, "Group", "PLAC", "Shapes", 3, "Verbose", false }}, ...
        'sigmaToMM',        0.248); ...
    struct('name', "Hickman_HALO", ...
        'constellationMat', fullfile(srcDir, "constellationHickmanHALO_v001.mat"), ...
        'noiseMat',         fullfile(srcDir, "noiseCharacterisation_hickmanHALO.mat"), ...
        'importerHandle',   @importDB_hickman_v002, ...
        'importerArgs',     {{ "Study", 2, "Group", "HALO", "Shapes", 3, "Verbose", false }}, ...
        'sigmaToMM',        0.248); ...
    };
nSets = numel(registry);

%% Per-dataset, per-trial
allLong = cell(nSets, 1);

for k = 1:nSets
    spec = registry{k};
    fprintf('--- Dataset %d/%d: %s ---\n', k, nSets, spec.name);

    if ~isfile(spec.constellationMat) || ~isfile(spec.noiseMat)
        warning('queryMonoSlope3:MissingFile', ...
            'Constellation or noise mat missing for %s, skipping.', spec.name);
        continue;
    end
    emp = load(spec.constellationMat, 'betaCanon', 'vgfAll', 'canonLabels');
    nz  = load(spec.noiseMat, 'bioResults');
    bio = nz.bioResults;
    nTri = size(emp.betaCanon, 1);

    if ~ismember('f0', bio.Properties.VariableNames)
        error('queryMonoSlope3:NoF0Column', ...
            'bioResults for %s has no f0 column.', spec.name);
    end
    if height(bio) ~= nTri
        error('queryMonoSlope3:RowMismatch', ...
            '%s: bioResults has %d rows but constellation has %d trials.', ...
            spec.name, height(bio), nTri);
    end

    trials = spec.importerHandle(spec.importerArgs{:});
    if numel(trials) ~= nTri
        error('queryMonoSlope3:LengthMismatch', ...
            '%s: importer returned %d trials but constellation has %d.', ...
            spec.name, numel(trials), nTri);
    end
    fsPerTri = arrayfun(@(t) double(t.fs), trials).';

    if ~ismember(opts.AlphaSource, bio.Properties.VariableNames)
        error('queryMonoSlope3:NoAlphaCol', ...
            'AlphaSource ''%s'' not in bioResults for %s.', opts.AlphaSource, spec.name);
    end
    alphaEmp   = double(bio.(opts.AlphaSource));
    sigmaEmpMM = double(bio.sigmaMean) * spec.sigmaToMM;
    f0Emp      = double(bio.f0);

    vgfCanon    = canonicaliseVGF(emp.vgfAll);
    vgfTrueEst  = median(vgfCanon, 2, 'omitnan');
    betaTrueEst = median(emp.betaCanon, 2, 'omitnan');
    vgfViaConversion = NaN(nTri, 1);
    for t = 1:nTri
        if ~isnan(vgfTrueEst(t)) && ~isnan(betaTrueEst(t))
            vgfViaConversion(t) = vgfTrueEst(t) * opts.VGFConversion(betaTrueEst(t), spec.sigmaToMM);
        end
    end
    vgfViaF0 = f0Emp * K_conv;
    nOutNew = sum(vgfViaF0 < min(monoBundle.VGFGrid) | vgfViaF0 > max(monoBundle.VGFGrid));
    fprintf('  alpha median %.3f, sigma median %.3fmm, f0-VGF out-of-grid %d/%d\n', ...
        median(alphaEmp, 'omitnan'), median(sigmaEmpMM, 'omitnan'), nOutNew, nTri);

    dummyBeta = NaN(1, nP);
    longRows  = cell(nTri * nP, 1);
    rowIdx = 0;
    for t = 1:nTri
        result = checkTrialInvertibility_v001( ...
            alphaEmp(t), sigmaEmpMM(t), fsPerTri(t), vgfViaF0(t), ...
            dummyBeta, monoBundle);
        for p = 1:nP
            rowIdx = rowIdx + 1;
            if isnan(result.aIdx)
                longRows{rowIdx} = struct('dataset', string(spec.name), ...
                    'trialIdx', t, 'pipeline', monoBundle.pipeOrder(p), ...
                    'alphaTrial', alphaEmp(t), 'sigmaMMTrial', sigmaEmpMM(t), ...
                    'snapWarn', "coord NaN, not snapped", ...
                    'gRise', NaN, 'gDesc', NaN, 'invRise', false, 'invDesc', false);
                continue;
            end
            a = result.aIdx; s = result.sIdx; f = result.fIdx; v = result.vIdx;
            longRows{rowIdx} = struct('dataset', string(spec.name), ...
                'trialIdx', t, 'pipeline', monoBundle.pipeOrder(p), ...
                'alphaTrial', alphaEmp(t), 'sigmaMMTrial', sigmaEmpMM(t), ...
                'snapWarn', result.snapWarnings, ...
                'gRise', monoBundle.betaSeg.rise.meanSlope(p, a, s, f, v), ...
                'gDesc', monoBundle.betaSeg.desc.meanSlope(p, a, s, f, v), ...
                'invRise', monoBundle.betaSeg.rise.invertible(p, a, s, f, v), ...
                'invDesc', monoBundle.betaSeg.desc.invertible(p, a, s, f, v));
        end
    end
    allLong{k} = struct2table([longRows{:}]);
    fprintf('\n');
end

perTrial = vertcat(allLong{~cellfun(@isempty, allLong)});

%% Aggregate per (dataset, pipeline) -- same shape as v002, plus alpha/sigma
[G, gKeys] = findgroups(perTrial(:, {'dataset','pipeline'}));
perDataPipe = gKeys;
perDataPipe.nTrials      = splitapply(@numel, perTrial.trialIdx, G);
perDataPipe.nInvRise     = splitapply(@(x) sum(x), perTrial.invRise, G);
perDataPipe.medianAlpha  = splitapply(@(x) median(x, 'omitnan'), perTrial.alphaTrial, G);
perDataPipe.medianSigmaMM= splitapply(@(x) median(x, 'omitnan'), perTrial.sigmaMMTrial, G);
perDataPipe.medianGRise  = splitapply(@(x) median(x, 'omitnan'), perTrial.gRise, G);
perDataPipe.medianGDesc  = splitapply(@(x) median(x, 'omitnan'), perTrial.gDesc, G);

fprintf('=== Per (dataset, pipeline): slope, alpha, sigma ===\n');
fprintf('%-14s %-10s %5s %8s %10s %8s %10s %10s\n', ...
    'dataset', 'pipeline', 'nTri', 'medAlpha', 'medSigMM', 'nInvR', 'medGRise', 'medGDesc');
for r = 1:height(perDataPipe)
    fprintf('%-14s %-10s %5d %8.3f %10.3f %8d %10.4f %10.4f\n', ...
        perDataPipe.dataset(r), perDataPipe.pipeline(r), perDataPipe.nTrials(r), ...
        perDataPipe.medianAlpha(r), perDataPipe.medianSigmaMM(r), perDataPipe.nInvRise(r), ...
        perDataPipe.medianGRise(r), perDataPipe.medianGDesc(r));
end

%% Dataset-level summary (mean across pipelines) x Finding #128 Part B ratio
% Part B ratios hardcoded from Finding #128's update note (2026-08-08) and
% original table -- disclosed as reused, not recomputed here.
ratioTable = table( ...
    ["Fraser";"Zarandi";"Cook_CTRL";"Cook_ASD";"Hickman_PLAC";"Hickman_HALO";"Dhieb"], ...
    [7.43; 3.13; 2.98; 3.11; 2.89; 3.96; 0.58], ...
    'VariableNames', {'dataset','semPartBRatio'});

[Gd, dKeys] = findgroups(perDataPipe.dataset);
dsSummary = table(dKeys, 'VariableNames', {'dataset'});
dsSummary.meanGRiseAcrossPipe = splitapply(@(x) mean(x, 'omitnan'), perDataPipe.medianGRise, Gd);
dsSummary.meanGDescAcrossPipe = splitapply(@(x) mean(x, 'omitnan'), perDataPipe.medianGDesc, Gd);
dsSummary.medianAlpha  = splitapply(@(x) mean(x, 'omitnan'), perDataPipe.medianAlpha, Gd);
dsSummary.medianSigmaMM= splitapply(@(x) mean(x, 'omitnan'), perDataPipe.medianSigmaMM, Gd);
dsSummary.invRiseFrac  = splitapply(@(x) mean(x, 'omitnan'), ...
    perDataPipe.nInvRise ./ perDataPipe.nTrials, Gd);

dsSummary = innerjoin(dsSummary, ratioTable, 'Keys', 'dataset');

fprintf('\n=== Dataset-level summary x Finding #128 Part B ratio ===\n');
fprintf('%-14s %8s %8s %8s %8s %10s %8s\n', ...
    'dataset', 'alpha', 'sigMM', 'gRise', 'gDesc', 'invRiseFr', 'ratioB');
for r = 1:height(dsSummary)
    fprintf('%-14s %8.3f %8.3f %8.4f %8.4f %10.3f %8.2f\n', ...
        dsSummary.dataset(r), dsSummary.medianAlpha(r), dsSummary.medianSigmaMM(r), ...
        dsSummary.meanGRiseAcrossPipe(r), dsSummary.meanGDescAcrossPipe(r), ...
        dsSummary.invRiseFrac(r), dsSummary.semPartBRatio(r));
end

if height(dsSummary) >= 4
    [rSlope, pSlope]   = corr(dsSummary.meanGRiseAcrossPipe, dsSummary.semPartBRatio, 'Type', 'Spearman');
    [rAlpha, pAlpha]   = corr(dsSummary.medianAlpha,        dsSummary.semPartBRatio, 'Type', 'Spearman');
    [rSigma, pSigma]   = corr(dsSummary.medianSigmaMM,      dsSummary.semPartBRatio, 'Type', 'Spearman');
    [rInvFr, pInvFr]   = corr(dsSummary.invRiseFrac,        dsSummary.semPartBRatio, 'Type', 'Spearman');
    fprintf('\nSpearman rank correlation with Finding #128 Part B ratio (n=%d datasets):\n', height(dsSummary));
    fprintf('  gRise (rise-branch slope) : r=%+.3f, p=%.4f\n', rSlope, pSlope);
    fprintf('  alpha (noise colour)      : r=%+.3f, p=%.4f\n', rAlpha, pAlpha);
    fprintf('  sigma (noise magnitude mm): r=%+.3f, p=%.4f\n', rSigma, pSigma);
    fprintf('  invRise fraction          : r=%+.3f, p=%.4f\n', rInvFr, pInvFr);
end

%% Save
matOut = fullfile(resDir, 'queryMonotonicSlopeAtEmpirical_v003.mat');
save(matOut, 'perTrial', 'perDataPipe', 'dsSummary', 'ratioTable', 'registry', 'opts', 'monoFile', 'K_conv', '-v7.3');
fprintf('\nSaved: %s\n', matOut);
fprintf('queryMonotonicSlopeAtEmpirical_v003 COMPLETE.\n');

end

% ======================================================================
function vgfCanon = canonicaliseVGF(vgfAll)
nTri   = size(vgfAll, 1);
nDeriv = size(vgfAll, 2);
nReg   = size(vgfAll, 3);
if nReg ~= 3
    error('queryMonoSlope3:VGFAllShape', ...
        'Expected vgfAll [nTri x nDeriv x 3]; got [%d x %d x %d].', nTri, nDeriv, nReg);
end
switch nDeriv
    case 2, canonIdx = [1, 2];
    case 3, canonIdx = [1, 3];
    otherwise
        error('queryMonoSlope3:VGFAllShape', ...
            'Unrecognised vgfAll nDeriv = %d. Expected 2 or 3.', nDeriv);
end
vgfCanon = NaN(nTri, 6);
col = 0;
for d = canonIdx
    for r = 1:3
        col = col + 1;
        vgfCanon(:, col) = vgfAll(:, d, r);
    end
end
end
