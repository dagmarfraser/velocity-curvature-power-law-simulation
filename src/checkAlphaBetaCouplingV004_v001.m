%% checkAlphaBetaCouplingV004_v001.m
%
% Re-verifies the "Link 2: alpha drives beta" finding from the 2026-04-10
% session ("Beta independence from f0 and sigma") against the current
% canonical v004 noise characterisation. That session used
% characteriseBiologicalNoise_v002.m; methods-findings.md's own canonical-
% methods note flags v003-era alpha as stale and v004 (IRASA on
% template-subtracted residuals) as canonical, so v002 -- two generations
% older -- needed re-checking before its correlation coefficients get
% cited anywhere. Only the alpha source changes here; betaCanon itself is
% reused unmodified from the existing constellation mats (beta is
% dimensionless/scale-invariant and was never computed from the noise
% characterisation, so it needs no re-run).
%
% Original finding (v002 alpha, for comparison):
%   Zarandi   SG-IRLS r=-0.03 p=0.71 (null) | BWFD-OLS r=-0.15 p=0.08 (null)
%   Cook(all) SG-IRLS r=+0.60 p=3.8e-20     | BWFD-OLS r=+0.50 p=5.3e-14
%   Dhieb     SG-IRLS r=-0.00 p=0.99 (null) | BWFD-OLS r=-0.16 p=0.22 (null)
%
% Uses ira_alphaMean (IRASA on template-subtracted residuals) as the
% canonical v004 alpha column -- verified present alongside the older
% alphaMean column in bioResults, not assumed.
%
% Created 2026-09-21. Dagmar Scott Fraser, d.s.fraser@bham.ac.uk

%% ===================== CONFIG ==============================
CONFIG.ProjectRoot = fileparts(fileparts(mfilename('fullpath'))); % assumes this file lives in src/
CONFIG.AlphaCol    = 'ira_alphaMean';
CONFIG.Datasets = struct('name', {}, 'constellationMat', {}, 'noiseMat', {});
CONFIG.Datasets(end+1) = struct('name', 'Zarandi', ...
    'constellationMat', 'constellationZarandi_v001.mat', 'noiseMat', 'noiseCharacterisation_zarandi.mat');
CONFIG.Datasets(end+1) = struct('name', 'Cook_CTRL', ...
    'constellationMat', 'constellationCook_v001.mat', 'noiseMat', 'noiseCharacterisation_cook.mat');
CONFIG.Datasets(end+1) = struct('name', 'Cook_ASD', ...
    'constellationMat', 'constellationCookASD_v001.mat', 'noiseMat', 'noiseCharacterisation_cookASD.mat');
CONFIG.Datasets(end+1) = struct('name', 'Dhieb', ...
    'constellationMat', 'constellationDhieb_v001.mat', 'noiseMat', 'noiseCharacterisation_dhieb.mat');
CONFIG.Datasets(end+1) = struct('name', 'Fraser', ...
    'constellationMat', 'constellationFraser_v001.mat', 'noiseMat', 'noiseCharacterisation_fraser.mat');
CONFIG.Datasets(end+1) = struct('name', 'Hickman_PLAC', ...
    'constellationMat', 'constellationHickmanPLAC_v001.mat', 'noiseMat', 'noiseCharacterisation_hickmanPLAC.mat');
CONFIG.Datasets(end+1) = struct('name', 'Hickman_HALO', ...
    'constellationMat', 'constellationHickmanHALO_v001.mat', 'noiseMat', 'noiseCharacterisation_hickmanHALO.mat');
%% ============================================================

srcDir = fullfile(CONFIG.ProjectRoot, 'src');
results = struct('dataset', {}, 'nTrials', {}, 'pipeline', {}, 'r', {}, 'p', {});

for i = 1:numel(CONFIG.Datasets)
    d = CONFIG.Datasets(i);
    cFile = fullfile(srcDir, d.constellationMat);
    nFile = fullfile(srcDir, d.noiseMat);
    if ~isfile(cFile)
        error('checkAlphaBetaCouplingV004_v001:NoConstellation', '%s', ...
            sprintf('Constellation mat not found: %s', cFile));
    end
    if ~isfile(nFile)
        error('checkAlphaBetaCouplingV004_v001:NoNoiseFile', '%s', ...
            sprintf('Noise file not found: %s', nFile));
    end

    C = load(cFile, 'betaCanon', 'canonLabels');
    N = load(nFile, 'bioResults');

    expectedLabels = ["BWFD-OLS","BWFD-LMLS","BWFD-IRLS","SG-OLS","SG-LMLS","SG-IRLS"];
    if ~isequal(string(C.canonLabels(:)), expectedLabels(:))
        error('checkAlphaBetaCouplingV004_v001:UnexpectedLabelOrder', '%s', ...
            sprintf('%s: canonLabels order does not match the assumed [%s]. Got [%s]. Do not proceed on a guess.', ...
            d.name, strjoin(expectedLabels, ','), strjoin(string(C.canonLabels(:)), ',')));
    end
    if ~ismember(CONFIG.AlphaCol, N.bioResults.Properties.VariableNames)
        error('checkAlphaBetaCouplingV004_v001:NoAlphaCol', '%s', ...
            sprintf('%s: bioResults lacks column "%s".', d.name, CONFIG.AlphaCol));
    end
    if height(N.bioResults) ~= size(C.betaCanon, 1)
        error('checkAlphaBetaCouplingV004_v001:RowMismatch', '%s', ...
            sprintf('%s: bioResults has %d rows but betaCanon has %d.', ...
            d.name, height(N.bioResults), size(C.betaCanon, 1)));
    end

    alphaV004 = N.bioResults.(CONFIG.AlphaCol);
    nTri = numel(alphaV004);

    for pipeIdx = [1, 6]  % BWFD-OLS, SG-IRLS -- matching the original two-pipeline comparison
        betaCol = C.betaCanon(:, pipeIdx);
        valid = isfinite(alphaV004) & isfinite(betaCol);
        [r, p] = corr(alphaV004(valid), betaCol(valid), 'Type', 'Pearson');
        results(end+1) = struct('dataset', d.name, 'nTrials', sum(valid), ...
            'pipeline', expectedLabels(pipeIdx), 'r', r, 'p', p); %#ok<AGROW>
    end
    fprintf('%-10s: %d trials, alpha (%s) median=%.3f [%.3f, %.3f]\n', ...
        d.name, nTri, CONFIG.AlphaCol, median(alphaV004, 'omitnan'), ...
        min(alphaV004), max(alphaV004));
end

resultsTable = struct2table(results);
fprintf('\n=== Alpha (v004, %s) vs Beta correlation, re-verified ===\n', CONFIG.AlphaCol);
fprintf('%-10s %-10s %6s %8s %10s\n', 'dataset', 'pipeline', 'n', 'r', 'p');
for i = 1:height(resultsTable)
    fprintf('%-10s %-10s %6d %8.3f %10.3g\n', string(resultsTable.dataset(i)), string(resultsTable.pipeline(i)), ...
        resultsTable.nTrials(i), resultsTable.r(i), resultsTable.p(i));
end

resDir = fullfile(srcDir, 'results');
if ~isfolder(resDir), mkdir(resDir); end
outFile = fullfile(resDir, 'checkAlphaBetaCouplingV004_v001.mat');
save(outFile, 'resultsTable', 'CONFIG', '-v7.3');
fprintf('\nSaved: %s\n', outFile);
