%% crossReferenceSimexMonotonicity_v001.m
%
% Joins a SIMEX corpus (runSimexBetaCorpus_v001/v002 output) against the
% pre-existing per-trial monotonicity classification in
% results/queryMonotonicSlopeAtEmpirical_v003.mat, and reports the same
% statistics Finding #219 used for Hickman: mean betaSimexCI width by
% invertibility class, Cohen's d, Welch's t-test, Wilcoxon rank-sum, and
% the two-way percentile overlap.
%
% GENERALISATION NOTE: queryMonotonicSlopeAtEmpirical_v003's own registry
% already covers Fraser, Zarandi, Cook_CTRL, Cook_ASD, Dhieb, Hickman_PLAC
% and Hickman_HALO -- so no new simulation is needed to extend Finding
% #219's cross-reference beyond Hickman, only this join.
%
% JOIN MECHANICS: perTrial (in queryMonotonicSlopeAtEmpirical_v003.mat)
% keys each row by (dataset, trialIdx, pipeline), where trialIdx is a
% plain 1:nTrial position, not a trialID string -- it was never carried
% through. This script rebuilds the trialIdx -> trialID mapping by
% calling each dataset's own importer with EXACTLY the args
% queryMonotonicSlopeAtEmpirical_v003's registry used (copied verbatim
% below), which is deterministic and therefore reproduces the same
% trial order. That mapping is then used to attach trialID to perTrial,
% which is joined against the SIMEX corpus's summaryTable on
% (trialID, pipeline).
%
% INVERTIBLE-DEFINITION CAVEAT: Finding #218/#219's own Hickman join was
% done ad hoc (evaluate_matlab_code, not saved as a script), so the exact
% rule used to collapse perTrial's invRise/invDesc into a single
% invertible/non-invertible split for that cross-reference is not
% independently documented anywhere. This script reports all three
% candidates (invRise alone, invDesc alone, invRise|invDesc) rather than
% silently picking one -- run this on Hickman FIRST and check which
% candidate reproduces Finding #219's own numbers (d=0.27, mean CI width
% 0.022 vs 0.017, 38%/66% overlap) before trusting any candidate's
% reading on a new dataset.
%
% Created 2026-09-21. Dagmar Scott Fraser, d.s.fraser@bham.ac.uk

%% ===================== CONFIG ==============================
CONFIG.ProjectRoot     = fileparts(fileparts(mfilename('fullpath'))); % assumes this file lives in src/
CONFIG.SimexCorpusFile = fullfile(CONFIG.ProjectRoot, 'results', 'simexBetaCorpus_hickmanPLAC-hickmanHALO_20260921_123014.mat');
CONFIG.MonotonicityFile = fullfile(CONFIG.ProjectRoot, 'src', 'results', ...
    'queryMonotonicSlopeAtEmpirical_v003.mat');
%% ============================================================

if isempty(CONFIG.SimexCorpusFile)
    error('crossReferenceSimexMonotonicity_v001:NoCorpusFile', '%s', ...
        'Set CONFIG.SimexCorpusFile to a runSimexBetaCorpus_v00x output mat before running.');
end

addpath(fullfile(CONFIG.ProjectRoot, 'src', 'functions'));

%% ----- Registry: dataset/group -> perTrial dataset name + noise mat -----
% NOTE (corrected from first draft): trialIdx in perTrial is keyed to
% bioResults' own row order in queryMonotonicSlopeAtEmpirical_v003 (that
% script indexes alphaEmp(t)/sigmaEmpMM(t) directly from bio row t) --
% NOT to a fresh call of the importer. Re-deriving trialIdx->trialID via
% a new importer call would silently assume importer-call order matches
% bioResults row order, which is an extra assumption the original script
% itself never makes. Using bioResults.trialID directly removes that
% assumption entirely.
registry = struct('simexDataset', {}, 'simexGroup', {}, 'monoName', {}, 'noiseFile', {});

registry(end+1) = struct('simexDataset', 'fraser', 'simexGroup', '', ...
    'monoName', "Fraser", 'noiseFile', 'noiseCharacterisation_fraser.mat');
registry(end+1) = struct('simexDataset', 'zarandi', 'simexGroup', '', ...
    'monoName', "Zarandi", 'noiseFile', 'noiseCharacterisation_zarandi.mat');
registry(end+1) = struct('simexDataset', 'cook', 'simexGroup', 'CTRL', ...
    'monoName', "Cook_CTRL", 'noiseFile', 'noiseCharacterisation_cook.mat');
registry(end+1) = struct('simexDataset', 'cook', 'simexGroup', 'ASD', ...
    'monoName', "Cook_ASD", 'noiseFile', 'noiseCharacterisation_cookASD.mat');
registry(end+1) = struct('simexDataset', 'dhieb', 'simexGroup', '', ...
    'monoName', "Dhieb", 'noiseFile', 'noiseCharacterisation_dhieb.mat');
registry(end+1) = struct('simexDataset', 'hickman', 'simexGroup', 'PLAC', ...
    'monoName', "Hickman_PLAC", 'noiseFile', 'noiseCharacterisation_hickmanPLAC.mat');
registry(end+1) = struct('simexDataset', 'hickman', 'simexGroup', 'HALO', ...
    'monoName', "Hickman_HALO", 'noiseFile', 'noiseCharacterisation_hickmanHALO.mat');

%% ----- Load both sources -----
if ~isfile(CONFIG.MonotonicityFile)
    error('crossReferenceSimexMonotonicity_v001:NoMonoFile', '%s', ...
        sprintf('Monotonicity file not found: %s. Run queryMonotonicSlopeAtEmpirical_v003 first.', ...
        CONFIG.MonotonicityFile));
end
M = load(CONFIG.MonotonicityFile, 'perTrial');
perTrial = M.perTrial;

if ~isfile(CONFIG.SimexCorpusFile)
    error('crossReferenceSimexMonotonicity_v001:NoSimexFile', '%s', ...
        sprintf('SIMEX corpus file not found: %s', CONFIG.SimexCorpusFile));
end
Sx = load(CONFIG.SimexCorpusFile, 'summaryTable');
simexTable = Sx.summaryTable;
if ~ismember('dataset', simexTable.Properties.VariableNames)
    % v001 (Hickman-only) predates the 'dataset'/'group' columns v002 added.
    simexTable.dataset = repmat({'hickman'}, height(simexTable), 1);
end
if ~ismember('group', simexTable.Properties.VariableNames)
    simexTable.group = simexTable.group; %#ok<NASGU> % placeholder; real v001 files already carry .group
end

%% ----- Which (dataset,group) pairs are actually present in this corpus -----
pairs = unique(simexTable(:, {'dataset','group'}), 'rows');
fprintf('SIMEX corpus %s carries %d (dataset,group) pair(s):\n', CONFIG.SimexCorpusFile, height(pairs));
for r = 1:height(pairs)
    fprintf('  %s / %s\n', string(pairs.dataset(r)), string(pairs.group(r)));
end

%% ----- Build trialIdx -> trialID map per pair, join, accumulate -----
joinedRows = {};

for r = 1:height(pairs)
    dsKey = lower(strtrim(string(pairs.dataset(r))));
    grKey = lower(strtrim(string(pairs.group(r))));

    regMatch = find(strcmpi({registry.simexDataset}, dsKey) & ...
        strcmpi(string({registry.simexGroup}), grKey));
    if isempty(regMatch)
        warning('crossReferenceSimexMonotonicity_v001:NoRegistryMatch', ...
            'No registry entry for dataset="%s" group="%s"; skipping (%d SIMEX rows unmatched).', ...
            dsKey, grKey, sum(strcmpi(string(simexTable.dataset), dsKey) & ...
                strcmpi(string(simexTable.group), grKey)));
        continue
    end
    spec = registry(regMatch);

    dsNoiseFile = fullfile(CONFIG.ProjectRoot, 'src', spec.noiseFile);
    if ~isfile(dsNoiseFile)
        error('crossReferenceSimexMonotonicity_v001:NoiseFileMissing', '%s', ...
            sprintf('Expected noise file not found: %s', dsNoiseFile));
    end
    Sn = load(dsNoiseFile, 'bioResults');
    nTri = height(Sn.bioResults);
    trialIdxToID = string(Sn.bioResults.trialID);  % row order IS trialIdx, per queryMonotonicSlopeAtEmpirical_v003's own indexing

    monoRows = perTrial(strcmp(string(perTrial.dataset), spec.monoName), :);
    if isempty(monoRows)
        warning('crossReferenceSimexMonotonicity_v001:NoMonoRows', ...
            'No perTrial rows for monoName="%s"; skipping.', spec.monoName);
        continue
    end
    if max(monoRows.trialIdx) > nTri
        error('crossReferenceSimexMonotonicity_v001:TrialCountMismatch', '%s', ...
            sprintf(['%s: perTrial max trialIdx=%d exceeds importer''s own trial count ' ...
                     '(%d). Registry importer args have drifted from what ' ...
                     'queryMonotonicSlopeAtEmpirical_v003 originally used -- do not ' ...
                     'proceed on a guess.'], spec.monoName, max(monoRows.trialIdx), nTri));
    end
    monoRows.trialID = trialIdxToID(monoRows.trialIdx);

    simexSub = simexTable(strcmpi(string(simexTable.dataset), dsKey) & ...
        strcmpi(string(simexTable.group), grKey) & simexTable.success, :);

    joined = innerjoin(monoRows, simexSub, ...
        'Keys', {'trialID', 'pipeline'}, ...
        'RightVariables', {'betaObserved', 'betaSimex', 'betaSimexCI_lo', 'betaSimexCI_hi'});
    joined.ciWidth = joined.betaSimexCI_hi - joined.betaSimexCI_lo;
    joined.simexDataset = repmat(dsKey, height(joined), 1);
    joined.simexGroup   = repmat(grKey, height(joined), 1);

    fprintf('%s/%s: %d trials imported, %d perTrial rows, %d joined to SIMEX successes\n', ...
        dsKey, grKey, nTri, height(monoRows), height(joined));

    joinedRows{end+1} = joined; %#ok<AGROW>
end

if isempty(joinedRows)
    error('crossReferenceSimexMonotonicity_v001:NothingJoined', '%s', ...
        'No (dataset,group) pair in the SIMEX corpus matched the monotonicity registry.');
end
allJoined = vertcat(joinedRows{:});

%% ----- Report, for all three candidate invertible-definitions -----
candidates = struct( ...
    'label',   {'invRise only', 'invDesc only', 'invRise | invDesc'}, ...
    'invMask', { ...
        @(t) t.invRise, ...
        @(t) t.invDesc, ...
        @(t) t.invRise | t.invDesc});

for c = 1:numel(candidates)
    fprintf('\n=== Candidate: %s ===\n', candidates(c).label);
    inv = candidates(c).invMask(allJoined);

    ciInv    = allJoined.ciWidth(inv);
    ciNonInv = allJoined.ciWidth(~inv);
    ciInv    = ciInv(isfinite(ciInv));
    ciNonInv = ciNonInv(isfinite(ciNonInv));

    if isempty(ciInv) || isempty(ciNonInv)
        fprintf('  One group is empty (invertible n=%d, non-invertible n=%d) -- skipping stats.\n', ...
            numel(ciInv), numel(ciNonInv));
        continue
    end

    meanInv    = mean(ciInv);
    meanNonInv = mean(ciNonInv);
    pooledSD   = sqrt(((numel(ciInv)-1)*var(ciInv) + (numel(ciNonInv)-1)*var(ciNonInv)) / ...
        (numel(ciInv) + numel(ciNonInv) - 2));
    cohensD    = (meanNonInv - meanInv) / pooledSD;

    [~, pT]  = ttest2(ciNonInv, ciInv, 'Vartype', 'unequal');
    pW       = ranksum(ciNonInv, ciInv);

    medNonInv = median(ciNonInv);
    medInv    = median(ciInv);
    pctInvExceedsNonInvMedian    = 100 * mean(ciInv    > medNonInv);
    pctNonInvExceedsInvMedian    = 100 * mean(ciNonInv > medInv);

    fprintf('  n invertible=%d, n non-invertible=%d\n', numel(ciInv), numel(ciNonInv));
    fprintf('  mean CI width: invertible=%.4f, non-invertible=%.4f\n', meanInv, meanNonInv);
    fprintf('  Cohen''s d = %.3f  (Welch t p=%.3g, Wilcoxon p=%.3g)\n', cohensD, pT, pW);
    fprintf('  overlap: %.0f%% of invertible exceed non-invertible median; %.0f%% of non-invertible exceed invertible median\n', ...
        pctInvExceedsNonInvMedian, pctNonInvExceedsInvMedian);
end

%% ----- Save -----
resDir = fullfile(CONFIG.ProjectRoot, 'src', 'results');
if ~isfolder(resDir), mkdir(resDir); end
[~, corpusStem] = fileparts(CONFIG.SimexCorpusFile);
outFile = fullfile(resDir, sprintf('crossRefSimexMonotonicity_%s.mat', corpusStem));
save(outFile, 'allJoined', 'CONFIG', '-v7.3');
fprintf('\nSaved: %s\n', outFile);
