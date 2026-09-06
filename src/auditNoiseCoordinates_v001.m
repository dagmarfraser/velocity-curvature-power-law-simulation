% AUDITNOISECOORDINATES_V001  Recompute every dataset's canonical (alpha, sigma)
% from source, and reconcile against the documented millimetre constants.
%
% WHY THIS EXISTS
%   Session 103's audit (2026-09-02) established the noise coordinates in
%   docs/EMPIRICAL_DATASETS.md from console code rather than a named script,
%   which is the exact "number with no traceable producer" problem that audit
%   was written to find. This script is that producer. Every figure in
%   EMPIRICAL_DATASETS.md Table 1 and in the paper's Section 2.3 should be
%   reproducible by running this.
%
% WHAT IT DOES
%   1. Loads src/noiseCharacterisation_<tag>.mat for each dataset.
%   2. Reports mean/median ira_alphaMean and sigmaMean over trials with both
%      finite (the same finite mask the downstream pipeline uses).
%   3. Applies each dataset's OWN SigmaToMM, read from its importDB_* function
%      rather than hardcoded here, and flags any dataset whose constant cannot
%      be found rather than substituting a default.
%
% KEY RESULT THIS GUARDS (Session 103)
%   Fraser and Pilot are distinct corpora on the same physical iPad Pro 13" M4
%   and use DELIBERATELY DIFFERENT constants: Fraser PX_PER_MM = 10.41793,
%   Pilot PX_PER_MM = 9.73. importDB_fraser_v001.m's own header says "Do NOT
%   use Pilot's 9.73 here". Their alphas agree to four significant figures
%   (4.2889 vs 4.2891) across independent corpora, which is a property of the
%   instrument and protocol, not a copying error. Their sigmas do not.
%
% USAGE
%   auditNoiseCoordinates_v001            % console table only, no files written
%
% Fraser, D.S. (2026)  v001

clearvars

%% CONFIG -----------------------------------------------------------------
CFG.srcDir       = fileparts(mfilename('fullpath'));

% {mat tag, importer stem}. Mat tags are the real filenames on disk, which do
% NOT all match their importer's name: cook.mat is Cook CTRL and cookASD.mat is
% the ASD arm, both imported by importDB_cook; the two Hickman arms share
% importDB_hickman. Listing the pair explicitly stops the importer lookup from
% silently failing and returning NaN for the arms.
CFG.datasets = { ...
    'fraser',       'fraser'; ...
    'pilot',        'pilot'; ...
    'zarandi',      'zarandi'; ...
    'cook',         'cook'; ...        % Cook CTRL
    'cookASD',      'cook'; ...        % Cook ASD, same rig and constant
    'hickmanPLAC',  'hickman'; ...
    'hickmanHALO',  'hickman'; ...
    'dhieb',        'dhieb'};

% Dagenais has no noiseCharacterisation mat: its coordinate (alpha 2.867,
% sigma 1.490 mm, Qualisys native mm so SigmaToMM = 1.0) comes from Finding #6's
% median-centroid table, not from this pipeline. It is deliberately absent here
% rather than faked, and hickmanPDM0/PDM1 are 1x1 stubs, also excluded.
% petersCONG is excluded too: Peters is out of corpus (EPP contamination,
% Section 2.6) and no constant is claimed for it anywhere.

% Disclosed overrides. The automatic lookup reads the constant out of the
% importer, which is the point: it cannot drift from the import path. But
% importDB_hickman_v003.m passes 0.248 POSITIONALLY at line 238 rather than
% assigning a named constant, so there is nothing for the parser to find.
% Rather than loosen the parser into something brittle, the value is declared
% here with its provenance and the script says out loud when it uses one.
% If you add a named constant to the Hickman importer, delete this entry.
CFG.constantOverride = { ...
    'hickman', 0.248, 'importDB_hickman_v003.m line 238, passed positionally (WACOM px, same rig as Cook)'};

CFG.alphaField   = 'ira_alphaMean';   % residual IRASA, the v004 canonical column
CFG.sigmaField   = 'sigmaMean';
CFG.verbose      = true;

if isempty(CFG.srcDir), CFG.srcDir = pwd; end
cd(CFG.srcDir);

fprintf('=== AUDIT: canonical noise coordinates from source ===\n');
fprintf('    run %s\n\n', datetime('now','Format','yyyy-MM-dd HH:mm'));

%% ------------------------------------------------------------------------
nDS = size(CFG.datasets,1);
rows = cell(nDS,1);

for i = 1:nDS
    tag    = CFG.datasets{i,1};
    impTag = CFG.datasets{i,2};
    matFile = fullfile(CFG.srcDir, sprintf('noiseCharacterisation_%s.mat', tag));

    if ~isfile(matFile)
        error('auditNoiseCoordinates:missingMat', '%s', sprintf( ...
            ['FAILED PATH: %s not found. This script''s dataset list is ' ...
             'explicit; a missing mat means the list and the disk have ' ...
             'diverged, which must be resolved rather than skipped.'], matFile));
    end

    S = load(matFile);
    fn = fieldnames(S);
    T  = S.(fn{1});
    if ~istable(T)
        error('auditNoiseCoordinates:notTable', '%s', ...
            sprintf('%s does not contain a table as its first variable.', matFile));
    end
    if ~all(ismember({CFG.alphaField, CFG.sigmaField}, T.Properties.VariableNames))
        error('auditNoiseCoordinates:missingColumn', '%s', sprintf( ...
            '%s lacks %s and/or %s.', matFile, CFG.alphaField, CFG.sigmaField));
    end

    a = T.(CFG.alphaField);
    s = T.(CFG.sigmaField);
    ok = isfinite(a) & isfinite(s);

    [k, kSrc] = localFindSigmaToMM(CFG.srcDir, impTag);

    if isnan(k)
        ov = find(strcmp(CFG.constantOverride(:,1), impTag), 1);
        if ~isempty(ov)
            k    = CFG.constantOverride{ov,2};
            kSrc = sprintf('OVERRIDE %.6f -- %s', k, CFG.constantOverride{ov,3});
            fprintf(['  NOTE: %s uses a declared override for SigmaToMM; the importer ' ...
                     'does not assign a named constant.\n'], tag);
        end
    end

    rows{i} = struct( ...
        'tag',        tag, ...
        'nTrials',    height(T), ...
        'nFinite',    sum(ok), ...
        'alphaMean',  mean(a(ok)), ...
        'alphaMedian',median(a(ok)), ...
        'sigmaNative',mean(s(ok)), ...
        'sigmaToMM',  k, ...
        'sigmaMM',    mean(s(ok)) * k, ...
        'constSource',kSrc);
end

rows = rows(~cellfun(@isempty, rows));

%% REPORT -----------------------------------------------------------------
fprintf('%-12s %7s %7s %10s %10s %12s %11s %9s\n', ...
    'dataset','nTrial','nFin','alphaMean','alphaMed','sigmaNative','SigmaToMM','sigma_mm');
fprintf('%s\n', repmat('-',1,86));
for i = 1:numel(rows)
    r = rows{i};
    if isnan(r.sigmaToMM)
        fprintf('%-12s %7d %7d %10.4f %10.4f %12.4f %11s %9s   <-- CONSTANT NOT FOUND\n', ...
            r.tag, r.nTrials, r.nFinite, r.alphaMean, r.alphaMedian, r.sigmaNative, 'NaN', 'NaN');
    else
        fprintf('%-12s %7d %7d %10.4f %10.4f %12.4f %11.6f %9.4f\n', ...
            r.tag, r.nTrials, r.nFinite, r.alphaMean, r.alphaMedian, ...
            r.sigmaNative, r.sigmaToMM, r.sigmaMM);
    end
end

fprintf('\nConstant provenance (read from importer source, never hardcoded here):\n');
for i = 1:numel(rows)
    fprintf('  %-12s %s\n', rows{i}.tag, rows{i}.constSource);
end

%% Fraser/Pilot guard -----------------------------------------------------
iF = find(strcmp(cellfun(@(r) r.tag, rows, 'uni', 0), 'fraser'), 1);
iP = find(strcmp(cellfun(@(r) r.tag, rows, 'uni', 0), 'pilot'),  1);
if ~isempty(iF) && ~isempty(iP)
    dA = abs(rows{iF}.alphaMean - rows{iP}.alphaMean);
    fprintf('\nFraser/Pilot alpha agreement : |%.4f - %.4f| = %.5f\n', ...
        rows{iF}.alphaMean, rows{iP}.alphaMean, dA);
    fprintf('Fraser/Pilot sigma (mm)      : %.4f vs %.4f  (constants %.6f vs %.6f)\n', ...
        rows{iF}.sigmaMM, rows{iP}.sigmaMM, rows{iF}.sigmaToMM, rows{iP}.sigmaToMM);
    if abs(rows{iF}.sigmaToMM - rows{iP}.sigmaToMM) < 1e-9
        warning('auditNoiseCoordinates:sharedConstant', '%s', ...
            ['Fraser and Pilot resolved to the SAME SigmaToMM. They are ' ...
             'documented as deliberately different (10.41793 vs 9.73 px/mm). ' ...
             'Check the importers before trusting these sigmas.']);
    end
end

fprintf('\nDone. No files written.\n');

%% ------------------------------------------------------------------------
function [k, src] = localFindSigmaToMM(srcDir, tag)
% Read the millimetre constant out of the dataset's own importer rather than
% duplicating it here, so this script cannot drift from the import path.
    k = NaN; src = 'not found';
    d = dir(fullfile(srcDir, 'functions', sprintf('importDB_%s*.m', tag)));
    if isempty(d)
        d = dir(fullfile(srcDir, sprintf('importDB_%s*.m', tag)));
    end
    if isempty(d), return, end

    [~, ord] = sort({d.name});          % highest version last
    d = d(ord);
    fp = fullfile(d(end).folder, d(end).name);
    L = splitlines(fileread(fp));

    for n = 1:numel(L)
        t = strtrim(L{n});
        if startsWith(t, '%'), continue, end

        tok = regexp(t, 'PX_PER_MM\s*=\s*([0-9.]+)', 'tokens', 'once');
        if ~isempty(tok)
            k = 1 / str2double(tok{1});
            src = sprintf('%s line %d: PX_PER_MM = %s', d(end).name, n, tok{1});
            return
        end
        tok = regexp(t, '(?:SIGMA_TO_MM|SigmaToMM)\s*=\s*([0-9.]+)', 'tokens', 'once');
        if ~isempty(tok)
            k = str2double(tok{1});
            src = sprintf('%s line %d: SigmaToMM = %s', d(end).name, n, tok{1});
            return
        end
    end
end
